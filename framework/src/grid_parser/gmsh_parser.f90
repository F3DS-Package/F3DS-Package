module class_gmsh_parser
    use typedef_module
    use face_type_module
    use stdio_module
    use string_utils_module
    use abstract_grid_parser
    use abstract_configuration

    implicit none

    private

    ! https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
    type, public, extends(grid_parser) :: gmsh_parser
        private

        ! params
        integer(int_kind) :: num_ghost_cells_ = 1
        integer(int_kind) :: sweep_section_limit_ = 2**31

        ! format info
        integer(int_kind) :: mesh_minor_version_ = 0
        integer(int_kind) :: mesh_major_version_ = 0
        integer(int_kind) :: mesh_filetype_      = 0
        integer(int_kind) :: mesh_datasize_      = 0

        ! convert boundary condition id
        integer(int_kind), allocatable :: gmsh_bc_to_f3ds_bc_(:)

        ! nodes
        real(real_kind), allocatable :: node_positions_(:,:)

        contains

        procedure, public, pass(self) :: parse
        procedure, public, pass(self) :: close
        procedure, public, pass(self) :: get_number_of_cells
        procedure, public, pass(self) :: get_number_of_faces
        procedure, public, pass(self) :: get_number_of_ghost_cells
        procedure, public, pass(self) :: get_number_of_boundary_faces
        procedure, public, pass(self) :: get_number_of_points
        procedure, public, pass(self) :: get_cells
        procedure, public, pass(self) :: get_faces
        procedure, public, pass(self) :: get_boundaries
        procedure, public, pass(self) :: get_cell_geometries

        procedure, private, pass(self) :: read_mesh_format
        procedure, private, pass(self) :: read_physical_names
        procedure, private, pass(self) :: read_nodes
        procedure, private, pass(self) :: read_elements
    end type gmsh_parser

    contains

    subroutine parse(self, config)
        class(gmsh_parser  ), intent(inout) :: self
        class(configuration), intent(inout) :: config

        ! for mesh input
        integer  (int_kind)              :: unit_number
        character(len=:   ), allocatable :: mesh_filename
        character(len=20  )              :: section
        integer  (int_kind)              :: io_statiment
        integer  (int_kind)              :: n_line

        ! for config
        character(len=:), allocatable :: error_msg
        logical :: found

        call config%get_char("Grid.Filepath", mesh_filename, found, "grid.msh")
        if(.not. found) call write_warring("'Grid.Filepath' is not found in configuration you set. To be set default value.")

        open(newunit=unit_number, iostat=io_statiment, file=mesh_filename, access='stream', form='formatted', status='old')
        if(io_statiment /= 0) call call_error("Can not read mesh file '"//mesh_filename//"'.")

        do while (.true.)
            read(unit_number, iostat=io_statiment) section

            if(io_statiment == -1) then
                close(unit_number)
                exit
            end if

            if (trim(section) == "$MeshFormat") then
                call self%read_mesh_format(unit_number)
            else if (trim(section) == "$Nodes")
                call self%read_nodes(unit_number)
            else if (trim(section) == "$Elements")
            else if (trim(section) == "$PhysicalNames")
                call self%read_physical_names(unit_number)
            else
                call write_debuginfo("Skip section '"//trim(section)//"'.")
            end if

            do n_line = 1, self%sweep_section_limit_, 1
                read(unit_number, iostat=io_statiment) section
                call write_debuginfo("Sweep line until end section. Now read '"//section//"'.")

                if (index(trim(section), "$End") == 1) exit

                if(io_statiment == -1) then
                    close(unit_number)
                    call call_error("EoF found. Section may be not closed.")
                end if

                if (n_line == self%sweep_section_limit_) then
                    close(unit_number)
                    call call_error("End section is not found.")
                end if
            end do
        end do
    end subroutine parse

    subroutine read_mesh_format(self, unit_number)
        class  (gmsh_parser), intent(inout) :: self
        integer(int_kind   ), intent(inout) :: unit_number

        character(len=20  ) :: mesh_ver_tmp
        integer  (int_kind) :: dot_index, mesh_ver_size

        ! mesh_ver_tmp will be write version information (e.g. '4.1')
        read(unit_number) mesh_ver_tmp, self%mesh_filetype_, self%mesh_datasize_
        dot_index     = index(trim(mesh_ver_tmp), ".")
        mesh_ver_size = size (trim(mesh_ver_tmp))
        read(trim(mesh_ver_tmp)[1          :dot_index-1  ], *) self%mesh_major_version_
        read(trim(mesh_ver_tmp)[dot_index+1:mesh_ver_size], *) self%mesh_minor_version_

        if (self%mesh_major_version_ /= 4 .or. self%mesh_minor_version_ /= 1) then
            call write_debuginfo("gmsh version = "//to_str(self%mesh_major_version_)//"."//to_str(self%mesh_minor_version_))
            close(unit_number)
            call call_error("Invalid mesh version.")
        end if
    end subroutine read_mesh_format

    subroutine read_physical_names(self, unit_number)
        class  (gmsh_parser), intent(inout) :: self
        integer(int_kind   ), intent(inout) :: unit_number

        integer  (int_kind               )              :: i, num_physical_names
        character(len=127                )              :: name
        integer  (int_kind               )              :: dim, tag, max_tag
        integer  (int_kind               ), allocatable :: tags(:)
        integer  (boundary_face_type_kind), allocatable :: bc_ids(:)
        integer  (int_kind               )              :: name_size

        read(unit_number) num_physical_names

        allocate(tags  (num_physical_names))
        allocate(bc_ids(num_physical_names))

        max_tag = 0
        do i = 1, num_physical_names, 1
            read(unit_number) dim, tag, name

            name_size = size(name)
            if (name[2:name_size-1] == "fluid") then
                tags  [i] = tag
                bc_ids[i] = -1
            else
                tags  [i] = tag
                bc_ids[i] = string_to_boundary_face_type(name[2:name_size-1])

                if (bc_ids[i] == unknown_face_type) then
                    close(unit_number)
                    call call_error("Unknown face type '"//name[2:name_size-1]//"' found.")
                end if
            end if

            if (tag > max_tag) max_tag = tag
        end do

        allocate(self%gmsh_bc_to_f3ds_bc_(max_tag))

        do i = 1, num_physical_names, 1
            self%gmsh_bc_to_f3ds_bc_(tags[i]) = bc_ids[i]
        end do
    end subroutine read_physical_names

    subroutine read_nodes(self, unit_number)
        class  (gmsh_parser), intent(inout) :: self
        integer(int_kind   ), intent(inout) :: unit_number

        integer(int_kind) :: entity_index, node_index
        integer(int_kind) :: num_entity_blocks
        integer(int_kind) :: num_nodes
        integer(int_kind) :: min_node_tag, max_node_tag
        integer(int_kind) :: entity_dim, entity_tag, parametric, num_nodes_entity
        integer(int_kind) :: node_tag

        read(unit_number) num_entity_blocks, num_nodes, min_node_tag, max_node_tag

        allocate(self%node_positions_(3, num_nodes))

        ! entity block (We don't use this infomation!!)
        do entity_index = 1, num_entity_blocks, 1
            read(unit_number) entity_dim, entity_tag, parametric, num_nodes_entity
            do node_index = 1, num_nodes_entity, 1
                read(unit_number) node_tag
            end do
        end do

        ! node block
        do node_index = 1, num_nodes, 1
            read(unit_number) self%node_positions_(1, node_index), self%node_positions_(2, node_index), self%node_positions_(3, node_index)
        end do
    end subroutine read_nodes

    subroutine read_elements(self, unit_number)
        class  (gmsh_parser), intent(inout) :: self
        integer(int_kind   ), intent(inout) :: unit_number

        integer(int_kind) :: entity_index, element_index
        integer(int_kind) :: num_entity_blocks, num_elements, min_element_tag, max_element_tag
        integer(int_kind) :: entity_dim, entity_tag, element_tag, num_elements_entity

        read(unit_number) num_entity_blocks, num_elements, min_element_tag, max_element_tag
    end subroutine read_elements
end module class_gmsh_parser