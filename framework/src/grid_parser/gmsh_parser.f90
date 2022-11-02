module class_gmsh_parser
    use typedef_module
    use face_type_module
    use stdio_module
    use string_utils_module
    use class_point_id_list
    use abstract_grid_parser
    use abstract_configuration

    implicit none

    private


    ! https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
    type, public, extends(grid_parser) :: gmsh_parser
        private

        ! params
        integer(int_kind) :: num_ghost_cells_ = 1
        integer(int_kind) :: sweep_section_limit_ = 2**16
        logical           :: parsed_ = .false.

        ! format info
        integer(int_kind) :: mesh_minor_version_ = 0
        integer(int_kind) :: mesh_major_version_ = 0
        integer(int_kind) :: mesh_filetype_      = 0
        integer(int_kind) :: mesh_datasize_      = 0

        ! convert boundary condition id
        integer(int_kind), allocatable :: gmsh_bc2f3ds_bc_(:)

        ! nodes
        integer(int_kind     )       :: number_of_nodes_
        real(real_kind), allocatable :: node_positions_(:,:)

        ! elements; TODO: Consider to convert elements to cells and boundary faces in parse() routine.
        integer(int_kind     )              :: number_of_elements_
        type   (point_id_list), allocatable :: elements_           (:)
        integer(int_kind     ), allocatable :: elements_entity_tag_(:)
        integer(int_kind     ), allocatable :: elements_type_      (:)

        contains

        procedure, public, pass(self) :: parse
        procedure, public, pass(self) :: close
        procedure, public, pass(self) :: get_number_of_cells
        procedure, public, pass(self) :: get_number_of_faces ! TODO: impl this!
        procedure, public, pass(self) :: get_number_of_ghost_cells
        procedure, public, pass(self) :: get_number_of_boundary_faces
        procedure, public, pass(self) :: get_number_of_points
        procedure, public, pass(self) :: get_cells ! TODO: impl this!
        procedure, public, pass(self) :: get_faces ! TODO: impl this!
        procedure, public, pass(self) :: get_boundaries ! TODO: impl this!
        procedure, public, pass(self) :: get_cell_geometries ! TODO: impl this!

        procedure, private, pass(self) :: read_mesh_format
        procedure, private, pass(self) :: read_physical_names
        procedure, private, pass(self) :: read_nodes
        procedure, private, pass(self) :: read_elements
        procedure, private, pass(self) :: supported_elemental_type
        procedure, private, pass(self) :: elemental_type2num_nodes
        procedure, private, pass(self) :: is_cell_element
        procedure, private, pass(self) :: is_face_element
    end type gmsh_parser

    contains

    function get_number_of_points(self) result(n)
        class  (gmsh_parser), intent(in) :: self
        integer(int_kind   )             :: n

        if(.not. self%parsed_)then
            call call_error("'parse' method of gmsh_parser is not called yet. But you call 'get_number_of_points' method.")
        end if

        n = self%number_of_nodes_
    end function get_number_of_ghost_cells

    function get_number_of_boundary_faces(self, type) result(n)
        class  (gmsh_parser            ), intent(in) :: self
        integer(boundary_face_type_kind), intent(in) :: type
        integer(int_kind   )             :: n

        integer(int_kind   )             :: i

        if(.not. self%parsed_)then
            call call_error("'parse' method of gmsh_parser is not called yet. But you call 'get_number_of_boundary_faces' method.")
        end if

        do i = 1, self%number_of_elements_, 1
            if (self%is_face_element(self%elements_type_(i))) then
                if (self%gmsh_bc2f3ds_bc_(self%elements_entity_tag_(i)) == type) n = n + 1
            end if
        end do
    end function get_number_of_boundary_faces

    function get_number_of_ghost_cells(self) result(n)
        class  (gmsh_parser), intent(in) :: self
        integer(int_kind   )             :: n

        integer(int_kind   )             :: i

        if(.not. self%parsed_)then
            call call_error("'parse' method of gmsh_parser is not called yet. But you call 'get_number_of_ghost_cells' method.")
        end if

        n = self%num_ghost_cells_
    end function get_number_of_ghost_cells

    function get_number_of_cells(self) result(n)
        class  (gmsh_parser), intent(in) :: self
        integer(int_kind   )             :: n

        integer(int_kind   )             :: i

        if(.not. self%parsed_)then
            call call_error("'parse' method of gmsh_parser is not called yet. But you call 'get_number_of_cells' method.")
        end if

        do i = 1, self%number_of_elements_, 1
            if (self%is_cell_element(self%elements_type_(i))) n = n + 1
        end do
    end function get_number_of_cells

    subroutine close(self)
        class(gmsh_parser  ), intent(inout) :: self

        if(.not. self%parsed_)then
            call call_error("'parse' method of gmsh_parser is not called yet. But you call 'close' method.")
        end if

        deallocate(self%gmsh_bc2f3ds_bc_    )
        deallocate(self%node_positions_     )
        deallocate(self%elements_           )
        deallocate(self%elements_entity_tag_)
        deallocate(self%elements_type_      )
    end subroutine

    subroutine parse(self, config)
        class(gmsh_parser  ), intent(inout) :: self
        class(configuration), intent(inout) :: config

        ! for mesh input
        integer  (int_kind)              :: unit_number
        character(len=:   ), allocatable :: mesh_filename
        character(len=20  )              :: section
        integer  (int_kind)              :: io_statiment
        integer  (int_kind)              :: n_line
        logical                          :: parsed_mesh_format    = .false.
        logical                          :: parsed_nodes          = .false.
        logical                          :: parsed_elements       = .false.
        logical                          :: parsed_physical_names = .false.

        ! for config
        character(len=:), allocatable :: error_msg
        logical :: found

        if(.not. self%parsed_)then
            call call_error("gmsh_parser is already parsed.")
        end if

        call config%get_char("Grid.Filepath", mesh_filename, found, "grid.msh")
        if(.not. found) call write_warring("'Grid.Filepath' is not found in configuration you set. To be set default value.")

        open(newunit=unit_number, iostat=io_statiment, file=mesh_filename, access='stream', form='formatted', status='old')
        if(io_statiment /= 0) call call_error("Can not read mesh file '"//mesh_filename//"'.")

        do while (.true.)
            do n_line = 1, self%sweep_section_limit_, 1
                read(unit_number, iostat=io_statiment) section
                call write_debuginfo("Sweep lines until start section. Now read '"//section//"'.")

                if (index(trim(section), "$") == 1) exit

                if(io_statiment == -1) then
                    close(unit_number)
                    exit
                end if

                if (n_line == self%sweep_section_limit_) then
                    close(unit_number)
                    exit
                end if
            end do

            if (trim(section) == "$MeshFormat") then
                call self%read_mesh_format(unit_number)
                parsed_mesh_format = .true.
            else if (trim(section) == "$Nodes") then
                call self%read_nodes(unit_number)
                parsed_nodes = .true.
            else if (trim(section) == "$Elements") then
                call self%read_elements(unit_number)
                parsed_elements = .true.
            else if (trim(section) == "$PhysicalNames") then
                call self%read_physical_names(unit_number)
                parsed_physical_names = .true.
            else
                call write_debuginfo("Skip unuse section '"//trim(section)//"'.")
            end if

            do n_line = 1, self%sweep_section_limit_, 1
                read(unit_number, iostat=io_statiment) section
                call write_debuginfo("Sweep lines until end section. Now read '"//section//"'.")

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

        if(                                     &
            (.not. parse_mesh_format    ) .and. &
            (.not. parsed_nodes         ) .and. &
            (.not. parsed_elements      ) .and. &
            (.not. parsed_physical_names)       &
        ) then
            call call_error("Not complite gmsh sections F3DS needed.")
        end if

        .not. self%parsed_ = .true.
    end subroutine parse

    subroutine read_mesh_format(self, unit_number)
        class  (gmsh_parser), intent(inout) :: self
        integer(int_kind   ), intent(inout) :: unit_number

        character(len=20  ) :: mesh_ver_tmp
        character(len=:   ), allocatable :: mesh_ver
        integer  (int_kind) :: dot_index, mesh_ver_size

        ! mesh_ver_tmp will be write version information (e.g. '4.1')
        read(unit_number) mesh_ver_tmp, self%mesh_filetype_, self%mesh_datasize_
        allocate(mesh_ver, source=trim(mesh_ver_tmp))
        dot_index     = index(mesh_ver, ".")
        mesh_ver_size = size (mesh_ver)
        read(mesh_ver(1          :dot_index-1  ), *) self%mesh_major_version_
        read(mesh_ver(dot_index+1:mesh_ver_size), *) self%mesh_minor_version_

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
            if (name(2:name_size-1) == "fluid") then
                tags  (i) = tag
                bc_ids(i) = -1
            else
                tags  (i) = tag
                bc_ids(i) = string_to_boundary_face_type(name(2:name_size-1))

                if (bc_ids(i) == unknown_face_type) then
                    close(unit_number)
                    call call_error("Unknown face type '"//name(2:name_size-1)//"' found.")
                end if
            end if

            if (tag > max_tag) max_tag = tag
        end do

        allocate(self%gmsh_bc2f3ds_bc_(max_tag))

        do i = 1, num_physical_names, 1
            self%gmsh_bc2f3ds_bc_(tags(i)) = bc_ids(i)
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
        integer(int_kind), allocatable :: node_tags(:)

        read(unit_number) num_entity_blocks, num_nodes, min_node_tag, max_node_tag

        allocate(self%node_positions_(3, num_nodes))
        self%number_of_nodes_ = num_nodes

        ! entity block
        do entity_index = 1, num_entity_blocks, 1
            read(unit_number) entity_dim, entity_tag, parametric, num_nodes_entity

            ! allocate
            allocate(node_tags(num_nodes_entity))

            ! tags
            do node_index = 1, num_nodes_entity, 1
                read(unit_number) node_tags(node_index)
                ! offset considering {@code min_node_tag}
                node_tags(node_index) = node_tags(node_index) - (min_node_tag - 1)
            end do

            ! node
            do node_index = 1, num_nodes_entity, 1
                read(unit_number) self%node_positions_(1, node_tags(node_index)), self%node_positions_(2, node_tags(node_index)), self%node_positions_(3, node_tags(node_index))
            end do

            ! deallocate
            deallocate(node_tags(num_nodes_entity))
        end do
    end subroutine read_nodes

    subroutine read_elements(self, unit_number)
        class  (gmsh_parser), intent(inout) :: self
        integer(int_kind   ), intent(inout) :: unit_number

        integer(int_kind) :: entity_index, element_index, node_index
        integer(int_kind) :: num_entity_blocks, num_elements, min_element_tag, max_element_tag
        integer(int_kind) :: entity_dim, entity_tag, element_type, num_elements_entity
        integer(int_kind) :: element_tag, num_nodes
        integer(int_kind), allocatable :: tmp_ids(:)

        read(unit_number) num_entity_blocks, num_elements, min_element_tag, max_element_tag

        allocate(self%elements_           (num_elements))
        allocate(self%elements_entity_tag_(num_elements))
        allocate(self%elements_type_      (num_elements))
        self%number_of_elements_ = num_elements

        do entity_index = 1, num_entity_blocks, 1
            read(unit_number) entity_dim, entity_tag, element_type, num_elements_entity


            if (.not. supported_elemental_type(element_type)) then
                call call_error("Unsupported elemental type"//to_str(element_type)//" found.")
            end if

            num_nodes = self%elemental_type2num_nodes(element_type)
            allocate(tmp_ids(num_nodes))

            do element_index = 1, num_elements_entity, 1
                read(unit_number) element_tag, tmp_ids

                ! offset considering {@code min_node_tag}
                element_tag = element_tag - (min_element_tag - 1)

                self%elements_(element_tag)%initialize(num_nodes)
                do node_index = 1, num_nodes, 1
                    self%elements_(element_tag)%set_point_id(node_index, tmp_ids(node_index))
                end do
                self%elements_type_      (element_tag) = element_type
                self%elements_entity_tag_(element_tag) = entity_tag
            end do

            deallocate(tmp_ids)
        end do
    end subroutine read_elements

    pure function is_face_element(num) result(face)
        integer(int_kind), intent(in) :: num
        logical           :: face
        face = (2 <= num) .and. (num <= 4)
    end function is_face_element

    pure function is_cell_element(num) result(cell)
        integer(int_kind), intent(in) :: num
        logical           :: cell
        cell = (5 <= num) .and. (num <= 7)
    end function is_cell_element

    pure function supported_elemental_type(num) result(supp)
        integer(int_kind), intent(in) :: num
        logical           :: supp
        ! 1st is a line. over 8th are higher order element.
        supp = (1 < num) .and. (num < 8)
    end function supported_elemental_type

    pure function elemental_type2num_nodes(num) result(n_node)
        integer(int_kind), intent(in) :: num
        integer(int_kind) :: n_node
        integer(int_kind), parameter :: map(33) = [integer(int_kind) :: &
            2, 3, 4, 4, 8, 6, 5, 3, 6, 9, 10, 27, 18, 14, 1, 8, 20, 15, 13, 9, 10, 12, 15, 15, 21, 4, 5, 6, 20, 35, 56, 64, 125]
        n_node = map(num)
    end function elemental_type2num_nodes
end module class_gmsh_parser