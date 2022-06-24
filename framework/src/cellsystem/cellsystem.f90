module class_cellsystem
    use json_module
    use typedef_module
    use stdio_module
    use boundary_type_module
    use class_point_id_list
    use abstract_grid_parser
    use abstract_result_writer

    implicit none

    private

    ! Data structure and rule:
    ! (local cell index)   -2        -1         0         1         2         3
    !                       |         |         |         |         |         |
    !
    !                  *---------*---------*---------*---------*---------*---------*
    !                  |         |         |         | n       |         |         | "n" is a normal vector.
    !                  |    o    |    o    |    o    x--> o    |    o    |    o    | Normal vector must be oriented to the right-hand cell.
    ! (cell index) ->  |    1    |    2    |    3    |    4    |    5    |    6    |
    !                  *---------*---------*---------*---------*---------*---------*
    !
    !   (inner cells) ______|_________|_________|    |    |_________|_________|________ (inner or ghost cells)
    !                                                |             *If fase is boundary, ghost cells are stored right-side.
    !                                      (boundary/inner face)
    type, public :: cellsystem
        private

        ! Real cells only (NOT include points of ghost/dummy cells)
        integer(int_kind) :: num_points

        ! Include ghost/dummy cells.
        integer(int_kind) :: num_cells

        ! Real cells only
        integer(int_kind) :: num_faces

        ! Also known as ghost cells
        integer(int_kind) :: num_local_cells

        ! Number of boundary faces
        integer(int_kind) :: num_outflow_faces
        integer(int_kind) :: num_slipwall_faces
        integer(int_kind) :: num_symmetric_faces

        ! ### Faces ###

        ! Reference for cell index
        ! Elm. 1) 1 : Number of faces                               , Face number
        ! Elm. 2) -Number of ghost cells + 1: Number of ghost cells , Local cell index
        integer(int_kind), allocatable :: face_to_cell_indexs(:,:)

        ! Normal vectors that is directed right-hand cell (local index is one)
        ! Elm. 1) 1 : Number of faces, face number
        ! Elm. 2) 1 : 3              , vector compornents
        real(real_kind), allocatable :: face_normal_vectors(:,:)

        ! Tangental vectors (1st)
        ! Elm. 1) 1 : Number of faces, face number
        ! Elm. 2) 1 : 3              , vector compornents
        real(real_kind), allocatable :: face_tangential1_vectors(:,:)

        ! Tangental vectors (2nd)
        ! Elm. 1) 1 : Number of faces, face number
        ! Elm. 2) 1 : 3              , vector compornents
        real(real_kind), allocatable :: face_tangential2_vectors(:,:)

        ! Face centor position
        ! Elm. 1) 1 : 3              , vector compornents
        ! Elm. 2) 1 : Number of faces, face number
        real(real_kind), allocatable :: face_positions(:,:)

        ! Face area
        ! Elm. 1) 1 : Number of faces, face number
        real(real_kind), allocatable :: face_areas(:)

        ! ### Cells ###

        ! Cell centor positions
        ! Elm. 1) 1 : 3                , vector compornents
        ! Elm. 2) 1 : {@code num_cells}, cell index
        real(real_kind), allocatable :: cell_centor_positions(:,:)

        ! Cell volumes
        ! Elm. 1) 1 : {@code num_cells}, cell index
        real(real_kind), allocatable :: cell_volumes(:)

        ! Elm. 1) 1 : {@code num_cells}, cell index
        logical, allocatable :: is_real_cell(:)

        ! Elm. 1) 1:3 = x, y, z
        ! Elm. 2) 1:{@code num_points}
        real(real_kind), allocatable :: points(:, :)

        ! Elm. 1) 1 : {@code num_cells}, cell index
        class(point_id_list), allocatable :: cell_geometries(:)

        ! Elm. 1) 1 : {@code num_cells}, cell index
        integer(type_kind), allocatable :: cell_types(:)

        ! ### Boundary References ###

        ! Elements are stored boundary face index.
        ! Elm. 1) 1 : num_{boundary condition}_faces
        integer(int_kind), public, allocatable :: outflow_face_indexs (:)
        integer(int_kind), public, allocatable :: slipwall_face_indexs(:)
        integer(int_kind), public, allocatable :: symmetric_face_indexs(:)

        contains

        ! ### Result writer ###
        procedure, public, pass(self) :: initialize_result_writer
        procedure, public, pass(self) :: write_result

        ! ### Cellsystem reader ###
        procedure, public, pass(self) :: read

        ! ### Getter ###
        procedure, public, pass(self) :: get_number_of_faces
        procedure, public, pass(self) :: get_number_of_local_cells
        procedure, public, pass(self) :: get_number_of_points
        procedure, public, pass(self) :: get_number_of_cells
        procedure, public, pass(self) :: get_number_of_outflow_faces
        procedure, public, pass(self) :: get_number_of_slipwall_faces
        procedure, public, pass(self) :: get_number_of_symmetric_faces

        ! ### Finalizer ###
        procedure, public, pass(self) :: finalize

        ! ### Utils ###
        procedure, private, pass(self) :: initialise_faces
        procedure, private, pass(self) :: initialise_cells
        procedure, private, pass(self) :: initialise_boundary_references
        procedure, private, pass(self) :: finalize_faces
        procedure, private, pass(self) :: finalize_cells
        procedure, private, pass(self) :: finalize_boundary_references
        procedure, private, pass(self) :: assign_boundary
    end type

    contains

    ! ###  Result writer ###
    subroutine initialize_result_writer(self, parser, json)
        class(cellsystem   ), intent(inout) :: self
        class(result_writer), intent(inout) :: parser
        class(json_file    ), intent(in   ) :: json

        call parser%initialize(self%is_real_cell, self%cell_geometries, self%cell_types, json)
    end subroutine initialize_result_writer

    subroutine write_result(self, parser, scolar_variables, vector_variables, time)
        class(cellsystem   ), intent(inout) :: self
        class(result_writer), intent(inout) :: parser
        real (real_kind    ), intent(in   ) :: scolar_variables(:,:)
        real (real_kind    ), intent(in   ) :: vector_variables(:,:)
        real (real_kind    ), intent(in   ) :: time

        call parser%write(scolar_variables, vector_variables, time)
    end subroutine write_result

    ! ### Cellsystem reader ###
    subroutine read(self, parser, filepath)
        class(cellsystem ), intent(inout) :: self
        class(grid_parser), intent(inout) :: parser
        character(len=*)  , intent(in   ) :: filepath

        integer(kind(boundary_type)), allocatable :: face_types(:)

        ! parse grid file
        call parser%parse(filepath)

        ! initialize
        allocate(face_types(parser%get_number_of_faces()))

        ! allocate grid
        call self%initialise_faces              (parser%get_number_of_faces (), parser%get_number_of_ghost_cells())
        call self%initialise_cells              (parser%get_number_of_points(), parser%get_number_of_cells      ())
        call self%initialise_boundary_references(parser%get_number_of_boundary_faces(outflow_boundary_type  ), &
                                                 parser%get_number_of_boundary_faces(slipwall_boundary_type ), &
                                                 parser%get_number_of_boundary_faces(symmetric_boundary_type))

        ! get grid data
        call parser%get_cells          (self%cell_centor_positions, self%cell_volumes, self%is_real_cell)
        call parser%get_faces          (self%face_to_cell_indexs, self%face_normal_vectors, self%face_tangential1_vectors, self%face_tangential2_vectors, self%face_positions, self%face_areas)
        call parser%get_cell_geometries(self%points, self%cell_geometries)
        call parser%get_boundaries     (face_types)

        ! assign boundary condition
        call self%assign_boundary(face_types)

        ! close file
        call parser%close()
    end subroutine read

    ! ### Getter ###
    pure function get_number_of_faces(self) result(n)
        class(cellsystem), intent(in) :: self
        integer(int_kind) :: n
        n = self%num_faces
    end function get_number_of_faces

    pure function get_number_of_local_cells(self) result(n)
        class(cellsystem), intent(in) :: self
        integer(int_kind) :: n
        n = self%num_local_cells
    end function get_number_of_local_cells

    pure function get_number_of_cells(self) result(n)
        class(cellsystem), intent(in) :: self
        integer(int_kind) :: n
        n = self%num_cells
    end function get_number_of_cells

    pure function get_number_of_points(self) result(n)
        class(cellsystem), intent(in) :: self
        integer(int_kind) :: n
        n = self%num_points
    end function

    pure function get_number_of_outflow_faces(self) result(n)
        class(cellsystem), intent(in) :: self
        integer(int_kind) :: n
        n = self%num_outflow_faces
    end function get_number_of_outflow_faces

    pure function get_number_of_slipwall_faces(self) result(n)
        class(cellsystem), intent(in) :: self
        integer(int_kind) :: n
        n = self%num_slipwall_faces
    end function get_number_of_slipwall_faces

    pure function get_number_of_symmetric_faces(self) result(n)
        class(cellsystem), intent(in) :: self
        integer(int_kind) :: n
        n = self%num_symmetric_faces
    end function get_number_of_symmetric_faces

    ! ### Finalizer ###
    subroutine finalize(self)
        class(cellsystem), intent(inout) :: self
        call self%finalize_cells()
        call self%finalize_faces()
        call self%finalize_boundary_references()
    end subroutine finalize

    ! ### Inner utils ###
    subroutine assign_boundary(self, face_types)
        class  (cellsystem         ), intent(inout) :: self
        integer(kind(boundary_type)), intent(in   ) :: face_types(:)

        integer(int_kind) :: index, outflow_index, slipwall_index, symmetric_index

        outflow_index   = 1
        slipwall_index  = 1
        symmetric_index = 1

        do index = 1, self%num_faces, 1
            if (face_types(index) == outflow_boundary_type) then
                self%outflow_face_indexs(outflow_index) = index
                outflow_index = outflow_index + 1
            else if (face_types(index) == slipwall_boundary_type) then
                self%slipwall_face_indexs(slipwall_index) = index
                slipwall_index = slipwall_index + 1
            else if (face_types(index) == symmetric_boundary_type) then
                self%symmetric_face_indexs(symmetric_index) = index
                symmetric_index = symmetric_index + 1
            else
                call call_error("Unknown boundary type found.")
            end if
        end do
    end subroutine assign_boundary

    subroutine initialise_faces(self, num_faces, num_local_cells)
        class(cellsystem), intent(inout) :: self
        integer(int_kind), intent(in) :: num_faces
        integer(int_kind), intent(in) :: num_local_cells

        if(num_faces < 1)then
            call call_error("Number of face must be set over zero.")
        end if
        self%num_faces = num_faces

        if(num_local_cells < 1)then
            call call_error("Number of ghost cells must be set over zero.")
        end if
        self%num_local_cells = num_local_cells

        if(allocated(self%face_to_cell_indexs))then
            call call_error("Array face_to_cell_indexs is already allocated. But you call the initialiser for faces.")
        end if
        allocate(self%face_to_cell_indexs(self%num_faces, -self%num_local_cells + 1:self%num_local_cells))

        if(allocated(self%face_normal_vectors))then
            call call_error("Array face_normal_vectors is already allocated. But you call the initialiser for faces.")
        end if
        allocate(self%face_normal_vectors(self%num_faces, 3))

        if(allocated(self%face_tangential1_vectors))then
            call call_error("Array face_tangential1_vectors is already allocated. But you call the initialiser for faces.")
        end if
        allocate(self%face_tangential1_vectors(self%num_faces, 3))

        if(allocated(self%face_tangential2_vectors))then
            call call_error("Array face_tangential2_vectors is already allocated. But you call the initialiser for faces.")
        end if
        allocate(self%face_tangential2_vectors(self%num_faces, 3))

        if(allocated(self%face_positions))then
            call call_error("Array face_positions is already allocated. But you call the initialiser for faces.")
        end if
        allocate(self%face_positions(self%num_faces, 3))

        if(allocated(self%face_areas))then
            call call_error("Array face_areas is already allocated. But you call the initialiser for faces.")
        end if
        allocate(self%face_areas(self%num_faces))
    end subroutine initialise_faces

    subroutine initialise_cells(self, num_points, num_cells)
        class(cellsystem), intent(inout) :: self
        integer(int_kind), intent(in) :: num_points, num_cells

        if(num_points < 1)then
            call call_error("Number of points must be set over zero.")
        end if
        self%num_points = num_points

        if(num_cells < 1)then
            call call_error("Number of cells must be set over zero.")
        end if
        self%num_cells = num_cells

        if(allocated(self%cell_centor_positions))then
            call call_error("Array cell_centor_positions is already allocated. But you call the initialiser for cells.")
        end if
        allocate(self%cell_centor_positions(self%num_cells, 3))

        if(allocated(self%cell_volumes))then
            call call_error("Array cell_volumes is already allocated. But you call the initialiser for cells.")
        end if
        allocate(self%cell_volumes(self%num_cells))

        if(allocated(self%is_real_cell))then
            call call_error("Array is_real_cell is already allocated. But you call the initialiser for cells.")
        end if
        allocate(self%is_real_cell(self%num_cells))

        if(allocated(self%points))then
            call call_error("Array points is already allocated. But you call the initialiser for cells.")
        end if
        allocate(self%points(self%num_points, 3))

        if(allocated(self%cell_geometries))then
            call call_error("Array cell_geometries is already allocated. But you call the initialiser for cells.")
        end if
        allocate(self%cell_geometries(self%num_cells))

        if(allocated(self%cell_types))then
            call call_error("Array cell_types is already allocated. But you call the initialiser for cells.")
        end if
        allocate(self%cell_types(self%num_cells))
    end subroutine initialise_cells

    subroutine initialise_boundary_references(self, num_outflow_faces, num_slipwall_faces, num_symetric_faces)
        class(cellsystem), intent(inout) :: self
        integer(int_kind), intent(in) :: num_outflow_faces
        integer(int_kind), intent(in) :: num_slipwall_faces
        integer(int_kind), intent(in) :: num_symetric_faces

        if(num_outflow_faces < 0)then
            call call_error("Number of outflow BC faces must be set over zero.")
        end if
        self%num_outflow_faces = num_outflow_faces

        if(num_slipwall_faces < 0)then
            call call_error("Number of slipwall BC faces must be set over zero.")
        end if
        self%num_slipwall_faces = num_slipwall_faces

        if(num_symetric_faces < 0)then
            call call_error("Number of symetric BC faces must be set over zero.")
        end if
        self%num_symmetric_faces = num_symetric_faces

        if(allocated(self%outflow_face_indexs))then
            call call_error("Array outflow_face_indexs is already allocated. But you call the initialiser for boundary_reference.")
        end if
        allocate(self%outflow_face_indexs(self%num_outflow_faces))

        if(allocated(self%slipwall_face_indexs))then
            call call_error("Array slipwall_face_indexs is already allocated. But you call the initialiser for boundary_reference.")
        end if
        allocate(self%slipwall_face_indexs(self%num_slipwall_faces))

        if(allocated(self%symmetric_face_indexs))then
            call call_error("Array symmetric_face_indexs is already allocated. But you call the initialiser for boundary_reference.")
        end if
        allocate(self%symmetric_face_indexs(self%num_symmetric_faces))
    end subroutine initialise_boundary_references

    subroutine finalize_faces(self)
        class(cellsystem), intent(inout) :: self

        if(.not. allocated(self%face_to_cell_indexs))then
            call call_error("Array face_to_cell_indexs are not allocated. But you call the finalizer for faces.")
        end if
        deallocate(self%face_to_cell_indexs)

        if(.not. allocated(self%face_normal_vectors))then
            call call_error("Array face_normal_vectors is not allocated. But you call the finalizer for faces.")
        end if
        deallocate(self%face_normal_vectors)

        if(.not. allocated(self%face_tangential1_vectors))then
            call call_error("Array face_tangential1_vectors is not allocated. But you call the finalizer for faces.")
        end if
        deallocate(self%face_tangential1_vectors)

        if(.not. allocated(self%face_tangential2_vectors))then
            call call_error("Array face_tangential2_vectors is not allocated. But you call the finalizer for faces.")
        end if
        deallocate(self%face_tangential2_vectors)

        if(.not. allocated(self%face_positions))then
            call call_error("Array face_positions is not allocated. But you call the finalizer for faces.")
        end if
        deallocate(self%face_positions)

        if(.not. allocated(self%face_areas))then
            call call_error("Array face_areas is not allocated. But you call the finalizer for faces.")
        end if
        deallocate(self%face_areas)

        self%num_faces = 0
        self%num_local_cells = 0
    end subroutine finalize_faces

    subroutine finalize_cells(self)
        class(cellsystem), intent(inout) :: self

        if(.not. allocated(self%cell_centor_positions))then
            call call_error("Array cell_centor_positions is not allocated. But you call the finalizer for cell_geometries_module.")
        end if
        deallocate(self%cell_centor_positions)

        if(.not. allocated(self%cell_volumes))then
            call call_error("Array cell_volumes is not allocated. But you call the finalizer for cell_geometries_module.")
        end if
        deallocate(self%cell_volumes)

        if(.not. allocated(self%is_real_cell))then
            call call_error("Array is_real_cell is not allocated. But you call the finalizer for cell_geometries_module.")
        end if
        deallocate(self%is_real_cell)

        if(.not. allocated(self%points))then
            call call_error("Array points is not allocated. But you call the finalizer for cell_geometries_module.")
        end if
        deallocate(self%points)

        if(.not. allocated(self%cell_geometries))then
            call call_error("Array cell_geometries is not allocated. But you call the finalizer for cell_geometries_module.")
        end if
        deallocate(self%cell_geometries)

        if(.not. allocated(self%cell_types))then
            call call_error("Array cell_types is not allocated. But you call the finalizer for cell_geometries_module.")
        end if
        deallocate(self%cell_types)

        self%num_points = 0
        self%num_cells  = 0
    end subroutine finalize_cells

    subroutine finalize_boundary_references(self)
        class(cellsystem), intent(inout) :: self

        if(.not. allocated(self%outflow_face_indexs))then
            call call_error("Array outflow_face_indexs is not allocated. But you call the finalizer for boundary_reference module.")
        end if
        deallocate(self%outflow_face_indexs)

        if(.not. allocated(self%slipwall_face_indexs))then
            call call_error("Array slipwall_face_indexs is not allocated. But you call the finalizer for boundary_reference module.")
        end if
        deallocate(self%slipwall_face_indexs)

        if(.not. allocated(self%symmetric_face_indexs))then
            call call_error("Array symmetric_face_indexs is not allocated. But you call the finalizer for boundary_reference module.")
        end if
        deallocate(self%symmetric_face_indexs)

        self%num_local_cells     = 0
        self%num_outflow_faces   = 0
        self%num_slipwall_faces  = 0
        self%num_symmetric_faces = 0
    end subroutine finalize_boundary_references
end module class_cellsystem