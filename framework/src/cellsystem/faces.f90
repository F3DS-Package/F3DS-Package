module faces_module
    use typedef_module
    use stdio_module

    implicit none

    private

    integer(int_kind) :: num_faces_
    integer(int_kind) :: num_local_cells_

    ! Reference for cells index
    ! (local cell index)   -2        -1         0         1         2         3
    !                       |         |         |         |         |         |
    !
    !                  *---------*---------*---------*---------*---------*---------*
    !                  |         |         |         |         |         |         |
    !                  |    o    |    o    |    o    x    o    |    o    |    o    |
    ! (cell index) ->  |    1    |    2    |    3    |    4    |    5    |    6    |
    !                  *---------*---------*---------*---------*---------*---------*
    !
    !   (inner cells) ______|_________|_________|    |    |_________|_________|________ (ghost cells)
    !                                                |             *If fase is boundary, ghost cells are stored right-side.
    !                                         (boundary face)
    !
    ! Elm. 2) -Number of ghost cells + 1: Number of ghost cells , Local cell index
    ! Elm. 3) 1 : Number of faces                               , Face number
    integer(int_kind), public, allocatable :: faces_reference_cell_index(:,:)

    ! Normal vectors that is directed right-hand cell (local index is one)
    ! Elm. 1) 1 : 3              , vector compornents
    ! Elm. 2) 1 : Number of faces, face number
    real(real_kind), public, allocatable :: faces_normal_vector(:,:)

    ! Tangental vectors (1st)
    ! Elm. 1) 1 : 3              , vector compornents
    ! Elm. 2) 1 : Number of faces, face number
    real(real_kind), public, allocatable :: faces_tangential1_vector(:,:)

    ! Tangental vectors (2nd)
    ! Elm. 1) 1 : 3              , vector compornents
    ! Elm. 2) 1 : Number of faces, face number
    real(real_kind), public, allocatable :: faces_tangential2_vector(:,:)

    ! Face centor position
    ! Elm. 1) 1 : 3              , vector compornents
    ! Elm. 2) 1 : Number of faces, face number
    real(real_kind), public, allocatable :: faces_position(:,:)

    ! Face area
    ! Elm. 1) 1 : Number of faces, face number
    real(real_kind), public, allocatable :: faces_area(:)

    public :: get_number_of_faces
    public :: get_number_of_local_cells
    public :: initialise_faces
    public :: finalize_faces

    contains

    pure function get_number_of_faces() result(n)
        integer(int_kind) :: n
        n = num_faces_
    end function get_number_of_faces

    pure function get_number_of_local_cells() result(n)
        integer(int_kind) :: n
        n = num_local_cells_
    end function get_number_of_local_cells

    subroutine initialise_faces(num_faces, num_local_cells)
        integer(int_kind), intent(in) :: num_faces
        integer(int_kind), intent(in) :: num_local_cells

        if(num_faces < 1)then
            call call_error("Number of face must be set over zero.")
        end if
        num_faces_ = num_faces

        if(num_local_cells < 1)then
            call call_error("Number of ghost cells must be set over zero.")
        end if
        num_local_cells_ = num_local_cells

        if(allocated(faces_reference_cell_index))then
            call call_error("Array faces_reference_cell_index is already allocated. But you call the initialiser for faces module.")
        end if
        allocate(faces_reference_cell_index(-num_local_cells + 1:num_local_cells, num_faces))

        if(allocated(faces_normal_vector))then
            call call_error("Array faces_normal_vector is already allocated. But you call the initialiser for faces module.")
        end if
        allocate(faces_normal_vector(3, num_faces))

        if(allocated(faces_tangential1_vector))then
            call call_error("Array faces_tangential1_vector is already allocated. But you call the initialiser for faces module.")
        end if
        allocate(faces_tangential1_vector(3, num_faces))

        if(allocated(faces_tangential2_vector))then
            call call_error("Array faces_tangential2_vector is already allocated. But you call the initialiser for faces module.")
        end if
        allocate(faces_tangential2_vector(3, num_faces))

        if(allocated(faces_position))then
            call call_error("Array faces_position is already allocated. But you call the initialiser for faces module.")
        end if
        allocate(faces_position(3, num_faces))

        if(allocated(faces_area))then
            call call_error("Array faces_area is already allocated. But you call the initialiser for faces module.")
        end if
        allocate(faces_area(num_faces))
    end subroutine initialise_faces

    subroutine finalize_faces()
        if(.not. allocated(faces_reference_cell_index))then
            call call_error("Array faces_reference_cell_index are not allocated. But you call the finalizer for faces module.")
        end if
        deallocate(faces_reference_cell_index)

        if(allocated(faces_normal_vector))then
            call call_error("Array faces_normal_vector is not allocated. But you call the finalizer for faces module.")
        end if
        deallocate(faces_normal_vector)

        if(allocated(faces_tangential1_vector))then
            call call_error("Array faces_tangential1_vector is not allocated. But you call the finalizer for faces module.")
        end if
        deallocate(faces_tangential1_vector)

        if(allocated(faces_tangential2_vector))then
            call call_error("Array faces_tangential2_vector is not allocated. But you call the finalizer for faces module.")
        end if
        deallocate(faces_tangential2_vector)

        if(allocated(faces_position))then
            call call_error("Array faces_position is not allocated. But you call the finalizer for faces module.")
        end if
        deallocate(faces_position)

        if(allocated(faces_area))then
            call call_error("Array faces_area is not allocated. But you call the finalizer for faces module.")
        end if
        deallocate(faces_area)

        num_faces_ = 0
        num_local_cells_ = 0
    end subroutine finalize_faces
end module faces_module