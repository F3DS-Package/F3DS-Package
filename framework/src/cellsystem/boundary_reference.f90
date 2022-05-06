module boundary_reference_module
    use typedef_module
    use stdio_module

    implicit none

    private

    integer(int_kind) :: num_ghost_cells_

    integer(int_kind) :: num_outflow_faces_
    integer(int_kind) :: num_slipwall_faces_
    integer(int_kind) :: num_symmetric_faces_

    ! This array is used by boundary faces to refer to ghost cells.
    ! Elements are stored cell index.
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
    !                                                |
    !                                         (boundary face)
    !
    ! Elm. 1) -num_ghost_cell + 1 : num_ghost_cell, local cell index
    ! Elm. 2) 1 : num_{boundary condition}_faces
    integer(int_kind), public, allocatable :: outflow_face_indexs (:, :)
    integer(int_kind), public, allocatable :: slipwall_face_indexs(:, :)
    integer(int_kind), public, allocatable :: symetric_face_indexs(:, :)

    public :: get_number_of_ghost_cells
    public :: get_number_of_outflow_faces
    public :: get_number_of_slipwall_faces
    public :: get_number_of_symmetric_faces
    public :: initialise_boundary_reference
    public :: finalize_boundary_reference

    contains

    pure function get_number_of_ghost_cells() result(n)
        integer(int_kind) :: n
        n = num_ghost_cells_
    end function get_number_of_ghost_cells

    pure function get_number_of_outflow_faces() result(n)
        integer(int_kind) :: n
        n = num_outflow_faces_
    end function get_number_of_outflow_faces

    pure function get_number_of_slipwall_faces() result(n)
        integer(int_kind) :: n
        n = num_slipwall_faces_
    end function get_number_of_slipwall_faces

    pure function get_number_of_symmetric_faces() result(n)
        integer(int_kind) :: n
        n = num_symmetric_faces_
    end function get_number_of_symmetric_faces

    subroutine initialise_boundary_reference(num_ghost_cells, num_outflow_faces, num_slipwall_faces, num_symetric_faces)
        integer(int_kind), intent(in) :: num_ghost_cells
        integer(int_kind), intent(in) :: num_outflow_faces
        integer(int_kind), intent(in) :: num_slipwall_faces
        integer(int_kind), intent(in) :: num_symetric_faces

        if(num_ghost_cells < 1)then
            call call_error("Number of ghost cells must be set over zero.")
        end if
        num_ghost_cells_ = num_ghost_cells

        if(num_outflow_faces < 0)then
            call call_error("Number of outflow BC faces must be set over zero.")
        end if
        num_outflow_faces_ = num_outflow_faces

        if(num_slipwall_faces < 0)then
            call call_error("Number of slipwall BC faces must be set over zero.")
        end if
        num_slipwall_faces_ = num_slipwall_faces

        if(num_symetric_faces < 1)then
            call call_error("Number of symetric BC faces must be set over zero.")
        end if
        num_symmetric_faces_ = num_symetric_faces

        if(allocated(outflow_face_indexs))then
            call call_error("Array outflow_face_indexs is already allocated. But you call the initialiser for boundary_reference module.")
        end if
        allocate(outflow_face_indexs(-num_ghost_cells_ + 1 : num_ghost_cells_, num_outflow_faces_))

        if(allocated(slipwall_face_indexs))then
            call call_error("Array slipwall_face_indexs is already allocated. But you call the initialiser for boundary_reference module.")
        end if
        allocate(slipwall_face_indexs(-num_ghost_cells_ + 1 : num_ghost_cells_, num_slipwall_faces_))

        if(allocated(symetric_face_indexs))then
            call call_error("Array symetric_face_indexs is already allocated. But you call the initialiser for boundary_reference module.")
        end if
        allocate(symetric_face_indexs(-num_ghost_cells_ + 1 : num_ghost_cells_, num_symmetric_faces_))
    end subroutine initialise_boundary_reference

    subroutine finalize_boundary_reference()
        if(.not. allocated(outflow_face_indexs))then
            call call_error("Array outflow_face_indexs is not allocated. But you call the finalizer for boundary_reference module.")
        end if
        deallocate(outflow_face_indexs)

        if(.not. allocated(slipwall_face_indexs))then
            call call_error("Array slipwall_face_indexs is not allocated. But you call the finalizer for boundary_reference module.")
        end if
        deallocate(slipwall_face_indexs)

        if(.not. allocated(symetric_face_indexs))then
            call call_error("Array symetric_face_indexs is not allocated. But you call the finalizer for boundary_reference module.")
        end if
        deallocate(symetric_face_indexs)

        num_ghost_cells_     = 0
        num_outflow_faces_   = 0
        num_slipwall_faces_  = 0
        num_symmetric_faces_ = 0
    end subroutine finalize_boundary_reference
end module boundary_reference_module