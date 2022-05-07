module cells_module
    use typedef_module
    use stdio_module

    implicit none

    private

    ! number of cells
    integer(int_kind):: num_cells_

    ! Cell centor positions
    ! Elm. 1) 1 : 3                , vector compornents
    ! Elm. 2) 1 : {@code num_cells}, cell index
    real(real_kind), public, allocatable :: cells_centor_position(:,:)

    ! Cell volumes
    ! Elm. 1) 1 : {@code num_cells}, cell index
    real(real_kind), public, allocatable :: cells_volume(:)

    ! Elm. 1) 1 : {@code num_cells}, cell index
    logical, public, allocatable :: cells_is_real_cell(:)

    public :: get_number_of_cells
    public :: initialise_cells
    public :: finalize_cells

    contains

    pure function get_number_of_cells() result(n)
        integer(int_kind) :: n
        n = num_cells_
    end function get_number_of_cells

    subroutine initialise_cells(num_cells)
        integer(int_kind), intent(in) :: num_cells

        if(num_cells < 1)then
            call call_error("Number of cells must be set over zero.")
        end if
        num_cells_ = num_cells

        if(allocated(cells_centor_position))then
            call call_error("Array cells_centor_position is already allocated. But you call the initialiser for cells module.")
        end if
        allocate(cells_centor_position(3, num_cells))

        if(allocated(cells_volume))then
            call call_error("Array cells_volume is already allocated. But you call the initialiser for cells module.")
        end if
        allocate(cells_volume(num_cells))

        if(allocated(cells_is_real_cell))then
            call call_error("Array cells_is_real_cell is already allocated. But you call the initialiser for cells module.")
        end if
        allocate(cells_is_real_cell(num_cells))
    end subroutine initialise_cells

    subroutine finalize_cells()
        if(.not. allocated(cells_centor_position))then
            call call_error("Array cells_centor_position is not allocated. But you call the finalizer for cells module.")
        end if
        deallocate(cells_centor_position)

        if(.not. allocated(cells_volume))then
            call call_error("Array cells_volume is not allocated. But you call the finalizer for cells module.")
        end if
        deallocate(cells_volume)

        if(.not. allocated(cells_is_real_cell))then
            call call_error("Array cells_is_real_cell is not allocated. But you call the finalizer for cells module.")
        end if
        deallocate(cells_is_real_cell)

        num_cells_ = 0
    end subroutine finalize_cells
end module cells_module