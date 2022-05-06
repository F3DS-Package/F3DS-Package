module cell_geometries_module
    use typedef_module
    use stdio_module
    use class_point_id_list

    implicit none

    private

    integer(int_kind) :: num_points_
    integer(int_kind) :: num_cells_

    ! Elm. 1) 1:3 = x, y, z
    ! Elm. 2) 1:{@code num_points}
    integer(real_kind), public, allocatable :: points(:, :)

    class(point_id_list), public, allocatable :: cell_geometries(:)

    public :: initialise_cell_geometries
    public :: finalize_cell_geometries
    public :: get_number_of_points

    contains

    pure function get_number_of_points() result(n)
        integer(int_kind) :: n
        n = num_points_
    end function

    subroutine initialise_cell_geometries(num_points, num_cells)
        integer(int_kind), intent(in) :: num_points, num_cells

        if(num_points < 1)then
            call call_error("Number of points must be set over zero.")
        end if
        num_points_ = num_points

        if(num_cells < 1)then
            call call_error("Number of cells must be set over zero.")
        end if
        num_cells_ = num_cells

        if(allocated(points))then
            call call_error("Array points is already allocated. But you call the initialiser for cell_geometries_module.")
        end if
        allocate(points(3, num_points_))

        if(allocated(cell_geometries))then
            call call_error("Array cell_geometries is already allocated. But you call the initialiser for cell_geometries_module.")
        end if
        allocate(cell_geometries(num_cells_))
    end subroutine initialise_cell_geometries

    subroutine finalize_cell_geometries()
        if(.not. allocated(points))then
            call call_error("Array points is not allocated. But you call the finalizer for cell_geometries_module.")
        end if
        deallocate(points)

        if(.not. allocated(cell_geometries))then
            call call_error("Array cell_geometries is not allocated. But you call the finalizer for cell_geometries_module.")
        end if
        deallocate(cell_geometries)

        num_points_ = 0
        num_cells_  = 0
    end subroutine finalize_cell_geometries
end module cell_geometries_module