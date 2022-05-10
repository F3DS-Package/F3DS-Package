module class_point_id_list
    use typedef_module
    use stdio_module

    implicit none

    private

    type, public :: point_id_list
        private

        integer(int_kind) :: num_points_
        integer(int_kind), allocatable :: point_ids_(:)

        contains

        procedure, public :: print_ids    => point_id_list_print_ids
        procedure, public :: print_points => point_id_list_print_points
        procedure, public :: initialize   => point_id_list_initialize
        procedure, public :: set_point_id => point_id_list_set_point_id
        procedure, public :: get_point_id => point_id_list_get_point_id
        procedure, public :: get_number_of_points => point_id_list_get_number_of_points
    end type point_id_list

    contains

    subroutine point_id_list_print_ids(self)
        class  (point_id_list), intent(in) :: self
        print *, self%num_points_, " points:", self%point_ids_
    end subroutine point_id_list_print_ids

    subroutine point_id_list_print_points(self, points)
        class(point_id_list), intent(in) :: self
        real (real_kind    ), intent(in) :: points(:,:)
        integer(int_kind) :: n
        do n = 1, self%num_points_, 1
            print *, "Point # ", n, ": ", points(:,self%point_ids_(n))
        end do
    end subroutine point_id_list_print_points

    pure function point_id_list_get_number_of_points(self) result(n)
        class  (point_id_list), intent(in) :: self
        integer(int_kind     )             :: n
        n = self%num_points_
    end function point_id_list_get_number_of_points

    subroutine point_id_list_initialize(self, num_points)
        class  (point_id_list), intent(inout) :: self
        integer(int_kind     ), intent(in   ) :: num_points

        if(num_points < 1)then
            call call_error("Number of points must be set over zero.")
        end if
        self%num_points_ = num_points
        allocate(self%point_ids_(self%num_points_))
    end subroutine point_id_list_initialize

    subroutine point_id_list_set_point_id(self, index, point_id)
        class  (point_id_list), intent(inout) :: self
        integer(int_kind     ), intent(in   ) :: index
        integer(int_kind     ), intent(in   ) :: point_id
        if(.not. allocated(self%point_ids_))then
            call call_error("'point_id_list' is not allocated. But you call 'set_point_id'.")
        end if
        if((index < 1) .or. (self%num_points_ < index))then
            call call_error("'index' you set to 'point_id_list' out of range.")
        end if
        self%point_ids_(index) = point_id
    end subroutine point_id_list_set_point_id

    pure function point_id_list_get_point_id(self, index) result(point_id)
        class  (point_id_list), intent(in) :: self
        integer(int_kind     ), intent(in) :: index
        integer(int_kind     )             :: point_id
        if((index < 1) .or. (self%num_points_ < index))then
            point_id = -1 ! error
        end if
        point_id = self%point_ids_(index)
    end function point_id_list_get_point_id
end module class_point_id_list