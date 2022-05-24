module sensor_class
    use typedef_module

    implicit none

    private

    type, public :: sensor
        private

        integer  (int_kind )              :: cell_id_
        character(:        ), allocatable :: output_filename_
        real     (real_kind)              :: output_timespan_
        real     (real_kind)              :: next_output_time_

        contains

        procedure, public, pass(self) :: initialize => initialize_by_id
    end type

    contains

    subroutine initialize_by_id(self, filename, frequency, cell_id)
        class    (sensor   ), intent(inout) :: self
        character(len=*    ), intent(in   ) :: filename
        real     (real_kind), intent(in   ) :: frequency
        integer  (int_kind ), intent(in   ) :: cell_id
        integer  (int_kind )                :: unit_number

        self%output_filename_  = filename
        self%output_timespan_  = 1.d0 / frequency
        self%cell_id_          = cell_id
        self%next_output_time_ = 0.d0

        open(newunit = unit_number, file = self%output_filename_, status = 'replace')
        write(unit_number, *) "# F3DS sensor data (cell id = ", self%cell_id_, ")"
        close(unit_number)
    end subroutine initialize_by_id

    subroutine write(self, time, primitive_values_set)
        class  (sensor   ), intent(inout) :: self
        real   (real_kind), intent(in   ) :: time
        real   (real_kind), intent(in   ) :: primitive_values_set(:,:)
        integer(int_kind )                :: unit_number

        if(time >= self%next_output_time_)then
            open(newunit = unit_number, file=self%output_filename_, status = 'old')
            write(unit_number, *) time, primitive_values_set(self%cell_id_, :)
            close(unit_number)
        end if

        self%next_output_time_ = self%next_output_time_ + self%output_timespan_
    end subroutine write
end module sensor_class