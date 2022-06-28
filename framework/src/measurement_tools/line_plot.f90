module line_plot_class
    use typedef_module
    use system_call_module

    implicit none

    private

    type, public :: line_plot
        private

        character(:        ), allocatable :: output_dir_
        integer  (int_kind ), allocatable :: cell_ids_  (:)
        real     (real_kind)              :: output_timespan_
        real     (real_kind)              :: next_output_time_
        integer  (int_kind )              :: n_output_
        integer  (int_kind )              :: n_cells_

        contains

        procedure, public, pass(self) :: initialize => initialize_by_ids
        procedure, public, pass(self) :: write
    end type line_plot

    contains

    subroutine initialize_by_ids(self, output_dir, frequency, cell_ids)
        class    (line_plot), intent(inout) :: self
        character(len=*    ), intent(in   ) :: output_dir
        real     (real_kind), intent(in   ) :: frequency
        integer  (int_kind ), intent(in   ) :: cell_ids(:)

        integer  (int_kind )                :: unit_number

        self%output_dir_           = output_dir
        self%output_timespan_      = 1.d0 / frequency
        self%next_output_time_     = 0.d0
        self%n_cells_              = size(cell_ids)
        self%n_output_             = 0
        allocate(self%cell_ids_, source=cell_ids)

        call make_dir("result/"//self%output_dir_)
    end subroutine

    subroutine write(self, time, values_set, cell_positions)
        class  (line_plot), intent(inout) :: self
        real   (real_kind), intent(in   ) :: time
        real   (real_kind), intent(in   ) :: values_set(:,:)
        real   (real_kind), intent(in   ) :: cell_positions      (:,:)

        integer  (int_kind )              :: unit_number, i
        character(7)                      :: cher_n_output

        if(time >= self%next_output_time_)then
            write(cher_n_output, "(i7.7)") self%n_output_
            open(newunit = unit_number, file="result/"//self%output_dir_//"/"//cher_n_output//".dat", status = 'replace')
            write(unit_number, *) "# time = ", time
            do i = 1, self%n_cells_, 1
                write(unit_number, *) cell_positions(self%cell_ids_(i), :), values_set(self%cell_ids_(i), :)
            end do
            close(unit_number)
            self%next_output_time_ = self%next_output_time_ + self%output_timespan_
            self%n_output_         = self%n_output_         + 1
        end if
    end subroutine write

end module line_plot_class