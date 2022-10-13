module class_control_volume
    use typedef_module
    use abstract_configuration
    use stdio_module
    use system_call_module

    implicit none

    private

    type, public :: control_volume
        private

        character(:        ), allocatable :: output_filename_
        integer  (int_kind ), allocatable :: cell_ids_  (:)
        real     (real_kind)              :: output_timespan_
        real     (real_kind)              :: next_output_time_
        integer  (int_kind )              :: n_output_cells_

        contains

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: write
        procedure, public, pass(self) :: is_writable
    end type control_volume

    contains

    subroutine initialize(self, config, cell_positions, is_real_cell, num_cells)
        class  (control_volume), intent(inout) :: self
        class  (configuration ), intent(inout) :: config
        real   (real_kind     ), intent(in   ) :: cell_positions(:,:)
        logical                , intent(in   ) :: is_real_cell  (:)
        integer(int_kind      ), intent(in   ) :: num_cells

        real (real_kind    ) :: frequency
        real (real_kind    ) :: min_point(3), max_point(3)
        integer(int_kind   ), allocatable :: tmp_ids(:)
        integer(int_kind   ) :: i, j, unit_number
        logical              :: found, in_x_range, in_y_range, in_z_range

        self%output_filename_      = "control_volume.dat"
        self%next_output_time_     = 0.d0

        call config%get_real("Control volume.Frequency", frequency, found, 1.d3)
        if(.not. found) call write_warring("'Control volume.Frequency' is not found in configuration file you set. To be set default value.")
        self%output_timespan_      = 1.d0 / frequency

        call config%get_real("Control volume.Min point.x", min_point(1), found)
        if(.not. found) call call_error("'Control volume.Min point.x' is not found in configuration file you set.")
        call config%get_real("Control volume.Min point.y", min_point(2), found)
        if(.not. found) call call_error("'Control volume.Min point.y' is not found in configuration file you set.")
        call config%get_real("Control volume.Min point.z", min_point(3), found)
        if(.not. found) call call_error("'Control volume.Min point.z' is not found in configuration file you set.")

        call config%get_real("Control volume.Max point.x", max_point(1), found)
        if(.not. found) call call_error("'Control volume.Max point.x' is not found in configuration file you set.")
        call config%get_real("Control volume.Max point.y", max_point(2), found)
        if(.not. found) call call_error("'Control volume.Max point.y' is not found in configuration file you set.")
        call config%get_real("Control volume.Max point.z", max_point(3), found)
        if(.not. found) call call_error("'Control volume.Max point.z' is not found in configuration file you set.")

        allocate(tmp_ids(num_cells))

        self%n_output_cells_ = 0
        do i = 1, num_cells, 1
            in_x_range = min_point(1) <= cell_positions(1,i) .and. cell_positions(1,i) <= max_point(1)
            in_y_range = min_point(2) <= cell_positions(2,i) .and. cell_positions(2,i) <= max_point(2)
            in_z_range = min_point(3) <= cell_positions(3,i) .and. cell_positions(3,i) <= max_point(3)
            if (in_x_range .and. in_y_range .and. in_z_range .and. is_real_cell(i)) then
                self%n_output_cells_ = self%n_output_cells_ + 1
                tmp_ids(self%n_output_cells_) = i
            end if
        end do

        allocate(self%cell_ids_, source=tmp_ids(1:self%n_output_cells_))
        deallocate(tmp_ids)

        call make_dir("result/")

        open(newunit = unit_number, file = "result/"//self%output_filename_, status = 'replace')
        write(unit_number, *) "# F3DS control volume data"
        close(unit_number)
    end subroutine initialize

    subroutine write(self, time, values_set, cell_positions, cell_volumes)
        class  (control_volume), intent(inout) :: self
        real   (real_kind     ), intent(in   ) :: time
        real   (real_kind     ), intent(in   ) :: values_set    (:,:)
        real   (real_kind     ), intent(in   ) :: cell_positions(:,:)
        real   (real_kind     ), intent(in   ) :: cell_volumes  (:)

        integer  (int_kind ) :: unit_number, i
        character(7        ) :: cher_n_output
        real     (real_kind) :: integrals(size(values_set(:,1)))

        if(time >= self%next_output_time_)then
            open(newunit = unit_number, file="result/"//self%output_filename_, status = 'old', position='append')
            do i = 1, self%n_output_cells_, 1
                integrals(:) = integrals(:) + values_set(:,self%cell_ids_(i)) * cell_volumes(self%cell_ids_(i))
            end do
            write(unit_number, *) time, integrals(:)
            close(unit_number)
            self%next_output_time_ = self%next_output_time_ + self%output_timespan_
        end if
    end subroutine write

    pure function is_writable(self, time) result(juge)
        class  (control_volume), intent(in) :: self
        real   (real_kind     ), intent(in) :: time
        logical                             :: juge
        juge = (time >= self%next_output_time_)
    end function is_writable
end module class_control_volume