module class_line_plotter
    use typedef_module
    use system_call_module
    use abstract_configuration
    use abstract_measurement_tool
    use vector_module
    use math_constant_module
    use stdio_module

    implicit none

    private

    type, public, extends(measurement_tool) :: line_plotter
        private

        character(:        ), allocatable :: output_dir_
        integer  (int_kind ), allocatable :: cell_ids_  (:)
        real     (real_kind)              :: output_timespan_
        real     (real_kind)              :: next_output_time_
        integer  (int_kind )              :: n_output_
        integer  (int_kind )              :: n_output_cells_
        logical                           :: is_enabled_

        contains

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: write
        procedure, public, pass(self) :: is_writable
    end type line_plotter

    contains

    subroutine initialize(     &
        self,                  &
        a_config,              &
        cell_positions,        &
        cell_volumes,          &
        is_real_cell,          &
        face_to_cell_index,    &
        face_centor_positions, &
        face_normal_vectors,   &
        face_areas,            &
        num_cells,             &
        num_faces,             &
        num_local_cells        &
    )
        class  (line_plotter ), intent(inout) :: self
        class  (configuration), intent(inout) :: a_config
        real   (real_kind    ), intent(in   ) :: cell_positions       (:, :)
        real   (real_kind    ), intent(in   ) :: cell_volumes         (:)
        logical               , intent(in   ) :: is_real_cell         (:)
        integer(int_kind     ), intent(in   ) :: face_to_cell_index   (:, :)
        real   (real_kind    ), intent(in   ) :: face_centor_positions(:, :)
        real   (real_kind    ), intent(in   ) :: face_normal_vectors  (:, :)
        real   (real_kind    ), intent(in   ) :: face_areas           (:)
        integer(int_kind     ), intent(in   ) :: num_cells
        integer(int_kind     ), intent(in   ) :: num_faces
        integer(int_kind     ), intent(in   ) :: num_local_cells

        real (real_kind    ) :: frequency
        real (real_kind    ) :: start_point(3), end_point(3), norm_vector(3)
        real (real_kind    ), allocatable :: points(:,:)
        integer(int_kind   ), allocatable :: tmp_ids(:)
        logical             , allocatable :: is_not_redundant_point(:)
        real (real_kind    ) :: distance, neighbor_distance
        integer(int_kind   ) :: i, j, neighbor_id, num_points
        logical              :: found

        self%output_dir_           = "line"
        self%next_output_time_     = 0.d0

        call a_config%get_bool("Line plotter.Enable", self%is_enabled_, found, .false.)
        if(.not. found) call write_warring("'Line plotter.Enable' is not found in configuration file you set. Disable the line plotter.")

        if (.not. self%is_enabled_) return

        call a_config%get_real("Line plotter.Frequency", frequency, found, 1.d3)
        if(.not. found) call write_warring("'Line plotter.Frequency' is not found in configuration file you set. To be set default value.")
        self%output_timespan_      = 1.d0 / frequency

        call a_config%get_real("Line plotter.Start point.x", start_point(1), found)
        if(.not. found) call call_error("'Line plotter.Start point.x' is not found in configuration file you set.")
        call a_config%get_real("Line plotter.Start point.y", start_point(2), found)
        if(.not. found) call call_error("'Line plotter.Start point.y' is not found in configuration file you set.")
        call a_config%get_real("Line plotter.Start point.z", start_point(3), found)
        if(.not. found) call call_error("'Line plotter.Start point.z' is not found in configuration file you set.")

        call a_config%get_real("Line plotter.End point.x", end_point(1), found)
        if(.not. found) call call_error("'Line plotter.End point.x' is not found in configuration file you set.")
        call a_config%get_real("Line plotter.End point.y", end_point(2), found)
        if(.not. found) call call_error("'Line plotter.End point.y' is not found in configuration file you set.")
        call a_config%get_real("Line plotter.End point.z", end_point(3), found)
        if(.not. found) call call_error("'Line plotter.End point.z' is not found in configuration file you set.")

        call a_config%get_int("Line plotter.Number of points", num_points, found, 1000)
        if(.not. found) call write_warring("'Line plotter.Number of points' is not found in configuration file you set. To be set default value.")

        allocate(points(3, num_points))
        allocate(tmp_ids(num_points))
        allocate(is_not_redundant_point(num_points),source=.true.)

        ! Make standardization vector
        norm_vector(:) = (end_point(:) - start_point(:)) / dble(num_points - 1)

        ! Generate points
        do i = 1, num_points, 1
            points(:,i) = start_point(:) + norm_vector(:) * dble(i - 1)
        end do

        ! Allocate id
        do i = 1, num_points, 1
            ! Nearest neighbor search
            neighbor_distance = 1.d0 / machine_epsilon
            do j = 1, num_cells, 1
                if(.not. is_real_cell(j)) cycle

                distance = vector_distance(points(:,i), cell_positions(:,j))
                if (distance < neighbor_distance) then
                    neighbor_distance = distance
                    neighbor_id = j
                end if
            end do
            tmp_ids(i) = neighbor_id
        end do

        ! Check redundant point ids
        do i = 1, num_points, 1
            if (is_not_redundant_point(i) .eqv. .false.) cycle

            do j = i + 1, num_points, 1
                if (tmp_ids(i) == tmp_ids(j)) then
                    is_not_redundant_point(j) = .false.
                end if
            end do
        end do

        ! Allocate
        allocate(self%cell_ids_, source = pack(tmp_ids(:), mask=is_not_redundant_point))

        self%n_output_cells_ = size(self%cell_ids_)
        self%n_output_       = 0

        call make_dir("result/"//self%output_dir_)
    end subroutine

    subroutine write(       &
        self,               &
        time,               &
        values_set,         &
        cell_positions,     &
        cell_volumes,       &
        face_to_cell_index, &
        face_areas,         &
        num_local_cells     &
    )
        class  (line_plotter), intent(inout) :: self
        real   (real_kind   ), intent(in   ) :: time
        real   (real_kind   ), intent(in   ) :: values_set        (:,:)
        real   (real_kind   ), intent(in   ) :: cell_positions    (:,:)
        real   (real_kind   ), intent(in   ) :: cell_volumes      (:)
        integer(int_kind    ), intent(in   ) :: face_to_cell_index(:,:)
        real   (real_kind   ), intent(in   ) :: face_areas        (:)
        integer(int_kind    ), intent(in   ) :: num_local_cells

        integer  (int_kind )              :: unit_number, i
        character(7)                      :: cher_n_output

        if (.not. self%is_enabled_) return

        if(time >= self%next_output_time_)then
            write(cher_n_output, '(i7.7)') self%n_output_
            open(newunit = unit_number, file="result/"//self%output_dir_//"/"//cher_n_output//".dat", status = 'replace')
            write(unit_number, '(a, g0)') "# time = ", time
            do i = 1, self%n_output_cells_, 1
                write(unit_number, '(*(g0, " "))') cell_positions(:,self%cell_ids_(i)), values_set(:,self%cell_ids_(i))
            end do
            close(unit_number)
            self%next_output_time_ = self%next_output_time_ + self%output_timespan_
            self%n_output_         = self%n_output_         + 1
        end if
    end subroutine write

    pure function is_writable(self, time) result(juge)
        class  (line_plotter), intent(in) :: self
        real   (real_kind   ), intent(in) :: time
        logical                           :: juge

        if (.not. self%is_enabled_) then
            juge = .false.
        else
            juge = (time >= self%next_output_time_)
        end if
    end function is_writable

end module class_line_plotter