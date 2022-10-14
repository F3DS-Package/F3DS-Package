module measurement_surface_class
    use typedef_module
    use vector_module
    use math_constant_module

    implicit none

    private

    type, public :: measurement_surface
        private

        integer  (int_kind ), allocatable :: face_ids_  (:)
        integer  (int_kind ), allocatable :: lhc_or_rhc_(:) ! 0 or 1
        character(:        ), allocatable :: output_filename_
        real     (real_kind)              :: output_timespan_
        real     (real_kind)              :: next_output_time_
        integer  (int_kind )              :: n_faces_

        contains

        procedure, public, pass(self) :: initialize => initialize_by_face_positions
        procedure, public, pass(self) :: write
        procedure,         pass(self) :: make_new_file
        procedure,         pass(self) :: is_point_in_region
    end type measurement_surface

    contains

    subroutine initialize_by_face_positions(self, filename, frequency, x_max, x_min, y_max, y_min, z_max, z_min, normal, face_to_cell_index, face_centor_positions, face_normal_vectors, is_real_cell)
        class    (measurement_surface), intent(inout) :: self
        character(len=*              ), intent(in   ) :: filename
        real     (real_kind          ), intent(in   ) :: frequency
        real     (real_kind          ), intent(in   ) :: x_max, x_min, y_max, y_min, z_max, z_min
        real     (real_kind          ), intent(in   ) :: normal               (3)
        integer  (int_kind           ), intent(in   ) :: face_to_cell_index   (:, :)
        real     (real_kind          ), intent(in   ) :: face_centor_positions(:, :)
        real     (real_kind          ), intent(in   ) :: face_normal_vectors  (:, :)
        logical                       , intent(in   ) :: is_real_cell         (:)

        integer  (int_kind           ), allocatable   :: tmp_face_ids  (:)
        integer  (int_kind           ), allocatable   :: tmp_lhc_or_rhc(:)
        integer  (int_kind           )                :: n_ghost_cells
        integer  (int_kind           )                :: n_faces, n_measurement_surface, face_index, rhc_index, lhc_index
        real     (real_kind          )                :: angle

        n_ghost_cells = size(face_to_cell_index   (1, :)) / 2
        n_faces       = size(face_centor_positions(:, 1))
        allocate(tmp_face_ids  (n_faces))
        allocate(tmp_lhc_or_rhc(n_faces))
        n_measurement_surface = 0
        do face_index = 1, n_faces, 1
            if(self%is_point_in_region(face_centor_positions(face_index, :), x_max, x_min, y_max, y_min, z_max, z_min))then
                angle = vector_angle(normal, face_normal_vectors(face_index, :))
                if(0.d0 - machine_epsilon <= angle .and. angle <= 0.d0 + machine_epsilon)then
                    lhc_index = face_to_cell_index(face_index, n_ghost_cells+0)
                    if(is_real_cell(lhc_index))then
                        n_measurement_surface = n_measurement_surface + 1
                        tmp_face_ids  (n_measurement_surface) = face_index
                        tmp_lhc_or_rhc(n_measurement_surface) = 0
                    end if
                else if(pi - machine_epsilon <= angle .and. angle <= pi + machine_epsilon)then
                    rhc_index = face_to_cell_index(face_index, n_ghost_cells+1)
                    if(is_real_cell(rhc_index))then
                        n_measurement_surface = n_measurement_surface + 1
                        tmp_face_ids  (n_measurement_surface) = face_index
                        tmp_lhc_or_rhc(n_measurement_surface) = 1
                    end if
                end if
            end if
        end do

        self%output_filename_  = filename
        self%output_timespan_  = 1.d0 / frequency
        allocate(self%face_ids_  , source = tmp_face_ids  (1:n_measurement_surface))
        allocate(self%lhc_or_rhc_, source = tmp_lhc_or_rhc(1:n_measurement_surface))
        self%next_output_time_ = 0.d0
        self%n_faces_          = n_measurement_surface

        call self%make_new_file(filename)
    end subroutine initialize_by_face_positions

    subroutine write(self, time, values_set, face_to_cell_index, face_areas)
        class  (measurement_surface), intent(inout) :: self
        real   (real_kind          ), intent(in   ) :: time
        real   (real_kind          ), intent(in   ) :: values_set        (:,:)
        integer(int_kind           ), intent(in   ) :: face_to_cell_index(:,:)
        real   (real_kind          ), intent(in   ) :: face_areas        (:)

        integer(int_kind           )                :: n_ghost_cells
        integer(int_kind           )                :: unit_number, id_index, cell_index, values_index, n_output_values
        real   (real_kind          )                :: output_values(size(values_set(1,:)))

        n_ghost_cells = size(face_to_cell_index(1, :)) / 2

        if(time >= self%next_output_time_)then
            n_output_values = size(values_set(1,:))
            do id_index = 1, self%n_faces_, 1
                cell_index = face_to_cell_index(self%face_ids_(id_index), n_ghost_cells+self%lhc_or_rhc_(id_index))
                do values_index = 1, n_output_values, 1
                    output_values(values_index) = output_values(values_index) + values_set(cell_index, values_index) * face_areas(self%face_ids_(id_index))
                end do
            end do
            open(newunit = unit_number, file="result/"//self%output_filename_, status = 'old', position = 'append')
            write(unit_number, *) time, output_values(:)
            close(unit_number)
            self%next_output_time_ = self%next_output_time_ + self%output_timespan_
        end if
    end subroutine write

    subroutine make_new_file(self, filename)
        class    (measurement_surface), intent(in) :: self
        character(len=*              ), intent(in) :: filename
        integer  (int_kind           )             :: unit_number

        open(newunit = unit_number, file = "result/"//self%output_filename_, status = 'replace')
        write(unit_number, *) "# F3DS measurement surface data"
        close(unit_number)
    end subroutine make_new_file

    pure function is_point_in_region(self, point, x_max, x_min, y_max, y_min, z_max, z_min) result(yes)
        class(measurement_surface), intent(in) :: self
        real (real_kind          ), intent(in) :: point(3)
        real (real_kind          ), intent(in) :: x_max, x_min, y_max, y_min, z_max, z_min
        logical                                :: yes
        yes = .false.
        if(x_min <= point(1) .and. point(1) <= x_max)then
            if(y_min <= point(2) .and. point(2) <= y_max)then
                if(z_min <= point(3) .and. point(3) <= z_max)then
                    yes = .true.
                end if
            end if
        end if
    end function is_point_in_region
end module measurement_surface_class