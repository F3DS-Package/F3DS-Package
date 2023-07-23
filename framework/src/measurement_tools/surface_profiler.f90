module class_surface_profiler
    use typedef_module
    use stdio_module
    use vector_module
    use math_constant_module
    use abstract_configuration

    implicit none

    private

    type, public :: surface_profiler
        private

        integer  (int_kind ), allocatable :: face_ids_  (:)
        integer  (int_kind ), allocatable :: lhc_or_rhc_(:) ! 0 or 1. 0 indicates helf-hand cell, and 1 indicates right-hand cell.
        character(:        ), allocatable :: output_filename_
        real     (real_kind)              :: output_timespan_
        real     (real_kind)              :: next_output_time_
        integer  (int_kind )              :: num_faces_
        logical                           :: is_enabled_

        contains

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: write
        procedure, public, pass(self) :: is_writable
    
        procedure,         pass(self) :: make_new_file
        procedure,         pass(self) :: surface_is_included_in_aabb
    end type surface_profiler

    contains

    subroutine initialize(self, config, face_to_cell_index, face_centor_positions, face_normal_vectors, is_real_cell, num_faces, num_local_cells)
        class    (surface_profiler), intent(inout) :: self
        class    (configuration   ), intent(inout) :: config
        integer  (int_kind        ), intent(in   ) :: face_to_cell_index   (:, :)
        real     (real_kind       ), intent(in   ) :: face_centor_positions(:, :)
        real     (real_kind       ), intent(in   ) :: face_normal_vectors  (:, :)
        logical                    , intent(in   ) :: is_real_cell         (:)
        integer  (int_kind        ), intent(in   ) :: num_faces
        integer  (int_kind        ), intent(in   ) :: num_local_cells

        ! from config
        real     (real_kind       ) :: frequency
        real     (real_kind       ) :: min_point(3), max_point(3)
        real     (real_kind       ) :: normal(3)

        integer  (int_kind        ), allocatable   :: tmp_face_ids  (:)
        integer  (int_kind        ), allocatable   :: tmp_lhc_or_rhc(:)
        integer  (int_kind        )                :: n_out_surface, face_index, rhc_index, lhc_index
        real     (real_kind       )                :: angle
        logical                                    :: found

        call config%get_bool("Surface profiler.Enable", self%is_enabled_, found, .false.)
        if(.not. found) call write_warring("'Surface profiler.Enable' is not found in configuration file you set. Disable this profiler.")

        if (.not. self%is_enabled_) return

        call config%get_real("Surface profiler.Frequency", frequency, found, 1.d3)
        if(.not. found) call write_warring("'Surface profiler.Frequency' is not found in configuration file you set. The default value is set.")
        self%output_timespan_      = 1.d0 / frequency

        call config%get_real("Surface profiler.Min point.x", min_point(1), found)
        if(.not. found) call call_error("'Surface profiler.Min point.x' is not found in configuration file you set.")
        call config%get_real("Surface profiler.Min point.y", min_point(2), found)
        if(.not. found) call call_error("'Surface profiler.Min point.y' is not found in configuration file you set.")
        call config%get_real("Surface profiler.Min point.z", min_point(3), found)
        if(.not. found) call call_error("'Surface profiler.Min point.z' is not found in configuration file you set.")

        call config%get_real("Surface profiler.Max point.x", max_point(1), found)
        if(.not. found) call call_error("'Surface profiler.Max point.x' is not found in configuration file you set.")
        call config%get_real("Surface profiler.Max point.y", max_point(2), found)
        if(.not. found) call call_error("'Surface profiler.Max point.y' is not found in configuration file you set.")
        call config%get_real("Surface profiler.Max point.z", max_point(3), found)
        if(.not. found) call call_error("'Surface profiler.Max point.z' is not found in configuration file you set.")

        call config%get_real("Surface profiler.Facial Direction.x", normal(1), found)
        if(.not. found) call call_error("'Surface profiler.Facial Direction' is not found in configuration file you set.")
        call config%get_real("Surface profiler.Facial Direction.y", normal(2), found)
        if(.not. found) call call_error("'Surface profiler.Facial Direction' is not found in configuration file you set.")
        call config%get_real("Surface profiler.Facial Direction.z", normal(3), found)
        if(.not. found) call call_error("'Surface profiler.Facial Direction' is not found in configuration file you set.")
        
        self%output_timespan_      = 1.d0 / frequency

        allocate(tmp_face_ids  (num_faces))
        allocate(tmp_lhc_or_rhc(num_faces))
        n_out_surface = 0
        do face_index = 1, num_faces, 1
            if(self%surface_is_included_in_aabb(face_centor_positions(:, face_index), max_point(1), min_point(1), max_point(2), min_point(2), max_point(3), min_point(3)))then
                angle = vector_angle(normal, face_normal_vectors(:, face_index))
                if(0.d0 - machine_epsilon <= angle .and. angle <= 0.d0 + machine_epsilon)then
                    lhc_index = face_to_cell_index(num_local_cells+0, face_index)
                    if(is_real_cell(lhc_index))then
                        n_out_surface = n_out_surface + 1
                        tmp_face_ids  (n_out_surface) = face_index
                        tmp_lhc_or_rhc(n_out_surface) = 0
                    end if
                else if(pi - machine_epsilon <= angle .and. angle <= pi + machine_epsilon)then
                    rhc_index = face_to_cell_index(num_local_cells+1, face_index)
                    if(is_real_cell(rhc_index))then
                        n_out_surface = n_out_surface + 1
                        tmp_face_ids  (n_out_surface) = face_index
                        tmp_lhc_or_rhc(n_out_surface) = 1
                    end if
                end if
            end if
        end do

        self%output_filename_  = "surface_profiler.dat"
        allocate(self%face_ids_  , source = tmp_face_ids  (1:n_out_surface))
        allocate(self%lhc_or_rhc_, source = tmp_lhc_or_rhc(1:n_out_surface))
        self%next_output_time_ = 0.d0
        self%num_faces_        = n_out_surface

        call self%make_new_file(self%output_filename_)
    end subroutine initialize

    subroutine write(self, time, values_set, face_to_cell_index, face_areas, num_local_cells)
        class  (surface_profiler), intent(inout) :: self
        real   (real_kind          ), intent(in   ) :: time
        real   (real_kind          ), intent(in   ) :: values_set        (:,:)
        integer(int_kind           ), intent(in   ) :: face_to_cell_index(:,:)
        real   (real_kind          ), intent(in   ) :: face_areas        (:)
        integer(int_kind           ), intent(in   ) :: num_local_cells

        integer(int_kind           )                :: unit_number, id_index, cell_index, values_index, n_output_values
        real   (real_kind          )                :: output_values(size(values_set(:,1)))

        if(time >= self%next_output_time_)then
            n_output_values = size(values_set(:,1))
            output_values(:) = 0.d0
            do id_index = 1, self%num_faces_, 1
                cell_index = face_to_cell_index(self%face_ids_(id_index), num_local_cells+self%lhc_or_rhc_(id_index))
                do values_index = 1, n_output_values, 1
                    output_values(values_index) = output_values(values_index) + values_set(values_index, cell_index) * face_areas(self%face_ids_(id_index))
                end do
            end do
            open(newunit = unit_number, file="result/"//self%output_filename_, status = 'old', position = 'append')
            write(unit_number, *) time, output_values(:)
            close(unit_number)
            self%next_output_time_ = self%next_output_time_ + self%output_timespan_
        end if
    end subroutine write

    pure function is_writable(self, time) result(juge)
        class  (surface_profiler), intent(in) :: self
        real   (real_kind       ), intent(in) :: time
        logical                           :: juge

        if (.not. self%is_enabled_) then
            juge = .false.
        else
            juge = (time >= self%next_output_time_)
        end if
    end function is_writable

    subroutine make_new_file(self, filename)
        class    (surface_profiler), intent(in) :: self
        character(len=*              ), intent(in) :: filename
        integer  (int_kind           )             :: unit_number

        open(newunit = unit_number, file = "result/"//self%output_filename_, status = 'replace')
        write(unit_number, *) "# F3DS measurement surface data"
        close(unit_number)
    end subroutine make_new_file

    pure function surface_is_included_in_aabb(self, point, x_max, x_min, y_max, y_min, z_max, z_min) result(yes)
        class(surface_profiler), intent(in) :: self
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
    end function surface_is_included_in_aabb
end module class_surface_profiler