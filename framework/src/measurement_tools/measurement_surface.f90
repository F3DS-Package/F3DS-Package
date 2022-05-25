module measurement_surface_class
    use typedef_module

    implicit none

    private

    type, public :: measurement_surface
        private

        integer  (int_kind ), allocatable :: face_ids_(:)
        character(:        ), allocatable :: output_filename_
        real     (real_kind)              :: output_timespan_
        real     (real_kind)              :: next_output_time_
        integer  (int_kind )              :: n_faces_

        contains

        procedure, public, pass(self) :: initialize => initialize_by_ids
        procedure, public, pass(self) :: write
    end type measurement_surface

    contains

    subroutine initialize_by_ids(self, filename, frequency, face_ids)
        class    (measurement_surface), intent(inout) :: self
        character(len=*              ), intent(in   ) :: filename
        real     (real_kind          ), intent(in   ) :: frequency
        integer  (int_kind           ), intent(in   ) :: face_ids(:)

        integer  (int_kind           )                :: unit_number

        self%output_filename_  = filename
        self%output_timespan_  = 1.d0 / frequency
        allocate(self%face_ids_, source = face_ids)
        self%next_output_time_ = 0.d0
        self%n_faces_          = size(face_ids)

        open(newunit = unit_number, file = self%output_filename_, status = 'replace')
        write(unit_number, *) "# F3DS measurement surface data"
        close(unit_number)
    end subroutine initialize_by_ids

    subroutine write(self, time, values_set, face_to_cell_index, face_areas, n_ghost_cells)
        class  (measurement_surface), intent(inout) :: self
        real   (real_kind          ), intent(in   ) :: time
        real   (real_kind          ), intent(in   ) :: values_set        (:,:)
        integer(int_kind           ), intent(in   ) :: face_to_cell_index(:,:)
        real   (real_kind          ), intent(in   ) :: face_areas        (:)
        integer(int_kind           ), intent(in   ) :: n_ghost_cells

        integer(int_kind           )                :: unit_number, face_index, rhc_index, lhc_index, values_index, n_output_values
        real   (real_kind          )                :: output_values(size(values_set(1,:)))

        if(time >= self%next_output_time_)then
            n_output_values = size(values_set(1,:))
            do face_index = 1, self%n_faces_, 1
                lhc_index = face_to_cell_index(face_index, n_ghost_cells+0)
                rhc_index = face_to_cell_index(face_index, n_ghost_cells+1)
                do values_index = 1, n_output_values, 1
                    output_values(values_index) = output_values(values_index) + values_set(rhc_index, values_index) * face_areas(face_index)
                end do
            end do
            open(newunit = unit_number, file=self%output_filename_, status = 'old', position='append')
            write(unit_number, *) time, output_values(:)
            close(unit_number)
            self%next_output_time_ = self%next_output_time_ + self%output_timespan_
        end if
    end subroutine write
end module measurement_surface_class