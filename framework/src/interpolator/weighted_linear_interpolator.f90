module class_weighted_linear_interpolator
    use typedef_module
    use abstract_interpolator
    use abstract_configuration
    use vector_module

    implicit none

    private

    type, public, extends(interpolator) :: weighted_linear_interpolator
        private

        contains

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: interpolate_face_variables
    end type weighted_linear_interpolator

    contains

    subroutine initialize(self, config)
        class(weighted_linear_interpolator), intent(inout) :: self
        class(configuration            ), intent(in   ) :: config
    end subroutine initialize

    function interpolate_face_variables(self, variables, face_to_cell_index, cell_positions, face_position, face_normal_vector, num_local_cells, num_variables, gradient_variables, velosity) result(face_variables)
        class  (weighted_linear_interpolator), intent(in)           :: self
        real   (real_kind                   ), intent(in)           :: variables         (:,:)
        integer(int_kind                    ), intent(in)           :: face_to_cell_index(:)
        real   (real_kind                   ), intent(in)           :: cell_positions    (:,:)
        real   (real_kind                   ), intent(in)           :: face_position     (3)
        real   (real_kind                   ), intent(in)           :: face_normal_vector(3)
        integer(int_kind                    ), intent(in)           :: num_local_cells
        integer(int_kind                    ), intent(in)           :: num_variables
        real   (real_kind                   ), intent(in), optional :: gradient_variables(:,:)
        real   (real_kind                   ), intent(in), optional :: velosity          (:,:)
        real   (real_kind                   )                       :: face_variables(num_variables)

        real   (real_kind) :: rhc_w, lhc_w, total_w

        associate(lhc_index => face_to_cell_index(num_local_cells+0), rhc_index => face_to_cell_index(num_local_cells+1))
            total_w = vector_magnitude(cell_positions(1:3, lhc_index) - cell_positions(1:3, rhc_index))
            lhc_w   = vector_magnitude(cell_positions(1:3, lhc_index) - face_position(1:3)) / total_w
            rhc_w   = vector_magnitude(cell_positions(1:3, rhc_index) - face_position(1:3)) / total_w
            face_variables(:) = rhc_w * variables(:, rhc_index) + lhc_w * variables(:, lhc_index)
        end associate
    end function interpolate_face_variables
end module class_weighted_linear_interpolator