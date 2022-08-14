module class_green_gasuu
    use typedef_module
    use abstract_gradient_calculator
    use abstract_configuration

    implicit none

    private

    type, public, extends(abstract_gradient_calculator) :: green_gasuu
        private

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: compute_gradient_residual_array
        procedure, public, pass(self) :: compute_gradient_residual
    end type

    contains

    subroutine initialize(self, a_configuration)
        class(green_gasuu  ), intent(inout) :: self
        class(configuration), intent(inout) :: a_configuration
    end subroutine initialize

    pure function compute_gradient_residual_array(self, rhc_variables, lhc_variables, rhc_position, lhc_position, face_normal_vector, face_area, num_variables) result(residual)
        class  (green_gasuu), intent(inout) :: self
        real   (real_kind  ), intent(in   ) :: rhc_variables           (:)
        real   (real_kind  ), intent(in   ) :: lhc_variables           (:)
        real   (real_kind  ), intent(in   ) :: lhc_position            (:)
        real   (real_kind  ), intent(in   ) :: rhc_position            (:)
        real   (real_kind  ), intent(in   ) :: face_normal_vector      (:)
        real   (real_kind  ), intent(in   ) :: face_position           (:)
        real   (real_kind  ), intent(in   ) :: face_area
        integer(int_kind   ), intent(in   ) :: num_variables
        real   (real_kind  )                :: residual                (num_variables*3)

        integer(int_kind) :: i

        do i = 0, num_variables-1, 1
            residual(3*i+1:3*i+3) = compure_gradient(rhc_variables(i+1), lhc_variables(i+1), cell_centor_position, cell_volume, face_normal_vector, face_area)
        end do
    end function compute_gradient_residual_array

    pure function compure_gradient_residual(self, rhc_variable, lhc_variable, rhc_position, lhc_position, face_normal_vector, face_area) result(residual)
        class  (green_gasuu), intent(inout) :: self
        real   (real_kind  ), intent(in   ) :: rhc_variable
        real   (real_kind  ), intent(in   ) :: lhc_variable
        real   (real_kind  ), intent(in   ) :: lhc_position            (:)
        real   (real_kind  ), intent(in   ) :: rhc_position            (:)
        real   (real_kind  ), intent(in   ) :: face_normal_vector      (:)
        real   (real_kind  ), intent(in   ) :: face_position           (:)
        real   (real_kind  ), intent(in   ) :: face_area
        real   (real_kind  )                :: residual                (3)

        residual(:) = (1.d0 / cell_volume) * 0.5d0 * (rhc_variable + lhc_variable) * face_area * face_normal_vector(:)
    end function compure_gradient_residual
end module class_green_gasuu