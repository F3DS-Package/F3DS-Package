module class_green_gasuu
    use typedef_module
    use abstract_gradient_calculator
    use abstract_configuration

    implicit none

    private

    type, public, extends(gradient_calculator) :: green_gasuu
        contains

        private

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: compute_residual
    end type

    contains

    subroutine initialize(self, a_configuration)
        class(green_gasuu  ), intent(inout) :: self
        class(configuration), intent(inout) :: a_configuration
    end subroutine initialize

    pure function compute_residual(self, lhc_variable, rhc_variable, lhc_position, rhc_position, face_normal_vector, face_position, face_area) result(residual)
        class  (green_gasuu), intent(in) :: self
        real   (real_kind  ), intent(in) :: lhc_variable
        real   (real_kind  ), intent(in) :: rhc_variable
        real   (real_kind  ), intent(in) :: lhc_position            (3)
        real   (real_kind  ), intent(in) :: rhc_position            (3)
        real   (real_kind  ), intent(in) :: face_normal_vector      (3)
        real   (real_kind  ), intent(in) :: face_position           (3)
        real   (real_kind  ), intent(in) :: face_area
        real   (real_kind  )             :: residual                (3)

        residual(:) = 0.5d0 * (rhc_variable + lhc_variable) * face_area * face_normal_vector(:)
    end function compute_residual
end module class_green_gasuu