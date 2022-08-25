module class_weighted_green_gasuu
    use typedef_module
    use vector_module
    use abstract_gradient_calculator
    use abstract_configuration

    implicit none

    private

    type, public, extends(gradient_calculator) :: weighted_green_gauss
        private

        contains

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: compute_residual
    end type

    contains

    subroutine initialize(self, a_configuration)
        class(weighted_green_gauss), intent(inout) :: self
        class(configuration       ), intent(inout) :: a_configuration
    end subroutine initialize

    pure function compute_residual(self, lhc_variable, rhc_variable, lhc_position, rhc_position, face_normal_vector, face_position, face_area) result(residual)
        class  (weighted_green_gauss), intent(in) :: self
        real   (real_kind           ), intent(in) :: lhc_variable
        real   (real_kind           ), intent(in) :: rhc_variable
        real   (real_kind           ), intent(in) :: lhc_position            (3)
        real   (real_kind           ), intent(in) :: rhc_position            (3)
        real   (real_kind           ), intent(in) :: face_normal_vector      (3)
        real   (real_kind           ), intent(in) :: face_position           (3)
        real   (real_kind           ), intent(in) :: face_area
        real   (real_kind           )             :: residual                (3)

        real   (real_kind) :: rhc_w, lhc_w, total_w

        total_w = vector_magnitude(lhc_position(:) - rhc_position (:))
        lhc_w   = vector_magnitude(lhc_position(:) - face_position(:)) / total_w
        rhc_w   = vector_magnitude(rhc_position(:) - face_position(:)) / total_w

        residual(:) = (rhc_w * rhc_variable + lhc_w * lhc_variable) * face_area * face_normal_vector(:)
    end function compute_residual
end module class_weighted_green_gasuu