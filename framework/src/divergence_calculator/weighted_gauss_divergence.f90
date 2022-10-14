module class_weighted_gauss_divergence
    use typedef_module
    use abstract_divergence_calculator
    use abstract_configuration
    use vector_module

    implicit none

    private

    type, public, extends(divergence_calculator) :: weighted_gauss_divergence
        private

        contains

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: compute_residual
    end type weighted_gauss_divergence

    contains

    subroutine initialize(self, config)
        class(weighted_gauss_divergence), intent(inout) :: self
        class(configuration            ), intent(in   ) :: config
    end subroutine initialize

    pure function compute_residual(self, lhc_vector, rhc_vector, lhc_position, rhc_position, face_normal_vector, face_position, face_area) result(residual)
        class(weighted_gauss_divergence), intent(in) :: self
        real (real_kind                ), intent(in) :: lhc_vector         (3)
        real (real_kind                ), intent(in) :: rhc_vector         (3)
        real (real_kind                ), intent(in) :: lhc_position       (3)
        real (real_kind                ), intent(in) :: rhc_position       (3)
        real (real_kind                ), intent(in) :: face_normal_vector (3)
        real (real_kind                ), intent(in) :: face_position      (3)
        real (real_kind                ), intent(in) :: face_area
        real (real_kind                )             :: residual

        real   (real_kind) :: rhc_w, lhc_w, total_w

        total_w = vector_magnitude(lhc_position(:) - rhc_position (:))
        lhc_w   = vector_magnitude(lhc_position(:) - face_position(:)) / total_w
        rhc_w   = vector_magnitude(rhc_position(:) - face_position(:)) / total_w

        residual = vector_multiply(                            &
                    (rhc_w * rhc_vector + lhc_w * lhc_vector), &
                    face_normal_vector * face_area             &
                )
    end function compute_residual
end module class_weighted_gauss_divergence