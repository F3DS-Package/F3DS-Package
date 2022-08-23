module abstract_divergence_calculator
    implicit none

    private

    type, public, abstract :: divergence_calculator
        contains

        procedure(initilaize_interface      ), pass(self), deferred :: initilaize
        procedure(compute_residual_interface), pass(self), deferred :: compute_residual
    end type divergence_calculator

    abstract interface
        subroutine initialize_interface(self, config)
            use abstract_configuration
            import divergence_calculator

            class(divergence_calculator), intent(inout) :: self
            class(configuration        ), intent(in   ) :: config
        end subroutine initialize_interface

        pure function compute_residual_interface(self, rhc_vector, lhc_vector, rhc_position, lhc_position, face_normal_vector, face_position, face_area) result(residual)
            use typedef_module
            import divergence_calculator

            class(divergence_calculator), intent(in) :: self
            real (real_kind            ), intent(in) :: rhc_vector         (3)
            real (real_kind            ), intent(in) :: lhc_vector         (3)
            real (real_kind            ), intent(in) :: rhc_position       (3)
            real (real_kind            ), intent(in) :: lhc_position       (3)
            real (real_kind            ), intent(in) :: face_normal_vector (3)
            real (real_kind            ), intent(in) :: face_position      (3)
            real (real_kind            ), intent(in) :: face_area
            real (real_kind            )             :: residual
        end function compute_residual_interface
    end interface
end module abstract_divergence_calculator