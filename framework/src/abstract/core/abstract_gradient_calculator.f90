module abstract_gradient_calculator
    implicit none

    private

    type, public, abstract :: gradient_calculator
        contains

        procedure(initialize_interface            ), pass(self), deferred :: initialize
        procedure(compute_residual_interface      ), pass(self), deferred :: compute_residual
    end type

    abstract interface

        subroutine initialize_interface(self, a_configuration)
            use abstract_configuration
            import gradient_calculator
            class(gradient_calculator), intent(inout) :: self
            class(configuration      ), intent(inout) :: a_configuration
        end subroutine initialize_interface

        pure function compute_residual_interface(self, lhc_variable, rhc_variable, lhc_position, rhc_position, face_normal_vector, face_position, face_area) result(residual)
            use typedef_module
            import gradient_calculator
            class  (gradient_calculator), intent(in   ) :: self
            real   (real_kind          ), intent(in   ) :: lhc_variable
            real   (real_kind          ), intent(in   ) :: rhc_variable
            real   (real_kind          ), intent(in   ) :: lhc_position            (3)
            real   (real_kind          ), intent(in   ) :: rhc_position            (3)
            real   (real_kind          ), intent(in   ) :: face_normal_vector      (3)
            real   (real_kind          ), intent(in   ) :: face_position           (3)
            real   (real_kind          ), intent(in   ) :: face_area
            real   (real_kind          )                :: residual                (3)
        end function compute_residual_interface
    end interface
end module abstract_gradient_calculator