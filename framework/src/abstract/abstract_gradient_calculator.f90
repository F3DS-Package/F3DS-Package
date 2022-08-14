module abstract_gradient_calculator
    implicit none

    private

    type, public, abstract :: gradient_calculator
        contains

        procedure(initialize_interface                     ), pass(self), deferred :: initialize
        procedure(compute_gradient_residual_array_interface), pass(self), deferred :: compute_gradient_residual_array
        procedure(compure_gradient_residual_interface      ), pass(self), deferred :: compure_gradient_residual
    end type

    abstract interface

        subroutine initialize_interface(self, a_configuration)
            use abstract_configuration
            import gradient_calculator
            class(gradient_calculator), intent(inout) :: self
            class(configuration      ), intent(inout) :: a_configuration
        end subroutine initialize_interface

        pure function compute_gradient_residual_array_interface(self, rhc_variables, lhc_variables, rhc_position, lhc_position, face_normal_vector, face_area, num_variables) result(residual)
            use typedef_module
            import gradient_calculator
            class  (gradient_calculator), intent(inout) :: self
            real   (real_kind          ), intent(in   ) :: lhc_variables           (:)
            real   (real_kind          ), intent(in   ) :: rhc_variables           (:)
            real   (real_kind          ), intent(in   ) :: lhc_position            (:)
            real   (real_kind          ), intent(in   ) :: rhc_position            (:)
            real   (real_kind          ), intent(in   ) :: face_normal_vector      (:)
            real   (real_kind          ), intent(in   ) :: face_position           (:)
            real   (real_kind          ), intent(in   ) :: face_area
            integer(int_kind           ), intent(in   ) :: num_variables
            real   (real_kind          )                :: residual                (num_variables)
        end function compute_gradient_residual_array_interface

        pure function compure_gradient_residual_interface(self, rhc_variable, lhc_variable, rhc_position, lhc_position, face_normal_vector, face_area) result(residual)
            use typedef_module
            import gradient_calculator
            class  (gradient_calculator), intent(inout) :: self
            real   (real_kind          ), intent(in   ) :: rhc_variable
            real   (real_kind          ), intent(in   ) :: lhc_variable
            real   (real_kind          ), intent(in   ) :: lhc_position            (:)
            real   (real_kind          ), intent(in   ) :: rhc_position            (:)
            real   (real_kind          ), intent(in   ) :: face_normal_vector      (:)
            real   (real_kind          ), intent(in   ) :: face_position           (:)
            real   (real_kind          ), intent(in   ) :: face_area
            real   (real_kind          )                :: residual
        end function compure_gradient_residual_interface
    end interface
end module abstract_gradient_calculator