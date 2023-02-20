module abstract_face_gradient_calculator
    implicit none

    private

    type, public, abstract :: face_gradient_calculator
        contains
        procedure(initialize_interface), pass(self), deferred :: initialize
        procedure(compute_interface   ), pass(self), deferred :: compute
    end type

    abstract interface
        subroutine initialize_interface(self, config)
            use abstract_configuration
            import face_gradient_calculator

            class(face_gradient_calculator), intent(inout) :: self
            class(configuration           ), intent(inout) :: config
        end subroutine initialize_interface

        pure function compute_interface(self, lhc_variables, rhc_variables, lhc_position, rhc_position, face_normal_vector, num_variables) result(gradient_variables)
            use typedef_module
            import face_gradient_calculator

            class  (face_gradient_calculator), intent(in) :: self
            real   (real_kind               ), intent(in) :: lhc_variables(:)
            real   (real_kind               ), intent(in) :: rhc_variables(:)
            real   (real_kind               ), intent(in) :: lhc_position (3)
            real   (real_kind               ), intent(in) :: rhc_position (3)
            real   (real_kind               ), intent(in) :: face_normal_vector(3)
            integer(int_kind                ), intent(in) :: num_variables
            real   (real_kind               )             :: gradient_variables(num_variables*3)
        end function compute_interface
    end interface
end module abstract_face_gradient_calculator