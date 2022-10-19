module abstract_face_gradient_interpolator
    implicit none

    private

    type, public, abstract :: face_gradient_interpolator
        contains
        procedure(initialize_interface ), pass(self), deferred :: initialize
        procedure(interpolate_interface), pass(self), deferred :: interpolate
    end type

    abstract interface
        subroutine initialize_interface(self, config)
            use abstract_configuration
            import face_gradient_interpolator

            class(face_gradient_interpolator), intent(inout) :: self
            class(configuration             ), intent(inout) :: config
        end subroutine initialize_interface

        pure function interpolate_interface(self, lhc_gradient_variables, rhc_gradient_variables, lhc_variables, rhc_variables, lhc_position, rhc_position, num_variables) result(face_gradient_variables)
            use typedef_module
            import face_gradient_interpolator

            class  (face_gradient_interpolator), intent(in) :: self
            real   (real_kind                 ), intent(in) :: lhc_gradient_variables(:)
            real   (real_kind                 ), intent(in) :: rhc_gradient_variables(:)
            real   (real_kind                 ), intent(in) :: lhc_variables         (:)
            real   (real_kind                 ), intent(in) :: rhc_variables         (:)
            real   (real_kind                 ), intent(in) :: lhc_position          (3)
            real   (real_kind                 ), intent(in) :: rhc_position          (3)
            integer(int_kind                  ), intent(in) :: num_variables

            real   (real_kind                 )             :: face_gradient_variables(num_variables*3)
        end function interpolate_interface
    end interface
end module abstract_face_gradient_interpolator