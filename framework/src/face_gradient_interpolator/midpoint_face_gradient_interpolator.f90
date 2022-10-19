module class_midpoint_face_gradient_interpolator
    use typedef_module
    use vector_module
    use abstract_configuration
    use abstract_face_gradient_interpolator

    implicit none

    type, public, extends(face_gradient_interpolator) :: midpoint_face_gradient_interpolator
        contains

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: interpolate
    end type

    contains

    subroutine initialize(self, config)
        class(midpoint_face_gradient_interpolator), intent(inout) :: self
        class(configuration                                ), intent(inout) :: config

        return
    end subroutine initialize

    pure function interpolate(self, lhc_gradient_variables, rhc_gradient_variables, lhc_variables, rhc_variables, lhc_position, rhc_position, num_variables) result(face_gradient_variables)
        class  (midpoint_face_gradient_interpolator), intent(in) :: self
        real   (real_kind                                    ), intent(in) :: lhc_gradient_variables(:)
        real   (real_kind                                    ), intent(in) :: rhc_gradient_variables(:)
        real   (real_kind                                    ), intent(in) :: lhc_variables         (:)
        real   (real_kind                                    ), intent(in) :: rhc_variables         (:)
        real   (real_kind                                    ), intent(in) :: lhc_position          (3)
        real   (real_kind                                    ), intent(in) :: rhc_position          (3)
        integer(int_kind                                     ), intent(in) :: num_variables
        real   (real_kind                                    )             :: face_gradient_variables(num_variables*3)

        integer(int_kind) :: i


        do i = 1, num_variables, 1
            associate(                                                                      &
                start_index   => 3*(i-1)+1                                                , &
                end_index     => 3*(i-1)+3                                                  &
            )
                face_gradient_variables(start_index:end_index) = 0.5d0 * (lhc_gradient_variables(start_index:end_index) + rhc_gradient_variables(start_index:end_index))
            end associate
        end do
    end function interpolate

end module class_midpoint_face_gradient_interpolator