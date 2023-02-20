module class_central_difference_face_gradient_calculator
    use typedef_module
    use abstract_face_gradient_calculator
    use abstract_configuration
    use vector_module

    implicit none

    private

    type, public, extends(face_gradient_calculator) :: central_difference_face_gradient_calculator
        contains
        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: compute
    end type

    contains

    subroutine initialize(self, config)
        class(central_difference_face_gradient_calculator), intent(inout) :: self
        class(configuration                              ), intent(inout) :: config
    end subroutine initialize

    pure function compute(self, lhc_variables, rhc_variables, lhc_position, rhc_position, face_normal_vector, num_variables) result(gradient_variables)
        class  (central_difference_face_gradient_calculator), intent(in) :: self
        real   (real_kind                                  ), intent(in) :: lhc_variables(:)
        real   (real_kind                                  ), intent(in) :: rhc_variables(:)
        real   (real_kind                                  ), intent(in) :: lhc_position (3)
        real   (real_kind                                  ), intent(in) :: rhc_position (3)
        real   (real_kind                                  ), intent(in) :: face_normal_vector(3)
        integer(int_kind                                   ), intent(in) :: num_variables
        real   (real_kind                                  )             :: gradient_variables(num_variables*3)

        integer(int_kind  ) :: i

        do i = 1, num_variables, 1
            associate(                                            &
                grad  => gradient_variables(3*(i-1)+1:3*(i-1)+3), &
                dphi  => lhc_variables(i) - rhc_variables(i)    , &
                dx    => lhc_position - rhc_position            , &
                theta => vector_angle(face_normal_vector, lhc_position - rhc_position) &
            )
                grad(1:3) = face_normal_vector * (1.0_real_kind / cos(theta)) * (dphi / vector_magnitude(dx))
            end associate
        end do
    end function compute
end module class_central_difference_face_gradient_calculator