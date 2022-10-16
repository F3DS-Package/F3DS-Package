module class_corrected_midpoint_face_gradient_interpolator
    ! See appendix C, H. Nishikawa, “Beyond Interface Gradient: A General Principle for Constructing Diffusion Schemes,” AIAA paper 2010-5093, 2010.

    use typedef_module
    use vector_module
    use abstract_configuration
    use abstract_face_gradient_interpolator

    implicit none

    type, public, extends(face_gradient_interpolator) :: corrected_midpoint_face_gradient_interpolator
        contains

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: interpolate
    end type

    contains

    subroutine initialize(self, config)
        class(corrected_midpoint_face_gradient_interpolator), intent(inout) :: self
        class(configuration                                ), intent(inout) :: config

        return
    end subroutine initialize

    pure function interpolate(self, lhc_gradient_variables, rhc_gradient_variables, lhc_variables, rhc_variables, lhc_position, rhc_position, num_variables) result(face_gradient_variables)
        class  (corrected_midpoint_face_gradient_interpolator), intent(in) :: self
        real   (real_kind                                    ), intent(in) :: lhc_gradient_variables(:)
        real   (real_kind                                    ), intent(in) :: rhc_gradient_variables(:)
        real   (real_kind                                    ), intent(in) :: lhc_variables         (:)
        real   (real_kind                                    ), intent(in) :: rhc_variables         (:)
        real   (real_kind                                    ), intent(in) :: lhc_position          (3)
        real   (real_kind                                    ), intent(in) :: rhc_position          (3)
        integer(int_kind                                     ), intent(in) :: num_variables
        real   (real_kind                                    )             :: face_gradient_variables(num_variables*3)

        integer(int_kind) :: i

        associate(                                                                      &
            midpoint_grad => 0.5d0 * (lhc_gradient_variables + rhc_gradient_variables), &
            magnitude     => vector_magnitude(lhc_position - rhc_position)            , &
            normalize     => vector_normalize(lhc_position - rhc_position)              &
        )
            do i = 1, num_variables, 1
                face_gradient_variables(i:i+2) = midpoint_grad(i:i+2)                                                                                               &
                                               + (vector_multiply(midpoint_grad(i:i+2), normalize) - (lhc_variables(i) - rhc_variables(i)) / magnitude) * normalize
            end do
        end associate
    end function interpolate

end module class_corrected_midpoint_face_gradient_interpolator