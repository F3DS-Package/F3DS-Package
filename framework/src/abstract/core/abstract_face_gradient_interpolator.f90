module abstract_face_gradient_interpolator
    ! """
    ! This module defines an abstract type `face_gradient_interpolator` and an abstract interface for initializing and interpolating
    ! face gradient variables.
    ! Module Attributes:
    !     - `face_gradient_interpolator`: An abstract type that contains deferred procedures for initializing and interpolating face
    !     gradient variables.
    ! Abstract Interface:
    !     - `initialize_interface`: A subroutine that initializes the face gradient interpolator with a given configuration.
    !         Parameters:
    !             - `self`: The face gradient interpolator object.
    !             - `config`: The configuration object.
    !     - `interpolate_interface`: A pure function that interpolates face gradient variables based on the given input.
    !         Parameters:
    !             - `self`: The face gradient interpolator object.
    !             - `lhc_gradient_variables`: The left-hand cell gradient variables.
    !             - `rhc_gradient_variables`: The right-hand cell gradient variables.
    !             - `lhc_variables`: The left-hand cell variables.
    !             - `rhc_variables`: The right-hand cell variables.
    !             - `lhc_position`: The position of the left-hand cell.
    !             - `rhc_position`: The position of the right-hand cell.
    !             - `num_variables`: The number of variables.
    !         Returns:
    !             - `face_gradient_variables`: The interpolated face gradient variables.
    ! """
    implicit none

    private

    type, public, abstract :: face_gradient_interpolator
        ! -----------------------------------------------------------------------
        ! Type: face_gradient_interpolator
        ! -----------------------------------------------------------------------
        ! This is an abstract type representing a face gradient interpolator.
        ! It contains two deferred procedures: initialize and interpolate.
        !
        ! Public Components:
        ! ------------------
        ! - initialize: A deferred procedure that initializes the interface of the interpolator.
        ! - interpolate: A deferred procedure that performs the interpolation.
        !
        ! -----------------------------------------------------------------------
        contains
        procedure(initialize_interface ), pass(self), deferred :: initialize
        procedure(interpolate_interface), pass(self), deferred :: interpolate
    end type

    abstract interface
        subroutine initialize_interface(self, config)
            ! > Initializes the interface of the face_gradient_interpolator object with the given configuration.
            ! > This subroutine initializes the interface of the face_gradient_interpolator object by taking in a configuration object. The
            ! configuration object is used to set the initial values and parameters of the face_gradient_interpolator.
            ! > @param self The face_gradient_interpolator object to be initialized.
            ! > @param config The configuration object containing the initial values and parameters.
            use abstract_configuration
            import face_gradient_interpolator

            class(face_gradient_interpolator), intent(inout) :: self
            class(configuration             ), intent(inout) :: config
        end subroutine initialize_interface

        pure function interpolate_interface(self, lhc_gradient_variables, rhc_gradient_variables, lhc_variables, rhc_variables, lhc_position, rhc_position, num_variables) result(face_gradient_variables)
            ! """
            ! Interpolates the gradient variables on the interface between two cells.
            ! Parameters:
            !     self (face_gradient_interpolator): The instance of the face_gradient_interpolator class.
            !     lhc_gradient_variables (real array): The gradient variables on the left-hand cell.
            !     rhc_gradient_variables (real array): The gradient variables on the right-hand cell.
            !     lhc_variables (real array): The variables on the left-hand cell.
            !     rhc_variables (real array): The variables on the right-hand cell.
            !     lhc_position (real array): The position of the left-hand cell.
            !     rhc_position (real array): The position of the right-hand cell.
            !     num_variables (integer): The number of variables.
            ! Returns:
            !     face_gradient_variables (real array): The interpolated gradient variables on the interface.
            ! """
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