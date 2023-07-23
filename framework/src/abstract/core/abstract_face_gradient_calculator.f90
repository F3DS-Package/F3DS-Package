module abstract_face_gradient_calculator
    ! """
    ! This module defines an abstract type `face_gradient_calculator` and an abstract interface for initializing and computing the
    ! gradient of variables across a face.
    ! Module Attributes:
    !     - `face_gradient_calculator`: An abstract type representing a face gradient calculator.
    !         - `initialize`: A deferred procedure for initializing the face gradient calculator.
    !         - `compute`: A deferred procedure for computing the gradient of variables across a face.
    ! Abstract Interface:
    !     - `initialize_interface`: A subroutine for initializing the face gradient calculator.
    !         - Parameters:
    !             - `self`: An inout argument of type `face_gradient_calculator` representing the instance of the calculator.
    !             - `config`: An inout argument of type `configuration` representing the configuration for the calculator.
    !     - `compute_interface`: A pure function for computing the gradient of variables across a face.
    !         - Parameters:
    !             - `self`: An input argument of type `face_gradient_calculator` representing the instance of the calculator.
    !             - `lhc_variables`: An input argument of type `real_kind` array representing the variables on the left-hand side of
    !             the face.
    !             - `rhc_variables`: An input argument of type `real_kind` array representing the variables on the right-hand side of
    !             the face.
    !             - `lhc_position`: An input argument of type `real_kind` array representing the position of the left-hand side of the
    !             face.
    !             - `rhc_position`: An input argument of type `real_kind` array representing the position of the right-hand side of
    !             the face.
    !             - `face_normal_vector`: An input argument of type `real_kind` array representing the normal vector of the face.
    !             - `num_variables`: An input argument of type `int_kind` representing the number of variables.
    !         - Returns:
    !             - `gradient_variables`: A `real_kind` array of size `num_variables*3` representing the gradient of variables across
    !             the face.
    ! """
    implicit none

    private

    type, public, abstract :: face_gradient_calculator
        ! -----------------------------------------------------------------------
        ! Type: face_gradient_calculator
        ! -----------------------------------------------------------------------
        ! This is an abstract type in Fortran used for calculating the gradient
        ! of a face in a computational fluid dynamics (CFD) simulation. It is
        ! intended to be inherited by concrete types that implement specific
        ! methods for initializing and computing the interface gradient.
        !
        ! Public Components:
        ! - initialize: Deferred procedure for initializing the interface.
        ! - compute: Deferred procedure for computing the interface gradient.
        !
        ! -----------------------------------------------------------------------
        contains
        procedure(initialize_interface), pass(self), deferred :: initialize
        procedure(compute_interface   ), pass(self), deferred :: compute
    end type

    abstract interface
        subroutine initialize_interface(self, config)
            ! > Initializes the interface of the face gradient calculator with the given configuration.
            ! > This subroutine initializes the interface of the face gradient calculator by taking in a configuration object.
            ! > The face gradient calculator object and the configuration object are passed as arguments to this subroutine.
            ! > @param self The face gradient calculator object to be initialized.
            ! > @param config The configuration object to be used for initialization.
            use abstract_configuration
            import face_gradient_calculator

            class(face_gradient_calculator), intent(inout) :: self
            class(configuration           ), intent(inout) :: config
        end subroutine initialize_interface

        pure function compute_interface(self, lhc_variables, rhc_variables, lhc_position, rhc_position, face_normal_vector, num_variables) result(gradient_variables)
            ! """
            ! Compute the interface gradient variables between two face gradient calculators.
            ! Parameters:
            !     self (face_gradient_calculator): The instance of the face gradient calculator.
            !     lhc_variables (array): The left-hand cell variables.
            !     rhc_variables (array): The right-hand cell variables.
            !     lhc_position (array): The position of the left-hand cell.
            !     rhc_position (array): The position of the right-hand cell.
            !     face_normal_vector (array): The normal vector of the face.
            !     num_variables (int): The number of variables.
            !
            ! Returns:
            !     gradient_variables (array): The computed gradient variables.
            ! """
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