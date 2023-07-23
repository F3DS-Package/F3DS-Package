module abstract_gradient_calculator
    ! """
    ! This module defines an abstract type `gradient_calculator` and an abstract interface for initializing the calculator and
    ! computing residuals.
    ! Module Attributes:
    !     - `gradient_calculator`: An abstract type that contains deferred procedures for initializing the calculator and computing
    !     residuals.
    ! Interface Functions:
    !     - `initialize_interface`: A subroutine that initializes the gradient calculator with a given configuration.
    !         Parameters:
    !             - `self`: The gradient calculator object.
    !             - `a_configuration`: The configuration object to be used for initialization.
    !     - `compute_residual_interface`: A pure function that computes the residual for a given set of variables and parameters.
    !         Parameters:
    !             - `self`: The gradient calculator object.
    !             - `lhc_variable`: The left-hand-side variable.
    !             - `rhc_variable`: The right-hand-side variable.
    !             - `lhc_position`: The left-hand-side position.
    !             - `rhc_position`: The right-hand-side position.
    !             - `face_normal_vector`: The normal vector of the face.
    !             - `face_position`: The position of the face.
    !             - `face_area`: The area of the face.
    !         Returns:
    !             - `residual`: The computed residual as a 3D vector.
    ! """
    implicit none

    private

    type, public, abstract :: gradient_calculator
        ! -----------------------------------------------------------------------
        ! Type: gradient_calculator
        ! -----------------------------------------------------------------------
        ! This is an abstract type in Fortran that represents a gradient calculator.
        ! It contains deferred procedures for initializing the interface and computing
        ! the residual. This type is intended to be inherited by other types that
        ! implement these procedures.
        !
        ! Public Components:
        ! - initialize: Deferred procedure for initializing the interface.
        ! - compute_residual: Deferred procedure for computing the residual.
        !
        ! -----------------------------------------------------------------------
        contains

        procedure(initialize_interface            ), pass(self), deferred :: initialize
        procedure(compute_residual_interface      ), pass(self), deferred :: compute_residual
    end type

    abstract interface

        subroutine initialize_interface(self, a_configuration)
            ! -----------------------------------------------------------------------
            ! Subroutine: initialize_interface
            ! -----------------------------------------------------------------------
            ! This subroutine initializes the interface between the gradient_calculator
            ! class and the configuration class. It takes in an instance of the
            ! gradient_calculator class and an instance of the configuration class as
            ! arguments.
            !
            ! Parameters:
            ! - self: An instance of the gradient_calculator class. It is passed as
            ! intent(inout), meaning it can be modified within the subroutine.
            ! - a_configuration: An instance of the configuration class. It is passed
            ! as intent(inout), meaning it can be modified within
            ! the subroutine.
            !
            ! -----------------------------------------------------------------------
            use abstract_configuration
            import gradient_calculator
            class(gradient_calculator), intent(inout) :: self
            class(configuration      ), intent(inout) :: a_configuration
        end subroutine initialize_interface

        pure function compute_residual_interface(self, lhc_variable, rhc_variable, lhc_position, rhc_position, face_normal_vector, face_position, face_area) result(residual)
            ! """
            ! Compute the residual for a given interface.
            ! Parameters:
            !     self (gradient_calculator): An instance of the gradient_calculator class.
            !     lhc_variable (real): The variable value on the left-hand cell.
            !     rhc_variable (real): The variable value on the right-hand cell.
            !     lhc_position (real array): The position of the left-hand cell (x, y, z coordinates).
            !     rhc_position (real array): The position of the right-hand cell (x, y, z coordinates).
            !     face_normal_vector (real array): The normal vector of the interface face (x, y, z components).
            !     face_position (real array): The position of the interface face (x, y, z coordinates).
            !     face_area (real): The area of the interface face.
            ! Returns:
            !     residual (real array): The computed residual (x, y, z components).
            ! """
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