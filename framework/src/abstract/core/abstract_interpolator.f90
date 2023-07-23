module abstract_interpolator
    ! """
    ! This module defines an abstract interpolator type and its associated interface for initializing and interpolating face
    ! variables.
    ! Module: abstract_interpolator
    ! Public Types:
    ! - interpolator: An abstract type representing an interpolator.
    ! Public Procedures:
    ! - initialize_interface: A deferred subroutine that initializes the interpolator with a given configuration.
    ! - interpolate_face_variables_interface: A deferred function that interpolates face variables based on the given input.
    ! Interface:
    ! - initialize_interface(self, config): Initializes the interpolator with the provided configuration.
    !     - Parameters:
    !         - self: The interpolator object to be initialized.
    !         - config: The configuration object used for initialization.
    ! - interpolate_face_variables_interface(self, variables, face_to_cell_index, cell_positions, face_position, face_normal_vector,
    ! num_local_cells, num_variables, gradient_variables, velosity): Interpolates face variables based on the given input.
    !     - Parameters:
    !         - self: The interpolator object used for interpolation.
    !         - variables: A 2D array of real numbers representing the variables.
    !         - face_to_cell_index: An array of integers representing the face-to-cell index.
    !         - cell_positions: A 2D array of real numbers representing the cell positions.
    !         - face_position: A 3D array of real numbers representing the face position.
    !         - face_normal_vector: A 3D array of real numbers representing the face normal vector.
    !         - num_local_cells: An integer representing the number of local cells.
    !         - num_variables: An integer representing the number of variables.
    !         - gradient_variables: An optional 2D array of real numbers representing the gradient variables.
    !         - velosity: An optional 2D array of real numbers representing the velocity.
    !     - Returns:
    !         - face_variables: A 1D array of real numbers representing the interpolated face variables.
    ! """
    implicit none

    private

    type, public, abstract :: interpolator
        ! -----------------------------------------------------------------------
        ! Type: interpolator
        ! -----------------------------------------------------------------------
        ! This is an abstract type representing an interpolator.
        !
        ! Public Components:
        ! - initialize: A deferred procedure that initializes the interpolator.
        ! - interpolate_face_variables: A deferred procedure that interpolates
        ! face variables.
        !
        ! -----------------------------------------------------------------------
        contains

        procedure(initialize_interface                ), pass(self), deferred :: initialize
        procedure(interpolate_face_variables_interface), pass(self), deferred :: interpolate_face_variables
    end type interpolator

    abstract interface
        subroutine initialize_interface(self, config)
            ! -----------------------------------------------------------------------
            ! Subroutine: initialize_interface
            ! -----------------------------------------------------------------------
            ! Description:
            ! This subroutine initializes the interface of the interpolator object
            ! by assigning the provided configuration object to the interpolator.
            !
            ! Parameters:
            ! self: interpolator (inout)
            ! - The interpolator object to be initialized.
            !
            ! config: configuration (in)
            ! - The configuration object to be assigned to the interpolator.
            !
            ! -----------------------------------------------------------------------
            use abstract_configuration
            import interpolator

            class(interpolator ), intent(inout) :: self
            class(configuration), intent(in   ) :: config
        end subroutine initialize_interface

        function interpolate_face_variables_interface(self, variables, face_to_cell_index, cell_positions, face_position, face_normal_vector, num_local_cells, num_variables, gradient_variables, velosity) result(face_variables)
            ! """
            ! Interpolates the variables from the cell centers to the face center using the provided interpolator.
            ! Parameters:
            !     self (interpolator): The interpolator object.
            !     variables (real(kind=real_kind), dimension(:,:)): The variables at the cell centers.
            !     face_to_cell_index (integer(kind=int_kind), dimension(:)): The mapping of face indices to cell indices.
            !     cell_positions (real(kind=real_kind), dimension(:,:)): The positions of the cell centers.
            !     face_position (real(kind=real_kind), dimension(3)): The position of the face center.
            !     face_normal_vector (real(kind=real_kind), dimension(3)): The normal vector of the face.
            !     num_local_cells (integer(kind=int_kind)): The number of local cells.
            !     num_variables (integer(kind=int_kind)): The number of variables.
            !     gradient_variables (real(kind=real_kind), dimension(:,:), optional): The gradient of the variables at the cell centers.
            !     velosity (real(kind=real_kind), dimension(:,:), optional): The velocity at the cell centers.
            ! Returns:
            !     face_variables (real(kind=real_kind)): The interpolated variables at the face center.
            ! """
            use typedef_module
            import interpolator

            class  (interpolator), intent(in)           :: self
            real   (real_kind   ), intent(in)           :: variables         (:,:)
            integer(int_kind    ), intent(in)           :: face_to_cell_index(:)
            real   (real_kind   ), intent(in)           :: cell_positions    (:,:)
            real   (real_kind   ), intent(in)           :: face_position     (3)
            real   (real_kind   ), intent(in)           :: face_normal_vector(3)
            integer(int_kind    ), intent(in)           :: num_local_cells
            integer(int_kind    ), intent(in)           :: num_variables
            real   (real_kind   ), intent(in), optional :: gradient_variables(:,:)
            real   (real_kind   ), intent(in), optional :: velosity          (:,:)
            real   (real_kind   )                       :: face_variables(num_variables)
        end function interpolate_face_variables_interface
    end interface
end module abstract_interpolator