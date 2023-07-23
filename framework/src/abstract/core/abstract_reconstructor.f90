module abstract_reconstructor
    ! """
    ! This module defines an abstract reconstructor and reconstructor generator for a specific application.
    ! Module Attributes:
    !     - reconstructor_generator: An abstract type that serves as a base for reconstructor generators.
    !         - generate_from_configuration: A deferred procedure that generates a reconstructor from a configuration object.
    !         - generate_from_name: A deferred procedure that generates a reconstructor from a given name.
    !     - reconstructor: An abstract type that serves as a base for reconstructors.
    !         - initialize: A deferred procedure that initializes a reconstructor.
    !         - reconstruct_lhc: A deferred procedure that reconstructs the left-hand cell.
    !         - reconstruct_rhc: A deferred procedure that reconstructs the right-hand cell.
    ! Interfaces:
    !     - generate_from_configuration_interface: An abstract interface for generating a reconstructor from a configuration object.
    !     - generate_from_name_interface: An abstract interface for generating a reconstructor from a given name.
    !     - initialize_interface: An abstract interface for initializing a reconstructor.
    !     - reconstruct_lhc_interface: An abstract interface for reconstructing the left-hand cell.
    !     - reconstruct_rhc_interface: An abstract interface for reconstructing the right-hand cell.
    ! """
    implicit none

    private

    type, public, abstract :: reconstructor_generator
        ! """
        ! Fortran type: reconstructor_generator
        ! This is an abstract type that serves as a base for reconstructor generator types.
        ! Attributes:
        ! - generate: A generic procedure that can be used to generate a reconstructor from either a configuration or a name.
        ! Methods:
        ! - generate_from_configuration: A deferred procedure that generates a reconstructor from a configuration.
        ! - generate_from_name: A deferred procedure that generates a reconstructor from a name.
        ! """
        contains
        generic, public :: generate => generate_from_configuration, generate_from_name
        procedure(generate_from_configuration_interface), pass(self), deferred :: generate_from_configuration
        procedure(generate_from_name_interface         ), pass(self), deferred :: generate_from_name
    end type reconstructor_generator

    type, public, abstract :: reconstructor
        ! -----------------------------------------------------------------------
        ! reconstructor
        ! -----------------------------------------------------------------------
        ! This abstract type represents a reconstructor object. It defines the
        ! interface for initializing and reconstructing left-hand and right-hand
        ! components.
        !
        ! Public Components:
        ! ------------------
        ! initialize : deferred procedure
        ! This procedure initializes the reconstructor object.
        !
        ! reconstruct_lhc : deferred procedure
        ! This procedure reconstructs the left-hand component.
        !
        ! reconstruct_rhc : deferred procedure
        ! This procedure reconstructs the right-hand component.
        !
        ! -----------------------------------------------------------------------
        contains

        procedure(initialize_interface     ), pass(self), deferred :: initialize
        procedure(reconstruct_lhc_interface), pass(self), deferred :: reconstruct_lhc
        procedure(reconstruct_rhc_interface), pass(self), deferred :: reconstruct_rhc
    end type

    abstract interface
        subroutine generate_from_configuration_interface(self, a_reconstructor, a_config)
            ! This subroutine generates a reconstructor from a given configuration.
            !
            ! Parameters:
            ! self: reconstructor_generator
            ! The reconstructor generator object.
            ! a_reconstructor: reconstructor pointer
            ! The pointer to the reconstructor object.
            ! a_config: configuration
            ! The configuration object.
            !
            ! Returns:
            ! None
            !
            ! Notes:
            ! - This subroutine requires the use of the abstract_configuration module.
            ! - It also requires the import of the reconstructor_generator and reconstructor modules.
            !
            use abstract_configuration
            import reconstructor_generator
            import reconstructor
            class(reconstructor_generator),          intent(inout) :: self
            class(reconstructor          ), pointer, intent(inout) :: a_reconstructor
            class(configuration          ),          intent(inout) :: a_config
        end subroutine generate_from_configuration_interface

        subroutine generate_from_name_interface(self, a_reconstructor, name)
            ! ---------------------------------------------------------------------
            ! This subroutine generates a reconstructor from a given name.
            !
            ! Parameters:
            ! - self: an inout object of type reconstructor_generator
            ! - a_reconstructor: an inout pointer to a reconstructor object
            ! - name: an input allocatable character string representing the name
            ! of the reconstructor to be generated
            !
            ! ---------------------------------------------------------------------
            use abstract_configuration
            import reconstructor_generator
            import reconstructor
            class    (reconstructor_generator),              intent(inout) :: self
            class    (reconstructor          ), pointer    , intent(inout) :: a_reconstructor
            character(len=:                  ), allocatable, intent(in   ) :: name
        end subroutine generate_from_name_interface
    end interface

    abstract interface
        subroutine initialize_interface(self, config, a_reconstructor_generator)
            ! -----------------------------------------------------------------------
            ! Subroutine: initialize_interface
            ! -----------------------------------------------------------------------
            ! Description:
            ! This subroutine initializes the interface of the reconstructor object.
            !
            ! Parameters:
            ! self: reconstructor (inout)
            ! The reconstructor object to be initialized.
            !
            ! config: configuration (inout)
            ! The configuration object to be used for initialization.
            !
            ! a_reconstructor_generator: reconstructor_generator (optional, inout)
            ! An optional reconstructor generator object to be used for initialization.
            !
            ! -----------------------------------------------------------------------
            use abstract_configuration
            import reconstructor
            import reconstructor_generator
            class(reconstructor          ),           intent(inout) :: self
            class(configuration          ),           intent(inout) :: config
            class(reconstructor_generator), optional, intent(inout) :: a_reconstructor_generator
        end subroutine initialize_interface

        pure function reconstruct_lhc_interface( &
            self                               , &
            primitive_variables_set            , &
            face_to_cell_index                 , &
            cell_centor_positions              , &
            face_centor_positions              , &
            face_index                         , &
            num_local_cells                    , &
            num_primitive_variables                ) result(reconstructed_primitive_variables)
            ! """
            ! Reconstructs the primitive variables at the interface of the Local Hydrodynamic Code (LHC).
            ! Parameters:
            !     self (reconstructor): An instance of the reconstructor class.
            !     primitive_variables_set (real array): Array containing the primitive variables set.
            !     face_to_cell_index (integer array): Array mapping face indices to cell indices.
            !     cell_centor_positions (real array): Array containing the positions of cell centers.
            !     face_centor_positions (real array): Array containing the positions of face centers.
            !     face_index (integer): Index of the face.
            !     num_local_cells (integer): Number of local cells.
            !     num_primitive_variables (integer): Number of primitive variables.
            ! Returns:
            !     reconstructed_primitive_variables (real array): Array containing the reconstructed primitive variables.
            ! Note:
            !     This function is a pure function, meaning it does not modify any of its input arguments.
            ! """

            use typedef_module
            import reconstructor

            class  (reconstructor), intent(in) :: self
            integer(int_kind     ), intent(in) :: face_index
            real   (real_kind    ), intent(in) :: primitive_variables_set          (:, :)
            integer(int_kind     ), intent(in) :: face_to_cell_index               (:, :)
            real   (real_kind    ), intent(in) :: cell_centor_positions            (:, :)
            real   (real_kind    ), intent(in) :: face_centor_positions            (:, :)
            integer(int_kind     ), intent(in) :: num_local_cells
            integer(int_kind     ), intent(in) :: num_primitive_variables
            real   (real_kind    )             :: reconstructed_primitive_variables(num_primitive_variables)
        end function reconstruct_lhc_interface

        pure function reconstruct_rhc_interface( &
            self                               , &
            primitive_variables_set            , &
            face_to_cell_index                 , &
            cell_centor_positions              , &
            face_centor_positions              , &
            face_index                         , &
            num_local_cells                    , &
            num_primitive_variables                ) result(reconstructed_primitive_variables)
            ! """
            ! Reconstructs the right-hand-side (RHS) interface of the primitive variables for a given face index.
            ! Parameters:
            !     self (reconstructor): The reconstructor object.
            !     primitive_variables_set (real, dimension(:,:)): Array of primitive variables.
            !     face_to_cell_index (integer, dimension(:,:)): Array mapping face indices to cell indices.
            !     cell_centor_positions (real, dimension(:,:)): Array of cell center positions.
            !     face_centor_positions (real, dimension(:,:)): Array of face center positions.
            !     face_index (integer): Index of the face for which the RHS interface is reconstructed.
            !     num_local_cells (integer): Number of local cells.
            !     num_primitive_variables (integer): Number of primitive variables.
            ! Returns:
            !     reconstructed_primitive_variables (real, dimension(:)): Array of reconstructed primitive variables for the RHS interface.
            ! """

            use typedef_module
            import reconstructor

            class  (reconstructor), intent(in) :: self
            integer(int_kind     ), intent(in) :: face_index
            real   (real_kind    ), intent(in) :: primitive_variables_set          (:, :)
            integer(int_kind     ), intent(in) :: face_to_cell_index               (:, :)
            real   (real_kind    ), intent(in) :: cell_centor_positions            (:, :)
            real   (real_kind    ), intent(in) :: face_centor_positions            (:, :)
            integer(int_kind     ), intent(in) :: num_local_cells
            integer(int_kind     ), intent(in) :: num_primitive_variables
            real   (real_kind    )             :: reconstructed_primitive_variables(num_primitive_variables)
        end function reconstruct_rhc_interface
    end interface
end module abstract_reconstructor