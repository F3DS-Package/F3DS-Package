module abstract_time_stepping
    ! """
    ! This module defines two abstract types: `time_stepping_generator` and `time_stepping`.
    ! The `time_stepping_generator` type is a public abstract type that contains two deferred procedures:
    ! `generate_from_configuration` and `generate_from_name`. These procedures are used to generate instances of the `time_stepping`
    ! type.
    ! The `time_stepping` type is also a public abstract type that contains several deferred procedures: `initialize`,
    ! `compute_next_stage`, `prepare_time_stepping`, `cleanup_stepping`, and `get_number_of_stages`. These procedures define the
    ! behavior of a time stepping algorithm.
    ! There are two abstract interfaces defined in this module. The first interface is used by the `time_stepping_generator` type to
    ! generate instances of the `time_stepping` type from a configuration object. The second interface is used by the
    ! `time_stepping_generator` type to generate instances of the `time_stepping` type from a name.
    ! The `generate_from_configuration_interface` subroutine takes a `time_stepping_generator` object, a pointer to a `time_stepping`
    ! object, and a configuration object as input. It generates an instance of the `time_stepping` type based on the given
    ! configuration.
    ! The `generate_from_name_interface` subroutine takes a `time_stepping_generator` object, a pointer to a `time_stepping` object,
    ! and a name as input. It generates an instance of the `time_stepping` type based on the given name.
    ! The `initialize_interface` subroutine is a deferred procedure that initializes a `time_stepping` object. It takes a
    ! configuration object, the number of cells, and the number of conservative variables as input.
    ! The `compute_next_stage_interface` function is a deferred function that computes the next stage of the time stepping algorithm.
    ! It takes a `time_stepping` object, a cell index, a stage number, a time increment, an array of conservative variables, and an
    ! array of residuals as input. It returns an updated array of conservative variables.
    ! The `prepare_time_stepping_interface` subroutine is a deferred procedure that prepares a `time_stepping` object for time
    ! stepping. It takes a `time_stepping` object, a cell index, an array of conservative variables, and an array of residuals as
    ! input.
    ! The `get_number_of_stages_interface` function is a deferred function that returns the number of stages in the time stepping
    ! algorithm. It takes a `time_stepping` object as input and returns an integer.
    ! Note: This module assumes the existence of the `abstract_configuration` module and the `typedef_module` module, which define the
    ! `configuration` type and various type definitions, respectively.
    ! """
    implicit none

    private

    type, public, abstract :: time_stepping_generator
        ! """
        ! Fortran type representing a time stepping generator.
        ! Attributes:
        ! - generate_from_configuration: A deferred procedure that generates a time stepping generator from a configuration.
        ! - generate_from_name: A deferred procedure that generates a time stepping generator from a name.
        ! Methods:
        ! - generate: A generic procedure that calls either `generate_from_configuration` or `generate_from_name` based on the input
        ! arguments.
        ! Note: This type is abstract and cannot be instantiated directly.
        ! """
        contains
        generic, public :: generate => generate_from_configuration, generate_from_name
        procedure(generate_from_configuration_interface), pass(self), deferred :: generate_from_configuration
        procedure(generate_from_name_interface         ), pass(self), deferred :: generate_from_name
    end type time_stepping_generator

    type, public, abstract :: time_stepping
        ! """
        ! Fortran type representing a time stepping algorithm.
        ! Attributes:
        !     initialize (procedure): Abstract method to initialize the time stepping algorithm.
        !     compute_next_stage (procedure): Abstract method to compute the next stage of the time stepping algorithm.
        !     prepare_time_stepping (procedure): Abstract method to prepare the time stepping algorithm.
        !     get_number_of_stages (procedure): Abstract method to get the number of stages in the time stepping algorithm.
        ! """
        contains

        procedure(initialize_interface           ), pass(self), deferred :: initialize
        procedure(compute_next_stage_interface   ), pass(self), deferred :: compute_next_stage
        procedure(prepare_time_stepping_interface), pass(self), deferred :: prepare_time_stepping
        !procedure(cleanup_stepping_interface    ), pass(self), deferred :: cleanup_stepping
        procedure(get_number_of_stages_interface ), pass(self), deferred :: get_number_of_stages
    end type time_stepping

    abstract interface
        subroutine generate_from_configuration_interface(self, a_time_stepping, a_config)
            ! This subroutine generates a time stepping object based on a given configuration object.
            ! Parameters:
            ! self: time_stepping_generator object
            ! The time stepping generator object that will generate the time stepping object.
            ! a_time_stepping: time_stepping object (pointer)
            ! The time stepping object that will be generated.
            ! a_config: configuration object
            ! The configuration object that will be used to generate the time stepping object.
            use abstract_configuration
            import time_stepping_generator
            import time_stepping
            class(time_stepping_generator),          intent(inout) :: self
            class(time_stepping          ), pointer, intent(inout) :: a_time_stepping
            class(configuration          ),          intent(inout) :: a_config
        end subroutine generate_from_configuration_interface

        subroutine generate_from_name_interface(self, a_time_stepping, name)
            ! """
            ! Generates a time stepping object based on the given name.
            ! Parameters:
            !     self (time_stepping_generator): The time stepping generator object.
            !     a_time_stepping (time_stepping): A pointer to the time stepping object.
            !     name (str): The name used to generate the time stepping object.
            ! """
            use abstract_configuration
            import time_stepping_generator
            import time_stepping
            class    (time_stepping_generator),              intent(inout) :: self
            class    (time_stepping          ), pointer    , intent(inout) :: a_time_stepping
            character(len=:                  ), allocatable, intent(in   ) :: name
        end subroutine generate_from_name_interface
    end interface

    abstract interface
        subroutine initialize_interface(self, config, num_cells, num_conservative_variables)
            ! """
            ! Initializes the interface of the time_stepping object.
            ! Parameters:
            !     self (time_stepping): The time_stepping object to be initialized.
            !     config (configuration): The configuration object containing the necessary parameters.
            !     num_cells (int): The number of cells in the simulation.
            !     num_conservative_variables (int): The number of conservative variables in the simulation.
            ! """
            use abstract_configuration
            use typedef_module
            import time_stepping
            class  (time_stepping), intent(inout) :: self
            class  (configuration), intent(in   ) :: config
            integer(int_kind     ), intent(in   ) :: num_cells
            integer(int_kind     ), intent(in   ) :: num_conservative_variables
        end subroutine initialize_interface

        pure function compute_next_stage_interface( &
            self                                  , &
            cell_index                            , &
            stage_num                             , &
            time_increment                        , &
            conservative_variables                , &
            residuals                                 ) result(updated_conservative_variables)
            ! """
            ! Compute the updated conservative variables at the next stage interface.
            ! Parameters:
            !     self (time_stepping): The time_stepping object.
            !     cell_index (int): The index of the cell.
            !     stage_num (int): The stage number.
            !     time_increment (float): The time increment.
            !     conservative_variables (array): The array of conservative variables.
            !     residuals (array): The array of residuals.
            ! Returns:
            !     updated_conservative_variables (array): The updated array of conservative variables.
            ! """

            use typedef_module
            import time_stepping

            class  (time_stepping      ), intent(in) :: self
            integer(int_kind           ), intent(in) :: cell_index
            integer(int_kind           ), intent(in) :: stage_num
            real   (real_kind          ), intent(in) :: time_increment
            real   (real_kind          ), intent(in) :: conservative_variables        (:)
            real   (real_kind          ), intent(in) :: residuals                     (:)
            real   (real_kind          )             :: updated_conservative_variables(1:size(conservative_variables))
        end function compute_next_stage_interface

        subroutine prepare_time_stepping_interface(   &
            self                                    , &
            cell_index                              , &
            conservative_variables                  , &
            residuals                                   )
            ! ----------------------------------------------------------------------------------------------
            ! This subroutine prepares the interface for time stepping.
            !
            ! Parameters:
            ! self: The time_stepping object.
            ! cell_index: The index of the cell.
            ! conservative_variables: Array of conservative variables.
            ! residuals: Array of residuals.
            ! ----------------------------------------------------------------------------------------------
            use typedef_module
            import time_stepping
            class  (time_stepping), intent(inout) :: self
            integer(int_kind     ), intent(in   ) :: cell_index
            real   (real_kind    ), intent(in   ) :: conservative_variables(:)
            real   (real_kind    ), intent(in   ) :: residuals             (:)
        end subroutine prepare_time_stepping_interface

        pure function get_number_of_stages_interface(self) result(n)
            ! -----------------------------------------------------------------------
            ! Function: get_number_of_stages_interface
            ! -----------------------------------------------------------------------
            ! Description:
            ! This function returns the number of stages in the time stepping
            ! algorithm.
            !
            ! Parameters:
            ! self: The time_stepping object for which the number of stages is
            ! to be determined.
            !
            ! Returns:
            ! n: The number of stages in the time stepping algorithm.
            !
            ! -----------------------------------------------------------------------
            use typedef_module
            import time_stepping

            class  (time_stepping), intent(in   ) :: self
            integer(int_kind     )                :: n
        end function get_number_of_stages_interface
    end interface
end module abstract_time_stepping