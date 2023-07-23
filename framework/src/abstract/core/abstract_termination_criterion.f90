module abstract_termination_criterion
    ! """
    ! This module defines two abstract types: `termination_criterion_generator` and `termination_criterion`.
    ! `termination_criterion_generator` is a public abstract type that contains two deferred procedures: `generate_from_configuration`
    ! and `generate_from_name`. It also has a generic procedure `generate` that can be used to call either
    ! `generate_from_configuration` or `generate_from_name`.
    ! `termination_criterion` is another public abstract type that contains two deferred procedures: `initialize` and `is_satisfied`.
    ! There are two abstract interfaces defined in this module:
    ! 1. `generate_from_configuration_interface`: This interface defines a subroutine that takes a `termination_criterion_generator`,
    ! a `termination_criterion`, and a `configuration` as input. It is used to generate a termination criterion from a configuration.
    ! 2. `generate_from_name_interface`: This interface defines a subroutine that takes a `termination_criterion_generator`, a
    ! `termination_criterion`, and a `name` as input. It is used to generate a termination criterion from a name.
    ! There are two more abstract interfaces defined in this module:
    ! 1. `initialize_interface`: This interface defines a subroutine that takes a `termination_criterion`, a `configuration`, and an
    ! optional `termination_criterion_generator` as input. It is used to initialize a termination criterion.
    ! 2. `is_satisfied_interface`: This interface defines a pure function that takes a `termination_criterion`, a `time`, and a
    ! `num_steps` as input and returns a logical value indicating whether the termination criterion is satisfied.
    ! Note: This module assumes the existence of other modules such as `abstract_configuration` and `typedef_module` which are not
    ! defined here.
    ! """
    implicit none

    type, public, abstract :: termination_criterion_generator
        ! """
        ! Fortran type representing a termination criterion generator.
        ! Attributes:
        !     generate (generic procedure): A generic procedure for generating a termination criterion.
        !     generate_from_configuration (deferred procedure): A deferred procedure for generating a termination criterion from a
        !     configuration.
        !     generate_from_name (deferred procedure): A deferred procedure for generating a termination criterion from a name.
        ! """
        contains
        generic, public :: generate => generate_from_configuration, generate_from_name
        procedure(generate_from_configuration_interface), pass(self), deferred :: generate_from_configuration
        procedure(generate_from_name_interface         ), pass(self), deferred :: generate_from_name
    end type termination_criterion_generator

    type, public, abstract :: termination_criterion
        ! -----------------------------------------------------------------------
        ! Type: termination_criterion
        ! -----------------------------------------------------------------------
        ! This is an abstract type representing a termination criterion for a
        ! program or algorithm. It contains two deferred procedures:
        !
        ! - initialize: This procedure is responsible for initializing the
        ! termination criterion.
        !
        ! - is_satisfied: This procedure is responsible for determining whether
        ! the termination criterion has been satisfied.
        !
        ! Note: This type is public and abstract, meaning that it can be
        ! inherited by other types but cannot be instantiated directly.
        ! -----------------------------------------------------------------------
        contains

        procedure(initialize_interface  ), pass(self), deferred :: initialize
        procedure(is_satisfied_interface), pass(self), deferred :: is_satisfied
    end type termination_criterion

    abstract interface
        subroutine generate_from_configuration_interface(self, a_termination_criterion, a_config)
            ! This subroutine generates a termination criterion from a given configuration.
            !
            ! Parameters:
            ! self: termination_criterion_generator
            ! The termination criterion generator object.
            ! a_termination_criterion: termination_criterion
            ! The termination criterion object to be generated.
            ! a_config: configuration
            ! The configuration object used to generate the termination criterion.
            !
            ! Returns:
            ! None
            !
            ! Notes:
            ! - The termination criterion is generated based on the given configuration.
            ! - The generated termination criterion is stored in the a_termination_criterion object.
            !
            ! Example usage:
            ! call generate_from_configuration_interface(self, a_termination_criterion, a_config)
            !
            use abstract_configuration
            import termination_criterion_generator
            import termination_criterion
            class(termination_criterion_generator), intent(inout) :: self
            class(termination_criterion          ), intent(inout) :: a_termination_criterion
            class(configuration                  ), intent(inout) :: a_config
        end subroutine generate_from_configuration_interface

        subroutine generate_from_name_interface(self, a_termination_criterion, name)
            ! ---------------------------------------------------------------------------
            ! This subroutine generates a termination criterion from a given name.
            !
            ! Parameters:
            ! self: termination_criterion_generator
            ! An instance of the termination_criterion_generator class.
            ! a_termination_criterion: termination_criterion
            ! An instance of the termination_criterion class.
            ! name: character(len=:), allocatable
            ! The name used to generate the termination criterion.
            !
            ! ---------------------------------------------------------------------------
            use abstract_configuration
            import termination_criterion_generator
            import termination_criterion
            class    (termination_criterion_generator),              intent(inout) :: self
            class    (termination_criterion          ),              intent(inout) :: a_termination_criterion
            character(len=:                          ), allocatable, intent(in   ) :: name
        end subroutine generate_from_name_interface
    end interface

    abstract interface
        subroutine initialize_interface(self, config, a_termination_criterion_generator)
            ! ---------------------------------------------------------------------------
            ! This subroutine initializes the interface for the termination criterion.
            !
            ! Parameters:
            ! self: termination_criterion (inout)
            ! The termination criterion object to be initialized.
            ! config: configuration (inout)
            ! The configuration object to be used for initialization.
            ! a_termination_criterion_generator: termination_criterion_generator (optional, inout)
            ! The termination criterion generator object to be used for initialization.
            !
            ! ---------------------------------------------------------------------------
            use abstract_configuration
            import termination_criterion
            import termination_criterion_generator

            class(termination_criterion          ),           intent(inout) :: self
            class(configuration                  ),           intent(inout) :: config
            class(termination_criterion_generator), optional, intent(inout) :: a_termination_criterion_generator
        end subroutine initialize_interface

        pure function is_satisfied_interface(self, time, num_steps) result(judge)
            ! """
            ! Check if the termination criterion is satisfied.
            ! Parameters:
            !     self (termination_criterion): The termination criterion object.
            !     time (real): The current time.
            !     num_steps (integer): The number of steps taken.
            ! Returns:
            !     judge (logical): True if the termination criterion is satisfied, False otherwise.
            ! """
            use typedef_module
            import termination_criterion

            class  (termination_criterion), intent(in) :: self
            real   (real_kind            ), intent(in) :: time
            integer(int_kind             ), intent(in) :: num_steps
            logical                                    :: judge
        end function is_satisfied_interface
    end interface
end module abstract_termination_criterion