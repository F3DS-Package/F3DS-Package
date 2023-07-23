module abstract_initial_condition_parser
    ! """
    ! This module defines an abstract type `initial_condition_parser` for parsing initial conditions.
    ! The `initial_condition_parser` type is an abstract type that contains deferred procedures for parsing initial conditions,
    ! closing the parser, and getting the set of conservative variables.
    ! Public Members:
    ! - `parse`: Deferred procedure for parsing initial conditions.
    ! - `close`: Deferred procedure for closing the parser.
    ! - `get_conservative_variables_set`: Deferred procedure for getting the set of conservative variables.
    ! Abstract Interface:
    ! - `parse_interface`: Subroutine for parsing initial conditions. Takes an `initial_condition_parser` object and a `configuration`
    ! object as arguments.
    ! - `close_interface`: Subroutine for closing the parser. Takes an `initial_condition_parser` object as an argument.
    ! - `get_conservative_variables_set_interface`: Subroutine for getting the set of conservative variables. Takes an
    ! `initial_condition_parser` object and a 2D array of `real` values as arguments.
    ! """
    implicit none

    private

    type, public, abstract :: initial_condition_parser
        ! """
        ! Fortran type representing an initial condition parser.
        ! Attributes:
        !     parse (procedure): Deferred procedure for parsing the initial condition.
        !     close (procedure): Deferred procedure for closing the initial condition parser.
        !     get_conservative_variables_set (procedure): Deferred procedure for getting the set of conservative variables.
        ! Note:
        !     This type is public and abstract.
        ! """
        contains
        procedure(parse_interface), pass(self), deferred :: parse
        procedure(close_interface), pass(self), deferred :: close

        procedure(get_conservative_variables_set_interface), pass(self), deferred :: get_conservative_variables_set
    end type initial_condition_parser

    abstract interface
        subroutine parse_interface(self, config)
            ! > Parses the interface between the `self` object and the `config` object.
            ! > This subroutine is responsible for parsing the interface between the `self` object, which is an instance of the
            ! `initial_condition_parser` class, and the `config` object, which is an instance of the `configuration` class.
            ! > @param self An instance of the `initial_condition_parser` class. This object will be modified during the parsing process.
            ! > @param config An instance of the `configuration` class. This object will be modified during the parsing process.
            ! > @note This subroutine assumes that the `self` and `config` objects have been properly initialized before calling this
            ! subroutine.
            ! > @note The parsing process involves extracting and processing information from the `self` and `config` objects to establish a
            ! connection between them.
            ! > @note This subroutine does not return any value.
            use abstract_configuration
            import initial_condition_parser
            class(initial_condition_parser), intent(inout) :: self
            class(configuration           ), intent(inout) :: config
        end subroutine parse_interface

        subroutine close_interface(self)
            ! -----------------------------------------------------------------------
            ! close_interface subroutine
            ! -----------------------------------------------------------------------
            ! Closes the interface of the initial_condition_parser class.
            !
            ! Parameters:
            ! self: class(initial_condition_parser)
            ! The instance of the initial_condition_parser class.
            ! -----------------------------------------------------------------------
            import initial_condition_parser
            class(initial_condition_parser), intent(inout) :: self
        end subroutine close_interface

        subroutine get_conservative_variables_set_interface(self, conservative_variables_set)
            ! -----------------------------------------------------------------------
            ! Subroutine: get_conservative_variables_set_interface
            ! -----------------------------------------------------------------------
            ! Description: This subroutine is used to get the conservative variables set
            ! from the initial condition parser object.
            !
            ! Parameters:
            ! - self: An object of the initial_condition_parser class. It is passed
            ! as an input/output argument.
            ! - conservative_variables_set: A 2D array of real numbers. It represents
            ! the set of conservative variables and is
            ! passed as an input/output argument.
            !
            ! -----------------------------------------------------------------------
            use typedef_module
            import initial_condition_parser
            class(initial_condition_parser), intent(inout) :: self
            real (real_kind               ), intent(inout) :: conservative_variables_set(:,:)
        end subroutine get_conservative_variables_set_interface
    end interface
end module abstract_initial_condition_parser