module abstract_configuration
    ! This module defines an abstract configuration type and its associated interface procedures.
    ! The `configuration` type is an abstract type that serves as a base for defining specific configuration types. It contains
    ! deferred procedures for parsing a configuration file and closing the configuration. It also contains deferred procedures for
    ! retrieving values of different types from the configuration.
    ! Interface Procedures:
    ! - `parse_interface(self, filepath)`: Parses the configuration file specified by `filepath` and updates the configuration object
    ! `self`.
    ! - `close_interface(self)`: Closes the configuration object `self`.
    ! - `get_real_interface(self, tag, val, found, default)`: Retrieves a real value associated with the given `tag` from the
    ! configuration object `self`. If the value is found, it is returned in `val`. If `found` is provided, it is set to `.true.` if
    ! the value is found, otherwise `.false.`. If `default` is provided, it is used as the default value if the value is not found.
    ! - `get_int_interface(self, tag, val, found, default)`: Retrieves an integer value associated with the given `tag` from the
    ! configuration object `self`. If the value is found, it is returned in `val`. If `found` is provided, it is set to `.true.` if
    ! the value is found, otherwise `.false.`. If `default` is provided, it is used as the default value if the value is not found.
    ! - `get_bool_interface(self, tag, val, found, default)`: Retrieves a boolean value associated with the given `tag` from the
    ! configuration object `self`. If the value is found, it is returned in `val`. If `found` is provided, it is set to `.true.` if
    ! the value is found, otherwise `.false.`. If `default` is provided, it is used as the default value if the value is not found.
    ! - `get_char_interface(self, tag, val, found, default)`: Retrieves a character value associated with the given `tag` from the
    ! configuration object `self`. If the value is found, it is returned in `val`. If `found` is provided, it is set to `.true.` if
    ! the value is found, otherwise `.false.`. If `default` is provided, it is used as the default value if the value is not found.
    ! - `get_real_array_interface(self, tag, val, found, default)`: Retrieves an array of real values associated with the given `tag`
    ! from the configuration object `self`. If the value is found, it is returned in `val`. If `found` is provided, it is set to
    ! `.true.` if the value is found, otherwise `.false.`. If `default` is provided, it is used as the default value if the value is
    ! not found.
    implicit none

    private

    type, public, abstract :: configuration
        ! Fortran type representing a configuration.
        ! Attributes:
        ! - parse: A deferred procedure that parses the configuration.
        ! - close: A deferred procedure that closes the configuration.
        ! - get_real: A deferred procedure that retrieves a real value from the configuration.
        ! - get_int: A deferred procedure that retrieves an integer value from the configuration.
        ! - get_bool: A deferred procedure that retrieves a boolean value from the configuration.
        ! - get_char: A deferred procedure that retrieves a character value from the configuration.
        ! - get_real_array: A deferred procedure that retrieves an array of real values from the configuration.
        contains

        procedure(parse_interface), pass(self), deferred :: parse
        procedure(close_interface), pass(self), deferred :: close

        procedure(get_real_interface), pass(self), deferred :: get_real
        procedure(get_int_interface ), pass(self), deferred :: get_int
        procedure(get_bool_interface), pass(self), deferred :: get_bool
        procedure(get_char_interface), pass(self), deferred :: get_char

        procedure(get_real_array_interface), pass(self), deferred :: get_real_array

        !procedure(add_real_interface), pass(self), deferred :: add_real
        !procedure(add_int_interface ), pass(self), deferred :: add_int
        !procedure(add_bool_interface), pass(self), deferred :: add_bool
        !procedure(add_char_interface), pass(self), deferred :: add_char
    end type

    abstract interface
        subroutine parse_interface(self, filepath)
            ! """
            ! Subroutine: parse_interface
            ! Description:
            ! This subroutine is used to parse an interface file and update the configuration object accordingly.
            ! Parameters:
            ! self: inout, class(configuration)
            ! The configuration object to be updated.
            ! filepath: in, character(len=*)
            ! The path to the interface file to be parsed.
            ! """
            import configuration
            class    (configuration), intent(inout) :: self
            character(len=*        ), intent(in   ) :: filepath
        end subroutine parse_interface

        subroutine close_interface(self)
            ! This subroutine is used to close the interface of the configuration object.
            ! Parameters:
            ! - self: The configuration object to close the interface of.
            import configuration
            class(configuration), intent(inout) :: self
        end subroutine close_interface

        subroutine get_real_interface(self, tag, val, found, default)
            ! """
            !     Retrieves the value of a real variable from the configuration
            !     object.
            !     Parameters:
            !         self (configuration): The configuration object.
            !         tag (str): The tag of the real variable to retrieve.
            !         val (float): The variable to store the retrieved value.
            !         found (bool, optional): Flag indicating if the variable was
            !         found in the configuration. Default is False.
            !         default (float, optional): Default value to use if the variable
            !         is not found. Default is None.
            !     Returns:
            !         None
            !     Raises:
            !         None
            ! """
            use typedef_module
            import configuration
            class    (configuration), intent(inout)            :: self
            character(len=*        ), intent(in   )            :: tag
            real     (real_kind    ), intent(inout)            :: val
            logical                 , intent(inout), optional  :: found
            real     (real_kind    ), intent(in   ), optional  :: default
        end subroutine get_real_interface

        subroutine get_int_interface(self, tag, val, found, default)
            ! """
            !     Retrieves an integer value from the configuration object.
            !     Parameters:
            !     -----------
            !     self : configuration
            !         The configuration object to retrieve the value from.
            !     tag : str
            !         The tag associated with the integer value.
            !     val : int
            !         The variable to store the retrieved integer value.
            !     found : bool, optional
            !         A logical variable indicating whether the value was found in the
            !         configuration object.
            !     default : int, optional
            !         The default value to use if the tag is not found in the
            !         configuration object.
            !     Notes:
            !     ------
            !     - The retrieved integer value is stored in the `val` variable.
            !     - If the `found` parameter is provided, it will be set to True if
            !     the value is found, and False otherwise.
            !     - If the `default` parameter is provided and the tag is not found,
            !     the `val` variable will be set to the default value.
            ! """
            use typedef_module
            import configuration
            class    (configuration), intent(inout)            :: self
            character(len=*        ), intent(in   )            :: tag
            integer  (int_kind     ), intent(inout)            :: val
            logical                 , intent(inout), optional  :: found
            integer  (int_kind     ), intent(in   ), optional  :: default
        end subroutine get_int_interface

        subroutine get_bool_interface(self, tag, val, found, default)
            ! """
            ! Subroutine: get_bool_interface
            ! Description:
            ! This subroutine is used to retrieve a boolean value from the configuration object. It searches for the specified tag and returns
            ! the corresponding value. If the tag is found, the value is updated in the 'val' variable. If the tag is not found, the 'found'
            ! variable is set to .false. If a default value is provided, it is returned when the tag is not found.
            ! Parameters:
            ! - self: The configuration object to retrieve the boolean value from.
            ! - tag: The tag to search for in the configuration object.
            ! - val: The variable to store the retrieved boolean value.
            ! - found (optional): A logical variable to indicate if the tag was found.
            ! - default (optional): The default value to return if the tag is not found.
            ! """
            use typedef_module
            import configuration
            class    (configuration), intent(inout)            :: self
            character(len=*        ), intent(in   )            :: tag
            logical                 , intent(inout)            :: val
            logical                 , intent(inout), optional  :: found
            logical                 , intent(in   ), optional  :: default
        end subroutine get_bool_interface

        subroutine get_char_interface(self, tag, val, found, default)
            ! """
            ! Subroutine: get_char_interface
            ! Description:
            ! This subroutine is used to retrieve a character value from the configuration object based on the provided tag. If the tag is
            ! found in the configuration object, the corresponding value is returned in the 'val' variable. If the tag is not found, the
            ! 'default' value is returned instead. The 'found' variable is set to .true. if the tag is found, and .false. otherwise.
            ! Parameters:
            ! - self: The configuration object to retrieve the value from.
            ! - tag: The tag to search for in the configuration object.
            ! - val: The character variable to store the retrieved value.
            ! - found (optional): A logical variable to indicate if the tag was found.
            ! - default (optional): The default value to return if the tag is not found.
            ! """
            use typedef_module
            import configuration
            class    (configuration), intent(inout)              :: self
            character(len=*        ), intent(in   )              :: tag
            character(len=:        ), intent(inout), allocatable :: val
            logical                 , intent(inout), optional    :: found
            character(len=*        ), intent(in   ), optional    :: default
        end subroutine get_char_interface

        subroutine get_real_array_interface(self, tag, val, found, default)
            ! """
            !     Retrieves a real array value from the configuration object.
            !     Parameters:
            !         self (configuration): The configuration object.
            !         tag (str): The tag of the real array value to retrieve.
            !         val (real array): The real array value to retrieve.
            !         found (bool, optional): Indicates whether the value was found in
            !         the configuration object.
            !         default (real array, optional): The default real array value to use
            !         if the tag is not found.
            !     Returns:
            !         None
            !     Raises:
            !         None
            ! """
            use typedef_module
            import configuration
            class    (configuration), intent(inout)            :: self
            character(len=*        ), intent(in   )            :: tag
            real     (real_kind    ), intent(inout)            :: val(:)
            logical                 , intent(inout), optional  :: found
            real     (real_kind    ), intent(in   ), optional  :: default(:)
        end subroutine get_real_array_interface

        !subroutine add_real_interface(self, tag, val)
        !    use typedef_module
        !    import configuration
        !    class    (configuration), intent(inout)            :: self
        !    character(len=*        ), intent(in   )            :: tag
        !    real     (real_kind    ), intent(inout)            :: val
        !end subroutine add_real_interface
        !subroutine add_int_interface(self, tag, val)
        !    use typedef_module
        !    import configuration
        !    class    (configuration), intent(inout)            :: self
        !    character(len=*        ), intent(in   )            :: tag
        !    integer  (int_kind     ), intent(inout)            :: val
        !end subroutine add_int_interface
        !subroutine add_bool_interface(self, tag, val)
        !    use typedef_module
        !    import configuration
        !    class    (configuration), intent(inout)            :: self
        !    character(len=*        ), intent(in   )            :: tag
        !    logical                 , intent(inout)            :: val
        !end subroutine add_bool_interface
        !subroutine add_char_interface(self, tag, val)
        !    use typedef_module
        !    import configuration
        !    class    (configuration), intent(inout)            :: self
        !    character(len=*        ), intent(in   )            :: tag
        !    character(len=*        ), intent(inout)            :: val
        !end subroutine add_char_interface
    end interface
end module abstract_configuration