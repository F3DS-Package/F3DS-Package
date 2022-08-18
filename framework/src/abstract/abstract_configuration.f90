module abstract_configuration
    implicit none

    private

    type, public, abstract :: configuration
        contains

        procedure(parse_interface), pass(self), deferred :: parse
        procedure(close_interface), pass(self), deferred :: close

        procedure(get_real_interface), pass(self), deferred :: get_real
        procedure(get_int_interface ), pass(self), deferred :: get_int
        procedure(get_bool_interface), pass(self), deferred :: get_bool
        procedure(get_char_interface), pass(self), deferred :: get_char

        procedure(add_real_interface), pass(self), deferred :: add_real
        procedure(add_int_interface ), pass(self), deferred :: add_int
        procedure(add_bool_interface), pass(self), deferred :: add_bool
        procedure(add_char_interface), pass(self), deferred :: add_char
    end type

    abstract interface
        subroutine parse_interface(self, filepath)
            import configuration
            class    (configuration), intent(inout) :: self
            character(len=*        ), intent(in   ) :: filepath
        end subroutine parse_interface

        subroutine close_interface(self)
            import configuration
            class(configuration), intent(inout) :: self
        end subroutine close_interface

        subroutine get_real_interface(self, tag, val, found, default)
            use typedef_module
            import configuration
            class    (configuration), intent(inout)            :: self
            character(len=*        ), intent(in   )            :: tag
            real     (real_kind    ), intent(inout)            :: val
            logical                 , intent(inout), optional  :: found
            real     (real_kind    ), intent(in   ), optional  :: default
        end subroutine get_real_interface

        subroutine get_int_interface(self, tag, val, found, default)
            use typedef_module
            import configuration
            class    (configuration), intent(inout)            :: self
            character(len=*        ), intent(in   )            :: tag
            integer  (int_kind     ), intent(inout)            :: val
            logical                 , intent(inout), optional  :: found
            integer  (int_kind     ), intent(in   ), optional  :: default
        end subroutine get_int_interface

        subroutine get_bool_interface(self, tag, val, found, default)
            use typedef_module
            import configuration
            class    (configuration), intent(inout)            :: self
            character(len=*        ), intent(in   )            :: tag
            logical                 , intent(inout)            :: val
            logical                 , intent(inout), optional  :: found
            logical                 , intent(in   ), optional  :: default
        end subroutine get_bool_interface

        subroutine get_char_interface(self, tag, val, found, default)
            use typedef_module
            import configuration
            class    (configuration), intent(inout)              :: self
            character(len=*        ), intent(in   )              :: tag
            character(len=:        ), intent(inout), allocatable :: val
            logical                 , intent(inout), optional    :: found
            character(len=*        ), intent(in   ), optional    :: default
        end subroutine get_char_interface

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