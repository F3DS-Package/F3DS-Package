module class_json_configuration
    use, intrinsic :: iso_fortran_env
    use json_module
    use typedef_module
    use stdio_module
    use abstract_configuration

    implicit none

    private

    type, public, extends(configuration) :: json_configuration
        private

        type     (json_file)              :: json
        character(len=:    ), allocatable :: error_msg
        character(len=:    ), allocatable :: string_val
        logical                           :: found, status_ok
        logical                           :: parsed = .false.

        contains

        procedure, public, pass(self) :: parse
        procedure, public, pass(self) :: close
        procedure, public, pass(self) :: get_real
        procedure, public, pass(self) :: get_int
        procedure, public, pass(self) :: get_bool
        procedure, public, pass(self) :: get_char
        !procedure, public, pass(self) :: add_real
        !procedure, public, pass(self) :: add_int
        !procedure, public, pass(self) :: add_bool
        !procedure, public, pass(self) :: add_char
    end type

    contains

    subroutine parse(self, filepath)
        class    (json_configuration), intent(inout) :: self
        character(len=*             ), intent(in   ) :: filepath

        if(self%parsed)then
            call call_error("'parse' method of json_configuration is already called.")
        end if

        call self%json%initialize()
        call self%json%load(filename=filepath)
        if (self%json%failed()) then
            call self%json%check_for_errors(self%status_ok, self%error_msg)
            call self%json%clear_exceptions()
            call self%json%destroy()
#ifdef _DEBUG
            print *, "DEBUG: JSONFortran error message:"
            print *, self%error_msg
            print *, "---------------------------------"
#endif
            call call_error("Can not read configration file '"//filepath//"'.")
        end if

        self%parsed = .true.
    end subroutine parse

    subroutine close(self)
        class(json_configuration), intent(inout) :: self

        call self%json%destroy()

        self%parsed = .false.
    end subroutine close

    subroutine get_real(self, tag, val, found, default)
        class    (json_configuration), intent(inout)           :: self
        character(len=*             ), intent(in   )           :: tag
        real     (real_kind         ), intent(inout)           :: val
        logical                      , intent(inout), optional :: found
        real     (real_kind         ), intent(in   ), optional :: default

        if(.not. self%parsed) call call_error("Configuration file is not parsed.")

        call self%json%get(tag, val, found)
        if(.not. found) val = default
    end subroutine get_real

    subroutine get_int(self, tag, val, found, default)
        class    (json_configuration), intent(inout)           :: self
        character(len=*             ), intent(in   )           :: tag
        integer  (int_kind          ), intent(inout)           :: val
        logical                      , intent(inout), optional :: found
        integer  (int_kind          ), intent(in   ), optional :: default

        if(.not. self%parsed) call call_error("Configuration file is not parsed.")

        call self%json%get(tag, val, found)
        if(.not. found) val = default
    end subroutine get_int

    subroutine get_bool(self, tag, val, found, default)
        class    (json_configuration), intent(inout)           :: self
        character(len=*             ), intent(in   )           :: tag
        logical                      , intent(inout)           :: val
        logical                      , intent(inout), optional :: found
        logical                      , intent(in   ), optional :: default

        if(.not. self%parsed) call call_error("Configuration file is not parsed.")

        call self%json%get(tag, val, found)
        if(.not. found) val = default
    end subroutine get_bool

    subroutine get_char(self, tag, val, found, default)
        class    (json_configuration), intent(inout)              :: self
        character(len=*             ), intent(in   )              :: tag
        character(len=:             ), intent(inout), allocatable :: val
        logical                      , intent(inout), optional    :: found
        character(len=*             ), intent(in   ), optional    :: default

        if(.not. self%parsed) call call_error("Configuration file is not parsed.")

        call self%json%get(tag, val, found)
        if(.not. found) val = default
    end subroutine get_char

    !subroutine add_real(self, tag, val)
    !    class    (json_configuration), intent(inout) :: self
    !    character(len=*             ), intent(in   ) :: tag
    !    real     (real_kind         ), intent(inout) :: val
    !end subroutine add_real
    !subroutine add_int(self, tag, val)
    !    class    (json_configuration), intent(inout) :: self
    !    character(len=*             ), intent(in   ) :: tag
    !    integer  (int_kind          ), intent(inout) :: val
    !end subroutine add_int
    !subroutine add_bool(self, tag, val)
    !    class    (json_configuration), intent(inout) :: self
    !    character(len=*             ), intent(in   ) :: tag
    !    logical                      , intent(inout) :: val
    !end subroutine add_bool
    !subroutine add_char(self, tag, val)
    !    class    (json_configuration), intent(inout) :: self
    !    character(len=*             ), intent(in   ) :: tag
    !    character(len=*             ), intent(inout) :: val
    !end subroutine add_char
end module class_json_configuration