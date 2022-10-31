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
        procedure, public, pass(self) :: get_real_array
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
            call write_debuginfo("Error occurred while parsing configration file by JSON Fortran. Error message is following:")
            print *, self%error_msg
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

        if (self%json%failed()) then
            call self%json%check_for_errors(self%status_ok, self%error_msg)
            call self%json%clear_exceptions()
            call self%json%destroy()
#ifdef _DEBUG
            call write_debuginfo("Error occurred while reading configration file by JSON Fortran. Error message is following:")
            print *, self%error_msg
#endif
        end if

        if(.not. found .and. present(default)) val = default
    end subroutine get_real

    subroutine get_int(self, tag, val, found, default)
        class    (json_configuration), intent(inout)           :: self
        character(len=*             ), intent(in   )           :: tag
        integer  (int_kind          ), intent(inout)           :: val
        logical                      , intent(inout), optional :: found
        integer  (int_kind          ), intent(in   ), optional :: default

        if(.not. self%parsed) call call_error("Configuration file is not parsed.")

        call self%json%get(tag, val, found)

        if (self%json%failed()) then
            call self%json%check_for_errors(self%status_ok, self%error_msg)
            call self%json%clear_exceptions()
            call self%json%destroy()
#ifdef _DEBUG
            call write_debuginfo("Error occurred while reading configration file by JSON Fortran. Error message is following:")
            print *, self%error_msg
#endif
        end if

        if(.not. found .and. present(default)) val = default
    end subroutine get_int

    subroutine get_bool(self, tag, val, found, default)
        class    (json_configuration), intent(inout)           :: self
        character(len=*             ), intent(in   )           :: tag
        logical                      , intent(inout)           :: val
        logical                      , intent(inout), optional :: found
        logical                      , intent(in   ), optional :: default

        if(.not. self%parsed) call call_error("Configuration file is not parsed.")

        call self%json%get(tag, val, found)

        if (self%json%failed()) then
            call self%json%check_for_errors(self%status_ok, self%error_msg)
            call self%json%clear_exceptions()
            call self%json%destroy()
#ifdef _DEBUG
            call write_debuginfo("Error occurred while reading configration file by JSON Fortran. Error message is following:")
            print *, self%error_msg
#endif
        end if

        if(.not. found .and. present(default)) val = default
    end subroutine get_bool

    subroutine get_char(self, tag, val, found, default)
        class    (json_configuration), intent(inout)              :: self
        character(len=*             ), intent(in   )              :: tag
        character(len=:             ), intent(inout), allocatable :: val
        logical                      , intent(inout), optional    :: found
        character(len=*             ), intent(in   ), optional    :: default

        if(.not. self%parsed) call call_error("Configuration file is not parsed.")

        call self%json%get(tag, val, found)

        if (self%json%failed()) then
            call self%json%check_for_errors(self%status_ok, self%error_msg)
            call self%json%clear_exceptions()
            call self%json%destroy()
#ifdef _DEBUG
            call write_debuginfo("Error occurred while reading configration file by JSON Fortran. Error message is following:")
            print *, self%error_msg
#endif
        end if

        if(.not. found .and. present(default)) val = default
    end subroutine get_char

    subroutine get_real_array(self, tag, val, found, default)
        class    (json_configuration), intent(inout)           :: self
        character(len=*             ), intent(in   )           :: tag
        real     (real_kind         ), intent(inout)           :: val(:)
        logical                      , intent(inout), optional :: found
        real     (real_kind         ), intent(in   ), optional :: default(:)

        real(json_RK), dimension(:), allocatable :: array_json_file
        logical                    :: match_array_range

        call call_error("'get_real_array' defined in 'json_configuration' can not be use in this version.")
        ! TODO: Fix segmentation fault

        if(.not. self%parsed) call call_error("Configuration file is not parsed.")

        if(present(default))then
            match_array_range = ((lbound(val(:), 1) == lbound(default(:), 1)) .and. (ubound(val(:), 1) == ubound(default(:), 1)))
            if(.not. match_array_range) call call_error("Array range of 'val' does not match 'default'.")
        end if

        call self%json%get(tag, array_json_file)

        if (self%json%failed()) then
            call self%json%check_for_errors(self%status_ok, self%error_msg)
            call self%json%clear_exceptions()
            call self%json%destroy()
#ifdef _DEBUG
            call write_debuginfo("Error occurred while reading configration file by JSON Fortran. Error message is following:")
            print *, self%error_msg
#endif
        end if

        if(.not. found .and. present(default))then
            val = default
            return
        end if

        match_array_range = ((lbound(val(:), 1) == lbound(array_json_file(:), 1)) .and. (ubound(val(:), 1) == ubound(array_json_file(:), 1)))
        if(.not. match_array_range) call call_error("Array range of 'val' dose not match array range of json path you set.")

        val = array_json_file
    end subroutine get_real_array

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