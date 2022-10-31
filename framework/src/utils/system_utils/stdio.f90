module stdio_module
    implicit none

    private

    character(len=5) :: color_green = char(27)//'[32m'
    character(len=5) :: color_yelow = char(27)//'[33m'
    character(len=5) :: color_red   = char(27)//'[31m'
    character(len=5) :: color_cyan  = char(27)//'[36m'
    character(len=4) :: esc_reset   = char(27)//'[0m'

    public :: call_error
    public :: write_message
    public :: write_warring
    public :: write_debuginfo

    contains

    subroutine call_error(message)
        character(len=*) :: message
        print *, color_red//"Error: "//esc_reset, message
        print *, "F3DS abend..."
        error stop
    end subroutine call_error

    subroutine write_message(message)
        character(len=*) :: message
        print *, color_green//"Message: "//esc_reset, message
    end subroutine write_message

    subroutine write_warring(message)
        character(len=*) :: message
        print *, color_yelow//"Warning: "//esc_reset, message
    end subroutine write_warring

    subroutine write_debuginfo(message)
        character(len=*) :: message
#ifdef _DEBUG
        print *, color_cyan//"Debug: "//esc_reset, message
#endif
    end subroutine write_debuginfo
end module stdio_module