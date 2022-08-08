module stdio_module
    implicit none

    public

    contains

    subroutine call_error(message)
        character(len=*) :: message
        print *, "\033[31m Error\033[m: ", message
        stop
    end subroutine call_error

    subroutine write_message(message)
        character(len=*) :: message
        print *, "\033[32m Message\033[m: ", message
    end subroutine write_message

    subroutine write_warring(message)
        character(len=*) :: message
        print *, "\033[33m Warning\033[m: ", message
    end subroutine write_warring
end module stdio_module