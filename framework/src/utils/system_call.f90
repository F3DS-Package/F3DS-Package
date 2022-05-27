module sytem_call_module
    implicit none

    private

    public :: make_dir

    contains

    subroutine make_dir(dir_name)
        character(len=*), intent(in) :: dir_name
        call system("if [ ! -d "//dir_name//" ]; then mkdir -p "//dir_name//"; fi")
    end subroutine
end module sytem_call_module