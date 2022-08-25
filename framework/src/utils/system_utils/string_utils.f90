module string_utils_module
    use typedef_module

    implicit none

    private

    public :: to_str
    interface to_str
        module procedure int_to_str ! add string conversion functions such as dble_to_str, ... in a future.
    end interface to_str

    contains

    pure function int_to_str(num, extra_digits) result(str)
        integer  (int_kind), intent(in)           :: num
        integer  (int_kind), intent(in), optional :: extra_digits
        character(:       ), allocatable          :: str

        integer  (int_kind)              :: fixed_num_digits
        integer  (int_kind)              :: num_digits
        integer  (int_kind)              :: dgt_digits
        character(:       ), allocatable :: str_digits
        character(:       ), allocatable :: fmt

        if(.not. present(extra_digits))then
            fixed_num_digits = 0
        else
            fixed_num_digits = extra_digits
        end if

        num_digits = get_int_digits_of(num)
        if (num_digits < fixed_num_digits) num_digits = fixed_num_digits
        dgt_digits = get_int_digits_of(num_digits)

        allocate(character(dgt_digits)::str_digits)
        write(str_digits,'(I0)') num_digits

        fmt = "(I"//str_digits//"."//str_digits//")"

        allocate(character(num_digits)::str)
        write(str,fmt) num
    end function

    pure function get_int_digits_of(num) result(num_digit)
        integer(int_kind), intent(in) :: num
        integer(int_kind)             :: num_digit

        num_digit = int(log10(dble(num)))+1
    end function get_int_digits_of

end module string_utils_module