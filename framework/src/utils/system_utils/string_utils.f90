module string_utils_module
    use typedef_module

    implicit none

    private

    public :: to_str => int_to_str

    contains

    pure function int_to_str(num, extra_digits) result(str)
        integer  (int_kind), intent(in)           :: num
        integer  (int_kind), intent(in), optional :: extra_digits
        character(:       ), allocatable          :: str

        integer  (int_kind)              :: num_digits
        integer  (int_kind)              :: dgt_digits
        character(:       ), allocatable :: str_digits
        character(:       ), allocatable :: fmt

        if(.not. present(extra_digits)) extra_digits = 0

        num_digits = get_int_digits_of(num)
        if (num_digits < extra_digits) num_digits = extra_digits
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