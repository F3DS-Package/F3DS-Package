module typedef_module
    use, intrinsic :: iso_fortran_env

    implicit none

    public

    integer, parameter :: real_kind  = real64
    integer, parameter :: int_kind   = int32
    integer, parameter :: type_kind  = int8
end module typedef_module