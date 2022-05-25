module math_constant_module
    use typedef_module

    implicit none

    real(real_kind), parameter :: pi = 3.14159265358979323846264338327950288d0
    real(real_kind), parameter :: machine_epsilon = 1.d-8
end module math_constant_module