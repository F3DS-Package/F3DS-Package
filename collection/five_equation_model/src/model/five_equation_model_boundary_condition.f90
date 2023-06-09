module five_equation_model_boundary_condition_module
    use typedef_module
    use vector_module

    implicit none

    private

    public :: outflow_bc, nonslip_wall_bc, slip_and_symmetric_bc, empty_bc

    contains

    pure function outflow_bc(inner_primitive_variables, num_primitive_variables) result(ghost_primitive_variables)
        real   (real_kind), intent(in) :: inner_primitive_variables(:)
        integer(int_kind ), intent(in) :: num_primitive_variables
        real   (real_kind)             :: ghost_primitive_variables(num_primitive_variables)
        ghost_primitive_variables(:) = inner_primitive_variables(:)
    end function outflow_bc

    pure function nonslip_wall_bc(inner_primitive_variables, num_primitive_variables) result(ghost_primitive_variables)
        real   (real_kind), intent(in) :: inner_primitive_variables(:)
        integer(int_kind ), intent(in) :: num_primitive_variables
        real   (real_kind)             :: ghost_primitive_variables(num_primitive_variables)
        ghost_primitive_variables(1:2) =         inner_primitive_variables(1:2)
        ghost_primitive_variables(3  ) = -1.d0 * inner_primitive_variables(3)  ! slipwall
        ghost_primitive_variables(4:7) =         inner_primitive_variables(4:7)
    end function nonslip_wall_bc

    pure function slip_and_symmetric_bc(inner_primitive_variables, num_primitive_variables) result(ghost_primitive_variables)
        real   (real_kind), intent(in) :: inner_primitive_variables(:)
        integer(int_kind ), intent(in) :: num_primitive_variables
        real   (real_kind)             :: ghost_primitive_variables(num_primitive_variables)
        ghost_primitive_variables(1:2) =         inner_primitive_variables(1:2)
        ghost_primitive_variables(3  ) = -1.d0 * inner_primitive_variables(3)
        ghost_primitive_variables(4:7) =         inner_primitive_variables(4:7)
    end function slip_and_symmetric_bc

    pure function empty_bc(inner_primitive_variables, num_primitive_variables) result(ghost_primitive_variables)
        real   (real_kind), intent(in) :: inner_primitive_variables(:)
        integer(int_kind ), intent(in) :: num_primitive_variables
        real   (real_kind)             :: ghost_primitive_variables(num_primitive_variables)
        ghost_primitive_variables(1:2) = inner_primitive_variables(1:2)
        ghost_primitive_variables(3  ) = inner_primitive_variables(3)
        ghost_primitive_variables(4:7) = inner_primitive_variables(4:7)
    end function empty_bc
end module five_equation_model_boundary_condition_module