module five_equation_model_boundary_condition_module
    use typedef_module
    use vector_module

    implicit none

    private

    public :: outflow_bc, slipwall_bc, symmetric_bc, empty_bc, gradient_volume_fraction_bc

    contains

    pure function outflow_bc(inner_primitive_variables, num_primitive_variables) result(ghost_primitive_variables)
        real   (real_kind), intent(in) :: inner_primitive_variables(:)
        integer(int_kind ), intent(in) :: num_primitive_variables
        real   (real_kind)             :: ghost_primitive_variables(num_primitive_variables)
        ghost_primitive_variables(:) = inner_primitive_variables(:)
    end function outflow_bc

    pure function slipwall_bc(inner_primitive_variables, num_primitive_variables) result(ghost_primitive_variables)
        real   (real_kind), intent(in) :: inner_primitive_variables(:)
        integer(int_kind ), intent(in) :: num_primitive_variables
        real   (real_kind)             :: ghost_primitive_variables(num_primitive_variables)
        ghost_primitive_variables(1:2) =         inner_primitive_variables(1:2)
        ghost_primitive_variables(3  ) = -1.d0 * inner_primitive_variables(3)
        ghost_primitive_variables(4:8) =         inner_primitive_variables(4:8)
    end function slipwall_bc

    pure function symmetric_bc(inner_primitive_variables, num_primitive_variables) result(ghost_primitive_variables)
        real   (real_kind), intent(in) :: inner_primitive_variables(:)
        integer(int_kind ), intent(in) :: num_primitive_variables
        real   (real_kind)             :: ghost_primitive_variables(num_primitive_variables)
        ghost_primitive_variables(1:2) =         inner_primitive_variables(1:2)
        ghost_primitive_variables(3  ) = -1.d0 * inner_primitive_variables(3)
        ghost_primitive_variables(4:8) =         inner_primitive_variables(4:8)
    end function symmetric_bc

    pure function empty_bc(inner_primitive_variables, num_primitive_variables) result(ghost_primitive_variables)
        real   (real_kind), intent(in) :: inner_primitive_variables(:)
        integer(int_kind ), intent(in) :: num_primitive_variables
        real   (real_kind)             :: ghost_primitive_variables(num_primitive_variables)
        ghost_primitive_variables(1:2) =         inner_primitive_variables(1:2)
        ghost_primitive_variables(3  ) = -1.d0 * inner_primitive_variables(3)
        ghost_primitive_variables(4:8) =         inner_primitive_variables(4:8)
    end function empty_bc

    pure function gradient_volume_fraction_bc(inner_gradient_values, num_gradient_values) result(ghost_gradient_values)
        real   (real_kind), intent(in) :: inner_gradient_values(:)
        integer(int_kind ), intent(in) :: num_gradient_values
        real   (real_kind)             :: ghost_gradient_values(num_gradient_values)
        ghost_gradient_values(1) = inner_gradient_values(1)
        ghost_gradient_values(2) = inner_gradient_values(2)
        ghost_gradient_values(3) = inner_gradient_values(3)
    end function gradient_volume_fraction_bc
end module five_equation_model_boundary_condition_module