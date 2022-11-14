module viscous_five_equation_model_boundary_condition_module
    use typedef_module
    use vector_module
    use math_constant_module

    implicit none

    private

    public :: outflow_bc, wall_bc, symmetric_bc, empty_bc
    public :: bc_for_normarized_gradient_volume_fraction, bc_for_normarized_gradient_volume_fraction_wall

    contains

    pure function outflow_bc(inner_primitive_variables, num_primitive_variables) result(ghost_primitive_variables)
        real   (real_kind), intent(in) :: inner_primitive_variables(:)
        integer(int_kind ), intent(in) :: num_primitive_variables
        real   (real_kind)             :: ghost_primitive_variables(num_primitive_variables)
        ghost_primitive_variables(:) = inner_primitive_variables(:)
    end function outflow_bc

    pure function wall_bc(inner_primitive_variables, num_primitive_variables) result(ghost_primitive_variables)
        real   (real_kind), intent(in) :: inner_primitive_variables(:)
        integer(int_kind ), intent(in) :: num_primitive_variables
        real   (real_kind)             :: ghost_primitive_variables(num_primitive_variables)
        ghost_primitive_variables(1:2) =         inner_primitive_variables(1:2)
        ghost_primitive_variables(3:5) = -1.d0 * inner_primitive_variables(3:5) ! nonslip-wall condition
        ghost_primitive_variables(6:8) =         inner_primitive_variables(6:8)
    end function wall_bc

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

    pure function bc_for_normarized_gradient_volume_fraction(inner_gradient_values, num_gradient_values) result(ghost_gradient_values)
        real   (real_kind), intent(in) :: inner_gradient_values(:)
        integer(int_kind ), intent(in) :: num_gradient_values
        real   (real_kind)             :: ghost_gradient_values(num_gradient_values)
        ghost_gradient_values(1  ) = -1.d0 * inner_gradient_values(1  )
        ghost_gradient_values(2:3) =         inner_gradient_values(2:3)
    end function bc_for_normarized_gradient_volume_fraction

    pure function bc_for_normarized_gradient_volume_fraction_wall(inner_gradient_values, num_gradient_values) result(ghost_gradient_values)
        real   (real_kind), intent(in) :: inner_gradient_values(:)
        integer(int_kind ), intent(in) :: num_gradient_values
        real   (real_kind)             :: ghost_gradient_values(num_gradient_values)
        real   (real_kind), parameter  :: contact_angle = (1.d0 / 6.d0) * pi ! 30 deg
        real   (real_kind), parameter  :: wall_normal  (1:3) = [real(real_kind) :: -1.d0, 0.0d0, 0.0d0]
        real   (real_kind)             :: wall_tangetal(1:3)

        wall_tangetal(1:3) = inner_gradient_values - (vector_multiply(wall_normal, inner_gradient_values) / vector_multiply(wall_normal, wall_normal)) * wall_normal
        ghost_gradient_values(1:3) = wall_normal(1:3) * cos(contact_angle) + wall_tangetal(1:3) * sin(contact_angle)
    end function bc_for_normarized_gradient_volume_fraction_wall
end module viscous_five_equation_model_boundary_condition_module