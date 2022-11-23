module viscous_five_equation_model_utils_module
    use typedef_module
    use stdio_module
    use abstract_eos
    use math_constant_module
    use vector_module
    use class_cellsystem
    use abstract_result_writer

    implicit none

    private

    public :: write_result
    public :: conservative_to_primitive
    public :: primitive_to_conservative
    public :: rotate_primitive
    public :: unrotate_primitive
    public :: rotate_gradient_value
    public :: unrotate_gradient_value

    contains

    subroutine write_result(a_cellsystem, a_result_writer, primitive_variables_set, surface_tension_variables_set)
        type (cellsystem   ), intent(inout) :: a_cellsystem
        class(result_writer), intent(inout) :: a_result_writer
        real (real_kind    ), intent(in   ) :: primitive_variables_set      (:,:)
        real (real_kind    ), intent(in   ) :: surface_tension_variables_set(:,:)

        real   (real_kind) :: density(a_cellsystem%get_number_of_cells())
        integer(int_kind ) :: cell_index

        call a_cellsystem%open_file   (a_result_writer)

        call write_message("Write a result "//a_cellsystem%get_filename(a_result_writer)//"...")

        do cell_index = 1, a_cellsystem%get_number_of_cells(), 1
            associate(                                          &
                rho1a1 => primitive_variables_set(1, cell_index), &
                rho2a2 => primitive_variables_set(2, cell_index), &
                z      => primitive_variables_set(7, cell_index)  &
            )
                density(cell_index) = rho1a1 + rho2a2
            end associate
        end do

        call a_cellsystem%write_scolar(a_result_writer, "Density"           , density                        )
        call a_cellsystem%write_scolar(a_result_writer, "Smoothed Density 1", primitive_variables_set(1  , :))
        call a_cellsystem%write_scolar(a_result_writer, "Smoothed Density 2", primitive_variables_set(2  , :))
        call a_cellsystem%write_vector(a_result_writer, "Velocity"          , primitive_variables_set(3:5, :))
        call a_cellsystem%write_scolar(a_result_writer, "Pressure"          , primitive_variables_set(6  , :))
        call a_cellsystem%write_scolar(a_result_writer, "Volume fraction"   , primitive_variables_set(7  , :))
        call a_cellsystem%write_vector(a_result_writer, "Interface normal vector" , surface_tension_variables_set(2:4, :))
        call a_cellsystem%write_scolar(a_result_writer, "Negative curvature"      , surface_tension_variables_set(  5, :))
        call a_cellsystem%close_file  (a_result_writer)
    end subroutine write_result

    pure function primitive_to_conservative(an_eos, primitive_variables, num_conservatives) result(conservative_variables)
        class  (eos         ), intent(in)  :: an_eos
        real   (real_kind   ), intent(in)  :: primitive_variables   (:)
        integer(int_kind    ), intent(in)  :: num_conservatives
        real   (real_kind   )              :: conservative_variables(num_conservatives)

        real(8) :: rho

        associate(                                 &
                rho1a1  => primitive_variables(1), &
                rho2a2  => primitive_variables(2), &
                u       => primitive_variables(3), &
                v       => primitive_variables(4), &
                w       => primitive_variables(5), &
                p       => primitive_variables(6), &
                alpha   => primitive_variables(7)  &
            )
            rho = rho1a1 + rho2a2

            conservative_variables(1) = rho1a1
            conservative_variables(2) = rho2a2
            conservative_variables(3) = u  * rho
            conservative_variables(4) = v  * rho
            conservative_variables(5) = w  * rho
            conservative_variables(6) = an_eos%compute_internal_energy_density(p, rho, alpha) + 0.5d0 * (u**2.d0 + v**2.d0 + w**2.d0) * rho
            conservative_variables(7) = alpha
        end associate
    end function primitive_to_conservative

    pure function conservative_to_primitive(an_eos, conservative_variables, num_primitives) result(primitive_variables)
        class  (eos      ), intent(in)  :: an_eos
        real   (real_kind), intent(in)  :: conservative_variables(:)
        integer(int_kind ), intent(in)  :: num_primitives
        real   (real_kind)              :: primitive_variables   (num_primitives)

        real(real_kind) :: rho, ie, mod_rho1a1, mod_rho2a2, mod_alpha

        associate(                                    &
                rho1a1  => conservative_variables(1), &
                rho2a2  => conservative_variables(2), &
                rho_u   => conservative_variables(3), &
                rho_v   => conservative_variables(4), &
                rho_w   => conservative_variables(5), &
                e       => conservative_variables(6), &
                alpha   => conservative_variables(7)  &
            )
            if(alpha < machine_epsilon)then
                mod_alpha  = 0.d0
                mod_rho1a1 = 0.d0
                mod_rho2a2 = rho2a2
            else if(alpha > 1.d0 - machine_epsilon)then
                mod_alpha  = 1.d0
                mod_rho1a1 = rho1a1
                mod_rho2a2 = 0.d0
            else
                mod_alpha  = alpha
                mod_rho1a1 = rho1a1
                mod_rho2a2 = rho2a2
            end if

            rho = mod_rho1a1 + mod_rho2a2
            ie  = e / rho - 0.5d0 * (rho_u**2.d0 + rho_v**2.d0 + rho_w**2.d0) / rho**2.d0

            primitive_variables(1) = mod_rho1a1
            primitive_variables(2) = mod_rho2a2
            if(rho < machine_epsilon)then
                primitive_variables(3:6) = 0.d0
            else
                primitive_variables(3) = rho_u / rho
                primitive_variables(4) = rho_v / rho
                primitive_variables(5) = rho_w / rho
                primitive_variables(6) = an_eos%compute_pressure(ie, rho, mod_alpha)
            end if
            primitive_variables(7) = mod_alpha
        end associate
    end function conservative_to_primitive

    pure function rotate_primitive( &
        primitives                , &
        face_normal_vector        , &
        face_tangential1_vector   , &
        face_tangential2_vector   , &
        num_primitives                ) result(p)

        real   (real_kind     ), intent(in)  :: primitives              (:)
        real   (real_kind     ), intent(in)  :: face_normal_vector      (3)
        real   (real_kind     ), intent(in)  :: face_tangential1_vector (3)
        real   (real_kind     ), intent(in)  :: face_tangential2_vector (3)
        integer(int_kind      ), intent(in)  :: num_primitives
        real   (real_kind     )              :: p(num_primitives)

        p(1  ) = primitives(1)
        p(2  ) = primitives(2)
        p(3:5) = vector_rotate(primitives(3:5), face_normal_vector, face_tangential1_vector, face_tangential2_vector)
        p(6  ) = primitives(6)
        p(7  ) = primitives(7)
        p(8  ) = primitives(8)
    end function rotate_primitive

    pure function unrotate_primitive( &
        primitives                  , &
        face_normal_vector          , &
        face_tangential1_vector     , &
        face_tangential2_vector     , &
        num_primitives                  ) result(p)

        real   (real_kind     ), intent(in)  :: primitives              (:)
        real   (real_kind     ), intent(in)  :: face_normal_vector      (3)
        real   (real_kind     ), intent(in)  :: face_tangential1_vector (3)
        real   (real_kind     ), intent(in)  :: face_tangential2_vector (3)
        integer(int_kind      ), intent(in)  :: num_primitives
        real   (real_kind     )              :: p(num_primitives)

        p(1  ) = primitives(1)
        p(2  ) = primitives(2)
        p(3:5) = vector_unrotate(primitives(3:5), face_normal_vector, face_tangential1_vector, face_tangential2_vector)
        p(6  ) = primitives(6)
        p(7  ) = primitives(7)
        p(8  ) = primitives(8)
    end function unrotate_primitive

    pure function rotate_gradient_value( &
        gradient_value                 , &
        face_normal_vector             , &
        face_tangential1_vector        , &
        face_tangential2_vector        , &
        num_values                         ) result(p)

        real   (real_kind     ), intent(in)  :: gradient_value          (:)
        real   (real_kind     ), intent(in)  :: face_normal_vector      (3)
        real   (real_kind     ), intent(in)  :: face_tangential1_vector (3)
        real   (real_kind     ), intent(in)  :: face_tangential2_vector (3)
        integer(int_kind      ), intent(in)  :: num_values
        real   (real_kind     )              :: p(num_values)

        p(1:3) = vector_rotate(gradient_value(1:3), face_normal_vector, face_tangential1_vector, face_tangential2_vector)
    end function rotate_gradient_value

    pure function unrotate_gradient_value( &
        gradient_value              , &
        face_normal_vector          , &
        face_tangential1_vector     , &
        face_tangential2_vector     , &
        num_values                  ) result(p)

        real   (real_kind     ), intent(in)  :: gradient_value          (:)
        real   (real_kind     ), intent(in)  :: face_normal_vector      (3)
        real   (real_kind     ), intent(in)  :: face_tangential1_vector (3)
        real   (real_kind     ), intent(in)  :: face_tangential2_vector (3)
        integer(int_kind      ), intent(in)  :: num_values
        real   (real_kind     )              :: p(num_values)

        p(1:3) = vector_unrotate(gradient_value(1:3), face_normal_vector, face_tangential1_vector, face_tangential2_vector)
    end function unrotate_gradient_value
end module viscous_five_equation_model_utils_module