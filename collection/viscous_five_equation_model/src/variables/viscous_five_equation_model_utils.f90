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
    public :: normalize_gradient_volume_fraction
    public :: rotate_gradient_value
    public :: unrotate_gradient_value
    public :: compute_smoothed_volume_fraction
    public :: curvature_preprocessing
    public :: interface_normal_smoothing_weight

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
                rho1 => primitive_variables_set(1, cell_index), &
                rho2 => primitive_variables_set(2, cell_index), &
                z    => primitive_variables_set(7, cell_index)  &
            )
                density(cell_index) = z * rho1 + (1.d0 - z) * rho2
            end associate
        end do

        call a_cellsystem%write_scolar(a_result_writer, "Density"        , density                        )
        call a_cellsystem%write_scolar(a_result_writer, "Density 1"      , primitive_variables_set(1  , :))
        call a_cellsystem%write_scolar(a_result_writer, "Density 2"      , primitive_variables_set(2  , :))
        call a_cellsystem%write_vector(a_result_writer, "Velocity"       , primitive_variables_set(3:5, :))
        call a_cellsystem%write_scolar(a_result_writer, "Pressure"       , primitive_variables_set(6  , :))
        call a_cellsystem%write_scolar(a_result_writer, "Volume fraction", primitive_variables_set(7  , :))
        call a_cellsystem%write_vector(a_result_writer, "Interface normal vector" , surface_tension_variables_set(2:4, :))
        call a_cellsystem%write_scolar(a_result_writer, "Curvature"               , surface_tension_variables_set(  5, :))
        call a_cellsystem%close_file  (a_result_writer)
    end subroutine write_result

    pure function primitive_to_conservative(an_eos, primitive_variables, num_conservatives) result(conservative_variables)
        class  (eos         ), intent(in)  :: an_eos
        real   (real_kind   ), intent(in)  :: primitive_variables   (:)
        integer(int_kind    ), intent(in)  :: num_conservatives
        real   (real_kind   )              :: conservative_variables(num_conservatives)

        real(8) :: rho

        associate(                                 &
                rho1    => primitive_variables(1), &
                rho2    => primitive_variables(2), &
                u       => primitive_variables(3), &
                v       => primitive_variables(4), &
                w       => primitive_variables(5), &
                p       => primitive_variables(6), &
                z1      => primitive_variables(7)  &
            )
            rho = rho1 * z1 + rho2 * (1.d0 - z1)
            conservative_variables(1) = rho1 * z1
            conservative_variables(2) = rho2 * (1.d0 - z1)
            conservative_variables(3) = u  * rho
            conservative_variables(4) = v  * rho
            conservative_variables(5) = w  * rho
            conservative_variables(6) = an_eos%compute_internal_energy_density(p, rho, z1) + 0.5d0 * (u**2.d0 + v**2.d0 + w**2.d0) * rho
            conservative_variables(7) = z1
        end associate
    end function primitive_to_conservative

    pure function conservative_to_primitive(an_eos, conservative_variables, num_primitives) result(primitive_variables)
        class  (eos      ), intent(in)  :: an_eos
        real   (real_kind), intent(in)  :: conservative_variables(:)
        integer(int_kind ), intent(in)  :: num_primitives
        real   (real_kind)              :: primitive_variables   (num_primitives)

        real(real_kind), parameter :: volume_fraction_limmit = 1.0d-6
        real(real_kind) :: rho, ie

        associate(                                    &
                rho1_z1 => conservative_variables(1), &
                rho2_z2 => conservative_variables(2), &
                rho_u   => conservative_variables(3), &
                rho_v   => conservative_variables(4), &
                rho_w   => conservative_variables(5), &
                e       => conservative_variables(6), &
                z1      => conservative_variables(7)  &
            )
            rho = rho1_z1 + rho2_z2
            if(z1 < volume_fraction_limmit)then
                primitive_variables(1) = 0.d0
                primitive_variables(2) = rho2_z2
                primitive_variables(7) = 0.d0
            else if(z1 > 1.d0 - volume_fraction_limmit)then
                primitive_variables(1) = rho1_z1
                primitive_variables(2) = 0.d0
                primitive_variables(7) = 1.d0
            else
                primitive_variables(1) = rho1_z1 / z1
                primitive_variables(2) = rho2_z2 / (1.d0 - z1)
                primitive_variables(7) = z1
            end if

            if(rho < machine_epsilon)then
                primitive_variables(3:6) = 0.d0
            else
                ie  = e / rho - 0.5d0 * (rho_u**2.d0 + rho_v**2.d0 + rho_w**2.d0) / rho**2.d0
                primitive_variables(3) = rho_u / rho
                primitive_variables(4) = rho_v / rho
                primitive_variables(5) = rho_w / rho
                primitive_variables(6) = an_eos%compute_pressure(ie, rho, primitive_variables(7))
            endif
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

    pure function compute_smoothed_volume_fraction(surface_tension_variables, primitive_variables, num_surface_tension_variables) result(dst_surface_tension_variables)
        real   (real_kind ), intent(in) :: surface_tension_variables(:), primitive_variables(:)
        integer(int_kind  ), intent(in) :: num_surface_tension_variables
        real   (real_kind )             :: dst_surface_tension_variables(num_surface_tension_variables)

        real   (real_kind )             :: heavest_volume_fraction
        real   (real_kind ), parameter  :: alpha = 0.1d0

        associate(                                     &
            rho1            => primitive_variables(1), &
            rho2            => primitive_variables(2), &
            volume_fraction => primitive_variables(7)  &
        )
            ! Select a heavest fluid.
            heavest_volume_fraction = 0.5d0 * (1.d0 + sign(1.d0,  rho1 - rho2)) * volume_fraction + 0.5d0 * (1.d0 + sign(1.d0,  rho2 - rho1)) * (1.d0 - volume_fraction)
            ! See [Garrick 2017, JCP]
            !dst_surface_tension_variables(1) = (heavest_volume_fraction**alpha) &
            !                                 / ((heavest_volume_fraction**alpha) + (1.d0 - heavest_volume_fraction)**alpha)
            dst_surface_tension_variables(1) = heavest_volume_fraction
        end associate
    end function compute_smoothed_volume_fraction

    pure function normalize_gradient_volume_fraction(surface_tension_variables, num_surface_tension_variables) result(dst_surface_tension_variables)
        real   (real_kind ), intent(in) :: surface_tension_variables(:)
        integer(int_kind  ), intent(in) :: num_surface_tension_variables
        real   (real_kind )             :: dst_surface_tension_variables(num_surface_tension_variables)

        associate(                                                          &
            mag        => vector_magnitude(surface_tension_variables(2:4)), &
            alpha      => surface_tension_variables(1)                    , &
            grad_alpha => surface_tension_variables(2:4)                    &
        )
            dst_surface_tension_variables(1) = surface_tension_variables(1)
            if((machine_epsilon < mag) .and. (machine_epsilon < alpha) .and. (alpha < 1.d0 - machine_epsilon))then
                ! Towerd a heavest fluid direction.
                dst_surface_tension_variables(2:4) = grad_alpha / mag
            else
                dst_surface_tension_variables(2:4) = 0.d0
            endif
        end associate
    end function normalize_gradient_volume_fraction

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

    pure function curvature_preprocessing(primitives, surface_tension_variables, num_variables) result(dst_primitives)
        real   (real_kind ), intent(in) :: surface_tension_variables(:), primitives(:)
        integer(int_kind  ), intent(in) :: num_variables
        real   (real_kind )             :: dst_primitives(num_variables)

        real(real_kind), parameter :: interface_threshold = 1e-2
        real(real_kind), parameter :: carvature_limit = 2.d0

        dst_primitives(1:7) = primitives(1:7)
        associate(                                 &
            rho1  => primitives(1)               , &
            rho2  => primitives(2)               , &
            z     => primitives(7)               , &
            kappa => surface_tension_variables(5)  &
        )
            if((interface_threshold < z) .and. (z < 1.d0 - interface_threshold))then
                dst_primitives(8) = -kappa!max(min(-kappa, carvature_limit), -carvature_limit)
            else
                dst_primitives(8) = 0.d0
            endif
        end associate
    end function curvature_preprocessing

    pure function interface_normal_smoothing_weight(own_cell_position, neighbor_cell_position, own_variables, neighbor_variables) result(weight)
        real(real_kind ), intent(in) :: own_cell_position     (3)
        real(real_kind ), intent(in) :: neighbor_cell_position(3)
        real(real_kind ), intent(in) :: own_variables         (:)
        real(real_kind ), intent(in) :: neighbor_variables    (:)
        real(real_kind )             :: weight

        real(real_kind), parameter :: smoothing_power = 2.d0 ! range is {@code smoothing_power} > 0.

        associate(alpha => neighbor_variables(1)) !{@code neighbor_variables} is {@code surface_tension_variables}
            weight = (alpha * (1.d0 - alpha))**smoothing_power
        end associate
    end function interface_normal_smoothing_weight
end module viscous_five_equation_model_utils_module