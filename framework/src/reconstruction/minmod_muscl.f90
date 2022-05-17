module minmod_muscl_module
    !**
    !* 3rd order TVD-MUSCL with minmod limmiter
    !*

    use typedef_module
    use vector_module

    implicit none

    private

    integer(int_kind ) :: num_ghost_cells_ = 3
    real   (real_kind), parameter :: kappa_ = 1.d0 / 3.d0

    public :: reconstruct_minmod_muscl

    contains

    pure function minmod(s) result(r)
        real(real_kind), intent(in) :: s
        real(real_kind)             :: r
        r = max(0.d0, min(1.d0, s))
    end function

    pure function reconstruct_leftside_minmod_muscl( &
        primitive_values_set             , &
        reference_cell_indexs_set        , &
        cell_positions                   , &
        face_positions                   , &
        face_index                          ) result(reconstructed_primitive)

        real   (real_kind), intent(in) :: primitive_values_set      (:, :)
        integer(int_kind ), intent(in) :: reference_cell_indexs_set (:, :)
        real   (real_kind), intent(in) :: cell_positions            (:, :)
        real   (real_kind), intent(in) :: face_positions            (:, :)
        integer(int_kind ), intent(in) :: face_index
        real   (real_kind)             :: reconstructed_primitive   (size(primitive_values_set(1, :)))

        integer(int_kind ) :: n_primitives, i
        real   (real_kind) :: w(3), p(3), cell_pos_l(3), cell_pos_r(3), face_pos(3)
        real   (real_kind) :: cell_cell_distance, cell_face_distanse, s

        n_primitives = size(primitive_values_set(1, :))

        do i = 1, n_primitives, 1
            associate(                                                                                      &
                v_m1 => primitive_values_set(reference_cell_indexs_set(face_index, num_ghost_cells_-1), i), &
                v    => primitive_values_set(reference_cell_indexs_set(face_index, num_ghost_cells_+0), i), &
                v_p1 => primitive_values_set(reference_cell_indexs_set(face_index, num_ghost_cells_+1), i)  &
            )
                s = (v - v_m1) / (v_p1 - v)
                reconstructed_primitive(i) = v &
                                           + 0.25d0 * (1.d0 - kappa_) * minmod(1.d0 / s) * (v    - v_m1) &
                                           + 0.25d0 * (1.d0 - kappa_) * minmod(       s) * (v_p1 - v   )
            end associate
        end do
    end function reconstruct_leftside_minmod_muscl

    pure function reconstruct_rightside_minmod_muscl( &
        primitive_values_set             , &
        reference_cell_indexs_set        , &
        cell_positions                   , &
        face_positions                   , &
        face_index                          ) result(reconstructed_primitive)

        real   (real_kind), intent(in) :: primitive_values_set      (:, :)
        integer(int_kind ), intent(in) :: reference_cell_indexs_set (:, :)
        real   (real_kind), intent(in) :: cell_positions            (:, :)
        real   (real_kind), intent(in) :: face_positions            (:, :)
        integer(int_kind ), intent(in) :: face_index
        real   (real_kind)             :: reconstructed_primitive   (size(primitive_values_set(1, :)))

        integer(int_kind ) :: n_primitives, i
        real   (real_kind) :: w(3), p(3), cell_pos_l(3), cell_pos_r(3), face_pos(3)
        real   (real_kind) :: cell_cell_distance, cell_face_distanse, s

        n_primitives = size(primitive_values_set(1, :))

        do i = 1, n_primitives, 1
            associate(                                                                                            &
                v_m1       => primitive_values_set(reference_cell_indexs_set(face_index, num_ghost_cells_+0), i), &
                v          => primitive_values_set(reference_cell_indexs_set(face_index, num_ghost_cells_+1), i), &
                v_p1       => primitive_values_set(reference_cell_indexs_set(face_index, num_ghost_cells_+2), i)  &
            )
                s = (v - v_m1) / (v_p1 - v)
                reconstructed_primitive(i) = v &
                                           - 0.25d0 * (1.d0 - kappa_) * minmod(1.d0 / s) * (v    - v_m1) &
                                           - 0.25d0 * (1.d0 - kappa_) * minmod(       s) * (v_p1 - v   )
            end associate
        end do
    end function reconstruct_rightside_minmod_muscl

    pure function reconstruct_minmod_muscl(        &
        primitive_values_set              , &
        cell_centor_positions             , &
        cell_volumes                      , &
        reference_cell_indexs_set         , &
        face_centor_positions             , &
        face_normal_vectors               , &
        face_tangential1_vectors          , &
        face_tangential2_vectors          , &
        face_areas                        , &
        face_index                        , &
        n_conservative_values             , &
        n_derivative_values              , &
        flux_function                     , &
        eos_pressure_function             , &
        eos_soundspeed_function           , &
        primitive_to_conservative_function, &
        integrated_element_function             ) result(element_lef_and_right_side)

        real   (real_kind), intent(in ) :: primitive_values_set     (:, :)
        real   (real_kind), intent(in ) :: cell_centor_positions    (:, :)
        real   (real_kind), intent(in ) :: cell_volumes             (:)
        integer(int_kind ), intent(in ) :: reference_cell_indexs_set(:,:)
        real   (real_kind), intent(in ) :: face_normal_vectors      (:,:)
        real   (real_kind), intent(in ) :: face_tangential1_vectors (:,:)
        real   (real_kind), intent(in ) :: face_tangential2_vectors (:,:)
        real   (real_kind), intent(in ) :: face_centor_positions    (:,:)
        real   (real_kind), intent(in ) :: face_areas               (:)
        integer(int_kind ), intent(in ) :: face_index
        integer(int_kind ), intent(in ) :: n_conservative_values
        integer(int_kind ), intent(in ) :: n_derivative_values

        real   (real_kind)              :: element_lef_and_right_side(2, n_conservative_values+n_derivative_values)

        interface
            pure function flux_function(       &
                left_conservative            , &
                left_main_velocity           , &
                left_density                 , &
                left_pressure                , &
                left_soundspeed              , &
                right_conservative           , &
                right_main_velocity          , &
                right_density                , &
                right_pressure               , &
                right_soundspeed             ) result(flux)

                use typedef_module
                real(real_kind), intent(in) :: left_conservative   (:)
                real(real_kind), intent(in) :: left_main_velocity
                real(real_kind), intent(in) :: left_density
                real(real_kind), intent(in) :: left_pressure
                real(real_kind), intent(in) :: left_soundspeed
                real(real_kind), intent(in) :: right_conservative   (:)
                real(real_kind), intent(in) :: right_main_velocity
                real(real_kind), intent(in) :: right_density
                real(real_kind), intent(in) :: right_pressure
                real(real_kind), intent(in) :: right_soundspeed
                real(real_kind)             :: flux(size(left_conservative))
            end function flux_function

            pure function eos_pressure_function(specific_internal_energy, density, volume_fruction) result(pressure)
                use typedef_module
                real(real_kind), intent(in) :: specific_internal_energy
                real(real_kind), intent(in) :: density
                real(real_kind), intent(in) :: volume_fruction
                real(real_kind)             :: pressure
            end function eos_pressure_function

            pure function eos_soundspeed_function(specific_internal_energy, density, volume_fruction) result(soundspeed)
                use typedef_module
                real(real_kind), intent(in) :: specific_internal_energy
                real(real_kind), intent(in) :: density
                real(real_kind), intent(in) :: volume_fruction
                real(real_kind)             :: soundspeed
            end function eos_soundspeed_function

            pure function primitive_to_conservative_function(primitive) result(conservative)
                use typedef_module
                real(real_kind), intent(in)  :: primitive   (:)
                real(real_kind), allocatable :: conservative(:)
            end function primitive_to_conservative_function

            pure function integrated_element_function( &
                reconstructed_leftside_primitive  , &
                reconstructed_rightside_primitive , &
                leftside_cell_volume              , &
                rightside_cell_volume             , &
                face_normal_vector                , &
                face_tangential1_vector           , &
                face_tangential2_vector           , &
                face_area                         , &
                n_conservative_values             , &
                n_derivative_values              , &
                flux_function                     , &
                eos_pressure_function             , &
                eos_soundspeed_function           , &
                primitive_to_conservative_function   ) result(element_lef_and_right_side)

                use typedef_module

                real   (real_kind), intent(in ) :: reconstructed_leftside_primitive  (:)
                real   (real_kind), intent(in ) :: reconstructed_rightside_primitive (:)
                real   (real_kind), intent(in ) :: leftside_cell_volume
                real   (real_kind), intent(in ) :: rightside_cell_volume
                real   (real_kind), intent(in ) :: face_normal_vector                (3)
                real   (real_kind), intent(in ) :: face_tangential1_vector           (3)
                real   (real_kind), intent(in ) :: face_tangential2_vector           (3)
                real   (real_kind), intent(in ) :: face_area
                integer(int_kind ), intent(in ) :: n_conservative_values
                integer(int_kind ), intent(in ) :: n_derivative_values
                real   (real_kind)              :: element_lef_and_right_side        (2, n_conservative_values+n_derivative_values)

                interface
                    pure function flux_function(       &
                        left_conservative            , &
                        left_main_velocity           , &
                        left_density                 , &
                        left_pressure                , &
                        left_soundspeed              , &
                        right_conservative           , &
                        right_main_velocity          , &
                        right_density                , &
                        right_pressure               , &
                        right_soundspeed             ) result(flux)

                        use typedef_module
                        real(real_kind), intent(in) :: left_conservative   (:)
                        real(real_kind), intent(in) :: left_main_velocity
                        real(real_kind), intent(in) :: left_density
                        real(real_kind), intent(in) :: left_pressure
                        real(real_kind), intent(in) :: left_soundspeed
                        real(real_kind), intent(in) :: right_conservative   (:)
                        real(real_kind), intent(in) :: right_main_velocity
                        real(real_kind), intent(in) :: right_density
                        real(real_kind), intent(in) :: right_pressure
                        real(real_kind), intent(in) :: right_soundspeed
                        real(real_kind)             :: flux(size(left_conservative))
                    end function flux_function

                    pure function eos_pressure_function(specific_internal_energy, density, volume_fruction) result(pressure)
                        use typedef_module
                        real(real_kind), intent(in) :: specific_internal_energy
                        real(real_kind), intent(in) :: density
                        real(real_kind), intent(in) :: volume_fruction
                        real(real_kind)             :: pressure
                    end function eos_pressure_function

                    pure function eos_soundspeed_function(specific_internal_energy, density, volume_fruction) result(soundspeed)
                        use typedef_module
                        real(real_kind), intent(in) :: specific_internal_energy
                        real(real_kind), intent(in) :: density
                        real(real_kind), intent(in) :: volume_fruction
                        real(real_kind)             :: soundspeed
                    end function eos_soundspeed_function

                    pure function primitive_to_conservative_function(primitive) result(conservative)
                        use typedef_module
                        real(real_kind), intent(in)  :: primitive   (:)
                        real(real_kind), allocatable :: conservative(:)
                    end function primitive_to_conservative_function
                end interface
            end function integrated_element_function
        end interface

        integer(int_kind )              :: lhc_index, rhc_index
        real   (real_kind), allocatable :: lhc_primitive   (:)
        real   (real_kind), allocatable :: rhc_primitive   (:)
        integer(int_kind )              :: n_primitives, i
        real   (real_kind)              :: delta

        lhc_index = reference_cell_indexs_set(face_index, num_ghost_cells_+0)
        rhc_index = reference_cell_indexs_set(face_index, num_ghost_cells_+1)

        lhc_primitive = reconstruct_leftside_minmod_muscl(      &
            primitive_values_set                       , &
            reference_cell_indexs_set                  , &
            cell_centor_positions                      , &
            face_centor_positions                      , &
            face_index                                   &
        )
        rhc_primitive = reconstruct_rightside_minmod_muscl(     &
            primitive_values_set                       , &
            reference_cell_indexs_set                  , &
            cell_centor_positions                      , &
            face_centor_positions                      , &
            face_index                                   &
        )

        element_lef_and_right_side = integrated_element_function(  &
            lhc_primitive                            , &
            rhc_primitive                            , &
            cell_volumes(lhc_index)                  , &
            cell_volumes(rhc_index)                  , &
            face_normal_vectors     (face_index, 1:3), &
            face_tangential1_vectors(face_index, 1:3), &
            face_tangential2_vectors(face_index, 1:3), &
            face_areas              (face_index)     , &
            n_conservative_values                    , &
            n_derivative_values                      , &
            flux_function                            , &
            eos_pressure_function                    , &
            eos_soundspeed_function                  , &
            primitive_to_conservative_function         &
        )
    end function reconstruct_minmod_muscl
end module minmod_muscl_module