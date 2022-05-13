module jiang_weno5_module
    !**
    !* NOTE: If you use a grid that has curvature, variables are not interpolated by 5th order.
    !* Efficient Implementation of Weighted ENO Schemes (https://www.sciencedirect.com/science/article/pii/S0021999196901308)
    !*

    use typedef_module
    use vector_module

    implicit none

    private

    integer(int_kind) :: num_ghost_cells_ = 3
    real(real_kind) :: acm_coef_ = 33.d0

    public :: reconstruct_jiang_weno5

    contains

    pure function minmod(v1, v2) result(s)
        real(real_kind), intent(in) :: v1, v2
        real(real_kind)             :: s
        s = 0.5d0 * (sign(1.d0, v1) + sign(1.d0, v2)) * min(abs(v1), abs(v2))
    end function

    pure function acm_parameter(v_m1, v, v_p1) result(a)
        real(real_kind), intent(in) :: v_m1, v, v_p1
        real(real_kind)             :: a
        a = acm_coef_ * (abs(v_p1 - 2.d0 * v + v_p1) / (abs(v_p1 - v) + abs(v - v_m1)))
    end function

    pure function compute_acm_delta(vl, vr, v_m1, v, v_p1) result(delta)
        real(real_kind), intent(in) :: vl, vr, v_m1, v, v_p1
        real(real_kind)             :: delta
        delta = minmod(0.5d0 * acm_parameter(v_m1, v, v_p1) * (vr - vl), (v_p1 - vr))
    end function

    pure function compute_delta(vl, vr, v_m1, v, v_p1) result(delta)
        real(real_kind), intent(in) :: vl, vr, v_m1, v, v_p1
        real(real_kind)             :: delta
        delta = minmod(vl - v, v_p1 - vr)
    end function

    pure function compute_weights(s, v_m2, v_m1, v, v_p1, v_p2) result(w)
        real(real_kind), intent(in) :: s, v_m2, v_m1, v, v_p1, v_p2
        real(real_kind)             :: w(3)
        real(real_kind)             :: b(3), g(3)
        real(real_kind)             :: total_w

        g(1) = (80.d0 * s**4.d0 - 160.d0 * s**3.d0 - 120.d0 * s**2.d0 + 200.d0 * s + 9.d0) &
                 / (80.d0 * (12.d0 * s**2.d0 + 12.d0 * s - 1.d0))
        g(2) = -1.d0 * (960.d0 * s**6.d0 - 5360 * s**4.d0 + 4548.d0 * s**2 - 49.d0) &
                 / (40.d0 * (12.d0 * s**2 - 12.d0 * s - 1.d0) * (12.d0 ** 2 + 12.d0 * s - 1.d0))
        g(3) = (80.d0 * s**4 + 160.d0 * s**3 - 120.d0 * s**2 - 200.d0 * s + 9.d0) &
                 / 80.d0 * (12.d0 * s**2 - 12.d0 * s - 1.d0)
        b(1) = (13.d0 / 12.d0) * (v_m2 - 2.d0 * v_m1 + v)**2.d0 &
             + (1.d0  / 4.d0 ) * (v_m2 - 4.d0 * v_m1 + 3.d0 * v)**2.d0
        b(2) = (13.d0 / 12.d0) * (v_m1 - 2.d0 * v + v_p1)**2.d0 &
             + (1.d0  / 4.d0 ) * (v_m1 - v_p1)**2.d0
        b(3) = (13.d0 / 12.d0) * (v - 2.d0 * v_p1 + v_p2)**2.d0 &
             + (1.d0  / 4.d0 ) * (3.d0 * v - 4.d0 * v_p1 + v_p2)**2.d0
        w(1) = g(1) / (b(1) + 1.d-8)**2.d0
        w(2) = g(2) / (b(2) + 1.d-8)**2.d0
        w(3) = g(3) / (b(3) + 1.d-8)**2.d0
        total_w = w(1) + w(2) + w(3)
        w(1) = w(1) / total_w
        w(2) = w(2) / total_w
        w(3) = w(3) / total_w
    end function compute_weights

    pure function compute_polynomials(s, v_m2, v_m1, v, v_p1, v_p2) result(p)
        real(real_kind), intent(in) :: s, v_m2, v_m1, v, v_p1, v_p2
        real(real_kind)             :: p(3)

        p(1) = (1.d0 / 24.d0) * (12.d0 * s**2.d0 + 36.d0 * s + 23.d0) * v    &
             - (1.d0 / 12.d0) * (12.d0 * s**2.d0 + 24.d0 * s - 1.d0 ) * v_m1 &
             + (1.d0 / 24.d0) * (12.d0 * s**2.d0 + 12.d0 * s - 1.d0 ) * v_m2
        p(2) = (1.d0 / 24.d0) * (12.d0 * s**2.d0 + 12.d0 * s - 1.d0 ) * v_p1 &
             - (1.d0 / 12.d0) * (12.d0 * s**2.d0             - 13.d0) * v    &
             + (1.d0 / 24.d0) * (12.d0 * s**2.d0 - 12.d0 * s - 1.d0 ) * v_m1
        p(3) = (1.d0 / 24.d0) * (12.d0 * s**2.d0 - 12.d0 * s - 1.d0 ) * v_p2 &
             - (1.d0 / 12.d0) * (12.d0 * s**2.d0 - 24.d0 * s - 1.d0 ) * v_p1 &
             + (1.d0 / 24.d0) * (12.d0 * s**2.d0 - 36.d0 * s + 23.d0) * v
    end function compute_polynomials

    pure function reconstruct_leftside_jiang_weno5( &
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
            cell_pos_l(1:3) = cell_positions(reference_cell_indexs_set(face_index, num_ghost_cells_+0), 1:3)
            cell_pos_r(1:3) = cell_positions(reference_cell_indexs_set(face_index, num_ghost_cells_+1), 1:3)
            face_pos  (1:3) = face_positions(face_index, 1:3)
            cell_cell_distance = vector_distance(cell_pos_l, cell_pos_r)
            cell_face_distanse = vector_distance(cell_pos_l, face_pos  )
            s = cell_face_distanse / cell_cell_distance
            associate(                                                                                              &
                v_m2 => primitive_values_set(reference_cell_indexs_set(face_index, num_ghost_cells_-2), i), &
                v_m1 => primitive_values_set(reference_cell_indexs_set(face_index, num_ghost_cells_-1), i), &
                v    => primitive_values_set(reference_cell_indexs_set(face_index, num_ghost_cells_+0), i), &
                v_p1 => primitive_values_set(reference_cell_indexs_set(face_index, num_ghost_cells_+1), i), &
                v_p2 => primitive_values_set(reference_cell_indexs_set(face_index, num_ghost_cells_+2), i)  &
            )
                w(1:3) = compute_weights    (s, v_m2, v_m1, v, v_p1, v_p2)
                p(1:3) = compute_polynomials(s, v_m2, v_m1, v, v_p1, v_p2)
                reconstructed_primitive(i) = w(1) * p(1) + w(2) * p(2) + w(3) * p(3)
            end associate
        end do
    end function reconstruct_leftside_jiang_weno5

    pure function reconstruct_rightside_jiang_weno5( &
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
            cell_pos_l(1:3) = cell_positions(reference_cell_indexs_set(face_index, num_ghost_cells_+0), 1:3)
            cell_pos_r(1:3) = cell_positions(reference_cell_indexs_set(face_index, num_ghost_cells_+1), 1:3)
            face_pos  (1:3) = face_positions(face_index, 1:3)
            cell_cell_distance = vector_distance(cell_pos_l, cell_pos_r)
            cell_face_distanse = vector_distance(cell_pos_l, face_pos  )
            s = - 1.d0 * cell_face_distanse / cell_cell_distance
            associate(                                                                                              &
                v_m2       => primitive_values_set(reference_cell_indexs_set(face_index, num_ghost_cells_-1), i), &
                v_m1       => primitive_values_set(reference_cell_indexs_set(face_index, num_ghost_cells_+0), i), &
                v          => primitive_values_set(reference_cell_indexs_set(face_index, num_ghost_cells_+1), i), &
                v_p1       => primitive_values_set(reference_cell_indexs_set(face_index, num_ghost_cells_+2), i), &
                v_p2       => primitive_values_set(reference_cell_indexs_set(face_index, num_ghost_cells_+3), i)  &
            )
                w(1:3) = compute_weights    (s, v_m2, v_m1, v, v_p1, v_p2)
                p(1:3) = compute_polynomials(s, v_m2, v_m1, v, v_p1, v_p2)
                reconstructed_primitive(i) = w(1) * p(1) + w(2) * p(2) + w(3) * p(3)
            end associate
        end do
    end function reconstruct_rightside_jiang_weno5

    pure function reconstruct_jiang_weno5(        &
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

        real   (real_kind)              :: element_lef_and_right_side(2, n_conservative_values)

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
                real   (real_kind)              :: element_lef_and_right_side        (2, n_conservative_values)

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

        lhc_primitive = reconstruct_leftside_jiang_weno5(      &
            primitive_values_set                       , &
            reference_cell_indexs_set                  , &
            cell_centor_positions                      , &
            face_centor_positions                      , &
            face_index                                   &
        )
        rhc_primitive = reconstruct_rightside_jiang_weno5(     &
            primitive_values_set                       , &
            reference_cell_indexs_set                  , &
            cell_centor_positions                      , &
            face_centor_positions                      , &
            face_index                                   &
        )

        n_primitives = size(primitive_values_set(1, :))
        do i = 1, n_primitives, 1
            associate(                                                                                            &
                v_m1 => primitive_values_set(reference_cell_indexs_set(face_index, num_ghost_cells_-1), i), &
                v    => primitive_values_set(reference_cell_indexs_set(face_index, num_ghost_cells_+0), i), &
                v_p1 => primitive_values_set(reference_cell_indexs_set(face_index, num_ghost_cells_+1), i), &
                vl   => lhc_primitive(i), &
                vr   => rhc_primitive(i)  &
            )
                delta = compute_delta(vl, vr, v_m1, v, v_p1)
                lhc_primitive(i) = lhc_primitive(i) + delta
                rhc_primitive(i) = rhc_primitive(i) - delta
            end associate
        end do

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
            flux_function                            , &
            eos_pressure_function                    , &
            eos_soundspeed_function                  , &
            primitive_to_conservative_function         &
        )
    end function reconstruct_jiang_weno5
end module jiang_weno5_module