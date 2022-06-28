module mp_weno5_js_module
    !**
    !* NOTE: If you use a grid that has curvature, variables are not interpolated by 5th order.
    !* Efficient Implementation of Weighted ENO Schemes (https://www.sciencedirect.com/science/article/pii/S0021999196901308)
    !*

    use typedef_module
    use vector_module
    use abstract_eos

    implicit none

    private

    integer(int_kind) :: num_ghost_cells_ = 3

    ! monotonicity-preserving (MP)
    real   (real_kind) :: alpha_ = 10.d0
    real   (real_kind) :: beta_  = 4.d0

    public :: reconstruct_mp_weno5_js
    public :: reconstruct_lhc_mp_weno5_js
    public :: reconstruct_rhc_mp_weno5_js

    contains

    pure function minmod(x, y) result(m)
        real(real_kind), intent(in) :: x, y
        real(real_kind)             :: m
        m = 0.5d0 * (sign(1.d0, x) + sign(1.d0, y)) * min(abs(x), abs(y))
    end function minmod

    pure function median(x, y, z) result(m)
        real(real_kind), intent(in) :: x, y, z
        real(real_kind)             :: m
        m = x + minmod(y - x, z - x)
    end function median

    pure function curvature_measure(v_m1, v, v_p1) result(c)
        real(real_kind), intent(in) :: v_m1, v, v_p1
        real(real_kind)             :: c
        c = v_p1 - 2.d0 * v + v_m1
    end function curvature_measure

    pure function compute_weights(s, v_m2, v_m1, v, v_p1, v_p2) result(w)
        real(real_kind), intent(in) :: s, v_m2, v_m1, v, v_p1, v_p2
        real(real_kind)             :: w(3)
        real(real_kind)             :: b(3), g(3)
        real(real_kind)             :: total_w

        if(s>0)then ! left-side
            g(1) = 0.1d0
            g(2) = 0.6d0
            g(3) = 0.3d0
        else
            g(1) = 0.3d0
            g(2) = 0.6d0
            g(3) = 0.1d0
        endif
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

        if(s>0)then ! left-side
            p(1) =  1.d0 / 3.d0 * v_m2 - 7.d0 / 6.d0 * v_m1 + 11.d0 / 6.d0 * v
            p(2) = -1.d0 / 6.d0 * v_m1 + 5.d0 / 6.d0 * v    + 1.d0  / 3.d0 * v_p1
            p(3) =  1.d0 / 3.d0 * v    + 5.d0 / 6.d0 * v_p1 - 1.d0  / 6.d0 * v_p2
        else
            p(1) = -1.d0  / 6.d0 * v_m2 + 5.d0 / 6.d0 * v_m1 + 1.d0 / 3.d0 * v
            p(2) =  1.d0  / 3.d0 * v_m1 + 5.d0 / 6.d0 * v    - 1.d0 / 6.d0 * v_p1
            p(3) =  11.d0 / 6.d0 * v    - 7.d0 / 6.d0 * v_p1 + 1.d0 / 3.d0 * v_p2
        endif
    end function compute_polynomials

    pure function reconstruct_lhc_mp_weno5_js( &
        primitive_values_set             , &
        face_to_cell_index        , &
        cell_positions                   , &
        face_positions                   , &
        face_index                          ) result(reconstructed_primitive)

        real   (real_kind), intent(in) :: primitive_values_set      (:, :)
        integer(int_kind ), intent(in) :: face_to_cell_index (:, :)
        real   (real_kind), intent(in) :: cell_positions            (:, :)
        real   (real_kind), intent(in) :: face_positions            (:, :)
        integer(int_kind ), intent(in) :: face_index
        real   (real_kind)             :: reconstructed_primitive   (size(primitive_values_set(1, :)))

        integer(int_kind ) :: n_primitives, i
        real   (real_kind) :: w(3), p(3), cell_pos_l(3), cell_pos_r(3), face_pos(3)
        real   (real_kind) :: cell_cell_distance, cell_face_distanse, s
        real   (real_kind) :: d_p1, d, d_m1, dm4
        real   (real_kind) :: v_ul, v_md, v_lc, v_min, v_max

        n_primitives = size(primitive_values_set(1, :))

        do i = 1, n_primitives, 1
            cell_pos_l(1:3) = cell_positions(face_to_cell_index(face_index, num_ghost_cells_+0), 1:3)
            cell_pos_r(1:3) = cell_positions(face_to_cell_index(face_index, num_ghost_cells_+1), 1:3)
            face_pos  (1:3) = face_positions(face_index, 1:3)
            cell_cell_distance = vector_distance(cell_pos_l, cell_pos_r)
            cell_face_distanse = vector_distance(cell_pos_l, face_pos  )
            s = 1.d0!cell_face_distanse / cell_cell_distance
            associate(                                                                                              &
                v_m2 => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_-2), i), &
                v_m1 => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_-1), i), &
                v    => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_+0), i), &
                v_p1 => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_+1), i), &
                v_p2 => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_+2), i)  &
            )
                ! WENO-JS
                w(1:3) = compute_weights    (s, v_m2, v_m1, v, v_p1, v_p2)
                p(1:3) = compute_polynomials(s, v_m2, v_m1, v, v_p1, v_p2)
                reconstructed_primitive(i) = w(1) * p(1) + w(2) * p(2) + w(3) * p(3)
                ! MP
                d_m1 = curvature_measure(v_m2, v_m1, v   )
                d    = curvature_measure(v_m1, v   , v_p1)
                d_p1 = curvature_measure(v   , v_p1, v_p2)
                dm4  = minmod(minmod(4.d0 * d - d_p1, 4.d0 * d_p1 - d), minmod(d, d_p1))
                v_ul = v + alpha_ * (v - v_m1)
                v_md = 0.5d0 * (v + v_p1 - dm4)
                v_lc = v + 0.5d0 * (v - v_m1) + (beta_ / 3.d0) * dm4
                v_min = max(min(v, v_p1, v_md), min(v, v_ul, v_lc))
                v_max = min(max(v, v_p1, v_md), max(v, v_ul, v_lc))
                reconstructed_primitive(i) = median(reconstructed_primitive(i), v_min, v_max)
            end associate
        end do
    end function reconstruct_lhc_mp_weno5_js

    pure function reconstruct_rhc_mp_weno5_js( &
        primitive_values_set             , &
        face_to_cell_index               , &
        cell_positions                   , &
        face_positions                   , &
        face_index                          ) result(reconstructed_primitive)

        real   (real_kind), intent(in) :: primitive_values_set      (:, :)
        integer(int_kind ), intent(in) :: face_to_cell_index (:, :)
        real   (real_kind), intent(in) :: cell_positions            (:, :)
        real   (real_kind), intent(in) :: face_positions            (:, :)
        integer(int_kind ), intent(in) :: face_index
        real   (real_kind)             :: reconstructed_primitive   (size(primitive_values_set(1, :)))

        integer(int_kind ) :: n_primitives, i
        real   (real_kind) :: w(3), p(3), cell_pos_l(3), cell_pos_r(3), face_pos(3)
        real   (real_kind) :: cell_cell_distance, cell_face_distanse, s
        real   (real_kind) :: d_p1, d, d_m1, dm4
        real   (real_kind) :: v_ul, v_md, v_lc, v_min, v_max

        n_primitives = size(primitive_values_set(1, :))

        do i = 1, n_primitives, 1
            cell_pos_l(1:3) = cell_positions(face_to_cell_index(face_index, num_ghost_cells_+0), 1:3)
            cell_pos_r(1:3) = cell_positions(face_to_cell_index(face_index, num_ghost_cells_+1), 1:3)
            face_pos  (1:3) = face_positions(face_index, 1:3)
            cell_cell_distance = vector_distance(cell_pos_l, cell_pos_r)
            cell_face_distanse = vector_distance(cell_pos_l, face_pos  )
            s = - 1.d0! * cell_face_distanse / cell_cell_distance
            associate(                                                                                              &
                v_m2       => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_-1), i), &
                v_m1       => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_+0), i), &
                v          => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_+1), i), &
                v_p1       => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_+2), i), &
                v_p2       => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_+3), i)  &
            )
                ! WENO5-JS
                w(1:3) = compute_weights    (s, v_m2, v_m1, v, v_p1, v_p2)
                p(1:3) = compute_polynomials(s, v_m2, v_m1, v, v_p1, v_p2)
                reconstructed_primitive(i) = w(1) * p(1) + w(2) * p(2) + w(3) * p(3)
                ! MP
                d_m1 = curvature_measure(v_m2, v_m1, v   )
                d    = curvature_measure(v_m1, v   , v_p1)
                d_p1 = curvature_measure(v   , v_p1, v_p2)
                dm4  = minmod(minmod(4.d0 * d - d_m1, 4.d0 * d_m1 - d), minmod(d, d_m1))
                v_ul = v + alpha_ * (v - v_p1)
                v_md = 0.5d0 * (v + v_m1 - dm4)
                v_lc = v + 0.5d0 * (v - v_p1) + (beta_ / 3.d0) * dm4
                v_min = max(min(v, v_m1, v_md), min(v, v_ul, v_lc))
                v_max = min(max(v, v_m1, v_md), max(v, v_ul, v_lc))
                reconstructed_primitive(i) = median(reconstructed_primitive(i), v_min, v_max)
            end associate
        end do
    end function reconstruct_rhc_mp_weno5_js

    pure function reconstruct_mp_weno5_js(  &
        primitive_values_set              , &
        cell_centor_positions             , &
        cell_volumes                      , &
        face_to_cell_index                , &
        face_centor_positions             , &
        face_normal_vectors               , &
        face_tangential1_vectors          , &
        face_tangential2_vectors          , &
        face_areas                        , &
        face_index                        , &
        n_conservative_values             , &
        n_derivative_values               , &
        an_eos                               , &
        flux_function                     , &
        primitive_to_conservative_function, &
        integrated_element_function             ) result(element)

        real   (real_kind  ), intent(in) :: primitive_values_set     (:, :)
        real   (real_kind  ), intent(in) :: cell_centor_positions    (:, :)
        real   (real_kind  ), intent(in) :: cell_volumes             (:)
        integer(int_kind   ), intent(in) :: face_to_cell_index       (:,:)
        real   (real_kind  ), intent(in) :: face_normal_vectors      (:,:)
        real   (real_kind  ), intent(in) :: face_tangential1_vectors (:,:)
        real   (real_kind  ), intent(in) :: face_tangential2_vectors (:,:)
        real   (real_kind  ), intent(in) :: face_centor_positions    (:,:)
        real   (real_kind  ), intent(in) :: face_areas               (:)
        integer(int_kind   ), intent(in) :: face_index
        integer(int_kind   ), intent(in) :: n_conservative_values
        integer(int_kind   ), intent(in) :: n_derivative_values
        class  (eos), intent(in) :: an_eos

        real   (real_kind)              :: element(2, n_conservative_values+n_derivative_values)

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

            pure function primitive_to_conservative_function(primitive, an_eos) result(conservative)
                use typedef_module
                use abstract_eos
                real (real_kind  ), intent(in)  :: primitive   (:)
                class(eos), intent(in)  :: an_eos
                real (real_kind  ), allocatable :: conservative(:)
            end function primitive_to_conservative_function

            pure function integrated_element_function( &
                reconstructed_lhc_primitive       , &
                reconstructed_rhc_primitive       , &
                lhc_primitive                     , &
                rhc_primitive                     , &
                lhc_cell_volume                   , &
                rhc_cell_volume                   , &
                face_normal_vector                , &
                face_tangential1_vector           , &
                face_tangential2_vector           , &
                face_area                         , &
                n_conservative_values             , &
                n_derivative_values               , &
                an_eos                               , &
                flux_function                     , &
                primitive_to_conservative_function   ) result(element)

                use typedef_module
                use abstract_eos

                real   (real_kind  ), intent(in) :: reconstructed_lhc_primitive (:)
                real   (real_kind  ), intent(in) :: reconstructed_rhc_primitive (:)
                real   (real_kind  ), intent(in) :: lhc_primitive               (:)
                real   (real_kind  ), intent(in) :: rhc_primitive               (:)
                real   (real_kind  ), intent(in) :: lhc_cell_volume
                real   (real_kind  ), intent(in) :: rhc_cell_volume
                real   (real_kind  ), intent(in) :: face_normal_vector                (3)
                real   (real_kind  ), intent(in) :: face_tangential1_vector           (3)
                real   (real_kind  ), intent(in) :: face_tangential2_vector           (3)
                real   (real_kind  ), intent(in) :: face_area
                integer(int_kind   ), intent(in) :: n_conservative_values
                integer(int_kind   ), intent(in) :: n_derivative_values
                class  (eos), intent(in) :: an_eos
                real   (real_kind)               :: element        (2, n_conservative_values+n_derivative_values)

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

                    pure function primitive_to_conservative_function(primitive, an_eos) result(conservative)
                        use typedef_module
                        use abstract_eos
                        real (real_kind  ), intent(in)  :: primitive   (:)
                        class(eos), intent(in)  :: an_eos
                        real (real_kind  ), allocatable :: conservative(:)
                    end function primitive_to_conservative_function
                end interface
            end function integrated_element_function
        end interface

        integer(int_kind )              :: lhc_index, rhc_index
        real   (real_kind), allocatable :: lhc_primitive   (:)
        real   (real_kind), allocatable :: rhc_primitive   (:)
        integer(int_kind )              :: n_primitives, i
        real   (real_kind)              :: delta

        lhc_index = face_to_cell_index(face_index, num_ghost_cells_+0)
        rhc_index = face_to_cell_index(face_index, num_ghost_cells_+1)

        lhc_primitive = reconstruct_lhc_mp_weno5_js(      &
            primitive_values_set                       , &
            face_to_cell_index                         , &
            cell_centor_positions                      , &
            face_centor_positions                      , &
            face_index                                   &
        )
        rhc_primitive = reconstruct_rhc_mp_weno5_js(     &
            primitive_values_set                       , &
            face_to_cell_index                         , &
            cell_centor_positions                      , &
            face_centor_positions                      , &
            face_index                                   &
        )

        element = integrated_element_function(  &
            lhc_primitive                            , &
            rhc_primitive                            , &
            primitive_values_set(lhc_index, :)       , &
            primitive_values_set(rhc_index, :)       , &
            cell_volumes(lhc_index)                  , &
            cell_volumes(rhc_index)                  , &
            face_normal_vectors     (face_index, 1:3), &
            face_tangential1_vectors(face_index, 1:3), &
            face_tangential2_vectors(face_index, 1:3), &
            face_areas              (face_index)     , &
            n_conservative_values                    , &
            n_derivative_values                      , &
            an_eos                                      , &
            flux_function                            , &
            primitive_to_conservative_function         &
        )
    end function reconstruct_mp_weno5_js
end module mp_weno5_js_module