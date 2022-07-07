module weno5_module
    !**
    !* NOTE: If you use a grid that has curvature, variables are not interpolated by 5th order.
    !*

    use typedef_module
    use vector_module
    use abstract_eos

    implicit none

    private

    integer(int_kind) :: num_ghost_cells_ = 3

    public :: reconstruct_weno5

    contains

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
        b(1) = (1.d0 / 3.d0) &
            * (10.d0 * v**2.d0 - 31.d0 * v_m1 * v + 11.d0 * v_m2 * v + 25.d0 * v_m1**2.d0 - 19.d0 * v_m2 * v_m1 + 4.d0 * v_m2**2.d0)
        b(2) = (1.d0 / 3.d0) &
            * (4.d0 * v_p1**2.d0 - 13.d0 * v * v_p1 + 5.d0 * v_m1 * v_p1 + 13.d0 * v**2.d0 - 13.d0 * v_m1 * v + 4.d0 * v_m1**2.d0)
        b(3) = (1.d0 / 3.d0) &
            * (4.d0 * v_p2**2.d0 - 19.d0 * v_p1 * v_p2 + 11.d0 * v * v_p2 + 25.d0 * v_p1**2.d0 - 31.d0 * v * v_p1 + 10.d0 * v**2)
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

    pure function reconstruct_lhc_weno5( &
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
        real   (real_kind)             :: reconstructed_primitive   (size(primitive_values_set(1,:)))

        integer(int_kind ) :: n_primitives, i
        real   (real_kind) :: w(3), p(3), cell_pos_l(3), cell_pos_r(3), face_pos(3)
        real   (real_kind) :: cell_cell_distance, cell_face_distanse, s

        n_primitives = size(primitive_values_set(1,:))

        do i = 1, n_primitives, 1
            cell_pos_l(1:3) = cell_positions(face_to_cell_index(face_index, num_ghost_cells_+0), 1:3)
            cell_pos_r(1:3) = cell_positions(face_to_cell_index(face_index, num_ghost_cells_+1), 1:3)
            face_pos  (1:3) = face_positions(face_index, 1:3)
            cell_cell_distance = vector_distance(cell_pos_l, cell_pos_r)
            cell_face_distanse = vector_distance(cell_pos_l, face_pos  )
            s = 0.5d0
            associate(                                                                                              &
                v_m2 => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_-2), i), &
                v_m1 => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_-1), i), &
                v    => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_+0), i), &
                v_p1 => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_+1), i), &
                v_p2 => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_+2), i)  &
            )
                w(1:3) = compute_weights    (s, v_m2, v_m1, v, v_p1, v_p2)
                p(1:3) = compute_polynomials(s, v_m2, v_m1, v, v_p1, v_p2)
                reconstructed_primitive(i) = w(1) * p(1) + w(2) * p(2) + w(3) * p(3)
            end associate
        end do
    end function reconstruct_lhc_weno5

    pure function reconstruct_rhc_weno5( &
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

        n_primitives = size(primitive_values_set(1,:))

        do i = 1, n_primitives, 1
            cell_pos_l(1:3) = cell_positions(face_to_cell_index(face_index, num_ghost_cells_+0), 1:3)
            cell_pos_r(1:3) = cell_positions(face_to_cell_index(face_index, num_ghost_cells_+1), 1:3)
            face_pos  (1:3) = face_positions(face_index, 1:3)
            cell_cell_distance = vector_distance(cell_pos_l, cell_pos_r)
            cell_face_distanse = vector_distance(cell_pos_l, face_pos  )
            s = - 0.5d0
            associate(                                                                                              &
                v_m2       => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_-1), i), &
                v_m1       => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_+0), i), &
                v          => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_+1), i), &
                v_p1       => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_+2), i), &
                v_p2       => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_+3), i)  &
            )
                w(1:3) = compute_weights    (s, v_m2, v_m1, v, v_p1, v_p2)
                p(1:3) = compute_polynomials(s, v_m2, v_m1, v, v_p1, v_p2)
                reconstructed_primitive(i) = w(1) * p(1) + w(2) * p(2) + w(3) * p(3)
            end associate
        end do
    end function reconstruct_rhc_weno5

    pure function reconstruct_weno5(        &
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
        real   (real_kind), allocatable :: reconstructed_lhc_primitive   (:)
        real   (real_kind), allocatable :: reconstructed_rhc_primitive   (:)

        lhc_index = face_to_cell_index(face_index, num_ghost_cells_+0)
        rhc_index = face_to_cell_index(face_index, num_ghost_cells_+1)

        reconstructed_lhc_primitive = reconstruct_lhc_weno5(      &
            primitive_values_set                       , &
            face_to_cell_index                         , &
            cell_centor_positions                      , &
            face_centor_positions                      , &
            face_index                                   &
        )
        reconstructed_rhc_primitive = reconstruct_rhc_weno5(     &
            primitive_values_set                       , &
            face_to_cell_index                         , &
            cell_centor_positions                      , &
            face_centor_positions                      , &
            face_index                                   &
        )

        element = integrated_element_function(  &
            reconstructed_lhc_primitive              , &
            reconstructed_rhc_primitive              , &
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
    end function reconstruct_weno5
end module weno5_module