module five_equation_model_rho_thinc_module
    use typedef_module
    use vector_module
    use math_constant_module
    use abstract_mixture_eos
    use minmod_muscl_module

    implicit none

    private

    real(real_kind)  , parameter :: specified_slope_parameter_ = 0.5d0
    integer(int_kind), parameter :: num_ghost_cells_ = 3

    public :: reconstruct_rho_thinc

    contains

    pure function reconstruct_lhc_rho_thinc( &
        primitive_values_set               , &
        face_to_cell_index                 , &
        cell_positions                     , &
        face_positions                     , &
        face_index                         ) result(reconstructed_primitive)

        real   (real_kind), intent(in) :: primitive_values_set      (:, :)
        integer(int_kind ), intent(in) :: face_to_cell_index        (:, :)
        real   (real_kind), intent(in) :: cell_positions            (:, :)
        real   (real_kind), intent(in) :: face_positions            (:, :)
        integer(int_kind ), intent(in) :: face_index
        real   (real_kind)             :: reconstructed_primitive   (size(primitive_values_set(1, :)))
    end function reconstruct_lhc_rho_thinc

    pure function reconstruct_rho_thinc(        &
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
        eos                               , &
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
        class  (mixture_eos), intent(in) :: eos

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

            pure function primitive_to_conservative_function(primitive, eos) result(conservative)
                use typedef_module
                use abstract_mixture_eos
                real (real_kind  ), intent(in)  :: primitive   (:)
                class(mixture_eos), intent(in)  :: eos
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
                eos                               , &
                flux_function                     , &
                primitive_to_conservative_function   ) result(element)

                use typedef_module
                use abstract_mixture_eos

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
                class  (mixture_eos), intent(in) :: eos
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

                    pure function primitive_to_conservative_function(primitive, eos) result(conservative)
                        use typedef_module
                        use abstract_mixture_eos
                        real (real_kind  ), intent(in)  :: primitive   (:)
                        class(mixture_eos), intent(in)  :: eos
                        real (real_kind  ), allocatable :: conservative(:)
                    end function primitive_to_conservative_function
                end interface
            end function integrated_element_function
        end interface

        integer(int_kind )              :: lhc_index, rhc_index
        real   (real_kind), allocatable :: lhc_primitive   (:)
        real   (real_kind), allocatable :: rhc_primitive   (:)

        real   (real_kind)              :: lhc_interface_location, lhc_sign, lhc_a, lhc_b, lhc_d, lhc_e
        real   (real_kind)              :: rhc_interface_location, rhc_sign, rhc_a, rhc_b, rhc_d, rhc_e

        lhc_index = face_to_cell_index(face_index, num_ghost_cells_+0)
        rhc_index = face_to_cell_index(face_index, num_ghost_cells_+1)

        ! u, v, w, and p are reconstructed by minmod MUSCL
        lhc_primitive = reconstruct_lhc_minmod_muscl(    &
            primitive_values_set                       , &
            face_to_cell_index                         , &
            cell_centor_positions                      , &
            face_centor_positions                      , &
            face_index                                   &
        )
        rhc_primitive = reconstruct_rhc_minmod_muscl(    &
            primitive_values_set                       , &
            face_to_cell_index                         , &
            cell_centor_positions                      , &
            face_centor_positions                      , &
            face_index                                   &
        )

        ! # rho1, rho2, and z1 are reconstructed by rho-THINC
        ! ## LHC
        lhc_sign = sign(1.d0, primitive_values_set(rhc_index, 7) - primitive_values_set(lhc_index, 7))
        associate(                                          &
            lhc_rho1 => primitive_values_set(lhc_index, 1), &
            lhc_rho2 => primitive_values_set(lhc_index, 2), &
            lhc_z1   => primitive_values_set(lhc_index, 7)  &
        )
            lhc_a    = exp (2.d0 * lhc_sign * specified_slope_parameter_)
            lhc_b    = exp (2.d0 * lhc_sign * specified_slope_parameter_ * lhc_z1)
            ! Following condition means "lhc_z1 = rhc_z1 = 1 or lhc_z1 = rhc_z1 = 0". If this condition is NOT satisfied, we use MUSCL.
            if ((.not. lhc_b == 1.d0) .and. (.not. lhc_a == lhc_b)) then
                lhc_interface_location = 1.d0 / (2.d0 * specified_slope_parameter_) * log((lhc_b - 1.d0) / (lhc_a - lhc_b))
                lhc_d = exp(2.d0 * specified_slope_parameter_ * lhc_interface_location)
                lhc_e = log((lhc_a * lhc_d + 1.d0)**2.d0 / (lhc_a * (lhc_d + 1.d0)**2.d0))
                lhc_primitive(1) =  (4.d0 * lhc_sign * specified_slope_parameter_ * lhc_rho1 * lhc_z1) &
                                 /  (lhc_e + 2.d0 * lhc_sign * specified_slope_parameter_)
                lhc_primitive(2) = -(4.d0 * lhc_sign * specified_slope_parameter_ * lhc_rho2 * (1.d0 - lhc_z1)) &
                                 /  (lhc_e - 2.d0 * lhc_sign * specified_slope_parameter_)
                lhc_primitive(7) = 0.5d0 * (1.d0 + tanh(specified_slope_parameter_ * (lhc_sign + lhc_interface_location)))
            end if
        end associate
        ! ## RHC
        rhc_sign = sign(1.d0, primitive_values_set(rhc_index, 7) - primitive_values_set(lhc_index, 7))
        associate(                                          &
            rhc_rho1 => primitive_values_set(rhc_index, 1), &
            rhc_rho2 => primitive_values_set(rhc_index, 2), &
            rhc_z1   => primitive_values_set(rhc_index, 7)  &
        )
            rhc_a    = exp (2.d0 * rhc_sign * specified_slope_parameter_)
            rhc_b    = exp (2.d0 * rhc_sign * specified_slope_parameter_ * rhc_z1)
            ! Following condition means "lhc_z1 = rhc_z1 = 1 or lhc_z1 = rhc_z1 = 0". If this condition is NOT satisfied, we use MUSCL.
            if ((.not. rhc_b == 1.d0) .and. (.not. rhc_a == rhc_b)) then
                rhc_interface_location = 1.d0 / (2.d0 * specified_slope_parameter_) * log((rhc_b - 1.d0) / (rhc_a - rhc_b))
                rhc_d = exp(2.d0 * specified_slope_parameter_ * rhc_interface_location)
                rhc_e = log((rhc_a * rhc_d + 1.d0)**2.d0 / (rhc_a * (rhc_d + 1.d0)**2.d0))
                rhc_primitive(1) =  (4.d0 * rhc_sign * specified_slope_parameter_ * rhc_rho1 * rhc_z1) &
                                 /  (rhc_e + 2.d0 * rhc_sign * specified_slope_parameter_)
                rhc_primitive(2) = -(4.d0 * rhc_sign * specified_slope_parameter_ * rhc_rho2 * (1.d0 - rhc_z1)) &
                                 /  (rhc_e - 2.d0 * rhc_sign * specified_slope_parameter_)
                rhc_primitive(7) = 0.5d0 * (1.d0 + tanh(specified_slope_parameter_ * rhc_interface_location))
            end if
        end associate

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
            eos                                      , &
            flux_function                            , &
            primitive_to_conservative_function         &
        )
    end function reconstruct_rho_thinc
end module five_equation_model_rho_thinc_module