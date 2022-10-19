module five_equation_model_module
    use vector_module
    use typedef_module
    use stdio_module
    use string_utils_module
    use abstract_configuration
    use abstract_eos
    use abstract_riemann_solver
    use five_equation_model_utils_module
    use matrix_module
    use vector_module

    implicit none

    private

    public :: compute_residual_element
    public :: spectral_radius

    contains

    pure function compute_mixture_density(primitive_variables) result(density)
        real (real_kind), intent(in) :: primitive_variables(:)
        real (real_kind)             :: density
        associate(                            &
            rho1   => primitive_variables(1), &
            rho2   => primitive_variables(2), &
            alpha1 => primitive_variables(7)  &
        )
            density = alpha1 * rho1 + (1.d0 - alpha1) * rho2
        end associate
    end function


    pure function compute_residual_element(       &
        an_eos                                  , &
        an_riemann_solver                       , &
        primitive_variables_lhc                 , &
        primitive_variables_rhc                 , &
        reconstructed_primitive_variables_lhc   , &
        reconstructed_primitive_variables_rhc   , &
        lhc_cell_volume                         , &
        rhc_cell_volume                         , &
        face_area                               , &
        face_normal_vector                      , &
        face_tangential1_vector                 , &
        face_tangential2_vector                 , &
        num_conservative_values                 , &
        num_primitive_values                      ) result(residual_element)

        class  (eos           ), intent(in) :: an_eos
        class  (riemann_solver), intent(in) :: an_riemann_solver
        real   (real_kind     ), intent(in) :: primitive_variables_lhc              (:)
        real   (real_kind     ), intent(in) :: primitive_variables_rhc              (:)
        real   (real_kind     ), intent(in) :: reconstructed_primitive_variables_lhc(:)
        real   (real_kind     ), intent(in) :: reconstructed_primitive_variables_rhc(:)
        real   (real_kind     ), intent(in) :: lhc_cell_volume
        real   (real_kind     ), intent(in) :: rhc_cell_volume
        real   (real_kind     ), intent(in) :: face_area
        real   (real_kind     ), intent(in) :: face_normal_vector     (3)
        real   (real_kind     ), intent(in) :: face_tangential1_vector(3)
        real   (real_kind     ), intent(in) :: face_tangential2_vector(3)
        integer(int_kind      ), intent(in) :: num_conservative_values
        integer(int_kind      ), intent(in) :: num_primitive_values

        real   (real_kind)                  :: residual_element(num_conservative_values, 1:2)

        ! ## face coordinate variables
        real(real_kind) :: local_coordinate_primitives_lhc  (num_primitive_values)
        real(real_kind) :: local_coordinate_primitives_rhc  (num_primitive_values)
        real(real_kind) :: local_coordinate_conservative_lhc(num_conservative_values)
        real(real_kind) :: local_coordinate_conservative_rhc(num_conservative_values)

        ! ## cell variables
        real(real_kind) :: lhc_soundspeed, lhc_pressure, lhc_density, lhc_main_velocity
        real(real_kind) :: rhc_soundspeed, rhc_pressure, rhc_density, rhc_main_velocity

        ! ## rieman solver
        real(real_kind) :: nonviscosity_flux      (num_conservative_values)
        real(real_kind) :: rieman_solver_features(7)

        ! ## face variables
        real(real_kind) :: numerical_velocity, interface_volume_fraction

        ! ## kapila Kdiv(u)
        real(real_kind) :: rhc_k, lhc_k


        ! # compute primitive-variables of face local coordinate
        ! ## left-side
        local_coordinate_primitives_lhc(:) = rotate_primitive( &
            reconstructed_primitive_variables_lhc           , &
            face_normal_vector                              , &
            face_tangential1_vector                         , &
            face_tangential2_vector                         , &
            num_primitive_values                              &
        )
        ! ## right-side
        local_coordinate_primitives_rhc(:) = rotate_primitive( &
            reconstructed_primitive_variables_rhc           , &
            face_normal_vector                              , &
            face_tangential1_vector                         , &
            face_tangential2_vector                         , &
            num_primitive_values                              &
        )

        ! # compute conservative-variables
        local_coordinate_conservative_lhc = primitive_to_conservative( &
            an_eos                                                   , &
            local_coordinate_primitives_lhc                          , &
            num_conservative_values                                    &
        )
        local_coordinate_conservative_rhc = primitive_to_conservative( &
            an_eos                                                   , &
            local_coordinate_primitives_rhc                          , &
            num_conservative_values                                    &
        )

        ! # compute EoS, main velosity, and fluxs
        associate(                                                 &
                rho1    => local_coordinate_primitives_lhc(1),     &
                rho2    => local_coordinate_primitives_lhc(2),     &
                u       => local_coordinate_primitives_lhc(3),     &
                v       => local_coordinate_primitives_lhc(4),     &
                w       => local_coordinate_primitives_lhc(5),     &
                p       => local_coordinate_primitives_lhc(6),     &
                z1      => local_coordinate_primitives_lhc(7)      &
            )
            lhc_density    = rho1 * z1 + rho2 * (1.d0 - z1)
            lhc_pressure   = p
            lhc_soundspeed = an_eos%compute_soundspeed(p, lhc_density, z1)
            lhc_main_velocity = u
        end associate
        associate(                                             &
                rho1    => local_coordinate_primitives_rhc(1), &
                rho2    => local_coordinate_primitives_rhc(2), &
                u       => local_coordinate_primitives_rhc(3), &
                v       => local_coordinate_primitives_rhc(4), &
                w       => local_coordinate_primitives_rhc(5), &
                p       => local_coordinate_primitives_rhc(6), &
                z1      => local_coordinate_primitives_rhc(7)  &
            )
            rhc_density    = rho1 * z1 + rho2 * (1.d0 - z1)
            rhc_pressure   = p
            rhc_soundspeed = an_eos%compute_soundspeed(p, rhc_density, z1)
            rhc_main_velocity = u
        end associate

        ! # compute flux
        ! ## local coordinate flux
        rieman_solver_features(:) = an_riemann_solver%compute_features( &
            lhc_main_velocity               , &
            lhc_density                     , &
            lhc_pressure                    , &
            lhc_soundspeed                  , &
            rhc_main_velocity               , &
            rhc_density                     , &
            rhc_pressure                    , &
            rhc_soundspeed                    &
        )
        nonviscosity_flux(1) = an_riemann_solver%compute_mass_flux( &
            local_coordinate_conservative_lhc(1), &
            lhc_main_velocity                   , &
            lhc_density                         , &
            lhc_pressure                        , &
            lhc_soundspeed                      , &
            local_coordinate_conservative_rhc(1), &
            rhc_main_velocity                   , &
            rhc_density                         , &
            rhc_pressure                        , &
            rhc_soundspeed                      , &
            rieman_solver_features                &
        )
        nonviscosity_flux(2) = an_riemann_solver%compute_mass_flux( &
            local_coordinate_conservative_lhc(2), &
            lhc_main_velocity                   , &
            lhc_density                         , &
            lhc_pressure                        , &
            lhc_soundspeed                      , &
            local_coordinate_conservative_rhc(2), &
            rhc_main_velocity                   , &
            rhc_density                         , &
            rhc_pressure                        , &
            rhc_soundspeed                      , &
            rieman_solver_features                &
        )
        nonviscosity_flux(3:5) = an_riemann_solver%compute_momentum_flux( &
            local_coordinate_conservative_lhc(3:5), &
            lhc_main_velocity                     , &
            lhc_density                           , &
            lhc_pressure                          , &
            lhc_soundspeed                        , &
            local_coordinate_conservative_rhc(3:5), &
            rhc_main_velocity                     , &
            rhc_density                           , &
            rhc_pressure                          , &
            rhc_soundspeed                        , &
            rieman_solver_features                  &
        )
        nonviscosity_flux(6) = an_riemann_solver%compute_energy_flux( &
            local_coordinate_conservative_lhc(6), &
            lhc_main_velocity                   , &
            lhc_density                         , &
            lhc_pressure                        , &
            lhc_soundspeed                      , &
            local_coordinate_conservative_rhc(6), &
            rhc_main_velocity                   , &
            rhc_density                         , &
            rhc_pressure                        , &
            rhc_soundspeed                      , &
            rieman_solver_features                &
        )
        nonviscosity_flux(7) = an_riemann_solver%compute_volume_fraction_flux( &
            local_coordinate_conservative_lhc(7), &
            lhc_main_velocity                   , &
            lhc_density                         , &
            lhc_pressure                        , &
            lhc_soundspeed                      , &
            local_coordinate_conservative_rhc(7), &
            rhc_main_velocity                   , &
            rhc_density                         , &
            rhc_pressure                        , &
            rhc_soundspeed                      , &
            rieman_solver_features                &
        )

        ! ## convert to global coordinate flux
        nonviscosity_flux(3:5) = vector_unrotate(  &
            nonviscosity_flux(3:5)              , &
            face_normal_vector                  , &
            face_tangential1_vector             , &
            face_tangential2_vector               &
        )

        ! # summation nonviscosity-flux
        residual_element(1:7, 1) = (-1.d0 / lhc_cell_volume) * nonviscosity_flux(:) * face_area
        residual_element(1:7, 2) = (+1.d0 / rhc_cell_volume) * nonviscosity_flux(:) * face_area

        ! # Compute interface values
        numerical_velocity = an_riemann_solver%compute_numerical_velocity( &
            lhc_main_velocity                   , &
            lhc_density                         , &
            lhc_pressure                        , &
            lhc_soundspeed                      , &
            rhc_main_velocity                   , &
            rhc_density                         , &
            rhc_pressure                        , &
            rhc_soundspeed                      , &
            rieman_solver_features                &
        )
        interface_volume_fraction = an_riemann_solver%compute_interface_value( &
            local_coordinate_primitives_lhc(7)  , &
            lhc_main_velocity                   , &
            lhc_density                         , &
            lhc_pressure                        , &
            lhc_soundspeed                      , &
            local_coordinate_primitives_rhc(7)  , &
            rhc_main_velocity                   , &
            rhc_density                         , &
            rhc_pressure                        , &
            rhc_soundspeed                      , &
            rieman_solver_features                &
        )

        ! # (-z1 - K) * div(u)
        associate(                                       &
            lhc_rhos    => primitive_variables_lhc(1:2), &
            lhc_p       => primitive_variables_lhc(6)  , &
            lhc_z1      => primitive_variables_lhc(7)  , &
            rhc_rhos    => primitive_variables_rhc(1:2), &
            rhc_p       => primitive_variables_rhc(6)  , &
            rhc_z1      => primitive_variables_rhc(7)    &
        )

            lhc_k = an_eos%compute_k(lhc_p, lhc_rhos, lhc_z1)
            rhc_k = an_eos%compute_k(rhc_p, rhc_rhos, rhc_z1)

            residual_element(7, 1) = residual_element(7, 1) &
                                   - (-lhc_z1 - lhc_k) * (1.d0 / lhc_cell_volume) * numerical_velocity * face_area
            residual_element(7, 2) = residual_element(7, 2) &
                                   + (-rhc_z1 - rhc_k) * (1.d0 / rhc_cell_volume) * numerical_velocity * face_area
        end associate
    end function compute_residual_element

    pure function spectral_radius(an_eos, primitive_variables, length) result(r)
        class(eos      ), intent(in) :: an_eos
        real (real_kind), intent(in) :: primitive_variables(:)
        real (real_kind), intent(in) :: length
        real (real_kind) :: r

        associate(                                                    &
            density  => compute_mixture_density(primitive_variables), &
            velocity => primitive_variables(3:5)                    , &
            pressure => primitive_variables(6)                      , &
            alpha1   => primitive_variables(7)                        &
        )
            r = an_eos%compute_soundspeed(pressure, density, alpha1) + vector_magnitude(velocity)
        end associate
    end function spectral_radius
end module five_equation_model_module