module five_equation_model_module
    use vector_module
    use typedef_module
    use math_constant_module
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

    ! for Kdiv(u) term
    logical :: ignore_kdivu_

    public :: flux_function
    public :: spectral_radius
    public :: initialize_model

    contains

    subroutine initialize_model(a_configuration)
        class  (configuration), intent(inout) :: a_configuration

        logical           :: found

        call a_configuration%get_bool      ("Model.Ignore kdivu", ignore_kdivu_, found, .false.)
        if(.not. found) call write_warring("'Model.Ignore kdivu' is not found in configuration you set. Apply this term.")
    end subroutine initialize_model

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

    pure function fix_reconstructed_primitive_variables(reconstructed_primitive_variables) result(fixed_reconstructed_primitive_variables)
        real   (real_kind), intent(in) :: reconstructed_primitive_variables(7)
        real   (real_kind)             :: fixed_reconstructed_primitive_variables(7)

        associate(                                         &
            rho1  => reconstructed_primitive_variables(1), &
            rho2  => reconstructed_primitive_variables(2), &
            u     => reconstructed_primitive_variables(3), &
            v     => reconstructed_primitive_variables(4), &
            w     => reconstructed_primitive_variables(5), &
            p     => reconstructed_primitive_variables(6), &
            alpha => reconstructed_primitive_variables(7)  &
        )
            fixed_reconstructed_primitive_variables(1) = max(rho1, 0.d0)
            fixed_reconstructed_primitive_variables(2) = max(rho2, 0.d0)
            fixed_reconstructed_primitive_variables(3) = u
            fixed_reconstructed_primitive_variables(4) = v
            fixed_reconstructed_primitive_variables(5) = w
            fixed_reconstructed_primitive_variables(6) = p
            fixed_reconstructed_primitive_variables(7) = max(0.d0, min(alpha, 1.d0))
        end associate
    end function fix_reconstructed_primitive_variables

    pure function compute_soundspeed(an_eos, pressure, densities, volume_fractions) result(c)
        class(eos      ), intent(in) :: an_eos
        real (real_kind), intent(in) :: pressure
        real (real_kind), intent(in) :: densities(:)
        real (real_kind), intent(in) :: volume_fractions(:)
        real (real_kind) :: c

        if(ignore_kdivu_)then
            c = an_eos%compute_mixture_soundspeed(pressure, densities, volume_fractions)
        else
            c = an_eos%compute_wood_soundspeed(pressure, densities, volume_fractions)
        endif
    end function compute_soundspeed

    pure function compute_k(an_eos, pressure, densities, volume_fraction) result(k)
        class(eos      ), intent(in) :: an_eos
        real (real_kind), intent(in) :: pressure
        real (real_kind), intent(in) :: densities(:)
        real (real_kind), intent(in) :: volume_fraction
        real (real_kind) :: k

        if((volume_fraction <= machine_epsilon) .or. (1.0_real_kind <= volume_fraction))then
            k = 0.0_real_kind
            return
        endif

        associate(                                                                                            &
            rho1_powc1 => densities(1) * an_eos%compute_soundspeed(pressure, densities(1), 1)**2.0_real_kind, &
            rho2_powc2 => densities(2) * an_eos%compute_soundspeed(pressure, densities(2), 2)**2.0_real_kind  &
        )
            k = (volume_fraction * (1.0_real_kind - volume_fraction) * (rho2_powc2 * rho1_powc1)) &
              / (volume_fraction * rho2_powc2 + (1.0_real_kind - volume_fraction) * rho1_powc1)
        end associate
    end function compute_k

    function flux_function(                       &
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
        num_primitive_values                      ) result(flux)

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

        real   (real_kind)                  :: flux(num_conservative_values, 1:2)

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
        local_coordinate_primitives_lhc(:) = rotate_primitive(                            &
            fix_reconstructed_primitive_variables(reconstructed_primitive_variables_lhc), &
            face_normal_vector                                                          , &
            face_tangential1_vector                                                     , &
            face_tangential2_vector                                                     , &
            num_primitive_values                                                          &
        )
        ! ## right-side
        local_coordinate_primitives_rhc(:) = rotate_primitive(                            &
            fix_reconstructed_primitive_variables(reconstructed_primitive_variables_rhc), &
            face_normal_vector                                                          , &
            face_tangential1_vector                                                     , &
            face_tangential2_vector                                                     , &
            num_primitive_values                                                          &
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
        associate(                                              &
                rhos => local_coordinate_primitives_lhc(1:2),   &
                u    => local_coordinate_primitives_lhc(3)  ,   &
                v    => local_coordinate_primitives_lhc(4)  ,   &
                w    => local_coordinate_primitives_lhc(5)  ,   &
                p    => local_coordinate_primitives_lhc(6)  ,   &
                z1   => local_coordinate_primitives_lhc(7)      &
            )
            lhc_density       = rhos(1) * z1 + rhos(2) * (1.d0 - z1)
            lhc_pressure      = p
            lhc_soundspeed    = compute_soundspeed(an_eos, p, rhos, [z1, 1.0_real_kind - z1])
            lhc_main_velocity = u
        end associate
        associate(                                               &
                rhos    => local_coordinate_primitives_rhc(1:2), &
                u       => local_coordinate_primitives_rhc(3)  , &
                v       => local_coordinate_primitives_rhc(4)  , &
                w       => local_coordinate_primitives_rhc(5)  , &
                p       => local_coordinate_primitives_rhc(6)  , &
                z1      => local_coordinate_primitives_rhc(7)    &
            )
            rhc_density       = rhos(1) * z1 + rhos(2) * (1.d0 - z1)
            rhc_pressure      = p
            rhc_soundspeed    = compute_soundspeed(an_eos, p, rhos, [z1, 1.0_real_kind - z1])
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
        flux(1:7, 1) = (-1.d0 / lhc_cell_volume) * nonviscosity_flux(:) * face_area
        flux(1:7, 2) = (+1.d0 / rhc_cell_volume) * nonviscosity_flux(:) * face_area

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
            lhc_alpha1  => primitive_variables_lhc(7)  , &
            rhc_rhos    => primitive_variables_rhc(1:2), &
            rhc_p       => primitive_variables_rhc(6)  , &
            rhc_alpha1  => primitive_variables_rhc(7)    &
        )

            if(ignore_kdivu_)then
                lhc_k = 0.d0
                rhc_k = 0.d0
            else
                lhc_k = compute_k(an_eos, lhc_p, lhc_rhos, lhc_alpha1)
                rhc_k = compute_k(an_eos, rhc_p, rhc_rhos, rhc_alpha1)
            end if

            flux(7, 1) = flux(7, 1) &
                                   + (lhc_alpha1 + lhc_k) * (1.d0 / lhc_cell_volume) * numerical_velocity * face_area
            flux(7, 2) = flux(7, 2) &
                                   - (rhc_alpha1 + rhc_k) * (1.d0 / rhc_cell_volume) * numerical_velocity * face_area
        end associate
    end function flux_function

    function spectral_radius(an_eos, primitive_variables, length) result(r)
        class(eos      ), intent(in) :: an_eos
        real (real_kind), intent(in) :: primitive_variables(:)
        real (real_kind), intent(in) :: length
        real (real_kind) :: r

        associate(                                 &
            densities => primitive_variables(1:2), &
            velocity  => primitive_variables(3:5), &
            pressure  => primitive_variables(6)  , &
            alpha1    => primitive_variables(7)    &
        )
            r = compute_soundspeed(an_eos, pressure, densities, [alpha1, 1.0_real_kind - alpha1]) + vector_magnitude(velocity)
        end associate
    end function spectral_radius
end module five_equation_model_module