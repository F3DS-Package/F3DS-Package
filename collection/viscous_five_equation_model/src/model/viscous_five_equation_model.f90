module viscous_five_equation_model_module
    use vector_module
    use typedef_module
    use stdio_module
    use string_utils_module
    use abstract_configuration
    use abstract_eos
    use abstract_riemann_solver
    use viscous_five_equation_model_utils_module
    use matrix_module
    use vector_module
    use math_constant_module

    implicit none

    private

    real(real_kind), parameter :: interface_threshold_ = 1e-8

    integer(int_kind )              :: num_phase_

    ! for Kdiv(u) term
    logical                         :: ignore_kdivu_

    ! for surface tension term
    logical                         :: ignore_surface_tension_
    real   (real_kind), allocatable :: surface_tension_(:)
    logical                         :: ignore_garrick_modification_
    real   (real_kind)              :: garrick_modification_coef_

    ! for gravity term
    real   (real_kind)              :: gravitational_acceleration_(3)

    ! for viscosity term (include thermal conduction)
    logical                         :: ignore_viscosity_
    real   (real_kind), allocatable :: dynamic_viscosity_(:)

    ! Average density is used to select a heaviest fluid.
    real   (real_kind), allocatable :: average_densities_(:)


    public :: initialize_model
    public :: compute_residual_element
    public :: compute_source_term
    public :: spectral_radius
    public :: compute_smoothed_volume_fraction
    public :: normalize_gradient_volume_fraction
    public :: curvature_smoothing_weight
    public :: curvature_preprocessing
    public :: apply_surface_tension

    contains

    pure function compute_pressure_jump(alpha1, alpha_heavy, interface_curvature) result(dp)
        real(real_kind), intent(in) :: alpha1, alpha_heavy, interface_curvature
        real(real_kind)             :: dp
        dp = mixture_surface_tension(alpha1) * interface_curvature * alpha_heavy * garrick_modification_coef_
    end function

    pure function apply_surface_tension() result(ret)
        logical :: ret
        if(ignore_surface_tension_)then
            ret = .false.
        else
            ret = .true.
        end if
    end function apply_surface_tension

    pure function apply_viscosity() result(ret)
        logical :: ret
        if(ignore_viscosity_)then
            ret = .false.
        else
            ret = .true.
        end if
    end function apply_viscosity

    pure function mixture_surface_tension(z1) result(sigma)
        real   (real_kind), intent(in) :: z1
        real   (real_kind)             :: sigma
        sigma = z1 * surface_tension_(1) + (1.d0 - z1) * surface_tension_(2)
    end function

    pure function mixture_dynamic_viscosity(z1) result(mu)
        real   (real_kind), intent(in) :: z1
        real   (real_kind)             :: mu
        mu = z1 * dynamic_viscosity_(1) + (1.d0 - z1) * dynamic_viscosity_(2)
    end function

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

    pure function get_heavest_volume_fraction(alpha1) result(alpha_heavy)
        real   (real_kind), intent(in) :: alpha1
        real   (real_kind)             :: alpha_heavy
        alpha_heavy = 0.5d0 * (1.d0 + sign(1.d0, average_densities_(1) - average_densities_(2))) * alpha1 + 0.5d0 * (1.d0 + sign(1.d0, average_densities_(2) - average_densities_(1))) * (1.d0 - alpha1)
    end function get_heavest_volume_fraction

    pure function fix_reconstructed_primitive_variables(reconstructed_primitive_variables) result(fixed_reconstructed_primitive_variables)
        real   (real_kind), intent(in) :: reconstructed_primitive_variables(8)
        real   (real_kind)             :: fixed_reconstructed_primitive_variables(8)

        associate(                                         &
            rho1  => reconstructed_primitive_variables(1), &
            rho2  => reconstructed_primitive_variables(2), &
            u     => reconstructed_primitive_variables(3), &
            v     => reconstructed_primitive_variables(4), &
            w     => reconstructed_primitive_variables(5), &
            p     => reconstructed_primitive_variables(6), &
            alpha => reconstructed_primitive_variables(7), &
            kappa => reconstructed_primitive_variables(8)  &
        )
            fixed_reconstructed_primitive_variables(1) = max(rho1, 0.d0)
            fixed_reconstructed_primitive_variables(2) = max(rho2, 0.d0)
            fixed_reconstructed_primitive_variables(3) = u
            fixed_reconstructed_primitive_variables(4) = v
            fixed_reconstructed_primitive_variables(5) = w
            fixed_reconstructed_primitive_variables(6) = p
            fixed_reconstructed_primitive_variables(7) = max(0.d0, min(alpha, 1.d0))
            fixed_reconstructed_primitive_variables(8) = kappa
        end associate
    end function fix_reconstructed_primitive_variables

    subroutine initialize_model(a_configuration, primitive_variables_set, num_cells)
        class  (configuration), intent(inout) :: a_configuration
        real   (real_kind    ), intent(in   ) :: primitive_variables_set(:,:)
        integer(int_kind     ), intent(in   ) :: num_cells

        logical           :: found
        logical           :: ignore_gravity = .false.
        integer(int_kind) :: i

        real   (real_kind), allocatable :: num_sums(:)

        call a_configuration%get_int    ("Phase.Number of phase", num_phase_, found)
        if(.not. found) call call_error("'Phase.Number of phase' is not found in configuration you set. Please check your configuration file.")

        if(.not. (num_phase_ == 2)) call call_error("This solver approve ONLY two phase flow.")

        allocate(surface_tension_  (num_phase_))
        allocate(dynamic_viscosity_(num_phase_))
        allocate(average_densities_(num_phase_))
        allocate(num_sums          (num_phase_))

        do i = 1, num_phase_, 1
            call a_configuration%get_real      ("Phase.Phase property "//to_str(i)//".Surface tension", surface_tension_(i), found, 0.d0)
            if(.not. found) call write_warring("'Phase.Phase property "//to_str(i)//".Surface tension' is not found in the configuration you set. Set to 0.")

            call a_configuration%get_real      ("Phase.Phase property "//to_str(i)//".Dynamic viscosity", dynamic_viscosity_(i), found, 0.d0)
            if(.not. found) call write_warring("'Phase.Phase property "//to_str(i)//".Dynamic viscosity' is not found in the configuration you set. Set to 0.")
        end do

        call a_configuration%get_bool      ("Model.Ignore surface tension", ignore_surface_tension_, found, .false.)
        if(.not. found) call write_warring("'Model.Ignore surface tension' is not found in configuration you set. Apply this term.")

        call a_configuration%get_bool      ("Model.Ignore Garrick's modification", ignore_garrick_modification_, found, .false.)
        if(.not. found) call write_warring("'Model.Ignore Garrick's modification' is not found in configuration you set. Apply this modification.")
        if(ignore_garrick_modification_)then
            garrick_modification_coef_ = 0.d0
        else
            garrick_modification_coef_ = 1.d0
        endif

        call a_configuration%get_bool      ("Model.Ignore viscosity", ignore_viscosity_, found, .false.)
        if(.not. found) call write_warring("'Model.Ignore viscosity' is not found in configuration you set. Apply this term.")

        call a_configuration%get_bool      ("Model.Ignore kdivu", ignore_kdivu_, found, .false.)
        if(.not. found) call write_warring("'Model.Ignore kdivu' is not found in configuration you set. Apply this term.")

        call a_configuration%get_real("Model.Gravitational acceleration.x", gravitational_acceleration_(1), found, 0.d0)
        if(.not. found) then
            call write_warring("'Model.Gravitational acceleration.x' is not found in configuration you set. Ignore gravity term.")
            ignore_gravity=.true.
        endif
        call a_configuration%get_real("Model.Gravitational acceleration.y", gravitational_acceleration_(2), found, 0.d0)
        if(.not. found) then
            call write_warring("'Model.Gravitational acceleration.y' is not found in configuration you set. Ignore gravity term.")
            ignore_gravity=.true.
        endif
        call a_configuration%get_real("Model.Gravitational acceleration.z", gravitational_acceleration_(3), found, 0.d0)
        if(.not. found) then
            call write_warring("'Model.Gravitational acceleration.z' is not found in configuration you set. Ignore gravity term.")
            ignore_gravity=.true.
        endif
        if (ignore_gravity) gravitational_acceleration_(:) = 0.d0

        ! compute {@code average_densities_}
        average_densities_(:) = 0.d0
        num_sums          (:) = 0.d0
        do i = 1, num_cells, 1
            associate(                                            &
                density1         => primitive_variables_set(1,i), &
                density2         => primitive_variables_set(2,i), &
                volume_fraction1 => primitive_variables_set(7,i)  &
            )
                if(volume_fraction1 > 1.d0 - machine_epsilon)then
                    average_densities_(1) = average_densities_(1) + density1
                    num_sums          (1) = num_sums(1) + 1.d0
                else if (volume_fraction1 < machine_epsilon)then
                    average_densities_(2) = average_densities_(2) + density2
                    num_sums          (2) = num_sums(2) + 1.d0
                end if
            end associate
        end do
        do i = 1, num_phase_, 1
            if(num_sums(i) > 0.d0)then
                average_densities_(i) = average_densities_(i) / num_sums(i)
            end if
        end do

        !sanity check
        if(ignore_surface_tension_)then
            surface_tension_(:) = 0.d0
        end if
        if(ignore_viscosity_)then
            dynamic_viscosity_(:) = 0.d0
        end if
    end subroutine initialize_model

    function compute_residual_element(            &
        an_eos                                  , &
        an_riemann_solver                       , &
        primitive_variables_lhc                 , &
        primitive_variables_rhc                 , &
        reconstructed_primitive_variables_lhc   , &
        reconstructed_primitive_variables_rhc   , &
        face_gradient_primitive_variables       , &
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
        real   (real_kind     ), intent(in) :: face_gradient_primitive_variables    (:)
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
        real(real_kind) :: lhc_soundspeed, lhc_pressure, lhc_density, lhc_main_velocity, lhc_pressure_jump
        real(real_kind) :: rhc_soundspeed, rhc_pressure, rhc_density, rhc_main_velocity, rhc_pressure_jump

        ! ## rieman solver
        real(real_kind) :: nonviscosity_flux      (num_conservative_values)
        real(real_kind) :: rieman_solver_features(7)

        ! ## face variables
        real(real_kind) :: numerical_velocity, interface_volume_fraction

        ! ## kapila Kdiv(u)
        real(real_kind) :: rhc_k, lhc_k

        ! ## surface tension
        real(real_kind) :: interface_alpha_heavy, rhc_alpha_heavy, lhc_alpha_heavy, interface_curvature

        ! ## viscosity variables
        real(real_kind) :: tau (3,3)
        real(real_kind) :: beta(3)
        real(real_kind) :: viscosity_flux(num_conservative_values)

        ! # compute primitive-variables of face local coordinate
        ! ## left-side
        local_coordinate_primitives_lhc(:) = rotate_primitive(                            &
            fix_reconstructed_primitive_variables(reconstructed_primitive_variables_lhc), &
            face_normal_vector                              , &
            face_tangential1_vector                         , &
            face_tangential2_vector                         , &
            num_primitive_values                              &
        )
        ! ## right-side
        local_coordinate_primitives_rhc(:) = rotate_primitive(                            &
            fix_reconstructed_primitive_variables(reconstructed_primitive_variables_rhc), &
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

        interface_curvature = 0.5d0 * (primitive_variables_lhc(8) + primitive_variables_rhc(8))

        ! # compute EoS, main velosity, and fluxs
        associate(                                                 &
                rho1    => local_coordinate_primitives_lhc(1),     &
                rho2    => local_coordinate_primitives_lhc(2),     &
                u       => local_coordinate_primitives_lhc(3),     &
                v       => local_coordinate_primitives_lhc(4),     &
                w       => local_coordinate_primitives_lhc(5),     &
                p       => local_coordinate_primitives_lhc(6),     &
                z1      => local_coordinate_primitives_lhc(7),     &
                curv    => local_coordinate_primitives_lhc(8)      &
            )
            lhc_density    = rho1 * z1 + rho2 * (1.d0 - z1)
            lhc_pressure   = p
            lhc_soundspeed = an_eos%compute_soundspeed(p, lhc_density, z1)
            lhc_main_velocity = u
            lhc_alpha_heavy   = get_heavest_volume_fraction(z1)
            lhc_pressure_jump = compute_pressure_jump(z1, lhc_alpha_heavy, interface_curvature)
        end associate
        associate(                                             &
                rho1    => local_coordinate_primitives_rhc(1), &
                rho2    => local_coordinate_primitives_rhc(2), &
                u       => local_coordinate_primitives_rhc(3), &
                v       => local_coordinate_primitives_rhc(4), &
                w       => local_coordinate_primitives_rhc(5), &
                p       => local_coordinate_primitives_rhc(6), &
                z1      => local_coordinate_primitives_rhc(7), &
                curv    => local_coordinate_primitives_rhc(8)  &
            )
            rhc_density    = rho1 * z1 + rho2 * (1.d0 - z1)
            rhc_pressure   = p
            rhc_soundspeed = an_eos%compute_soundspeed(p, rhc_density, z1)
            rhc_main_velocity = u
            rhc_alpha_heavy   = get_heavest_volume_fraction(z1)
            rhc_pressure_jump = compute_pressure_jump(z1, rhc_alpha_heavy, interface_curvature)
        end associate

        ! # compute flux
        ! ## local coordinate flux
        rieman_solver_features(:) = an_riemann_solver%compute_features( &
            lhc_main_velocity               , &
            lhc_density                     , &
            lhc_pressure - lhc_pressure_jump, &
            lhc_soundspeed                  , &
            rhc_main_velocity               , &
            rhc_density                     , &
            rhc_pressure - rhc_pressure_jump, &
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
            if(ignore_kdivu_)then
                lhc_k = 0.d0
                rhc_k = 0.d0
            else
                lhc_k = an_eos%compute_k(lhc_p, lhc_rhos, lhc_z1)
                rhc_k = an_eos%compute_k(rhc_p, rhc_rhos, rhc_z1)
            end if
            residual_element(7, 1) = residual_element(7, 1) &
                                   + (lhc_z1 + lhc_k) * (1.d0 / lhc_cell_volume) * numerical_velocity * face_area
            residual_element(7, 2) = residual_element(7, 2) &
                                   - (rhc_z1 + rhc_k) * (1.d0 / rhc_cell_volume) * numerical_velocity * face_area
        end associate

        ! # surface tension term
        if(apply_surface_tension())then
            associate(                                               &
                lhc_rho1  => primitive_variables_lhc(1),             &
                lhc_rho2  => primitive_variables_lhc(2),             &
                lhc_z1    => primitive_variables_lhc(7),             &
                lhc_curv  => primitive_variables_lhc(8),             &
                rhc_rho1  => primitive_variables_rhc(1),             &
                rhc_rho2  => primitive_variables_rhc(2),             &
                rhc_z1    => primitive_variables_rhc(7),             &
                rhc_curv  => primitive_variables_rhc(8)              &
            )
                interface_alpha_heavy = get_heavest_volume_fraction(lhc_z1)
                lhc_alpha_heavy       = get_heavest_volume_fraction(lhc_z1)
                rhc_alpha_heavy       = get_heavest_volume_fraction(rhc_z1)
                residual_element(3:5, 1) = residual_element(3:5, 1) &
                                       + mixture_surface_tension(lhc_z1) * lhc_curv * (1.d0 / lhc_cell_volume) * interface_alpha_heavy * face_area * face_normal_vector(1:3)
                residual_element(3:5, 2) = residual_element(3:5, 2) &
                                       - mixture_surface_tension(rhc_z1) * rhc_curv * (1.d0 / rhc_cell_volume) * interface_alpha_heavy * face_area * face_normal_vector(1:3)
                residual_element(6, 1) = residual_element(6, 1) &
                                       + mixture_surface_tension(lhc_z1) * lhc_curv * (1.d0 / lhc_cell_volume) * (interface_alpha_heavy * numerical_velocity - lhc_alpha_heavy * numerical_velocity) * face_area
                residual_element(6, 2) = residual_element(6, 2) &
                                       - mixture_surface_tension(rhc_z1) * rhc_curv * (1.d0 / rhc_cell_volume) * (interface_alpha_heavy * numerical_velocity - rhc_alpha_heavy * numerical_velocity) * face_area
            end associate
        end if

        ! # viscosity term (Stokes hypothesis)
        if(apply_viscosity())then
            associate(                                                                                                         &
                dudx => face_gradient_primitive_variables(7 )                                                                , &
                dudy => face_gradient_primitive_variables(8 )                                                                , &
                dudz => face_gradient_primitive_variables(9 )                                                                , &
                dvdx => face_gradient_primitive_variables(10)                                                                , &
                dvdy => face_gradient_primitive_variables(11)                                                                , &
                dvdz => face_gradient_primitive_variables(12)                                                                , &
                dwdx => face_gradient_primitive_variables(13)                                                                , &
                dwdy => face_gradient_primitive_variables(14)                                                                , &
                dwdz => face_gradient_primitive_variables(15)                                                                , &
                u    =>                                0.5d0 * (primitive_variables_lhc(3:5) + primitive_variables_rhc(3:5)) , &
                mu   => mixture_dynamic_viscosity(0.5d0 * (primitive_variables_lhc(7  ) + primitive_variables_rhc(7  ))), &
                n    => face_normal_vector                                                                                     &
            )
                tau (1,1) = 2.d0 * mu * dudx - (2.d0 / 3.d0) * mu * (dudx + dvdy + dwdz)
                tau (2,2) = 2.d0 * mu * dvdy - (2.d0 / 3.d0) * mu * (dudx + dvdy + dwdz)
                tau (3,3) = 2.d0 * mu * dwdz - (2.d0 / 3.d0) * mu * (dudx + dvdy + dwdz)
                tau (1,2) = mu * (dudy + dvdx)
                tau (1,3) = mu * (dudz + dwdx)
                tau (2,3) = mu * (dwdy + dvdz)
                tau (2,1) = tau(1,2)
                tau (3,1) = tau(1,3)
                tau (3,2) = tau(2,3)
                beta(1)   = vector_multiply(tau(1,:), u) ! not inculde heat diffusion !
                beta(2)   = vector_multiply(tau(2,:), u)
                beta(3)   = vector_multiply(tau(3,:), u)
                viscosity_flux(1:2) = 0.d0
                viscosity_flux(3:5) = matrix_multiply(tau , n)
                viscosity_flux(6  ) = vector_multiply(beta, n)
                viscosity_flux(7  ) = 0.d0
                viscosity_flux(3:5) = vector_unrotate( &
                    viscosity_flux(3:5)              , &
                    face_normal_vector               , &
                    face_tangential1_vector          , &
                    face_tangential2_vector            &
                )
                residual_element(1:7, 1) = residual_element(1:7, 1) + (1.d0 / lhc_cell_volume) * viscosity_flux(1:7) * face_area
                residual_element(1:7, 2) = residual_element(1:7, 2) - (1.d0 / rhc_cell_volume) * viscosity_flux(1:7) * face_area
            end associate
        end if
    end function compute_residual_element

    function compute_source_term(variables, num_conservative_values) result(source)
        real   (real_kind), intent(in) :: variables(:)
        integer(int_kind ), intent(in) :: num_conservative_values
        real   (real_kind)             :: source(num_conservative_values)

        associate(                                               &
            density  => compute_mixture_density(variables), &
            g        => gravitational_acceleration_(:)    , &
            velocity => variables(3:5)                           &
        )
            source(1:2) = 0.d0
            source(3:5) = density * g
            source(6  ) = vector_multiply(density * g, velocity)
            source(7  ) = 0.d0
        end associate
    end function compute_source_term

    function spectral_radius(an_eos, primitive_variables, length) result(r)
        class(eos      ), intent(in) :: an_eos
        real (real_kind), intent(in) :: primitive_variables(:)
        real (real_kind), intent(in) :: length
        real (real_kind) :: r

        associate(                                                         &
            density  => compute_mixture_density(primitive_variables), &
            velocity => primitive_variables(3:5)                         , &
            pressure => primitive_variables(6)                           , &
            alpha1   => primitive_variables(7)                             &
        )
            r = an_eos%compute_soundspeed(pressure, density, alpha1) + vector_magnitude(velocity) + mixture_dynamic_viscosity(alpha1) / (density * length)
        end associate
    end function spectral_radius

    function compute_smoothed_volume_fraction(surface_tension_variables, primitive_variables, num_surface_tension_variables) result(dst_surface_tension_variables)
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
            heavest_volume_fraction = get_heavest_volume_fraction(volume_fraction)
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
            if(machine_epsilon < mag)then
                ! Towerd a heavest fluid direction.
                dst_surface_tension_variables(2:4) = grad_alpha / mag
            else
                dst_surface_tension_variables(2:4) = 0.d0
            endif
        end associate
    end function normalize_gradient_volume_fraction

    pure function curvature_smoothing_weight(own_cell_position, neighbor_cell_position, own_variables, neighbor_variables) result(weight)
        real(real_kind ), intent(in) :: own_cell_position     (3)
        real(real_kind ), intent(in) :: neighbor_cell_position(3)
        real(real_kind ), intent(in) :: own_variables         (:)
        real(real_kind ), intent(in) :: neighbor_variables    (:)
        real(real_kind )             :: weight

        real(real_kind), parameter :: smoothing_power = 2.d0 ! range is {@code smoothing_power} > 0.

        associate(alpha => neighbor_variables(5))
            weight = (alpha * (1.d0 - alpha))**smoothing_power
        end associate
    end function curvature_smoothing_weight

    pure function curvature_preprocessing(primitives, surface_tension_variables, num_variables) result(dst_primitives)
        real   (real_kind ), intent(in) :: surface_tension_variables(:), primitives(:)
        integer(int_kind  ), intent(in) :: num_variables
        real   (real_kind )             :: dst_primitives(num_variables)

        dst_primitives(1:7) = primitives(1:7)
        associate(                                 &
            rho1  => primitives(1)               , &
            rho2  => primitives(2)               , &
            z     => primitives(7)               , &
            kappa => surface_tension_variables(5)  &
        )
            if((interface_threshold_ < z) .and. (z < 1.d0 - interface_threshold_))then
                dst_primitives(8) = -kappa
            else
                dst_primitives(8) = 0.d0
            endif
        end associate
    end function curvature_preprocessing
end module viscous_five_equation_model_module