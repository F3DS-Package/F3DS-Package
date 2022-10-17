module class_five_equation_model_space_discretization
    use vector_module
    use typedef_module
    use stdio_module
    use string_utils_module
    use abstract_configuration
    use abstract_eos
    use abstract_riemann_solver
    use abstract_model
    use five_equation_model_variables_module
    use matrix_module
    use vector_module

    implicit none

    private

    type, public, extends(model) :: five_equation_model_space_discretization
        private

        integer(int_kind )              :: num_phase_

        ! for Kdiv(u) term
        logical                         :: ignore_kdivu_

        ! for surface tension term
        real   (real_kind), allocatable :: surface_tension_(:)

        ! for gravity term
        real   (real_kind)              :: gravitational_acceleration_(3)

        ! for viscosity term (include thermal conduction)
        real   (real_kind), allocatable :: dynamic_viscosity_(:)

        contains

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: compute_residual_element
        procedure, public, pass(self) :: compute_source_term
        procedure, public, pass(self) :: spectral_radius

        procedure :: mixture_surface_tension
        procedure :: mixture_dynamic_viscosity
        procedure :: compute_mixture_density
    end type five_equation_model_space_discretization


    contains

    pure function mixture_surface_tension(self, z1) result(sigma)
        class  (five_equation_model_space_discretization), intent(in) :: self
        real   (real_kind), intent(in) :: z1
        real   (real_kind)             :: sigma
        sigma = z1 * self%surface_tension_(1) + (1.d0 - z1) * self%surface_tension_(2)
    end function

    pure function mixture_dynamic_viscosity(self, z1) result(mu)
        class  (five_equation_model_space_discretization), intent(in) :: self
        real   (real_kind), intent(in) :: z1
        real   (real_kind)             :: mu
        mu = z1 * self%dynamic_viscosity_(1) + (1.d0 - z1) * self%dynamic_viscosity_(2)
    end function

    pure function compute_mixture_density(self, primitive_variables) result(density)
        class(five_equation_model_space_discretization), intent(in) :: self
        real (real_kind                               ), intent(in) :: primitive_variables(:)
        real (real_kind                               )             :: density
        associate(                            &
            rho1   => primitive_variables(1), &
            rho2   => primitive_variables(2), &
            alpha1 => primitive_variables(7)  &
        )
            density = alpha1 * rho1 + (1.d0 - alpha1) * rho2
        end associate
    end function

    subroutine initialize(self, a_configuration)
        class  (five_equation_model_space_discretization), intent(inout) :: self
        class  (configuration                           ), intent(inout) :: a_configuration

        logical           :: found
        integer(int_kind) :: i

        call a_configuration%get_int    ("Phase.Number of phase", self%num_phase_, found)
        if(.not. found) call call_error("'Phase.Number of phase' is not found in configuration you set. Please check your configuration file.")

        allocate(self%surface_tension_  (self%num_phase_))
        allocate(self%dynamic_viscosity_(self%num_phase_))

        do i = 1, self%num_phase_, 1
            call a_configuration%get_real      ("Phase.Phase property "//to_str(i)//".Surface tension", self%surface_tension_(i), found, 0.d0)
            if(.not. found) call write_warring("'Phase.Phase property "//to_str(i)//".Surface tension' is not found in the configuration you set. Set to 0.")

            call a_configuration%get_real      ("Phase.Phase property "//to_str(i)//".Dynamic viscosity", self%dynamic_viscosity_(i), found, 0.d0)
            if(.not. found) call write_warring("'Phase.Phase property "//to_str(i)//".Dynamic viscosity' is not found in the configuration you set. Set to 0.")
        end do

        call a_configuration%get_bool      ("Model.Ignore kdivu", self%ignore_kdivu_, found, .false.)
        if(.not. found) call write_warring("'Model.Ignore kdivu' is not found in configuration you set. Apply this term.")

        call a_configuration%get_real_array("Model.Gravitational acceleration", self%gravitational_acceleration_, found, [0.d0, 0.d0, 0.d0])
        if(.not. found) call write_warring("'Model.Gravitational acceleration' is not found in configuration you set. Ignore gravity term.")
    end subroutine initialize

    pure function compute_residual_element(       &
        self                                    , &
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

        class  (five_equation_model_space_discretization), intent(in) :: self

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

        ! ## viscosity variables
        real(real_kind) :: tau (3,3)
        real(real_kind) :: beta(3)
        real(real_kind) :: viscosity_flux(num_conservative_values)

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
                z1      => local_coordinate_primitives_lhc(7),     &
                curv    => local_coordinate_primitives_lhc(8)      &
            )
            lhc_density    = rho1 * z1 + rho2 * (1.d0 - z1)
            lhc_pressure   = p
            lhc_soundspeed = an_eos%compute_soundspeed(p, lhc_density, z1)
            lhc_main_velocity = u
            lhc_pressure_jump = self%mixture_surface_tension(z1) * curv * (1.d0 - z1)
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
            rhc_pressure_jump = self%mixture_surface_tension(z1) * curv * (1.d0 - z1)
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
            if(self%ignore_kdivu_)then
                lhc_k = 0.d0
                rhc_k = 0.d0
            else
                lhc_k = an_eos%compute_k(lhc_p, lhc_rhos, lhc_z1)
                rhc_k = an_eos%compute_k(rhc_p, rhc_rhos, rhc_z1)
            end if
            residual_element(7, 1) = residual_element(7, 1) &
                                   - (-lhc_z1 - lhc_k) * (1.d0 / lhc_cell_volume) * numerical_velocity * face_area
            residual_element(7, 2) = residual_element(7, 2) &
                                   + (-rhc_z1 - rhc_k) * (1.d0 / rhc_cell_volume) * numerical_velocity * face_area
        end associate

        ! # surface tension term
        associate(                                               &
            lhc_z1    => primitive_variables_lhc(7),             &
            lhc_curv  => primitive_variables_lhc(8),             &
            rhc_z1    => primitive_variables_rhc(7),             &
            rhc_curv  => primitive_variables_rhc(8),             &
            alpha_l   => (1.d0 - interface_volume_fraction)      & ! In this solver, primary volume fraction z1 is a gas volume fraction! thus we need 1.d0 - {@code interface_volume_fraction}.
        )
            residual_element(3:5, 1) = residual_element(3:5, 1) &
                                   + self%mixture_surface_tension(lhc_z1) * lhc_curv * (1.d0 / lhc_cell_volume) * alpha_l * face_area * face_normal_vector(1:3)
            residual_element(3:5, 2) = residual_element(3:5, 2) &
                                   - self%mixture_surface_tension(rhc_z1) * rhc_curv * (1.d0 / rhc_cell_volume) * alpha_l * face_area * face_normal_vector(1:3)
            ! numerical velocity is defined in the direction toward the left cell
            residual_element(6, 1) = residual_element(6, 1) &
                                   - self%mixture_surface_tension(lhc_z1) * lhc_curv * (1.d0 / lhc_cell_volume) * (alpha_l * numerical_velocity - numerical_velocity) * face_area
            residual_element(6, 2) = residual_element(6, 2) &
                                   + self%mixture_surface_tension(rhc_z1) * rhc_curv * (1.d0 / rhc_cell_volume) * (alpha_l * numerical_velocity - numerical_velocity) * face_area
        end associate

        ! # viscosity term (Stokes hypothesis)
        associate(                                                                         &
            dudx => face_gradient_primitive_variables(7 )                                , &
            dudy => face_gradient_primitive_variables(8 )                                , &
            dudz => face_gradient_primitive_variables(9 )                                , &
            dvdx => face_gradient_primitive_variables(10)                                , &
            dvdy => face_gradient_primitive_variables(11)                                , &
            dvdz => face_gradient_primitive_variables(12)                                , &
            dwdx => face_gradient_primitive_variables(13)                                , &
            dwdy => face_gradient_primitive_variables(14)                                , &
            dwdz => face_gradient_primitive_variables(15)                                , &
            u    => 0.5d0 * (primitive_variables_lhc(3:5) + primitive_variables_rhc(3:5)), &
            mu   => self%mixture_dynamic_viscosity(interface_volume_fraction)            , &
            n    => face_normal_vector                                                     &
        )
            tau(1,1) = 2.d0 * mu * dudx - (2.d0 / 3.d0) * mu * (dudx + dvdy + dwdz)
            tau(2,2) = 2.d0 * mu * dvdy - (2.d0 / 3.d0) * mu * (dudx + dvdy + dwdz)
            tau(3,3) = 2.d0 * mu * dwdz - (2.d0 / 3.d0) * mu * (dudx + dvdy + dwdz)
            tau(1,2) = mu * (dudy + dvdx)
            tau(1,3) = mu * (dudz + dwdx)
            tau(2,3) = mu * (dwdy + dvdz)
            tau(2,1) = tau(1,2)
            tau(3,1) = tau(1,3)
            tau(3,2) = tau(2,3)
            beta(1) = vector_multiply(tau(1,:), u) ! not inculde heat diffusion !
            beta(2) = vector_multiply(tau(2,:), u)
            beta(3) = vector_multiply(tau(3,:), u)

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

            residual_element(1:7, 1) = residual_element(1:7, 1) - (1.d0 / lhc_cell_volume) * viscosity_flux(:) * face_area
            residual_element(1:7, 2) = residual_element(1:7, 2) + (1.d0 / rhc_cell_volume) * viscosity_flux(:) * face_area
        end associate
    end function compute_residual_element

    pure function compute_source_term(self, primitive_variables, gradient_primitive_variables, cell_volume, num_primitive_values) result(source)
        class  (five_equation_model_space_discretization), intent(in) :: self
        real   (real_kind), intent(in) :: primitive_variables          (:)
        real   (real_kind), intent(in) :: gradient_primitive_variables (:)
        real   (real_kind), intent(in) :: cell_volume
        integer(int_kind ), intent(in) :: num_primitive_values
        real   (real_kind)             :: source(num_primitive_values)

        associate(                                                         &
            density  => self%compute_mixture_density(primitive_variables), &
            g        => self%gravitational_acceleration_(:)              , &
            velocity => primitive_variables(3:5)                           &
        )
            source(3:5) = density * g
            source(6  ) = vector_multiply(density * g, velocity)
        end associate
    end function compute_source_term

    pure function spectral_radius(self, an_eos, primitive_variables, length) result(r)
        class  (five_equation_model_space_discretization), intent(in) :: self
        class(eos      ), intent(in) :: an_eos
        real (real_kind), intent(in) :: primitive_variables(:)
        real (real_kind), intent(in) :: length
        real (real_kind) :: r

        associate(                                                         &
            density  => self%compute_mixture_density(primitive_variables), &
            velocity => primitive_variables(3:5)                         , &
            pressure => primitive_variables(6)                           , &
            alpha1   => primitive_variables(7)                             &
        )
            r = an_eos%compute_soundspeed(pressure, density, alpha1) + vector_magnitude(velocity) + self%mixture_dynamic_viscosity(alpha1) / (density * length)
        end associate
    end function spectral_radius
end module class_five_equation_model_space_discretization