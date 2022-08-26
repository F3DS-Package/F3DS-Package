module class_hllc
    use typedef_module
    use abstract_riemann_solver
    use abstract_configuration

    implicit none

    private

    type, public, extends(riemann_solver) :: hllc
        contains

        private

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: compute_features
        procedure, public, pass(self) :: compute_mass_flux
        procedure, public, pass(self) :: compute_momentum_flux
        procedure, public, pass(self) :: compute_energy_flux
        procedure, public, pass(self) :: compute_volume_fraction_flux
        procedure, public, pass(self) :: compute_numerical_velocity
        procedure, public, pass(self) :: compute_interface_value
    end type

    contains

    subroutine initialize(self, config)
        class(hllc         ), intent(inout) :: self
        class(configuration), intent(inout) :: config
    end subroutine initialize

    pure function compute_features(            &
        self                                 , &
        left_main_velocity                   , &
        left_density                         , &
        left_pressure                        , &
        left_soundspeed                      , &
        right_main_velocity                  , &
        right_density                        , &
        right_pressure                       , &
        right_soundspeed                         ) result(features)

        class(hllc          ), intent(in   ) :: self
        real (real_kind     ), intent(in   ) :: left_main_velocity
        real (real_kind     ), intent(in   ) :: left_density
        real (real_kind     ), intent(in   ) :: left_pressure
        real (real_kind     ), intent(in   ) :: left_soundspeed
        real (real_kind     ), intent(in   ) :: right_main_velocity
        real (real_kind     ), intent(in   ) :: right_density
        real (real_kind     ), intent(in   ) :: right_pressure
        real (real_kind     ), intent(in   ) :: right_soundspeed
        real (real_kind     ), allocatable   :: features(:)

        real (real_kind     )                :: ave_vel, ave_c

        allocate(features(7)) ! 1: s_left, 2: s_right, 3: s_muinus, 4: s_puls, 5: s_mid, 6: c_star_left, 7: c_star_right

        associate(                                 &
            lhc_u        => left_main_velocity   , &
            lhc_rho      => left_density         , &
            lhc_p        => left_pressure        , &
            lhc_c        => left_soundspeed      , &
            rhc_u        => right_main_velocity  , &
            rhc_rho      => right_density        , &
            rhc_p        => right_pressure       , &
            rhc_c        => right_soundspeed     , &
            s_left       => features(1)          , &
            s_right      => features(2)          , &
            s_muinus     => features(3)          , &
            s_puls       => features(4)          , &
            s_mid        => features(5)          , &
            c_star_left  => features(6)          , &
            c_star_right => features(7)            &
        )
            ave_vel      = 0.5d0 * (lhc_u + rhc_u)
            ave_c        = 0.5d0 * (lhc_c + rhc_c)
            s_left       = min(ave_vel - ave_c, lhc_u - lhc_c)
            s_right      = max(ave_vel + ave_c, rhc_u + rhc_c)
            s_muinus     = min(0.d0, s_left )
            s_puls       = max(0.d0, s_right)
            s_mid        = (rhc_p - lhc_p + lhc_rho * lhc_u * (s_left  - lhc_u)        &
                                          - rhc_rho * rhc_u * (s_right - rhc_u))       &
                         / (lhc_rho * (s_left  - lhc_u) - rhc_rho * (s_right - rhc_u))
            c_star_left  = (s_left  - lhc_u) / (s_left  - s_mid)
            c_star_right = (s_right - rhc_u) / (s_right - s_mid)
        end associate
    end function compute_features

    pure function compute_mass_flux( &
        self                       , &
        left_conservative_mass     , &
        left_main_velocity         , &
        left_density               , &
        left_pressure              , &
        left_soundspeed            , &
        right_conservative_mass    , &
        right_main_velocity        , &
        right_density              , &
        right_pressure             , &
        right_soundspeed           , &
        features                       ) result(mass_flux)

        class(hllc          ), intent(in   ) :: self
        real (real_kind     ), intent(in   ) :: left_conservative_mass
        real (real_kind     ), intent(in   ) :: left_main_velocity
        real (real_kind     ), intent(in   ) :: left_density
        real (real_kind     ), intent(in   ) :: left_pressure
        real (real_kind     ), intent(in   ) :: left_soundspeed
        real (real_kind     ), intent(in   ) :: right_conservative_mass
        real (real_kind     ), intent(in   ) :: right_main_velocity
        real (real_kind     ), intent(in   ) :: right_density
        real (real_kind     ), intent(in   ) :: right_pressure
        real (real_kind     ), intent(in   ) :: right_soundspeed
        real (real_kind     ), intent(in   ) :: features(:)
        real (real_kind     )                :: mass_flux

        real (real_kind) :: f_left, f_right, q_star_left, q_star_right

        associate(                                 &
            lhc_u        => left_main_velocity   , &
            lhc_rho      => left_density         , &
            lhc_p        => left_pressure        , &
            lhc_c        => left_soundspeed      , &
            rhc_u        => right_main_velocity  , &
            rhc_rho      => right_density        , &
            rhc_p        => right_pressure       , &
            rhc_c        => right_soundspeed     , &
            s_left       => features(1)          , &
            s_right      => features(2)          , &
            s_muinus     => features(3)          , &
            s_puls       => features(4)          , &
            s_mid        => features(5)          , &
            c_star_left  => features(6)          , &
            c_star_right => features(7)            &
        )
            f_left       = left_conservative_mass  * lhc_u
            f_right      = right_conservative_mass * rhc_u
            q_star_left  = left_conservative_mass  * c_star_left
            q_star_right = right_conservative_mass * c_star_right
            mass_flux = 0.5d0 * (1.d0 + sign(1.d0, s_mid))                           &
                      * (f_left + s_muinus * (q_star_left - left_conservative_mass)) &
                      + 0.5d0 * (1.d0 - sign(1.d0, s_mid))                           &
                      * (f_right + s_puls * (q_star_right - right_conservative_mass))
        end associate
    end function compute_mass_flux

    pure function compute_momentum_flux( &
        self                           , &
        left_conservative_momentum     , &
        left_main_velocity             , &
        left_density                   , &
        left_pressure                  , &
        left_soundspeed                , &
        right_conservative_momentum    , &
        right_main_velocity            , &
        right_density                  , &
        right_pressure                 , &
        right_soundspeed               , &
        features                           ) result(momentum_flux)

        class(hllc          ), intent(in   ) :: self
        real (real_kind     ), intent(in   ) :: left_conservative_momentum(3)
        real (real_kind     ), intent(in   ) :: left_main_velocity
        real (real_kind     ), intent(in   ) :: left_density
        real (real_kind     ), intent(in   ) :: left_pressure
        real (real_kind     ), intent(in   ) :: left_soundspeed
        real (real_kind     ), intent(in   ) :: right_conservative_momentum(3)
        real (real_kind     ), intent(in   ) :: right_main_velocity
        real (real_kind     ), intent(in   ) :: right_density
        real (real_kind     ), intent(in   ) :: right_pressure
        real (real_kind     ), intent(in   ) :: right_soundspeed
        real (real_kind     ), intent(in   ) :: features(:)
        real (real_kind     )                :: momentum_flux(3)

        real (real_kind) :: f_left(3), f_right(3), q_star_left(3), q_star_right(3)

        associate(                                 &
            lhc_u        => left_main_velocity   , &
            lhc_rho      => left_density         , &
            lhc_p        => left_pressure        , &
            lhc_c        => left_soundspeed      , &
            rhc_u        => right_main_velocity  , &
            rhc_rho      => right_density        , &
            rhc_p        => right_pressure       , &
            rhc_c        => right_soundspeed     , &
            s_left       => features(1)          , &
            s_right      => features(2)          , &
            s_muinus     => features(3)          , &
            s_puls       => features(4)          , &
            s_mid        => features(5)          , &
            c_star_left  => features(6)          , &
            c_star_right => features(7)            &
        )
            f_left      (1) = left_conservative_momentum (1) * lhc_u + lhc_p
            f_left      (2) = left_conservative_momentum (2) * lhc_u
            f_left      (3) = left_conservative_momentum (3) * lhc_u
            f_right     (1) = right_conservative_momentum(1) * rhc_u + rhc_p
            f_right     (2) = right_conservative_momentum(2) * rhc_u
            f_right     (3) = right_conservative_momentum(3) * rhc_u
            q_star_left (1) = c_star_left  * lhc_rho * s_mid
            q_star_left (2) = c_star_left  * left_conservative_momentum(2)
            q_star_left (3) = c_star_left  * left_conservative_momentum(3)
            q_star_right(1) = c_star_right * rhc_rho * s_mid
            q_star_right(2) = c_star_right * right_conservative_momentum(2)
            q_star_right(3) = c_star_right * right_conservative_momentum(3)
            momentum_flux(:) = 0.5d0 * (1.d0 + sign(1.d0, s_mid))                                        &
                             * (f_left(:) + s_muinus * (q_star_left(:) - left_conservative_momentum(:))) &
                             + 0.5d0 * (1.d0 - sign(1.d0, s_mid))                                        &
                             * (f_right(:) + s_puls * (q_star_right(:) - right_conservative_momentum(:)))
        end associate
    end function compute_momentum_flux

    pure function compute_energy_flux( &
        self                         , &
        left_conservative_energy     , &
        left_main_velocity           , &
        left_density                 , &
        left_pressure                , &
        left_soundspeed              , &
        right_conservative_energy    , &
        right_main_velocity          , &
        right_density                , &
        right_pressure               , &
        right_soundspeed             , &
        features                         ) result(energy_flux)

        class(hllc          ), intent(in   ) :: self
        real (real_kind     ), intent(in   ) :: left_conservative_energy
        real (real_kind     ), intent(in   ) :: left_main_velocity
        real (real_kind     ), intent(in   ) :: left_density
        real (real_kind     ), intent(in   ) :: left_pressure
        real (real_kind     ), intent(in   ) :: left_soundspeed
        real (real_kind     ), intent(in   ) :: right_conservative_energy
        real (real_kind     ), intent(in   ) :: right_main_velocity
        real (real_kind     ), intent(in   ) :: right_density
        real (real_kind     ), intent(in   ) :: right_pressure
        real (real_kind     ), intent(in   ) :: right_soundspeed
        real (real_kind     ), intent(in   ) :: features(:)
        real (real_kind     )                :: energy_flux

        real (real_kind) :: f_left, f_right, q_star_left, q_star_right

        associate(                                 &
            lhc_u        => left_main_velocity   , &
            lhc_rho      => left_density         , &
            lhc_p        => left_pressure        , &
            lhc_c        => left_soundspeed      , &
            rhc_u        => right_main_velocity  , &
            rhc_rho      => right_density        , &
            rhc_p        => right_pressure       , &
            rhc_c        => right_soundspeed     , &
            s_left       => features(1)          , &
            s_right      => features(2)          , &
            s_muinus     => features(3)          , &
            s_puls       => features(4)          , &
            s_mid        => features(5)          , &
            c_star_left  => features(6)          , &
            c_star_right => features(7)            &
        )
            f_left       = left_conservative_energy  * lhc_u
            f_right      = right_conservative_energy * rhc_u
            q_star_left  = c_star_left  * (left_conservative_energy  + (s_mid - lhc_u) * (lhc_rho * s_mid + lhc_p / (s_left  - lhc_u)))
            q_star_right = c_star_right * (right_conservative_energy + (s_mid - rhc_u) * (rhc_rho * s_mid + rhc_p / (s_right - rhc_u)))
            energy_flux = 0.5d0 * (1.d0 + sign(1.d0, s_mid))                           &
                      * (f_left + s_muinus * (q_star_left - left_conservative_energy)) &
                      + 0.5d0 * (1.d0 - sign(1.d0, s_mid))                           &
                      * (f_right + s_puls * (q_star_right - right_conservative_energy))
        end associate
    end function compute_energy_flux

    pure function compute_volume_fraction_flux( &
        self                                  , &
        left_conservative_volume_fraction     , &
        left_main_velocity                    , &
        left_density                          , &
        left_pressure                         , &
        left_soundspeed                       , &
        right_conservative_volume_fraction    , &
        right_main_velocity                   , &
        right_density                         , &
        right_pressure                        , &
        right_soundspeed                      , &
        features                                  ) result(volume_fraction_flux)

        class(hllc          ), intent(in   ) :: self
        real (real_kind     ), intent(in   ) :: left_conservative_volume_fraction
        real (real_kind     ), intent(in   ) :: left_main_velocity
        real (real_kind     ), intent(in   ) :: left_density
        real (real_kind     ), intent(in   ) :: left_pressure
        real (real_kind     ), intent(in   ) :: left_soundspeed
        real (real_kind     ), intent(in   ) :: right_conservative_volume_fraction
        real (real_kind     ), intent(in   ) :: right_main_velocity
        real (real_kind     ), intent(in   ) :: right_density
        real (real_kind     ), intent(in   ) :: right_pressure
        real (real_kind     ), intent(in   ) :: right_soundspeed
        real (real_kind     ), intent(in   ) :: features(:)
        real (real_kind     )                :: volume_fraction_flux

        real (real_kind) :: f_left, f_right, q_star_left, q_star_right

        associate(                                 &
            lhc_u        => left_main_velocity   , &
            lhc_rho      => left_density         , &
            lhc_p        => left_pressure        , &
            lhc_c        => left_soundspeed      , &
            rhc_u        => right_main_velocity  , &
            rhc_rho      => right_density        , &
            rhc_p        => right_pressure       , &
            rhc_c        => right_soundspeed     , &
            s_left       => features(1)          , &
            s_right      => features(2)          , &
            s_muinus     => features(3)          , &
            s_puls       => features(4)          , &
            s_mid        => features(5)          , &
            c_star_left  => features(6)          , &
            c_star_right => features(7)            &
        )
            f_left       = left_conservative_volume_fraction  * lhc_u
            f_right      = right_conservative_volume_fraction * rhc_u
            q_star_left  = c_star_left  * left_conservative_volume_fraction
            q_star_right = c_star_right * right_conservative_volume_fraction
            volume_fraction_flux = 0.5d0 * (1.d0 + sign(1.d0, s_mid))                           &
                      * (f_left + s_muinus * (q_star_left - left_conservative_volume_fraction)) &
                      + 0.5d0 * (1.d0 - sign(1.d0, s_mid))                                      &
                      * (f_right + s_puls * (q_star_right - right_conservative_volume_fraction))
        end associate
    end function compute_volume_fraction_flux

    pure function compute_numerical_velocity(  &
        self                                 , &
        left_main_velocity                   , &
        left_density                         , &
        left_pressure                        , &
        left_soundspeed                      , &
        right_main_velocity                  , &
        right_density                        , &
        right_pressure                       , &
        right_soundspeed                     , &
        features                                  ) result(numerical_velocity)

        class(hllc     ), intent(in)  :: self
        real (real_kind), intent(in)  :: left_main_velocity
        real (real_kind), intent(in)  :: left_density
        real (real_kind), intent(in)  :: left_pressure
        real (real_kind), intent(in)  :: left_soundspeed
        real (real_kind), intent(in)  :: right_main_velocity
        real (real_kind), intent(in)  :: right_density
        real (real_kind), intent(in)  :: right_pressure
        real (real_kind), intent(in)  :: right_soundspeed
        real (real_kind), intent(in)  :: features(:)
        real (real_kind)              :: numerical_velocity

        associate(                                 &
            lhc_u        => left_main_velocity   , &
            lhc_rho      => left_density         , &
            lhc_p        => left_pressure        , &
            lhc_c        => left_soundspeed      , &
            rhc_u        => right_main_velocity  , &
            rhc_rho      => right_density        , &
            rhc_p        => right_pressure       , &
            rhc_c        => right_soundspeed     , &
            s_left       => features(1)          , &
            s_right      => features(2)          , &
            s_muinus     => features(3)          , &
            s_puls       => features(4)          , &
            s_mid        => features(5)          , &
            c_star_left  => features(6)          , &
            c_star_right => features(7)            &
        )
            numerical_velocity &
                = 0.5d0 * (1.d0 + sign(1.d0, s_mid)) * (lhc_u + s_muinus * (c_star_left  - 1.d0)) &
                + 0.5d0 * (1.d0 - sign(1.d0, s_mid)) * (rhc_u + s_puls   * (c_star_right - 1.d0))
        end associate
    end function compute_numerical_velocity

    pure function compute_interface_value(  &
        self                              , &
        left_value                        , &
        left_main_velocity                , &
        left_density                      , &
        left_pressure                     , &
        left_soundspeed                   , &
        right_value                       , &
        right_main_velocity               , &
        right_density                     , &
        right_pressure                    , &
        right_soundspeed                  , &
        features                               ) result(interface_value)
        class(hllc     ), intent(in)  :: self
        real (real_kind), intent(in)  :: left_value
        real (real_kind), intent(in)  :: left_main_velocity
        real (real_kind), intent(in)  :: left_density
        real (real_kind), intent(in)  :: left_pressure
        real (real_kind), intent(in)  :: left_soundspeed
        real (real_kind), intent(in)  :: right_value
        real (real_kind), intent(in)  :: right_main_velocity
        real (real_kind), intent(in)  :: right_density
        real (real_kind), intent(in)  :: right_pressure
        real (real_kind), intent(in)  :: right_soundspeed
        real (real_kind), intent(in)  :: features(:)
        real (real_kind)              :: interface_value

        associate(                                 &
            lhc_u        => left_main_velocity   , &
            lhc_rho      => left_density         , &
            lhc_p        => left_pressure        , &
            lhc_c        => left_soundspeed      , &
            rhc_u        => right_main_velocity  , &
            rhc_rho      => right_density        , &
            rhc_p        => right_pressure       , &
            rhc_c        => right_soundspeed     , &
            s_left       => features(1)          , &
            s_right      => features(2)          , &
            s_muinus     => features(3)          , &
            s_puls       => features(4)          , &
            s_mid        => features(5)          , &
            c_star_left  => features(6)          , &
            c_star_right => features(7)            &
        )
            interface_value &
                = 0.5d0 * (1.d0 + sign(1.d0, s_mid)) * left_value &
                + 0.5d0 * (1.d0 - sign(1.d0, s_mid)) * right_value
        end associate
    end function compute_interface_value
end module class_hllc