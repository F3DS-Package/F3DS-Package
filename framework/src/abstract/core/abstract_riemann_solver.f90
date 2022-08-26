module abstract_riemann_solver
    implicit none

    private

    type, public, abstract :: riemann_solver
        contains

        procedure(initialize_interface                  ), pass(self), deferred :: initialize
        procedure(compute_features_interface            ), pass(self), deferred :: compute_features
        procedure(compute_mass_flux_interface           ), pass(self), deferred :: compute_mass_flux
        procedure(compute_momentum_flux_interface       ), pass(self), deferred :: compute_momentum_flux
        procedure(compute_energy_flux_interface         ), pass(self), deferred :: compute_energy_flux
        procedure(compute_volume_fraction_flux_interface), pass(self), deferred :: compute_volume_fraction_flux
        procedure(compute_numerical_velocity_interface  ), pass(self), deferred :: compute_numerical_velocity
        procedure(compute_interface_value_interface     ), pass(self), deferred :: compute_interface_value
    end type

    abstract interface
        subroutine initialize_interface(self, config)
            use abstract_configuration
            import riemann_solver
            class(riemann_solver), intent(inout) :: self
            class(configuration ), intent(inout) :: config
        end subroutine initialize_interface

        pure function compute_features_interface(  &
            self                                 , &
            left_main_velocity                   , &
            left_density                         , &
            left_pressure                        , &
            left_soundspeed                      , &
            right_main_velocity                  , &
            right_density                        , &
            right_pressure                       , &
            right_soundspeed                         ) result(features)

            use typedef_module
            import riemann_solver

            class(riemann_solver), intent(in)  :: self
            real (real_kind     ), intent(in)  :: left_main_velocity
            real (real_kind     ), intent(in)  :: left_density
            real (real_kind     ), intent(in)  :: left_pressure
            real (real_kind     ), intent(in)  :: left_soundspeed
            real (real_kind     ), intent(in)  :: right_main_velocity
            real (real_kind     ), intent(in)  :: right_density
            real (real_kind     ), intent(in)  :: right_pressure
            real (real_kind     ), intent(in)  :: right_soundspeed
            real (real_kind     ), allocatable :: features(:)
        end function compute_features_interface

        pure function compute_mass_flux_interface( &
            self                                 , &
            left_conservative_mass               , &
            left_main_velocity                   , &
            left_density                         , &
            left_pressure                        , &
            left_soundspeed                      , &
            right_conservative_mass              , &
            right_main_velocity                  , &
            right_density                        , &
            right_pressure                       , &
            right_soundspeed                     , &
            features                                 ) result(mass_flux)

            use typedef_module
            import riemann_solver

            class(riemann_solver), intent(in) :: self
            real (real_kind     ), intent(in) :: left_conservative_mass
            real (real_kind     ), intent(in) :: left_main_velocity
            real (real_kind     ), intent(in) :: left_density
            real (real_kind     ), intent(in) :: left_pressure
            real (real_kind     ), intent(in) :: left_soundspeed
            real (real_kind     ), intent(in) :: right_conservative_mass
            real (real_kind     ), intent(in) :: right_main_velocity
            real (real_kind     ), intent(in) :: right_density
            real (real_kind     ), intent(in) :: right_pressure
            real (real_kind     ), intent(in) :: right_soundspeed
            real (real_kind     ), intent(in) :: features(:)
            real (real_kind     )             :: mass_flux
        end function compute_mass_flux_interface

        pure function compute_momentum_flux_interface( &
            self                                 , &
            left_conservative_momentum           , &
            left_main_velocity                   , &
            left_density                         , &
            left_pressure                        , &
            left_soundspeed                      , &
            right_conservative_momentum          , &
            right_main_velocity                  , &
            right_density                        , &
            right_pressure                       , &
            right_soundspeed                     , &
            features                                 ) result(momentum_flux)

            use typedef_module
            import riemann_solver

            class(riemann_solver), intent(in) :: self
            real (real_kind     ), intent(in) :: left_conservative_momentum(3)
            real (real_kind     ), intent(in) :: left_main_velocity
            real (real_kind     ), intent(in) :: left_density
            real (real_kind     ), intent(in) :: left_pressure
            real (real_kind     ), intent(in) :: left_soundspeed
            real (real_kind     ), intent(in) :: right_conservative_momentum(3)
            real (real_kind     ), intent(in) :: right_main_velocity
            real (real_kind     ), intent(in) :: right_density
            real (real_kind     ), intent(in) :: right_pressure
            real (real_kind     ), intent(in) :: right_soundspeed
            real (real_kind     ), intent(in) :: features(:)
            real (real_kind     )             :: momentum_flux(3)
        end function compute_momentum_flux_interface

        pure function compute_energy_flux_interface( &
            self                                   , &
            left_conservative_energy               , &
            left_main_velocity                     , &
            left_density                           , &
            left_pressure                          , &
            left_soundspeed                        , &
            right_conservative_energy              , &
            right_main_velocity                    , &
            right_density                          , &
            right_pressure                         , &
            right_soundspeed                       , &
            features                                   ) result(energy_flux)

            use typedef_module
            import riemann_solver

            class(riemann_solver), intent(in) :: self
            real (real_kind     ), intent(in) :: left_conservative_energy
            real (real_kind     ), intent(in) :: left_main_velocity
            real (real_kind     ), intent(in) :: left_density
            real (real_kind     ), intent(in) :: left_pressure
            real (real_kind     ), intent(in) :: left_soundspeed
            real (real_kind     ), intent(in) :: right_conservative_energy
            real (real_kind     ), intent(in) :: right_main_velocity
            real (real_kind     ), intent(in) :: right_density
            real (real_kind     ), intent(in) :: right_pressure
            real (real_kind     ), intent(in) :: right_soundspeed
            real (real_kind     ), intent(in) :: features(:)
            real (real_kind     )             :: energy_flux
        end function compute_energy_flux_interface

        pure function compute_volume_fraction_flux_interface( &
            self                                            , &
            left_conservative_volume_fraction               , &
            left_main_velocity                              , &
            left_density                                    , &
            left_pressure                                   , &
            left_soundspeed                                 , &
            right_conservative_volume_fraction              , &
            right_main_velocity                             , &
            right_density                                   , &
            right_pressure                                  , &
            right_soundspeed                                , &
            features                                            ) result(volume_fraction_flux)

            use typedef_module
            import riemann_solver

            class(riemann_solver), intent(in) :: self
            real (real_kind     ), intent(in) :: left_conservative_volume_fraction
            real (real_kind     ), intent(in) :: left_main_velocity
            real (real_kind     ), intent(in) :: left_density
            real (real_kind     ), intent(in) :: left_pressure
            real (real_kind     ), intent(in) :: left_soundspeed
            real (real_kind     ), intent(in) :: right_conservative_volume_fraction
            real (real_kind     ), intent(in) :: right_main_velocity
            real (real_kind     ), intent(in) :: right_density
            real (real_kind     ), intent(in) :: right_pressure
            real (real_kind     ), intent(in) :: right_soundspeed
            real (real_kind     ), intent(in) :: features(:)
            real (real_kind     )             :: volume_fraction_flux
        end function compute_volume_fraction_flux_interface

        pure function compute_numerical_velocity_interface(  &
            self                                           , &
            left_main_velocity                             , &
            left_density                                   , &
            left_pressure                                  , &
            left_soundspeed                                , &
            right_main_velocity                            , &
            right_density                                  , &
            right_pressure                                 , &
            right_soundspeed                               , &
            features                                            ) result(numerical_velocity)

            use typedef_module
            import riemann_solver

            class(riemann_solver), intent(in)  :: self
            real (real_kind     ), intent(in)  :: left_main_velocity
            real (real_kind     ), intent(in)  :: left_density
            real (real_kind     ), intent(in)  :: left_pressure
            real (real_kind     ), intent(in)  :: left_soundspeed
            real (real_kind     ), intent(in)  :: right_main_velocity
            real (real_kind     ), intent(in)  :: right_density
            real (real_kind     ), intent(in)  :: right_pressure
            real (real_kind     ), intent(in)  :: right_soundspeed
            real (real_kind     ), intent(in)  :: features(:)
            real (real_kind     )              :: numerical_velocity
        end function compute_numerical_velocity_interface

        pure function compute_interface_value_interface(  &
            self                                        , &
            left_value                                  , &
            left_main_velocity                          , &
            left_density                                , &
            left_pressure                               , &
            left_soundspeed                             , &
            right_value                                 , &
            right_main_velocity                         , &
            right_density                               , &
            right_pressure                              , &
            right_soundspeed                            , &
            features                                         ) result(interface_value)

            use typedef_module
            import riemann_solver

            class(riemann_solver), intent(in)  :: self
            real (real_kind     ), intent(in)  :: left_value
            real (real_kind     ), intent(in)  :: left_main_velocity
            real (real_kind     ), intent(in)  :: left_density
            real (real_kind     ), intent(in)  :: left_pressure
            real (real_kind     ), intent(in)  :: left_soundspeed
            real (real_kind     ), intent(in)  :: right_value
            real (real_kind     ), intent(in)  :: right_main_velocity
            real (real_kind     ), intent(in)  :: right_density
            real (real_kind     ), intent(in)  :: right_pressure
            real (real_kind     ), intent(in)  :: right_soundspeed
            real (real_kind     ), intent(in)  :: features(:)
            real (real_kind     )              :: interface_value
        end function compute_interface_value_interface
    end interface
end module abstract_riemann_solver