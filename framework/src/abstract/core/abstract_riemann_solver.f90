module abstract_riemann_solver
    ! """
    ! This module defines an abstract type `riemann_solver` and an abstract interface for various methods related to solving Riemann
    ! problems.
    ! The `riemann_solver` type contains the following deferred procedures:
    ! - `initialize`: Initializes the Riemann solver with a given configuration.
    ! - `compute_features`: Computes the features of the Riemann problem, such as main velocities, densities, pressures, and
    ! soundspeeds.
    ! - `compute_mass_flux`: Computes the mass flux across the interface of the Riemann problem.
    ! - `compute_momentum_flux`: Computes the momentum flux across the interface of the Riemann problem.
    ! - `compute_energy_flux`: Computes the energy flux across the interface of the Riemann problem.
    ! - `compute_volume_fraction_flux`: Computes the volume fraction flux across the interface of the Riemann problem.
    ! - `compute_numerical_velocity`: Computes the numerical velocity across the interface of the Riemann problem.
    ! - `compute_interface_value`: Computes the interface value across the interface of the Riemann problem.
    ! The abstract interface defines the following methods:
    ! - `initialize_interface`: Initializes the Riemann solver with a given configuration.
    ! - `compute_features_interface`: Computes the features of the Riemann problem.
    ! - `compute_mass_flux_interface`: Computes the mass flux across the interface of the Riemann problem.
    ! - `compute_momentum_flux_interface`: Computes the momentum flux across the interface of the Riemann problem.
    ! - `compute_energy_flux_interface`: Computes the energy flux across the interface of the Riemann problem.
    ! - `compute_volume_fraction_flux_interface`: Computes the volume fraction flux across the interface of the Riemann problem.
    ! - `compute_numerical_velocity_interface`: Computes the numerical velocity across the interface of the Riemann problem.
    ! - `compute_interface_value_interface`: Computes the interface value across the interface of the Riemann problem.
    ! """
    implicit none

    private

    type, public, abstract :: riemann_solver
        ! -----------------------------------------------------------------------
        ! Type: riemann_solver
        ! -----------------------------------------------------------------------
        ! This is an abstract type representing a Riemann solver. It contains
        ! deferred procedures that need to be implemented by derived types.
        !
        ! Public Components:
        ! ------------------
        ! - initialize: Initializes the Riemann solver.
        ! - compute_features: Computes the features of the interface.
        ! - compute_mass_flux: Computes the mass flux across the interface.
        ! - compute_momentum_flux: Computes the momentum flux across the interface.
        ! - compute_energy_flux: Computes the energy flux across the interface.
        ! - compute_volume_fraction_flux: Computes the volume fraction flux across the interface.
        ! - compute_numerical_velocity: Computes the numerical velocity across the interface.
        ! - compute_interface_value: Computes the interface value.
        !
        ! -----------------------------------------------------------------------
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
            ! > Initializes the interface between the Riemann solver and the configuration.
            ! This subroutine sets up the necessary connections between the Riemann solver
            ! and the configuration object.
            !
            ! @param self The Riemann solver object.
            ! @param config The configuration object.
            !
            ! @return None
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
            ! """
            ! Compute features using the Riemann solver.
            ! Parameters:
            !     self (riemann_solver): An instance of the Riemann solver class.
            !     left_main_velocity (real): The main velocity on the left side.
            !     left_density (real): The density on the left side.
            !     left_pressure (real): The pressure on the left side.
            !     left_soundspeed (real): The soundspeed on the left side.
            !     right_main_velocity (real): The main velocity on the right side.
            !     right_density (real): The density on the right side.
            !     right_pressure (real): The pressure on the right side.
            !     right_soundspeed (real): The soundspeed on the right side.
            ! Returns:
            !     features (real): An array of computed features.
            ! This function computes features using the Riemann solver. It takes in the necessary parameters and returns an array of computed
            ! features.
            ! """

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
            ! """
            ! Compute the mass flux across the interface using the Riemann solver.
            ! Parameters:
            !     self (riemann_solver): An instance of the Riemann solver class.
            !     left_conservative_mass (real): The conservative mass on the left side of the interface.
            !     left_main_velocity (real): The main velocity on the left side of the interface.
            !     left_density (real): The density on the left side of the interface.
            !     left_pressure (real): The pressure on the left side of the interface.
            !     left_soundspeed (real): The soundspeed on the left side of the interface.
            !     right_conservative_mass (real): The conservative mass on the right side of the interface.
            !     right_main_velocity (real): The main velocity on the right side of the interface.
            !     right_density (real): The density on the right side of the interface.
            !     right_pressure (real): The pressure on the right side of the interface.
            !     right_soundspeed (real): The soundspeed on the right side of the interface.
            !     features (real array): An array of additional features.
            ! Returns:
            !     mass_flux (real): The computed mass flux across the interface.
            ! """

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
            ! """
            ! Compute the momentum flux at the interface between two regions.
            ! Parameters:
            !     self (riemann_solver): The Riemann solver object.
            !     left_conservative_momentum (array): Array of size 3 containing the conservative momentum on the left side of the interface.
            !     left_main_velocity (float): The main velocity on the left side of the interface.
            !     left_density (float): The density on the left side of the interface.
            !     left_pressure (float): The pressure on the left side of the interface.
            !     left_soundspeed (float): The soundspeed on the left side of the interface.
            !     right_conservative_momentum (array): Array of size 3 containing the conservative momentum on the right side of the
            !     interface.
            !     right_main_velocity (float): The main velocity on the right side of the interface.
            !     right_density (float): The density on the right side of the interface.
            !     right_pressure (float): The pressure on the right side of the interface.
            !     right_soundspeed (float): The soundspeed on the right side of the interface.
            !     features (array): Array containing additional features.
            ! Returns:
            !     momentum_flux (array): Array of size 3 containing the computed momentum flux at the interface.
            ! """

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
            ! """
            ! Compute the energy flux across the interface using the Riemann solver.
            ! Parameters:
            !     self (riemann_solver): An instance of the Riemann solver class.
            !     left_conservative_energy (real): The conservative energy on the left side of the interface.
            !     left_main_velocity (real): The main velocity on the left side of the interface.
            !     left_density (real): The density on the left side of the interface.
            !     left_pressure (real): The pressure on the left side of the interface.
            !     left_soundspeed (real): The soundspeed on the left side of the interface.
            !     right_conservative_energy (real): The conservative energy on the right side of the interface.
            !     right_main_velocity (real): The main velocity on the right side of the interface.
            !     right_density (real): The density on the right side of the interface.
            !     right_pressure (real): The pressure on the right side of the interface.
            !     right_soundspeed (real): The soundspeed on the right side of the interface.
            !     features (real array): An array of additional features.
            ! Returns:
            !     energy_flux (real): The computed energy flux across the interface.
            ! """

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
            ! """
            ! Compute the volume fraction flux at the interface using the Riemann solver.
            ! Parameters:
            !     self (riemann_solver): An instance of the Riemann solver class.
            !     left_conservative_volume_fraction (real): The conservative volume fraction on the left side of the interface.
            !     left_main_velocity (real): The main velocity on the left side of the interface.
            !     left_density (real): The density on the left side of the interface.
            !     left_pressure (real): The pressure on the left side of the interface.
            !     left_soundspeed (real): The soundspeed on the left side of the interface.
            !     right_conservative_volume_fraction (real): The conservative volume fraction on the right side of the interface.
            !     right_main_velocity (real): The main velocity on the right side of the interface.
            !     right_density (real): The density on the right side of the interface.
            !     right_pressure (real): The pressure on the right side of the interface.
            !     right_soundspeed (real): The soundspeed on the right side of the interface.
            !     features (real array): An array of additional features.
            !
            ! Returns:
            !     volume_fraction_flux (real): The computed volume fraction flux at the interface.
            ! """

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
            ! """
            ! Compute the numerical velocity at the interface using the given Riemann solver.
            ! Parameters:
            !     self (riemann_solver): An instance of the Riemann solver class.
            !     left_main_velocity (real): The main velocity on the left side of the interface.
            !     left_density (real): The density on the left side of the interface.
            !     left_pressure (real): The pressure on the left side of the interface.
            !     left_soundspeed (real): The soundspeed on the left side of the interface.
            !     right_main_velocity (real): The main velocity on the right side of the interface.
            !     right_density (real): The density on the right side of the interface.
            !     right_pressure (real): The pressure on the right side of the interface.
            !     right_soundspeed (real): The soundspeed on the right side of the interface.
            !     features (real array): An array of additional features.
            ! Returns:
            !     numerical_velocity (real): The computed numerical velocity at the interface.
            ! """

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
            ! """
            ! Compute the interface value using the Riemann solver.
            ! Parameters:
            !     self (riemann_solver): An instance of the Riemann solver class.
            !     left_value (real): The value on the left side of the interface.
            !     left_main_velocity (real): The main velocity on the left side of the interface.
            !     left_density (real): The density on the left side of the interface.
            !     left_pressure (real): The pressure on the left side of the interface.
            !     left_soundspeed (real): The soundspeed on the left side of the interface.
            !     right_value (real): The value on the right side of the interface.
            !     right_main_velocity (real): The main velocity on the right side of the interface.
            !     right_density (real): The density on the right side of the interface.
            !     right_pressure (real): The pressure on the right side of the interface.
            !     right_soundspeed (real): The soundspeed on the right side of the interface.
            !     features (real array): An array of additional features.
            ! Returns:
            !     interface_value (real): The computed interface value.
            ! """

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