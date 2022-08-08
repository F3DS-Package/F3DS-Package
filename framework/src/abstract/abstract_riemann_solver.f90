module abstract_riemann_solver
    implicit none

    private

    type, public, abstract :: riemann_solver
        contains

        procedure(initialize_interface  ), pass(self), deferred :: initialize
        procedure(compute_flux_interface), pass(self), deferred :: compute_flux
    end type

    abstract interface
        subroutine initialize_interface(self, config)
            use abstract_configuration
            import riemann_solver
            class(riemann_solver), intent(inout) :: self
            class(configuration ), intent(inout) :: config
        end subroutine initialize_interface

        pure function compute_flux_interface( &
            self                            , &
            left_conservative               , &
            left_main_velocity              , &
            left_density                    , &
            left_pressure                   , &
            left_soundspeed                 , &
            right_conservative              , &
            right_main_velocity             , &
            right_density                   , &
            right_pressure                  , &
            right_soundspeed                ) result(flux)

            use typedef_module
            import riemann_solver

            class(riemann_solver), intent(in) :: self
            real (real_kind     ), intent(in) :: left_conservative   (:)
            real (real_kind     ), intent(in) :: left_main_velocity
            real (real_kind     ), intent(in) :: left_density
            real (real_kind     ), intent(in) :: left_pressure
            real (real_kind     ), intent(in) :: left_soundspeed
            real (real_kind     ), intent(in) :: right_conservative  (:)
            real (real_kind     ), intent(in) :: right_main_velocity
            real (real_kind     ), intent(in) :: right_density
            real (real_kind     ), intent(in) :: right_pressure
            real (real_kind     ), intent(in) :: right_soundspeed
            real (real_kind     )             :: flux(size(left_conservative))
        end function compute_flux_interface
    end interface
end module abstract_riemann_solver