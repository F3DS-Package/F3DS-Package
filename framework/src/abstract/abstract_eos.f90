module abstract_eos
    ! TODO: Support to simulation of more than 2 phases flow.
    implicit none

    private

    type, public, abstract :: eos
        contains
        procedure(initialize_interface                      ), pass(self), deferred :: initialize
        procedure(compute_pressure_interface                ), pass(self), deferred :: compute_pressure
        procedure(compute_internal_energy_density_interface ), pass(self), deferred :: compute_internal_energy_density
        procedure(compute_soundspeed_interface              ), pass(self), deferred :: compute_soundspeed
        !procedure(compute_mixture_soundspeed_interface      ), pass(self), deferred :: compute_mixture_soundspeed
        !procedure(compute_frozen_soundspeed_interface       ), pass(self), deferred :: compute_frozen_soundspeed
        procedure(compute_wood_soundspeed_interface         ), pass(self), deferred :: compute_wood_soundspeed
        procedure(compute_k_interface                       ), pass(self), deferred :: compute_k
    end type eos

    abstract interface
        subroutine initialize_interface(self, a_configuration)
            use abstract_configuration
            import eos
            class(eos          ), intent(inout) :: self
            class(configuration), intent(in   ) :: a_configuration
        end subroutine initialize_interface

        pure function compute_pressure_interface(self, specific_internal_energy, density, volume_fraction) result(pressure)
            use typedef_module
            import eos
            class(eos), intent(in) :: self
            real (real_kind  ), intent(in) :: specific_internal_energy
            real (real_kind  ), intent(in) :: density
            real (real_kind  ), intent(in) :: volume_fraction
            real (real_kind  )             :: pressure
        end function compute_pressure_interface

        pure function compute_internal_energy_density_interface(self, pressure, density, volume_fraction) result(specific_internal_energy)
            use typedef_module
            import eos
            class(eos), intent(in) :: self
            real (real_kind  ), intent(in) :: pressure
            real (real_kind  ), intent(in) :: density
            real (real_kind  ), intent(in) :: volume_fraction
            real (real_kind  )             :: specific_internal_energy
        end function compute_internal_energy_density_interface

        pure function compute_mixture_soundspeed_interface(self, pressure, density, volume_fraction) result(soundspeed)
            use typedef_module
            import eos
            class(eos), intent(in) :: self
            real (real_kind  ), intent(in) :: pressure
            real (real_kind  ), intent(in) :: density
            real (real_kind  ), intent(in) :: volume_fraction
            real (real_kind  )             :: soundspeed
        end function compute_mixture_soundspeed_interface

        pure function compute_soundspeed_interface(self, pressure, density, volume_fraction) result(soundspeed)
            use typedef_module
            import eos
            class(eos), intent(in) :: self
            real (real_kind  ), intent(in) :: pressure
            real (real_kind  ), intent(in) :: density
            real (real_kind  ), intent(in) :: volume_fraction
            real (real_kind  )             :: soundspeed
        end function compute_soundspeed_interface

        pure function compute_frozen_soundspeed_interface(self, pressure, densities, volume_fraction) result(soundspeed)
            use typedef_module
            import eos
            class(eos), intent(in) :: self
            real (real_kind  ), intent(in) :: pressure
            real (real_kind  ), intent(in) :: densities(:)
            real (real_kind  ), intent(in) :: volume_fraction
            real (real_kind  )             :: soundspeed
        end function compute_frozen_soundspeed_interface

        pure function compute_wood_soundspeed_interface(self, pressure, densities, volume_fraction) result(soundspeed)
            use typedef_module
            import eos
            class(eos), intent(in) :: self
            real (real_kind  ), intent(in) :: pressure
            real (real_kind  ), intent(in) :: densities(:)
            real (real_kind  ), intent(in) :: volume_fraction
            real (real_kind  )             :: soundspeed
        end function compute_wood_soundspeed_interface

        pure function compute_k_interface(self, pressure, densities, volume_fraction) result(k)
            use typedef_module
            import eos
            class(eos), intent(in) :: self
            real (real_kind  ), intent(in) :: pressure
            real (real_kind  ), intent(in) :: densities(:)
            real (real_kind  ), intent(in) :: volume_fraction
            real (real_kind  )             :: k
        end function compute_k_interface
    end interface
end module abstract_eos