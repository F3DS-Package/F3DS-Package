module abstract_mixture_eos
    implicit none

    private

    type, public, abstract :: mixture_eos
        contains
        procedure(initialize_interface                      ), pass(self), deferred :: initialize
        procedure(compute_pressure_interface                ), pass(self), deferred :: compute_pressure
        procedure(compute_soundspeed_interface              ), pass(self), deferred :: compute_soundspeed
        procedure(compute_internal_energy_density_interface ), pass(self), deferred :: compute_internal_energy_density
        procedure(compute_wood_soundspeed_interface         ), pass(self), deferred :: compute_wood_soundspeed
        procedure(compute_k_interface                       ), pass(self), deferred :: compute_k
    end type mixture_eos

    abstract interface
        subroutine initialize_interface(self, specific_heat_ratio_fluid1, specific_heat_ratio_fluid2, reference_pressure_fluid1, reference_pressure_fluid2)
            use typedef_module
            import mixture_eos
            class(mixture_eos), intent(inout) :: self
            real (real_kind  ), intent(in)    :: specific_heat_ratio_fluid1
            real (real_kind  ), intent(in)    :: specific_heat_ratio_fluid2
            real (real_kind  ), intent(in)    :: reference_pressure_fluid1
            real (real_kind  ), intent(in)    :: reference_pressure_fluid2
        end subroutine initialize_interface

        pure function compute_pressure_interface(self, specific_internal_energy, density, volume_fruction) result(pressure)
            use typedef_module
            import mixture_eos
            class(mixture_eos), intent(in) :: self
            real (real_kind  ), intent(in) :: specific_internal_energy
            real (real_kind  ), intent(in) :: density
            real (real_kind  ), intent(in) :: volume_fruction
            real (real_kind  )             :: pressure
        end function compute_pressure_interface

        pure function compute_soundspeed_interface(self, pressure, density, volume_fruction) result(soundspeed)
            use typedef_module
            import mixture_eos
            class(mixture_eos), intent(in) :: self
            real (real_kind  ), intent(in) :: pressure
            real (real_kind  ), intent(in) :: density
            real (real_kind  ), intent(in) :: volume_fruction
            real (real_kind  )             :: soundspeed
        end function compute_soundspeed_interface

        pure function compute_internal_energy_density_interface(self, pressure, density, volume_fruction) result(specific_internal_energy)
            use typedef_module
            import mixture_eos
            class(mixture_eos), intent(in) :: self
            real (real_kind  ), intent(in) :: pressure
            real (real_kind  ), intent(in) :: density
            real (real_kind  ), intent(in) :: volume_fruction
            real (real_kind  )             :: specific_internal_energy
        end function compute_internal_energy_density_interface

        pure function compute_wood_soundspeed_interface(self, pressure, density1, density2, volume_fruction) result(soundspeed)
            use typedef_module
            import mixture_eos
            class(mixture_eos), intent(in) :: self
            real (real_kind  ), intent(in) :: pressure
            real (real_kind  ), intent(in) :: density1
            real (real_kind  ), intent(in) :: density2
            real (real_kind  ), intent(in) :: volume_fruction
            real (real_kind  )             :: soundspeed
        end function compute_wood_soundspeed_interface

        pure function compute_k_interface(self, pressure, density1, density2, volume_fruction) result(k)
            use typedef_module
            import mixture_eos
            class(mixture_eos), intent(in) :: self
            real (real_kind  ), intent(in) :: pressure
            real (real_kind  ), intent(in) :: density1
            real (real_kind  ), intent(in) :: density2
            real (real_kind  ), intent(in) :: volume_fruction
            real (real_kind  )             :: k
        end function compute_k_interface
    end interface
end module abstract_mixture_eos