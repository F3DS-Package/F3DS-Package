module class_stiffened_gas_eos
    use typedef_module
    use abstract_eos
    use abstract_configuration
    use math_constant_module
    use stdio_module
    use string_utils_module

    implicit none

    private

    type, public, extends(eos) :: stiffened_gas_eos
        private

        integer(int_kind )              :: num_phase_
        real   (real_kind), allocatable :: specific_heat_ratio_(:)
        real   (real_kind), allocatable :: reference_pressure_ (:)

        contains

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: compute_pressure
        procedure, public, pass(self) :: compute_internal_energy_density
        procedure, public, pass(self) :: compute_soundspeed
        !procedure, public, pass(self) :: compute_mixture_soundspeed
        !procedure, public, pass(self) :: compute_frozen_soundspeed
        procedure, public, pass(self) :: compute_wood_soundspeed
        procedure, public, pass(self) :: compute_k

        procedure, pass(self) :: compute_mixture_specific_heat_ratio
        procedure, pass(self) :: compute_mixture_reference_pressure
        procedure, pass(self) :: compute_phi
    end type stiffened_gas_eos

    contains

    subroutine initialize(self, a_configuration)
        class(stiffened_gas_eos), intent(inout) :: self
        class(configuration    ), intent(inout) :: a_configuration

        logical           :: found
        integer(int_kind) :: i

        call a_configuration%get_int    ("Phase.Number of Phase", self%num_phase_, found)
        if(.not. found) call call_error("'Phase.Number of Phase' is not found in configuration you set. Please check your configuration file.")

        if(self%num_phase_ > 2) call call_error("This solver does not support simulation of more than 2 phases flow.")

        allocate(self%specific_heat_ratio_(self%num_phase_))
        allocate(self%reference_pressure_ (self%num_phase_))

        do i = 1, self%num_phase_, 1
            call a_configuration%get_real   ("Phase.Phase property "//to_str(i)//".Specific heat ratio", self%specific_heat_ratio_(i), found)
            if(.not. found) call call_error("'Phase.Phase property "//to_str(i)//".Specific heat ratio' is not found in configuration you set. Please check your configuration file.")

            call a_configuration%get_real   ("Phase.Phase property "//to_str(i)//".Reference pressure", self%reference_pressure_(i), found)
            if(.not. found) call call_error("'Phase.Phase property "//to_str(i)//".Reference pressure' is not found in configuration you set. Please check your configuration file.")
        end do
    end subroutine

    pure function compute_pressure(self, specific_internal_energy, density, volume_fraction) result(pressure)
        class(stiffened_gas_eos), intent(in) :: self
        real(real_kind), intent(in) :: specific_internal_energy
        real(real_kind), intent(in) :: density
        real(real_kind), intent(in) :: volume_fraction
        real(real_kind)             :: pressure
        real(real_kind)             :: g, pref

        g    = self%compute_mixture_specific_heat_ratio(volume_fraction)
        pref = self%compute_mixture_reference_pressure (volume_fraction, g)

        pressure = (g - 1.d0) * specific_internal_energy * density - g * pref
    end function compute_pressure

    pure function compute_soundspeed(self, pressure, density, volume_fraction) result(soundspeed)
        class(stiffened_gas_eos), intent(in) :: self
        real(real_kind), intent(in) :: pressure
        real(real_kind), intent(in) :: density
        real(real_kind), intent(in) :: volume_fraction
        real(real_kind)             :: soundspeed
        real(real_kind)             :: g, pref

        g    = self%compute_mixture_specific_heat_ratio(volume_fraction)
        pref = self%compute_mixture_reference_pressure (volume_fraction, g)
        soundspeed = sqrt(g * (pressure + pref) / density)
    end function compute_soundspeed

    pure function compute_internal_energy_density(self, pressure, density, volume_fraction) result(specific_internal_energy)
        class(stiffened_gas_eos), intent(in) :: self
        real(real_kind), intent(in) :: pressure
        real(real_kind), intent(in) :: density
        real(real_kind), intent(in) :: volume_fraction
        real(real_kind)             :: specific_internal_energy
        real(real_kind)             :: g, pref

        g    = self%compute_mixture_specific_heat_ratio(volume_fraction)
        pref = self%compute_mixture_reference_pressure (volume_fraction, g)
        specific_internal_energy = (pressure + g * pref) / (g - 1.d0)
    end function compute_internal_energy_density

    pure function compute_wood_soundspeed(self, pressure, densities, volume_fraction) result(soundspeed)
        class(stiffened_gas_eos), intent(in) :: self
        real(real_kind), intent(in) :: pressure
        real(real_kind), intent(in) :: densities(:)
        real(real_kind), intent(in) :: volume_fraction
        real(real_kind)             :: soundspeed
        real(real_kind)             :: g1, pref1, g2, pref2
        real(real_kind)             :: squared_c1, squared_c2

        if (volume_fraction < machine_epsilon) then
            squared_c2 = self%specific_heat_ratio_(2) * (pressure + self%reference_pressure_(2)) / densities(2)
            soundspeed = 1.d0 / densities(2) * squared_c2
        else if (1.d0 - machine_epsilon < volume_fraction) then
            squared_c1 = self%specific_heat_ratio_(1) * (pressure + self%reference_pressure_(2)) / densities(2)
            soundspeed = 1.d0 / densities(1) * squared_c1
        else
            squared_c1 = self%specific_heat_ratio_(1) * (pressure + self%reference_pressure_(1)) / densities(1)
            squared_c2 = self%specific_heat_ratio_(2) * (pressure + self%reference_pressure_(2)) / densities(2)
            soundspeed = volume_fraction / densities(1) * squared_c1 + (1.d0 - volume_fraction)  / densities(2) * squared_c2
        end if
    end function compute_wood_soundspeed

    pure function compute_k(self, pressure, densities, volume_fraction) result(k)
        class(stiffened_gas_eos), intent(in) :: self
        real(real_kind), intent(in) :: pressure
        real(real_kind), intent(in) :: densities(:)
        real(real_kind), intent(in) :: volume_fraction
        real(real_kind)             :: k
        real(real_kind)             :: g1, pref1, g2, pref2
        real(real_kind)             :: squared_c1, squared_c2

        if (volume_fraction < machine_epsilon) then
            k = 0.d0
        else if (1.d0 - machine_epsilon < volume_fraction) then
            k = 0.d0
        else
            squared_c1 = self%specific_heat_ratio_(1) * (pressure + self%reference_pressure_(1)) / densities(1)
            squared_c2 = self%specific_heat_ratio_(2) * (pressure + self%reference_pressure_(2)) / densities(2)
            k = volume_fraction * (1.d0 - volume_fraction) * (densities(1) * squared_c1 - densities(2) * squared_c2) &
              / (volume_fraction * densities(2) * squared_c2 + (1.d0 - volume_fraction) * densities(1) * squared_c1)
        end if
    end function compute_k

    pure function compute_mixture_specific_heat_ratio(self, volume_fraction) result(g)
        class(stiffened_gas_eos), intent(in) :: self
        real(real_kind), intent(in) :: volume_fraction
        real(real_kind)             :: g
        g = 1.d0 + ((self%specific_heat_ratio_(1) - 1.d0) * (self%specific_heat_ratio_(2) - 1.d0)) &
                 / (volume_fraction * (self%specific_heat_ratio_(2) - 1.d0) + (1.d0 - volume_fraction) * (self%specific_heat_ratio_(1) - 1.d0))
    end function compute_mixture_specific_heat_ratio

    pure function compute_mixture_reference_pressure(self, volume_fraction, mixture_specific_heat_ratio) result(p)
        class(stiffened_gas_eos), intent(in) :: self
        real(real_kind), intent(in) :: volume_fraction, mixture_specific_heat_ratio
        real(real_kind)             :: p
        p = (mixture_specific_heat_ratio - 1.d0) / mixture_specific_heat_ratio &
          * self%compute_phi(volume_fraction)
    end function compute_mixture_reference_pressure

    pure function compute_phi(self, volume_fraction) result(phi)
        class(stiffened_gas_eos), intent(in) :: self
        real(real_kind), intent(in) :: volume_fraction
        real(real_kind)             :: phi
        phi =         volume_fraction  * self%specific_heat_ratio_(1) * self%reference_pressure_(1) / (self%specific_heat_ratio_(1) - 1.d0) &
            + (1.d0 - volume_fraction) * self%specific_heat_ratio_(2) * self%reference_pressure_(2) / (self%specific_heat_ratio_(2) - 1.d0)
    end function compute_phi
end module class_stiffened_gas_eos