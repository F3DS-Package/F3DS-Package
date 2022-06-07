module class_stiffened_gas_mixture_eos
    use typedef_module
    use abstract_mixture_eos
    use math_constant_module

    implicit none

    private

    type, public, extends(mixture_eos) :: stiffened_gas_mixture_eos
        private

        real(real_kind) :: specific_heat_ratio_fluid1_
        real(real_kind) :: specific_heat_ratio_fluid2_
        real(real_kind) :: reference_pressure_fluid1_
        real(real_kind) :: reference_pressure_fluid2_

        contains

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: compute_pressure
        procedure, public, pass(self) :: compute_soundspeed
        procedure, public, pass(self) :: compute_internal_energy_density
        procedure, public, pass(self) :: compute_wood_soundspeed
        procedure, public, pass(self) :: compute_k

        procedure, pass(self) :: compute_mixture_specific_heat_ratio
        procedure, pass(self) :: compute_mixture_reference_pressure
        procedure, pass(self) :: compute_pi
    end type stiffened_gas_mixture_eos

    contains

    subroutine initialize(self, specific_heat_ratio_fluid1, specific_heat_ratio_fluid2, reference_pressure_fluid1, reference_pressure_fluid2)
        class(stiffened_gas_mixture_eos), intent(inout) :: self
        real(real_kind), intent(in) :: specific_heat_ratio_fluid1
        real(real_kind), intent(in) :: specific_heat_ratio_fluid2
        real(real_kind), intent(in) :: reference_pressure_fluid1
        real(real_kind), intent(in) :: reference_pressure_fluid2

        self%specific_heat_ratio_fluid1_ = specific_heat_ratio_fluid1
        self%specific_heat_ratio_fluid2_ = specific_heat_ratio_fluid2
        self%reference_pressure_fluid1_  = reference_pressure_fluid1
        self%reference_pressure_fluid2_  = reference_pressure_fluid2
    end subroutine

    pure function compute_pressure(self, specific_internal_energy, density, volume_fruction) result(pressure)
        class(stiffened_gas_mixture_eos), intent(in) :: self
        real(real_kind), intent(in) :: specific_internal_energy
        real(real_kind), intent(in) :: density
        real(real_kind), intent(in) :: volume_fruction
        real(real_kind)             :: pressure
        real(real_kind)             :: g, pref

        g    = self%compute_mixture_specific_heat_ratio(volume_fruction)
        pref = self%compute_mixture_reference_pressure (volume_fruction, g)

        pressure = (g - 1.d0) * specific_internal_energy * density - g * pref
    end function compute_pressure

    pure function compute_soundspeed(self, pressure, density, volume_fruction) result(soundspeed)
        class(stiffened_gas_mixture_eos), intent(in) :: self
        real(real_kind), intent(in) :: pressure
        real(real_kind), intent(in) :: density
        real(real_kind), intent(in) :: volume_fruction
        real(real_kind)             :: soundspeed
        real(real_kind)             :: g, pref

        g    = self%compute_mixture_specific_heat_ratio(volume_fruction)
        pref = self%compute_mixture_reference_pressure (volume_fruction, g)
        soundspeed = sqrt(g * (pressure + pref) / density)
    end function compute_soundspeed

    pure function compute_internal_energy_density(self, pressure, density, volume_fruction) result(specific_internal_energy)
        class(stiffened_gas_mixture_eos), intent(in) :: self
        real(real_kind), intent(in) :: pressure
        real(real_kind), intent(in) :: density
        real(real_kind), intent(in) :: volume_fruction
        real(real_kind)             :: specific_internal_energy
        real(real_kind)             :: g, pref

        g    = self%compute_mixture_specific_heat_ratio(volume_fruction)
        pref = self%compute_mixture_reference_pressure (volume_fruction, g)
        specific_internal_energy = (pressure + g * pref) / (g - 1.d0)
    end function compute_internal_energy_density

    pure function compute_wood_soundspeed(self, pressure, density1, density2, volume_fruction) result(soundspeed)
        class(stiffened_gas_mixture_eos), intent(in) :: self
        real(real_kind), intent(in) :: pressure
        real(real_kind), intent(in) :: density1
        real(real_kind), intent(in) :: density2
        real(real_kind), intent(in) :: volume_fruction
        real(real_kind)             :: soundspeed
        real(real_kind)             :: g1, pref1, g2, pref2
        real(real_kind)             :: squared_c1, squared_c2

        if (volume_fruction < machine_epsilon) then
            squared_c2 = self%specific_heat_ratio_fluid2_ * (pressure + self%reference_pressure_fluid2_) / density2
            soundspeed = 1.d0 / density2 * squared_c2
        else if (1.d0 - machine_epsilon < volume_fruction) then
            squared_c1 = self%specific_heat_ratio_fluid1_ * (pressure + self%reference_pressure_fluid1_) / density1
            soundspeed = 1.d0 / density1 * squared_c1
        else
            squared_c1 = self%specific_heat_ratio_fluid1_ * (pressure + self%reference_pressure_fluid1_) / density1
            squared_c2 = self%specific_heat_ratio_fluid2_ * (pressure + self%reference_pressure_fluid2_) / density2
            soundspeed = volume_fruction / density1 * squared_c1 + (1.d0 - volume_fruction) / density2 * squared_c2
        end if
    end function compute_wood_soundspeed

    pure function compute_k(self, pressure, density1, density2, volume_fruction) result(k)
        class(stiffened_gas_mixture_eos), intent(in) :: self
        real(real_kind), intent(in) :: pressure
        real(real_kind), intent(in) :: density1
        real(real_kind), intent(in) :: density2
        real(real_kind), intent(in) :: volume_fruction
        real(real_kind)             :: k
        real(real_kind)             :: g1, pref1, g2, pref2
        real(real_kind)             :: squared_c1, squared_c2

        if (volume_fruction < machine_epsilon) then
            k = 0.d0
        else if (1.d0 - machine_epsilon < volume_fruction) then
            k = 0.d0
        else
            squared_c1 = self%specific_heat_ratio_fluid1_ * (pressure + self%reference_pressure_fluid1_) / density1
            squared_c2 = self%specific_heat_ratio_fluid2_ * (pressure + self%reference_pressure_fluid2_) / density2
            k = volume_fruction * (1.d0 - volume_fruction) * (density1 * squared_c1 - density2 * squared_c2) &
              / (volume_fruction * density2 * squared_c2 + (1.d0 - volume_fruction) * density1 * squared_c1)
        end if
    end function compute_k

    pure function compute_mixture_specific_heat_ratio(self, volume_fruction) result(g)
        class(stiffened_gas_mixture_eos), intent(in) :: self
        real(real_kind), intent(in) :: volume_fruction
        real(real_kind)             :: g
        g = 1.d0 + ((self%specific_heat_ratio_fluid1_ - 1.d0) * (self%specific_heat_ratio_fluid2_ - 1.d0)) &
                 / (volume_fruction * (self%specific_heat_ratio_fluid2_ - 1.d0) + (1.d0 - volume_fruction) * (self%specific_heat_ratio_fluid1_ - 1.d0))
    end function compute_mixture_specific_heat_ratio

    pure function compute_mixture_reference_pressure(self, volume_fruction, mixture_specific_heat_ratio) result(p)
        class(stiffened_gas_mixture_eos), intent(in) :: self
        real(real_kind), intent(in) :: volume_fruction, mixture_specific_heat_ratio
        real(real_kind)             :: p
        p = (mixture_specific_heat_ratio - 1.d0) / mixture_specific_heat_ratio &
          * self%compute_pi(volume_fruction)
    end function compute_mixture_reference_pressure

    pure function compute_pi(self, volume_fruction) result(pi)
        class(stiffened_gas_mixture_eos), intent(in) :: self
        real(real_kind), intent(in) :: volume_fruction
        real(real_kind)             :: pi
        pi =         volume_fruction  * self%specific_heat_ratio_fluid1_ * self%reference_pressure_fluid1_ / (self%specific_heat_ratio_fluid1_ - 1.d0) &
           + (1.d0 - volume_fruction) * self%specific_heat_ratio_fluid2_ * self%reference_pressure_fluid2_ / (self%specific_heat_ratio_fluid2_ - 1.d0)
    end function compute_pi
end module class_stiffened_gas_mixture_eos