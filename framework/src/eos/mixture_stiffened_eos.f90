module mixture_stiffened_eos_module
    use typedef_module

    implicit none

    private

    real(real_kind) :: specific_heat_ratio_fluid1_
    real(real_kind) :: specific_heat_ratio_fluid2_
    real(real_kind) :: reference_pressure_fluid1_
    real(real_kind) :: reference_pressure_fluid2_

    public :: initialize_mixture_stiffened_eos
    public :: compute_pressure_mixture_stiffened_eos
    public :: compute_soundspeed_mixture_stiffened_eos

    contains

    subroutine initialize_mixture_stiffened_eos(specific_heat_ratio_fluid1, specific_heat_ratio_fluid2, reference_pressure_fluid1, reference_pressure_fluid2)
        real(real_kind), intent(in) :: specific_heat_ratio_fluid1
        real(real_kind), intent(in) :: specific_heat_ratio_fluid2
        real(real_kind), intent(in) :: reference_pressure_fluid1
        real(real_kind), intent(in) :: reference_pressure_fluid2

        specific_heat_ratio_fluid1_ = specific_heat_ratio_fluid1
        specific_heat_ratio_fluid2_ = specific_heat_ratio_fluid2
        reference_pressure_fluid1_  = reference_pressure_fluid1
        reference_pressure_fluid2_  = reference_pressure_fluid2
    end subroutine

    pure function compute_pressure_mixture_stiffened_eos(specific_internal_energy, density, volume_fruction) result(pressure)
        use typedef_module
        real(real_kind), intent(in) :: specific_internal_energy
        real(real_kind), intent(in) :: density
        real(real_kind), intent(in) :: volume_fruction
        real(real_kind)             :: pressure

        pressure = (specific_internal_energy * density - compute_pi(volume_fruction)) / compute_gamma(volume_fruction)
    end function compute_pressure_mixture_stiffened_eos

    pure function compute_soundspeed_mixture_stiffened_eos(specific_internal_energy, density, volume_fruction) result(soundspeed)
        use typedef_module
        real(real_kind), intent(in) :: specific_internal_energy
        real(real_kind), intent(in) :: density
        real(real_kind), intent(in) :: volume_fruction
        real(real_kind)             :: soundspeed
        real(real_kind)             :: g, p, pref

        g    = compute_mixture_specific_heat_ratio(volume_fruction)
        pref = compute_mixture_reference_pressure (volume_fruction, g)
        p    = compute_pressure_mixture_stiffened_eos(specific_internal_energy, density, volume_fruction)
        soundspeed = sqrt(g * (p + pref) / density)
    end function compute_soundspeed_mixture_stiffened_eos

    pure function compute_mixture_specific_heat_ratio(volume_fruction) result(g)
        real(real_kind), intent(in) :: volume_fruction
        real(real_kind)             :: g
        g = 1.d0 + ((specific_heat_ratio_fluid1_ - 1.d0) * (specific_heat_ratio_fluid2_ - 1.d0)) &
                 / (volume_fruction * (specific_heat_ratio_fluid2_ - 1.d0) + (1.d0 - volume_fruction) * (specific_heat_ratio_fluid1_ - 1.d0))
    end function compute_mixture_specific_heat_ratio

    pure function compute_mixture_reference_pressure(volume_fruction, mixture_specific_heat_ratio) result(p)
        real(real_kind), intent(in) :: volume_fruction, mixture_specific_heat_ratio
        real(real_kind)             :: p
        p = (mixture_specific_heat_ratio - 1.d0) / mixture_specific_heat_ratio &
          * compute_pi(volume_fruction)
    end function compute_mixture_reference_pressure

    pure function compute_gamma(volume_fruction) result(g)
        real(real_kind), intent(in) :: volume_fruction
        real(real_kind)             :: g
        g =         volume_fruction  / (specific_heat_ratio_fluid1_ - 1.d0) &
          + (1.d0 - volume_fruction) / (specific_heat_ratio_fluid2_ - 1.d0)
    end function compute_gamma

    pure function compute_pi(volume_fruction) result(pi)
        real(real_kind), intent(in) :: volume_fruction
        real(real_kind)             :: pi
        pi =         volume_fruction  * specific_heat_ratio_fluid1_ * reference_pressure_fluid1_ / (specific_heat_ratio_fluid1_ - 1.d0) &
           + (1.d0 - volume_fruction) * specific_heat_ratio_fluid2_ * reference_pressure_fluid2_ / (specific_heat_ratio_fluid2_ - 1.d0)
    end function compute_pi
end module mixture_stiffened_eos_module