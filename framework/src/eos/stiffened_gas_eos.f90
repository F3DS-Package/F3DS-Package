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
        procedure, public, pass(self) :: compute_isobaric_pressure
        procedure, public, pass(self) :: compute_internal_energy_density
        procedure, public, pass(self) :: compute_isobaric_internal_energy_density
        procedure, public, pass(self) :: compute_soundspeed
        procedure, public, pass(self) :: compute_mixture_soundspeed
        procedure, public, pass(self) :: compute_frozen_soundspeed
        procedure, public, pass(self) :: compute_wood_soundspeed

        procedure, pass(self) :: compute_pi
        procedure, pass(self) :: compute_gamma
        procedure, pass(self) :: compute_mixture_density
    end type stiffened_gas_eos

    contains

    subroutine initialize(self, a_configuration, an_eos_generator)
        class(stiffened_gas_eos),           intent(inout) :: self
        class(configuration    ),           intent(inout) :: a_configuration
        class(eos_generator    ), optional, intent(inout) :: an_eos_generator

        logical           :: found
        integer(int_kind) :: i

        call a_configuration%get_int    ("Phase.Number of phase", self%num_phase_, found)
        if(.not. found) call call_error("'Phase.Number of phase' is not found in configuration you set. Please check your configuration file.")

        allocate(self%specific_heat_ratio_(self%num_phase_))
        allocate(self%reference_pressure_ (self%num_phase_))

        do i = 1, self%num_phase_, 1
            call a_configuration%get_real   ("Phase.Phase property "//to_str(i)//".Specific heat ratio", self%specific_heat_ratio_(i), found)
            if(.not. found) call call_error("'Phase.Phase property "//to_str(i)//".Specific heat ratio' is not found in configuration you set. Please check your configuration file.")

            call a_configuration%get_real   ("Phase.Phase property "//to_str(i)//".Reference pressure", self%reference_pressure_(i), found)
            if(.not. found) call call_error("'Phase.Phase property "//to_str(i)//".Reference pressure' is not found in configuration you set. Please check your configuration file.")
        end do
    end subroutine

    pure function compute_pressure(self, specific_internal_energy, density, phase_num) result(pressure)
        class  (stiffened_gas_eos), intent(in)           :: self
        real   (real_kind        ), intent(in)           :: specific_internal_energy
        real   (real_kind        ), intent(in)           :: density
        integer(int_kind         ), intent(in), optional :: phase_num
        real   (real_kind        )                       :: pressure

        integer :: n
        if(present(phase_num))then
            n = phase_num
        else
            n = 1
        endif

        pressure = (self%specific_heat_ratio_(n) - 1.0_real_kind) * specific_internal_energy * density &
                 - self%specific_heat_ratio_ (n) * self%reference_pressure_(n)
    end function compute_pressure

    pure function compute_isobaric_pressure(self, specific_internal_energy, densities, volume_fractions) result(pressure)
        class(stiffened_gas_eos), intent(in) :: self
        real (real_kind        ), intent(in) :: specific_internal_energy
        real (real_kind        ), intent(in) :: densities(:)
        real (real_kind        ), intent(in) :: volume_fractions(:)
        real (real_kind        )             :: pressure

        pressure = (specific_internal_energy * self%compute_mixture_density(densities, volume_fractions) - self%compute_pi(volume_fractions)) / self%compute_gamma(volume_fractions)
    end function compute_isobaric_pressure

    pure function compute_internal_energy_density(self, pressure, density, phase_num) result(energy_density)
        class  (stiffened_gas_eos), intent(in)           :: self
        real   (real_kind        ), intent(in)           :: pressure
        real   (real_kind        ), intent(in)           :: density
        integer(int_kind         ), intent(in), optional :: phase_num
        real   (real_kind        )                       :: energy_density

        integer :: n
        if(present(phase_num))then
            n = phase_num
        else
            n = 1
        endif

        energy_density = (pressure + self%specific_heat_ratio_(n) * self%reference_pressure_(n)) &
                       / (self%specific_heat_ratio_(n) - 1.0_real_kind)
    end function compute_internal_energy_density

    pure function compute_isobaric_internal_energy_density(self, pressure, densities, volume_fractions) result(energy_density)
        class  (stiffened_gas_eos), intent(in) :: self
        real   (real_kind        ), intent(in) :: pressure
        real   (real_kind        ), intent(in) :: densities(:)
        real   (real_kind        ), intent(in) :: volume_fractions(:)
        real   (real_kind        )             :: energy_density

        energy_density = self%compute_gamma(volume_fractions) * pressure + self%compute_pi(volume_fractions)
    end function compute_isobaric_internal_energy_density

    pure function compute_soundspeed(self, pressure, density, phase_num) result(soundspeed)
        class  (stiffened_gas_eos), intent(in)           :: self
        real   (real_kind        ), intent(in)           :: pressure
        real   (real_kind        ), intent(in)           :: density
        integer(int_kind         ), intent(in), optional :: phase_num
        real   (real_kind        )                       :: soundspeed

        integer :: n
        if(present(phase_num))then
            n = phase_num
        else
            n = 1
        endif

        ! Protect square root of a negative number {@code max(pressure + self%reference_pressure_(n), 0.0_real_kind)} and zero division {@code max(density, machine_epsilon)}.
        soundspeed = sqrt(self%specific_heat_ratio_(n) * max(pressure + self%reference_pressure_(n), 0.0_real_kind) / max(density, machine_epsilon))
    end function compute_soundspeed

    pure function compute_mixture_soundspeed(self, pressure, densities, volume_fractions) result(soundspeed)
        class(stiffened_gas_eos), intent(in) :: self
        real (real_kind        ), intent(in) :: pressure
        real (real_kind        ), intent(in) :: densities(:)
        real (real_kind        ), intent(in) :: volume_fractions(:)
        real (real_kind        )             :: soundspeed

        real   (real_kind) :: sigma
        integer(int_kind ) :: n

        sigma = 0.0_real_kind
        do n = 1, self%num_phase_, 1
            sigma = sigma &
                  + volume_fractions(n) * densities(n) * self%compute_soundspeed(pressure, densities(n), n) ** 2.0_real_kind &
                  / (self%specific_heat_ratio_(n) - 1.0_real_kind)
        end do
        soundspeed = sqrt(sigma / (self%compute_mixture_density(densities, volume_fractions) * self%compute_gamma(volume_fractions)))
    end function compute_mixture_soundspeed

    pure function compute_frozen_soundspeed(self, pressures, densities, volume_fractions) result(soundspeed)
        class(stiffened_gas_eos), intent(in) :: self
        real (real_kind        ), intent(in) :: pressures(:)
        real (real_kind        ), intent(in) :: densities(:)
        real (real_kind        ), intent(in) :: volume_fractions(:)
        real (real_kind        )             :: soundspeed

        real   (real_kind) :: sigma, mixture_density
        integer(int_kind ) :: n

        mixture_density = self%compute_mixture_density(densities, volume_fractions)

        sigma = 0.0_real_kind
        do n = 1, self%num_phase_, 1
            associate(mass_fraction => (volume_fractions(n) * densities(n) / mixture_density))
                sigma = sigma + mass_fraction * self%compute_soundspeed(pressures(n), densities(n), n) ** 2.0_real_kind
            end associate
        end do
        soundspeed = sqrt(sigma)
    end function compute_frozen_soundspeed

    pure function compute_wood_soundspeed(self, pressure, densities, volume_fractions) result(soundspeed)
        class(stiffened_gas_eos), intent(in) :: self
        real (real_kind        ), intent(in) :: pressure
        real (real_kind        ), intent(in) :: densities(:)
        real (real_kind        ), intent(in) :: volume_fractions(:)
        real (real_kind        )             :: soundspeed

        real   (real_kind) :: sigma
        integer(int_kind ) :: n

        sigma = 0.0_real_kind
        do n = 1, self%num_phase_, 1
            if(densities(n) <= 0.0_real_kind)cycle
            sigma = sigma + volume_fractions(n) / (densities(n) * self%compute_soundspeed(pressure, densities(n), n))
        end do
        soundspeed = (self%compute_mixture_density(densities, volume_fractions) * sigma)**(-0.5_real_kind)
    end function compute_wood_soundspeed

    pure function compute_pi(self, volume_fractions) result(pi)
        class(stiffened_gas_eos), intent(in) :: self
        real (real_kind        ), intent(in) :: volume_fractions(:)
        real (real_kind        )             :: pi

        integer :: n

        pi = 0.0_real_kind
        do n = 1, self%num_phase_, 1
            pi = pi + volume_fractions(n) * self%specific_heat_ratio_(n) * self%reference_pressure_(n) / (self%specific_heat_ratio_(n) - 1.0_real_kind)
        end do
    end function compute_pi

    pure function compute_gamma(self, volume_fractions) result(gamma)
        class(stiffened_gas_eos), intent(in) :: self
        real (real_kind        ), intent(in) :: volume_fractions(:)
        real (real_kind        )             :: gamma

        integer :: n

        gamma = 0.0_real_kind
        do n = 1, self%num_phase_, 1
            gamma = gamma + volume_fractions(n) / (self%specific_heat_ratio_(n) - 1.0_real_kind)
        end do
    end function compute_gamma

    pure function compute_mixture_density(self, densities, volume_fractions) result(density)
        class(stiffened_gas_eos), intent(in) :: self
        real (real_kind        ), intent(in) :: densities(:)
        real (real_kind        ), intent(in) :: volume_fractions(:)
        real (real_kind        )             :: density

        integer :: n

        density = 0.0_real_kind
        do n = 1, self%num_phase_, 1
            density = density + densities(n) * volume_fractions(n)
        end do
    end function compute_mixture_density
end module class_stiffened_gas_eos