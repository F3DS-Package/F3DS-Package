module abstract_eos
    implicit none

    private

    type, public, abstract :: eos_generator
        contains
        generic, public :: generate => generate_from_configuration, generate_from_name
        procedure(generate_from_configuration_interface), pass(self), deferred :: generate_from_configuration
        procedure(generate_from_name_interface         ), pass(self), deferred :: generate_from_name
    end type eos_generator

    type, public, abstract :: eos
        contains
        procedure(initialize_interface                              ), pass(self), deferred :: initialize
        procedure(compute_pressure_interface                        ), pass(self), deferred :: compute_pressure
        procedure(compute_isobaric_pressure_interface               ), pass(self), deferred :: compute_isobaric_pressure
        procedure(compute_internal_energy_density_interface         ), pass(self), deferred :: compute_internal_energy_density
        procedure(compute_isobaric_internal_energy_density_interface), pass(self), deferred :: compute_isobaric_internal_energy_density
        procedure(compute_soundspeed_interface                      ), pass(self), deferred :: compute_soundspeed
        procedure(compute_mixture_soundspeed_interface              ), pass(self), deferred :: compute_mixture_soundspeed
        procedure(compute_frozen_soundspeed_interface               ), pass(self), deferred :: compute_frozen_soundspeed
        procedure(compute_wood_soundspeed_interface                 ), pass(self), deferred :: compute_wood_soundspeed
    end type eos

    abstract interface
        subroutine generate_from_configuration_interface(self, an_eos, a_config)
            use abstract_configuration
            import eos_generator
            import eos
            class(eos_generator),          intent(inout) :: self
            class(eos          ), pointer, intent(inout) :: an_eos
            class(configuration),          intent(inout) :: a_config
        end subroutine generate_from_configuration_interface

        subroutine generate_from_name_interface(self, an_eos, name)
            use abstract_configuration
            import eos_generator
            import eos
            class    (eos_generator),              intent(inout) :: self
            class    (eos          ), pointer    , intent(inout) :: an_eos
            character(len=:        ), allocatable, intent(in   ) :: name
        end subroutine generate_from_name_interface
    end interface

    abstract interface
        subroutine initialize_interface(self, a_configuration, an_eos_generator)
            use abstract_configuration
            import eos
            import eos_generator
            class(eos          ),           intent(inout) :: self
            class(configuration),           intent(inout) :: a_configuration
            class(eos_generator), optional, intent(inout) :: an_eos_generator
        end subroutine initialize_interface

        pure function compute_pressure_interface(self, specific_internal_energy, density, phase_num) result(pressure)
            use typedef_module
            import eos
            class  (eos      ), intent(in)           :: self
            real   (real_kind), intent(in)           :: specific_internal_energy
            real   (real_kind), intent(in)           :: density
            integer(int_kind ), intent(in), optional :: phase_num
            real   (real_kind)                       :: pressure
        end function compute_pressure_interface

        pure function compute_isobaric_pressure_interface(self, specific_internal_energy, densities, volume_fractions) result(pressure)
            use typedef_module
            import eos
            class(eos      ), intent(in) :: self
            real (real_kind), intent(in) :: specific_internal_energy
            real (real_kind), intent(in) :: densities(:)
            real (real_kind), intent(in) :: volume_fractions(:)
            real (real_kind)             :: pressure
        end function compute_isobaric_pressure_interface

        pure function compute_internal_energy_density_interface(self, pressure, density, phase_num) result(energy_density)
            use typedef_module
            import eos
            class  (eos      ), intent(in)           :: self
            real   (real_kind), intent(in)           :: pressure
            real   (real_kind), intent(in)           :: density
            integer(int_kind ), intent(in), optional :: phase_num
            real   (real_kind)                       :: energy_density
        end function compute_internal_energy_density_interface

        pure function compute_isobaric_internal_energy_density_interface(self, pressure, densities, volume_fractions) result(energy_density)
            use typedef_module
            import eos
            class(eos      ), intent(in) :: self
            real (real_kind), intent(in) :: pressure
            real (real_kind), intent(in) :: densities(:)
            real (real_kind), intent(in) :: volume_fractions(:)
            real (real_kind)             :: energy_density
        end function compute_isobaric_internal_energy_density_interface

        pure function compute_soundspeed_interface(self, pressure, density, phase_num) result(soundspeed)
            use typedef_module
            import eos
            class  (eos      ), intent(in)           :: self
            real   (real_kind), intent(in)           :: pressure
            real   (real_kind), intent(in)           :: density
            integer(int_kind ), intent(in), optional :: phase_num
            real   (real_kind)                       :: soundspeed
        end function compute_soundspeed_interface

        pure function compute_mixture_soundspeed_interface(self, pressure, densities, volume_fractions) result(soundspeed)
            use typedef_module
            import eos
            class(eos      ), intent(in) :: self
            real (real_kind), intent(in) :: pressure
            real (real_kind), intent(in) :: densities(:)
            real (real_kind), intent(in) :: volume_fractions(:)
            real (real_kind)             :: soundspeed
        end function compute_mixture_soundspeed_interface

        pure function compute_wood_soundspeed_interface(self, pressure, densities, volume_fractions) result(soundspeed)
            use typedef_module
            import eos
            class(eos      ), intent(in) :: self
            real (real_kind), intent(in) :: pressure
            real (real_kind), intent(in) :: densities(:)
            real (real_kind), intent(in) :: volume_fractions(:)
            real (real_kind)             :: soundspeed
        end function compute_wood_soundspeed_interface

        pure function compute_frozen_soundspeed_interface(self, pressures, densities, volume_fractions) result(soundspeed)
            use typedef_module
            import eos
            class(eos      ), intent(in) :: self
            real (real_kind), intent(in) :: pressures(:)
            real (real_kind), intent(in) :: densities(:)
            real (real_kind), intent(in) :: volume_fractions(:)
            real (real_kind)             :: soundspeed
        end function compute_frozen_soundspeed_interface
    end interface
end module abstract_eos