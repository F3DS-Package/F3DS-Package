module abstract_eos
    ! > @module abstract_eos
    ! This module defines abstract types and interfaces for equations of state (EOS).
    ! > @type eos_generator
    ! Abstract type for generating EOS objects.
    ! @abstract
    ! > @type eos
    ! Abstract type for equations of state.
    ! @abstract
    ! > @interface
    ! Abstract interface for generating EOS objects from a configuration.
    ! > @interface
    ! Abstract interface for initializing an EOS object.
    implicit none

    private

    type, public, abstract :: eos_generator
        ! """
        ! Fortran type representing an equation of state (EOS) generator.
        ! Attributes:
        ! - generate: A generic procedure for generating an EOS. It can be called with either the `generate_from_configuration` or
        ! `generate_from_name` procedure.
        ! - generate_from_configuration: A deferred procedure for generating an EOS from a configuration.
        ! - generate_from_name: A deferred procedure for generating an EOS from a name.
        ! Note: This type is abstract and cannot be instantiated directly.
        ! """
        contains
        generic, public :: generate => generate_from_configuration, generate_from_name
        procedure(generate_from_configuration_interface), pass(self), deferred :: generate_from_configuration
        procedure(generate_from_name_interface         ), pass(self), deferred :: generate_from_name
    end type eos_generator

    type, public, abstract :: eos
        ! This type represents an equation of state (EOS) for a material.
        ! It provides methods for initializing the EOS and computing various properties.
        ! Initializes the EOS.
        ! Computes the pressure of the material.
        ! Computes the isobaric pressure of the material.
        ! Computes the internal energy density of the material.
        ! Computes the isobaric internal energy density of the material.
        ! Computes the sound speed of the material.
        ! Computes the mixture sound speed of the material.
        ! Computes the frozen sound speed of the material.
        ! Computes the Wood sound speed of the material.
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
            ! -----------------------------------------------------------------------
            ! Subroutine: generate_from_configuration_interface
            ! -----------------------------------------------------------------------
            ! This subroutine generates an equation of state (EOS) using the provided
            ! configuration object. It is a part of the eos_generator class.
            !
            ! Parameters:
            ! - self: An instance of the eos_generator class. It is passed as an
            ! inout argument.
            ! - an_eos: A pointer to an instance of the eos class. It is passed as
            ! an inout argument.
            ! - a_config: An instance of the configuration class. It is passed as
            ! an inout argument.
            !
            ! -----------------------------------------------------------------------
            use abstract_configuration
            import eos_generator
            import eos
            class(eos_generator),          intent(inout) :: self
            class(eos          ), pointer, intent(inout) :: an_eos
            class(configuration),          intent(inout) :: a_config
        end subroutine generate_from_configuration_interface

        subroutine generate_from_name_interface(self, an_eos, name)
            ! -----------------------------------------------------------------------
            ! This subroutine generates an equation of state (EOS) object based on
            ! the given name.
            !
            ! Parameters:
            ! self:           An instance of the eos_generator class.
            ! an_eos:         A pointer to an instance of the eos class.
            ! name:           The name of the EOS to generate.
            !
            ! -----------------------------------------------------------------------
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
            ! -----------------------------------------------------------------------
            ! This subroutine initializes the interface of the `eos` object.
            !
            ! Parameters:
            ! self:               The `eos` object to be initialized.
            ! a_configuration:    The configuration object to be used.
            ! an_eos_generator:   (Optional) The EOS generator object to be used.
            !
            ! -----------------------------------------------------------------------
            use abstract_configuration
            import eos
            import eos_generator
            class(eos          ),           intent(inout) :: self
            class(configuration),           intent(inout) :: a_configuration
            class(eos_generator), optional, intent(inout) :: an_eos_generator
        end subroutine initialize_interface

        pure function compute_pressure_interface(self, specific_internal_energy, density, phase_num) result(pressure)
            ! """
            ! Compute the pressure at the interface.
            ! Parameters:
            !     self (eos): An instance of the equation of state class.
            !     specific_internal_energy (real): The specific internal energy.
            !     density (real): The density.
            !     phase_num (integer, optional): The phase number. Default is None.
            ! Returns:
            !     pressure (real): The computed pressure at the interface.
            ! """
            use typedef_module
            import eos
            class  (eos      ), intent(in)           :: self
            real   (real_kind), intent(in)           :: specific_internal_energy
            real   (real_kind), intent(in)           :: density
            integer(int_kind ), intent(in), optional :: phase_num
            real   (real_kind)                       :: pressure
        end function compute_pressure_interface

        pure function compute_isobaric_pressure_interface(self, specific_internal_energy, densities, volume_fractions) result(pressure)
            ! """
            ! Compute the isobaric pressure at the interface.
            ! Parameters:
            !     self (eos): The equation of state object.
            !     specific_internal_energy (real): The specific internal energy.
            !     densities (real array): The array of densities.
            !     volume_fractions (real array): The array of volume fractions.
            ! Returns:
            !     pressure (real): The computed isobaric pressure.
            ! """
            use typedef_module
            import eos
            class(eos      ), intent(in) :: self
            real (real_kind), intent(in) :: specific_internal_energy
            real (real_kind), intent(in) :: densities(:)
            real (real_kind), intent(in) :: volume_fractions(:)
            real (real_kind)             :: pressure
        end function compute_isobaric_pressure_interface

        pure function compute_internal_energy_density_interface(self, pressure, density, phase_num) result(energy_density)
            ! """
            ! Compute the internal energy density.
            ! Parameters:
            !     self (eos): The equation of state object.
            !     pressure (real): The pressure.
            !     density (real): The density.
            !     phase_num (optional, int): The phase number. Default is None.
            ! Returns:
            !     energy_density (real): The computed internal energy density.
            ! """
            use typedef_module
            import eos
            class  (eos      ), intent(in)           :: self
            real   (real_kind), intent(in)           :: pressure
            real   (real_kind), intent(in)           :: density
            integer(int_kind ), intent(in), optional :: phase_num
            real   (real_kind)                       :: energy_density
        end function compute_internal_energy_density_interface

        pure function compute_isobaric_internal_energy_density_interface(self, pressure, densities, volume_fractions) result(energy_density)
            ! """
            ! Compute the isobaric internal energy density.
            ! Parameters:
            !     self (eos): The equation of state object.
            !     pressure (real): The pressure.
            !     densities (real array): The array of densities.
            !     volume_fractions (real array): The array of volume fractions.
            ! Returns:
            !     energy_density (real): The computed isobaric internal energy density.
            ! """
            use typedef_module
            import eos
            class(eos      ), intent(in) :: self
            real (real_kind), intent(in) :: pressure
            real (real_kind), intent(in) :: densities(:)
            real (real_kind), intent(in) :: volume_fractions(:)
            real (real_kind)             :: energy_density
        end function compute_isobaric_internal_energy_density_interface

        pure function compute_soundspeed_interface(self, pressure, density, phase_num) result(soundspeed)
            ! """
            ! Compute the speed of sound at the interface.
            ! Parameters:
            !     self (eos): The equation of state object.
            !     pressure (real): The pressure at the interface.
            !     density (real): The density at the interface.
            !     phase_num (int, optional): The phase number. Default is None.
            ! Returns:
            !     soundspeed (real): The speed of sound at the interface.
            ! """
            use typedef_module
            import eos
            class  (eos      ), intent(in)           :: self
            real   (real_kind), intent(in)           :: pressure
            real   (real_kind), intent(in)           :: density
            integer(int_kind ), intent(in), optional :: phase_num
            real   (real_kind)                       :: soundspeed
        end function compute_soundspeed_interface

        pure function compute_mixture_soundspeed_interface(self, pressure, densities, volume_fractions) result(soundspeed)
            ! """
            ! Compute the speed of sound in a mixture of substances.
            ! Parameters:
            !     self (eos): An instance of the equation of state class.
            !     pressure (real): The pressure of the mixture.
            !     densities (real array): An array of densities of the substances in the mixture.
            !     volume_fractions (real array): An array of volume fractions of the substances in the mixture.
            ! Returns:
            !     soundspeed (real): The speed of sound in the mixture.
            ! """
            use typedef_module
            import eos
            class(eos      ), intent(in) :: self
            real (real_kind), intent(in) :: pressure
            real (real_kind), intent(in) :: densities(:)
            real (real_kind), intent(in) :: volume_fractions(:)
            real (real_kind)             :: soundspeed
        end function compute_mixture_soundspeed_interface

        pure function compute_wood_soundspeed_interface(self, pressure, densities, volume_fractions) result(soundspeed)
            ! """
            ! Compute the speed of sound in a wood material.
            ! Parameters:
            !     self (eos): The equation of state object.
            !     pressure (real): The pressure of the wood material.
            !     densities (real array): The densities of the wood material components.
            !     volume_fractions (real array): The volume fractions of the wood material components.
            ! Returns:
            !     soundspeed (real): The speed of sound in the wood material.
            ! """
            use typedef_module
            import eos
            class(eos      ), intent(in) :: self
            real (real_kind), intent(in) :: pressure
            real (real_kind), intent(in) :: densities(:)
            real (real_kind), intent(in) :: volume_fractions(:)
            real (real_kind)             :: soundspeed
        end function compute_wood_soundspeed_interface

        pure function compute_frozen_soundspeed_interface(self, pressures, densities, volume_fractions) result(soundspeed)
            ! """
            ! Compute the frozen sound speed at the interface.
            ! Parameters:
            !     self (eos): The equation of state object.
            !     pressures (array): Array of pressures.
            !     densities (array): Array of densities.
            !     volume_fractions (array): Array of volume fractions.
            ! Returns:
            !     soundspeed (float): The frozen sound speed at the interface.
            ! """
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