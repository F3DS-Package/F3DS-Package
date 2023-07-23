module abstract_time_increment_controller
    ! """
    ! This module defines an abstract interface for time increment controllers.
    ! The module contains two abstract types: `time_increment_controller_generator` and `time_increment_controller`.
    ! The `time_increment_controller_generator` type is a public abstract type that contains two deferred procedures:
    ! `generate_from_configuration` and `generate_from_name`. These procedures are used to generate instances of
    ! `time_increment_controller` based on a configuration object or a name, respectively.
    ! The `time_increment_controller` type is also a public abstract type that contains four deferred procedures: `initialize`,
    ! `get_constant_dt`, `compute_local_dt`, and `returns_constant`. These procedures define the behavior of a time increment
    ! controller.
    ! The module also defines two abstract interfaces: `generate_from_configuration_interface` and `generate_from_name_interface`.
    ! These interfaces are used to define the signatures of the `generate_from_configuration` and `generate_from_name` procedures,
    ! respectively.
    ! The module also defines four abstract interfaces: `initialize_interface`, `get_constant_dt_interface`,
    ! `compute_local_dt_interface`, and `returns_constant_interface`. These interfaces are used to define the signatures of the
    ! `initialize`, `get_constant_dt`, `compute_local_dt`, and `returns_constant` procedures, respectively.
    ! """
    implicit none

    private

    type, public, abstract :: time_increment_controller_generator
        ! """
        ! Fortran type representing a time increment controller generator.
        ! Attributes:
        !     generate (generic procedure): A generic procedure for generating a time increment controller.
        !     generate_from_configuration (deferred procedure): A deferred procedure for generating a time increment controller from a
        !     configuration.
        !     generate_from_name (deferred procedure): A deferred procedure for generating a time increment controller from a name.
        ! """
        contains
        generic, public :: generate => generate_from_configuration, generate_from_name
        procedure(generate_from_configuration_interface), pass(self), deferred :: generate_from_configuration
        procedure(generate_from_name_interface         ), pass(self), deferred :: generate_from_name
    end type time_increment_controller_generator

    type, public, abstract :: time_increment_controller
        ! """
        ! Fortran type representing a time increment controller.
        ! Attributes:
        !     initialize (procedure): Abstract method to initialize the interface.
        !     get_constant_dt (procedure): Abstract method to get the constant time increment.
        !     compute_local_dt (procedure): Abstract method to compute the local time increment.
        !     returns_constant (procedure): Abstract method to check if the time increment is constant.
        ! """
        contains

        procedure(initialize_interface      ), pass(self), deferred :: initialize
        procedure(get_constant_dt_interface ), pass(self), deferred :: get_constant_dt
        procedure(compute_local_dt_interface), pass(self), deferred :: compute_local_dt
        procedure(returns_constant_interface), pass(self), deferred :: returns_constant
    end type time_increment_controller

    abstract interface
        subroutine generate_from_configuration_interface(self, a_time_increment_controller, a_config)
            ! This subroutine generates a time increment controller from a given configuration.
            !
            ! Parameters:
            ! self: time_increment_controller_generator
            ! The time increment controller generator object.
            ! a_time_increment_controller: time_increment_controller, pointer
            ! The pointer to the time increment controller object.
            ! a_config: configuration
            ! The configuration object.
            !
            ! Returns:
            ! None
            !
            ! Notes:
            ! - This subroutine requires the use of the abstract_configuration module.
            ! - It imports the time_increment_controller_generator and time_increment_controller modules.
            !
            use abstract_configuration
            import time_increment_controller_generator
            import time_increment_controller
            class(time_increment_controller_generator),          intent(inout) :: self
            class(time_increment_controller          ), pointer, intent(inout) :: a_time_increment_controller
            class(configuration                      ),          intent(inout) :: a_config
        end subroutine generate_from_configuration_interface

        subroutine generate_from_name_interface(self, a_time_increment_controller, name)
            ! ---------------------------------------------------------------------------
            ! Generates a time increment controller from a given name.
            !
            ! Parameters:
            ! self: time_increment_controller_generator
            ! The generator object used to generate the time increment controller.
            ! a_time_increment_controller: time_increment_controller, pointer
            ! The pointer to the generated time increment controller.
            ! name: character(len=:), allocatable
            ! The name used to generate the time increment controller.
            ! ---------------------------------------------------------------------------
            use abstract_configuration
            import time_increment_controller_generator
            import time_increment_controller
            class    (time_increment_controller_generator),              intent(inout) :: self
            class    (time_increment_controller          ), pointer    , intent(inout) :: a_time_increment_controller
            character(len=:                              ), allocatable, intent(in   ) :: name
        end subroutine generate_from_name_interface
    end interface

    abstract interface
        subroutine initialize_interface(self, config)
            ! > Initializes the interface of the time increment controller with the given configuration.
            ! > This subroutine initializes the interface of the time increment controller by taking in the self object and the config object.
            ! The self object is of type time_increment_controller and is passed as an inout argument, meaning that it can be modified within
            ! the subroutine. The config object is of type configuration and is also passed as an inout argument.
            ! > Parameters:
            ! >   - self: The time increment controller object to be initialized.
            ! >   - config: The configuration object to be used for initialization.
            use abstract_configuration
            import time_increment_controller
            import time_increment_controller_generator

            class(time_increment_controller          ),           intent(inout) :: self
            class(configuration                      ),           intent(inout) :: config
        end subroutine initialize_interface

        function get_constant_dt_interface(self) result(dt)
            ! -----------------------------------------------------------------------
            ! Function: get_constant_dt_interface
            ! -----------------------------------------------------------------------
            ! Description:
            ! This function returns the constant time increment (dt) from the
            ! time increment controller object.
            !
            ! Parameters:
            ! self: time_increment_controller
            ! The time increment controller object.
            !
            ! Returns:
            ! dt: real
            ! The constant time increment.
            ! -----------------------------------------------------------------------
            use typedef_module
            import time_increment_controller

            class  (time_increment_controller), intent(in) :: self
            real   (real_kind                )             :: dt
        end function get_constant_dt_interface

        pure function compute_local_dt_interface(self, cell_volume, surface_area, spectral_radius) result(dt)
            ! """
            ! Compute the local time increment for the interface.
            ! Parameters:
            !     self (time_increment_controller): The time increment controller object.
            !     cell_volume (real): The volume of the cell.
            !     surface_area (real): The surface area of the cell.
            !     spectral_radius (real): The spectral radius.
            ! Returns:
            !     dt (real): The computed local time increment.
            ! """
            use typedef_module
            import time_increment_controller

            class  (time_increment_controller), intent(in) :: self
            real   (real_kind                ), intent(in) :: cell_volume
            real   (real_kind                ), intent(in) :: surface_area
            real   (real_kind                ), intent(in) :: spectral_radius
            real   (real_kind                )             :: dt
        end function compute_local_dt_interface

        pure function returns_constant_interface(self) result(ret)
            ! -----------------------------------------------------------------------
            ! Function Name: returns_constant_interface
            ! -----------------------------------------------------------------------
            ! Description:
            ! This function returns a logical value indicating whether the input
            ! object is a constant interface. The function is a pure function,
            ! meaning it does not modify any variables outside its scope.
            !
            ! Parameters:
            ! self: The input object of type time_increment_controller. It is
            ! passed by reference and is intended to be read-only.
            !
            ! Returns:
            ! ret: A logical value indicating whether the input object is a
            ! constant interface. If the object is a constant interface,
            ! ret is set to .TRUE., otherwise it is set to .FALSE.
            !
            ! -----------------------------------------------------------------------
            import time_increment_controller

            class  (time_increment_controller), intent(in) :: self
            logical                                        :: ret
        end function returns_constant_interface
    end interface
end module abstract_time_increment_controller