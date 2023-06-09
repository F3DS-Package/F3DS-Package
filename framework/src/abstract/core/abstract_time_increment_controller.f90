module abstract_time_increment_controller
    implicit none

    private

    type, public, abstract :: time_increment_controller_generator
        contains
        generic, public :: generate => generate_from_configuration, generate_from_name
        procedure(generate_from_configuration_interface), pass(self), deferred :: generate_from_configuration
        procedure(generate_from_name_interface         ), pass(self), deferred :: generate_from_name
    end type time_increment_controller_generator

    type, public, abstract :: time_increment_controller
        contains

        procedure(initialize_interface      ), pass(self), deferred :: initialize
        procedure(get_constant_dt_interface ), pass(self), deferred :: get_constant_dt
        procedure(compute_local_dt_interface), pass(self), deferred :: compute_local_dt
        procedure(returns_constant_interface), pass(self), deferred :: returns_constant
    end type time_increment_controller

    abstract interface
        subroutine generate_from_configuration_interface(self, a_time_increment_controller, a_config)
            use abstract_configuration
            import time_increment_controller_generator
            import time_increment_controller
            class(time_increment_controller_generator),          intent(inout) :: self
            class(time_increment_controller          ), pointer, intent(inout) :: a_time_increment_controller
            class(configuration                      ),          intent(inout) :: a_config
        end subroutine generate_from_configuration_interface

        subroutine generate_from_name_interface(self, a_time_increment_controller, name)
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
            use abstract_configuration
            import time_increment_controller
            import time_increment_controller_generator

            class(time_increment_controller          ),           intent(inout) :: self
            class(configuration                      ),           intent(inout) :: config
        end subroutine initialize_interface

        function get_constant_dt_interface(self) result(dt)
            use typedef_module
            import time_increment_controller

            class  (time_increment_controller), intent(in) :: self
            real   (real_kind                )             :: dt
        end function get_constant_dt_interface

        pure function compute_local_dt_interface(self, cell_volume, surface_area, spectral_radius) result(dt)
            use typedef_module
            import time_increment_controller

            class  (time_increment_controller), intent(in) :: self
            real   (real_kind                ), intent(in) :: cell_volume
            real   (real_kind                ), intent(in) :: surface_area
            real   (real_kind                ), intent(in) :: spectral_radius
            real   (real_kind                )             :: dt
        end function compute_local_dt_interface

        pure function returns_constant_interface(self) result(ret)
            import time_increment_controller

            class  (time_increment_controller), intent(in) :: self
            logical                                        :: ret
        end function returns_constant_interface
    end interface
end module abstract_time_increment_controller