module default_time_increment_controller_generator_module
    use abstract_configuration
    use stdio_module
    use abstract_time_increment_controller
    use class_constant_time_increment_controller
    use class_adaptive_time_increment_controller

    implicit none

    private

    type, public, extends(time_increment_controller_generator) :: default_time_increment_controller_generator
        contains
        procedure, pass(self) :: generate_from_configuration
        procedure, pass(self) :: generate_from_name
    end type default_time_increment_controller_generator

    contains

    subroutine generate_from_configuration(self, a_time_increment_controller, a_config)
        class(default_time_increment_controller_generator),          intent(inout) :: self
        class(time_increment_controller                  ), pointer, intent(inout) :: a_time_increment_controller
        class(configuration                              )         , intent(inout) :: a_config

        logical          :: found
        character(len=:), allocatable :: name

        call a_config%get_char("Time increment control.Controller name", name, found)
        if(.not. found) call call_error("'Time increment control.Controller name' is not found in configuration file you set.")

        call self%generate_from_name(a_time_increment_controller, name)
    end subroutine generate_from_configuration

    subroutine generate_from_name(self, a_time_increment_controller, name)
        class    (default_time_increment_controller_generator),              intent(inout) :: self
        class    (time_increment_controller                  ), pointer    , intent(inout) :: a_time_increment_controller
        character(len=:                                      ), allocatable, intent(in   ) :: name

        if(name == "Adaptive controller") then
            allocate(adaptive_time_increment_controller :: a_time_increment_controller)
        else if(name == "Constant controller") then
            allocate(constant_time_increment_controller :: a_time_increment_controller)
        else
            call call_error("Unknow time increment controller name '"//name//"' is found.")
        end if
    end subroutine generate_from_name
end module default_time_increment_controller_generator_module