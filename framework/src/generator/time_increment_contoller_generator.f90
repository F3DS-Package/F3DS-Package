module time_increment_controller_generator_module
    use abstract_configuration
    use stdio_module
    use abstract_time_increment_controller
    use class_constant_time_increment_controller
    use class_adaptive_time_increment_controller

    implicit none

    private

    public :: f3ds_time_increment_controller_generator

    interface f3ds_time_increment_controller_generator
        module procedure with_config, with_name
    end interface f3ds_time_increment_controller_generator

    contains

    subroutine with_config(config, controller)
        class(configuration            )         , intent(inout) :: config
        class(time_increment_controller), pointer, intent(inout) :: controller

        logical          :: found
        character(len=:), allocatable :: name

        call config%get_char("Time incriment control.Controller name", name, found)
        if(.not. found) call call_error("'Time incriment control.Controller name' is not found in configuration file you set.")

        call with_name(name, controller)
    end subroutine with_config

    subroutine with_name(name, controller)
        character(len=:                    ), allocatable, intent(in   ) :: name
        class    (time_increment_controller), pointer    , intent(inout) :: controller

        if(name == "Adaptive controller") then
            allocate(adaptive_time_increment_controller :: controller)
        else if(name == "Constant controller") then
            allocate(constant_time_increment_controller :: controller)
        else
            call call_error("Unknow time increment controller name '"//name//"' is found.")
        end if
    end subroutine with_name
end module time_increment_controller_generator_module