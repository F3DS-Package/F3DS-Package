module time_stepping_generator_module
    use abstract_configuration
    use stdio_module
    use abstract_time_stepping
    use class_second_order_tvd_rk
    use class_third_order_tvd_rk
    use class_forward_euler

    implicit none

    private

    public :: f3ds_time_stepping_generator

    interface f3ds_time_stepping_generator
        module procedure with_config, with_name
    end interface f3ds_time_stepping_generator

    contains

    subroutine with_config(config, method)
        class(configuration), intent(inout) :: config
        class(time_stepping), pointer, intent(inout) :: method

        logical          :: found
        character(len=:), allocatable :: name

        call config%get_char("Time stepping.Name", name, found)
        if(.not. found) call call_error("'Time stepping.Name' is not found in configuration file you set.")

        call with_name(name, method)
    end subroutine with_config

    subroutine with_name(name, method)
        character(len=:), allocatable, intent(in) :: name
        class(time_stepping), pointer, intent(inout) :: method

        if(name == "2nd order TVD RK") then
            allocate(second_order_tvd_rk :: method)
        else if(name == "3rd order TVD RK") then
            allocate(third_order_tvd_rk :: method)
        else if(name == "Forward Euler") then
            allocate(forward_euler :: method)
        else
            call call_error("Unknow time stepping name '"//name//"' is found.")
        end if
    end subroutine with_name
end module time_stepping_generator_module