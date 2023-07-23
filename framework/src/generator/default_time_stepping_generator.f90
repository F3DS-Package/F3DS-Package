module default_time_stepping_generator_module
    use abstract_configuration
    use stdio_module
    use abstract_time_stepping
    use class_second_order_tvd_rk
    use class_third_order_tvd_rk
    use class_forward_euler

    implicit none

    private

    type, public, extends(time_stepping_generator) :: default_time_stepping_generator
        contains
        procedure, pass(self) :: generate_from_configuration
        procedure, pass(self) :: generate_from_name
    end type default_time_stepping_generator

    contains

    subroutine generate_from_configuration(self, a_time_stepping, a_config)
        class(default_time_stepping_generator),          intent(inout) :: self
        class(time_stepping                  ), pointer, intent(inout) :: a_time_stepping
        class(configuration                  )         , intent(inout) :: a_config

        logical          :: found
        character(len=:), allocatable :: name

        call a_config%get_char("Time stepping.Name", name, found)
        if(.not. found) call call_error("'Time stepping.Name' is not found in configuration file you set.")

        call self%generate_from_name(a_time_stepping, name)
    end subroutine generate_from_configuration

    subroutine generate_from_name(self, a_time_stepping, name)
        class    (default_time_stepping_generator),              intent(inout) :: self
        class    (time_stepping                  ), pointer    , intent(inout) :: a_time_stepping
        character(len=:                          ), allocatable, intent(in   ) :: name

        if(name == "2nd order TVD RK") then
            allocate(second_order_tvd_rk :: a_time_stepping)
        else if(name == "3rd order TVD RK") then
            allocate(third_order_tvd_rk :: a_time_stepping)
        else if(name == "Forward Euler") then
            allocate(forward_euler :: a_time_stepping)
        else
            call call_error("Unknow time stepping name '"//name//"' is found.")
        end if
    end subroutine generate_from_name
end module default_time_stepping_generator_module