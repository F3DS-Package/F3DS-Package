module abstract_termination_criterion
    implicit none

    type, public, abstract :: termination_criterion_generator
        contains
        generic, public :: generate => generate_from_configuration, generate_from_name
        procedure(generate_from_configuration_interface), pass(self), deferred :: generate_from_configuration
        procedure(generate_from_name_interface         ), pass(self), deferred :: generate_from_name
    end type termination_criterion_generator

    type, public, abstract :: termination_criterion
        contains

        procedure(initialize_interface  ), pass(self), deferred :: initialize
        procedure(is_satisfied_interface), pass(self), deferred :: is_satisfied
    end type termination_criterion

    abstract interface
        subroutine generate_from_configuration_interface(self, a_termination_criterion, a_config)
            use abstract_configuration
            import termination_criterion_generator
            import termination_criterion
            class(termination_criterion_generator), intent(inout) :: self
            class(termination_criterion          ), intent(inout) :: a_termination_criterion
            class(configuration                  ), intent(inout) :: a_config
        end subroutine generate_from_configuration_interface

        subroutine generate_from_name_interface(self, a_termination_criterion, name)
            use abstract_configuration
            import termination_criterion_generator
            import termination_criterion
            class    (termination_criterion_generator),              intent(inout) :: self
            class    (termination_criterion          ),              intent(inout) :: a_termination_criterion
            character(len=:                          ), allocatable, intent(in   ) :: name
        end subroutine generate_from_name_interface
    end interface

    abstract interface
        subroutine initialize_interface(self, config, a_termination_criterion_generator)
            use abstract_configuration
            import termination_criterion
            import termination_criterion_generator

            class(termination_criterion          ),           intent(inout) :: self
            class(configuration                  ),           intent(inout) :: config
            class(termination_criterion_generator), optional, intent(inout) :: a_termination_criterion_generator
        end subroutine initialize_interface

        pure function is_satisfied_interface(self, time, num_steps) result(judge)
            use typedef_module
            import termination_criterion

            class  (termination_criterion), intent(in) :: self
            real   (real_kind            ), intent(in) :: time
            integer(int_kind             ), intent(in) :: num_steps
            logical                                    :: judge
        end function is_satisfied_interface
    end interface
end module abstract_termination_criterion