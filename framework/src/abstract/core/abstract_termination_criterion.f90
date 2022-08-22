module abstract_termination_criterion
    implicit none

    type, public, abstract :: termination_criterion
        contains

        procedure(initialize_interface  ), pass(self), deferred :: initialize
        procedure(is_satisfied_interface), pass(self), deferred :: is_satisfied
    end type termination_criterion

    abstract interface
        subroutine initialize_interface(self, config)
            use acstract_configuration
            import termination_criterion

            class(termination_criterion) :: self
            class(configuration        ) :: config
        end subroutine initialize_interface

        pure function is_satisfied_interface(self, time, num_steps) result(judge)
            use typedef_module
            import termination_criterion

            class  (termination_criterion) :: self
            real   (real_kind            ) :: time
            integer(int_kind             ) :: num_steps
            logical                        :: judge
        end function is_satisfied_interface
    end interface
end module abstract_termination_criterion