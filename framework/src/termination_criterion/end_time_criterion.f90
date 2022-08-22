module class_end_time_criterion
    use typedef_module
    use abstract_termination_criterion
    use abstract_configuration

    implicit none

    private

    type, public, extends(termination_criterion) :: end_time_criterion
        private

        integer(int_kind) :: end_time

        contains

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: is_satisfied
    end type end_time_criterion

    contains

    subroutine initialize(self, config)
        class(end_time_criterion) :: self
        class(configuration     ) :: config
        logical :: found

        call config%get_real("Termination criterion.End time", self%end_time, found)
        if(.not. found) call call_error("'Termination criterion.End time' is not found in configuration file you set.")
    end subroutine initialize

    pure function is_satisfied(self, time, num_steps) result(judge)
        class  (termination_criterion) :: self
        real   (real_kind            ) :: time
        integer(int_kind             ) :: num_steps
        logical                        :: judge

        judge = (time <= self%end_time)
    end function is_satisfied
end module class_end_time_criterion