submodule(class_cellsystem) termination_criterion_impl
    implicit none

    contains

    module subroutine initialize_termination_ctiterion(self, criterion, config, num_conservative_variables)
        class  (cellsystem           ), intent(inout) :: self
        class  (termination_criterion), intent(inout) :: criterion
        class  (configuration        ), intent(inout) :: config
        integer(int_kind             ), intent(in   ) :: num_conservative_variables

#ifdef _DEBUG
        call write_debuginfo("In initialize_termination_ctiterion(), cellsystem.")
#endif
        call criterion%initialize(config)
    end subroutine initialize_termination_ctiterion

    module pure function satisfy_termination_criterion(self, criterion) result(judge)
        class(cellsystem           ), intent(in) :: self
        class(termination_criterion), intent(in) :: criterion
        logical :: judge
        judge = criterion%is_satisfied(self%time, self%num_steps)
    end function satisfy_termination_criterion
end submodule