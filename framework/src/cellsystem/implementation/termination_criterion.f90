submodule(class_cellsystem) termination_criterion_impl
    implicit none

    contains

    module pure function satisfy_termination_criterion(self, criterion) result(judge)
        class(cellsystem           ), intent(in) :: self
        class(termination_criterion), intent(in) :: criterion
        logical :: judge
        judge = criterion%is_satisfied(self%time, self%num_steps)
    end function satisfy_termination_criterion
end submodule