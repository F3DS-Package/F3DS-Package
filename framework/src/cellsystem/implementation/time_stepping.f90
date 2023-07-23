submodule(class_cellsystem) time_stepping_impl
    implicit none

    contains

    module subroutine compute_next_stage_primitive_rank2(self, a_time_stepping, an_eos, stage_num, conservative_variables_set, primitive_variables_set, residual_set, num_primitive_variables, &
        conservative_to_primitive_function)

        class  (cellsystem         ), intent(inout) :: self
        class  (time_stepping      ), intent(inout) :: a_time_stepping
        class  (eos                ), intent(in   ) :: an_eos
        integer(int_kind           ), intent(in   ) :: stage_num
        real   (real_kind          ), intent(inout) :: conservative_variables_set(:,:)
        real   (real_kind          ), intent(inout) :: primitive_variables_set   (:,:)
        real   (real_kind          ), intent(inout) :: residual_set              (:,:)
        integer(int_kind           ), intent(in   ) :: num_primitive_variables

        procedure(conservative_to_primitive_function_interface) :: conservative_to_primitive_function

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In compute_next_stage_primitive_rank2(), cellsystem.")
#endif

!$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            conservative_variables_set(:,i) = a_time_stepping%compute_next_stage(i, stage_num, self%time_increment, conservative_variables_set(:,i), residual_set(:,i))
            residual_set              (:,i) = 0.d0
            primitive_variables_set   (:,i) = conservative_to_primitive_function(an_eos, conservative_variables_set(:,i), num_primitive_variables)
        end do
    end subroutine compute_next_stage_primitive_rank2

    module subroutine compute_next_stage_rank2(self, a_time_stepping, stage_num, conservative_variables_set, residual_set)
        class  (cellsystem         ), intent(inout) :: self
        class  (time_stepping      ), intent(inout) :: a_time_stepping
        integer(int_kind           ), intent(in   ) :: stage_num
        real   (real_kind          ), intent(inout) :: conservative_variables_set(:,:)
        real   (real_kind          ), intent(inout) :: residual_set              (:,:)

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In compute_next_stage_rank2(), cellsystem.")
#endif

!$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            conservative_variables_set(:,i) = a_time_stepping%compute_next_stage(i, stage_num, self%time_increment, conservative_variables_set(:,i), residual_set(:,i))
            residual_set              (:,i) = 0.d0
        end do
    end subroutine compute_next_stage_rank2

    module subroutine compute_next_stage_rank1(self, a_time_stepping, stage_num, conservative_variable_set, residual_set)

        class  (cellsystem         ), intent(inout) :: self
        class  (time_stepping      ), intent(inout) :: a_time_stepping
        integer(int_kind           ), intent(in   ) :: stage_num
        real   (real_kind          ), intent(inout) :: conservative_variable_set(:)
        real   (real_kind          ), intent(inout) :: residual_set             (:)

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In compute_next_stage_rank1(), cellsystem.")
#endif

!$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            ! Fortran cannot implicitly convert a 1-dimensional array of 1 elements to a scalar.
            ! In the following, the built-in function sum is used to convert to a scalar.
            conservative_variable_set(i) = sum( a_time_stepping%compute_next_stage(i, stage_num, self%time_increment, [conservative_variable_set(i)], [residual_set(i)]) )
            residual_set             (i) = 0.d0
        end do
    end subroutine compute_next_stage_rank1

    module subroutine prepare_time_stepping_rank2(    &
        self                      , &
        a_time_stepping           , &
        conservative_variables_set, &
        residual_set                  )
        class(cellsystem   ), intent(inout) :: self
        class(time_stepping), intent(inout) :: a_time_stepping
        real (real_kind    ), intent(inout) :: conservative_variables_set(:,:)
        real (real_kind    ), intent(inout) :: residual_set              (:,:)

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In prepare_time_stepping_rank2(), cellsystem.")
#endif

!$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            call a_time_stepping%prepare_time_stepping(i, conservative_variables_set(:,i), residual_set(:,i))
        end do
    end subroutine

    module subroutine prepare_time_stepping_rank1(    &
        self                      , &
        a_time_stepping           , &
        conservative_variables_set, &
        residual_set                  )
        class(cellsystem   ), intent(inout) :: self
        class(time_stepping), intent(inout) :: a_time_stepping
        real (real_kind    ), intent(inout) :: conservative_variables_set(:)
        real (real_kind    ), intent(inout) :: residual_set              (:)

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In prepare_time_stepping_rank1(), cellsystem.")
#endif

!$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            call a_time_stepping%prepare_time_stepping(i, [conservative_variables_set(i)], [residual_set(i)])
        end do
    end subroutine

    module pure function get_number_of_stages(self, a_time_stepping) result(n)
        class  (cellsystem   ), intent(in) :: self
        class  (time_stepping), intent(in) :: a_time_stepping
        integer(int_kind     )                :: n
        n = a_time_stepping%get_number_of_stages()
    end function get_number_of_stages
end submodule