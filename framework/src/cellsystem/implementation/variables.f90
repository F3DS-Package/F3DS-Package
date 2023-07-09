submodule (class_cellsystem) variables_impl
    implicit none

    contains

    module subroutine initialize_variables_rank2(self, variables_set, num_variables)
        class  (cellsystem), intent(inout)              :: self
        real   (real_kind ), intent(inout), allocatable :: variables_set(:,:)
        integer(int_kind  ), intent(in   )              :: num_variables

#ifdef _DEBUG
        call write_debuginfo("In initialize_variables_rank2(), cellsystem.")
#endif

        if(.not. self%read_cellsystem) call call_error("'read' subroutine is not called. You should call with following steps: first you call 'read' subroutine, next you initialze variables with 'initialze' subroutine. Please check your code.")

        if(allocated(variables_set))then
            call call_error("Array variables_set is allocated. But you call 'initialize' subroutine.")
        end if
        allocate(variables_set(1:num_variables, 1:self%num_cells))
        variables_set(:,:) = 0._real_kind
    end subroutine initialize_variables_rank2

    module subroutine initialize_variables_rank1(self, variables_set)
        class  (cellsystem), intent(inout)              :: self
        real   (real_kind ), intent(inout), allocatable :: variables_set(:)

#ifdef _DEBUG
        call write_debuginfo("In initialize_variables_rank1(), cellsystem.")
#endif

        if(.not. self%read_cellsystem) call call_error("'read' subroutine is not called. You should call with following steps: first you call 'read' subroutine, next you initialze variables with 'initialze' subroutine. Please check your code.")

        if(allocated(variables_set))then
            call call_error("Array variables_set is allocated. But you call 'initialize' subroutine.")
        end if
        allocate(variables_set(1:self%num_cells))
        variables_set(:) = 0.d0
    end subroutine initialize_variables_rank1

    module subroutine read_initial_condition_rank2(self, an_initial_condition_parser, config, conservative_variables_set)
        class  (cellsystem              ), intent(inout) :: self
        class  (initial_condition_parser), intent(inout) :: an_initial_condition_parser
        class  (configuration           ), intent(inout) :: config
        real   (real_kind               ), intent(inout) :: conservative_variables_set(:,:)

        character(len=:), allocatable :: filepath
        logical          :: found
        integer          :: i

#ifdef _DEBUG
        call write_debuginfo("In read_initial_condition_rank2(), cellsystem.")
#endif

        call an_initial_condition_parser%parse(config)
        call an_initial_condition_parser%get_conservative_variables_set(conservative_variables_set)
        call an_initial_condition_parser%close()
    end subroutine read_initial_condition_rank2

    module subroutine read_initial_condition_rank1(self, an_initial_condition_parser, config, conservative_variables_set)
        class  (cellsystem              ), intent(inout) :: self
        class  (initial_condition_parser), intent(inout) :: an_initial_condition_parser
        class  (configuration           ), intent(inout) :: config
        real   (real_kind               ), intent(inout) :: conservative_variables_set(:)

        character(len=:), allocatable :: filepath
        logical          :: found
        integer          :: i
        real(real_kind), allocatable :: temp_variables(:,:)

#ifdef _DEBUG
        call write_debuginfo("In read_initial_condition_rank1(), cellsystem.")
#endif

        allocate(temp_variables(1, size(conservative_variables_set)))

        call an_initial_condition_parser%parse(config)
        call an_initial_condition_parser%get_conservative_variables_set(temp_variables)
        call an_initial_condition_parser%close()

        conservative_variables_set(:) = temp_variables(1,:)
    end subroutine read_initial_condition_rank1

    module subroutine conservative_to_primitive_variables_all(self, an_eos, conservative_variables_set, primitive_variables_set, num_primitive_variables, conservative_to_primitive_function)
        class  (cellsystem              ), intent(inout) :: self
        class  (eos                     ), intent(in   ) :: an_eos
        real   (real_kind               ), intent(in   ) :: conservative_variables_set(:,:)
        real   (real_kind               ), intent(inout) :: primitive_variables_set(:,:)
        integer(int_kind                ), intent(in   ) :: num_primitive_variables

        procedure(conservative_to_primitive_function_interface) :: conservative_to_primitive_function

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In conservative_to_primitive_variables_all(), cellsystem.")
#endif

!$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            primitive_variables_set(:,i) = conservative_to_primitive_function(an_eos, conservative_variables_set(:,i), num_primitive_variables)
        end do
    end subroutine conservative_to_primitive_variables_all

    module subroutine operate_cellwise_rank2(self, variables_set, operator_subroutine)
        class  (cellsystem  ), intent(inout) :: self
        real   (real_kind   ), intent(inout) :: variables_set(:,:)

        procedure(operator_subroutine_rank1_interface) :: operator_subroutine

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In operate_cellwise_rank2(), cellsystem.")
#endif

!$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            call operator_subroutine(variables_set(:,i))
        end do
    end subroutine operate_cellwise_rank2

    module subroutine operate_cellwise_rank2_rank2(self, primary_variables_set, secondary_variables_set, operator_subroutine)
        class  (cellsystem  ), intent(inout) :: self
        real   (real_kind   ), intent(inout) :: primary_variables_set  (:,:)
        real   (real_kind   ), intent(in   ) :: secondary_variables_set(:,:)

        procedure(operator_subroutine_rank1_rank1_interface) :: operator_subroutine

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In operate_cellwise_rank2_rank2(), cellsystem.")
#endif

!$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            call operator_subroutine(primary_variables_set(:,i), secondary_variables_set(:,i))
        end do
    end subroutine operate_cellwise_rank2_rank2

    module subroutine operate_cellwise_rank2_rank1(self, primary_variables_set, secondary_variable_set, operator_subroutine)
        class  (cellsystem  ), intent(inout) :: self
        real   (real_kind   ), intent(inout) :: primary_variables_set  (:,:)
        real   (real_kind   ), intent(in   ) :: secondary_variable_set (:)

        procedure(operator_subroutine_rank1_rank0_interface) :: operator_subroutine

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In operate_cellwise_rank2_rank1(), cellsystem.")
#endif

!$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            call operator_subroutine(primary_variables_set(:,i), secondary_variable_set(i))
        end do
    end subroutine operate_cellwise_rank2_rank1

    module subroutine smooth_variables(self, variables_set, weight_function)
        class  (cellsystem  ), intent(inout) :: self
        real   (real_kind   ), intent(inout) :: variables_set(:,:)
        procedure(weight_function_interface) :: weight_function

        integer(int_kind ) :: face_index, cell_index, rhc_index, lhc_index
        real   (real_kind) :: smoothed_variables_set(size(variables_set(:,1)), size(variables_set(1,:))), total_weight_set(size(variables_set(1,:)))
        real   (real_kind) :: lhc_w, rhc_w

#ifdef _DEBUG
        call write_debuginfo("In smooth_variables(), cellsystem.")
#endif

!$omp parallel do private(cell_index)
        do cell_index = 1, self%num_cells, 1
            smoothed_variables_set(:,cell_index) = 0._real_kind
            total_weight_set        (cell_index) = 0._real_kind
        end do

!$omp parallel do private(face_index, rhc_index, lhc_index, lhc_w, rhc_w)
        do face_index = 1, self%num_faces, 1
            lhc_index = self%face_to_cell_indexes(self%num_local_cells - 0, face_index)
            rhc_index = self%face_to_cell_indexes(self%num_local_cells + 1, face_index)

            associate(                                               &
                x_lhc => self%cell_centor_positions(1:3, lhc_index), &
                x_rhc => self%cell_centor_positions(1:3, rhc_index), &
                v_lhc => variables_set             ( : , lhc_index), &
                v_rhc => variables_set             ( : , rhc_index)  &
            )
                lhc_w = weight_function(x_lhc, x_rhc, v_lhc, v_rhc)
                if(self%is_real_cell(rhc_index))then
                    rhc_w = weight_function(x_lhc, x_rhc, v_lhc, v_rhc)
                else
                    rhc_w = 0.d0
                endif

                smoothed_variables_set(:, lhc_index) = smoothed_variables_set(:, lhc_index) + lhc_w * v_lhc
                smoothed_variables_set(:, rhc_index) = smoothed_variables_set(:, rhc_index) + rhc_w * v_rhc

                total_weight_set(lhc_index) = total_weight_set(lhc_index) + lhc_w
                total_weight_set(rhc_index) = total_weight_set(rhc_index) + rhc_w
            end associate
        end do

!$omp parallel do private(cell_index)
        do cell_index = 1, self%num_cells, 1
            if(total_weight_set(cell_index) > 0._real_kind)then
                variables_set(:,cell_index) = smoothed_variables_set(:,cell_index) / total_weight_set(cell_index)
            else
                variables_set(:,cell_index) = 0._real_kind
            end if
        end do
    end subroutine smooth_variables

    module subroutine substitute_rank2(self, variables, val)
        class(cellsystem  ), intent(inout) :: self
        real (real_kind   ), intent(inout) :: variables(:,:)
        real (real_kind   ), intent(in   ) :: val

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In substitute_rank2(), cellsystem.")
#endif
!$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            variables(:,i) = val
        end do
    end subroutine substitute_rank2

    module subroutine substitute_zeros_rank2(self, variables)
        class(cellsystem  ), intent(inout) :: self
        real (real_kind   ), intent(inout) :: variables(:,:)

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In substitute_zeros_rank2(), cellsystem.")
#endif
        call self%substitute_rank2(variables, 0._real_kind)
    end subroutine substitute_zeros_rank2

    module subroutine substitute_rank1(self, variables, val)
        class(cellsystem  ), intent(inout) :: self
        real (real_kind   ), intent(inout) :: variables(:)
        real (real_kind   ), intent(in   ) :: val

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In substitute_rank1(), cellsystem.")
#endif
!$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            variables(i) = val
        end do
    end subroutine substitute_rank1

    module subroutine substitute_zeros_rank1(self, variables)
        class(cellsystem  ), intent(inout) :: self
        real (real_kind   ), intent(inout) :: variables(:)

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In substitute_zeros_rank1(), cellsystem.")
#endif
        call self%substitute_rank1(variables, 0._real_kind)
    end subroutine substitute_zeros_rank1
end submodule