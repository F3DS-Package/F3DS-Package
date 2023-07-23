submodule(class_cellsystem) initialization_impl
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

    module subroutine initialize_eos(self, an_eos, config, num_conservative_variables)
        class  (cellsystem    ), intent(inout) :: self
        class  (eos           ), intent(inout) :: an_eos
        class  (configuration ), intent(inout) :: config
        integer(int_kind      ), intent(in   ) :: num_conservative_variables

#ifdef _DEBUG
        call write_debuginfo("In initialize_eos(), cellsystem.")
#endif
        call an_eos%initialize(config)
    end subroutine initialize_eos

    module subroutine initialize_time_increment_controller(self, controller, config, num_conservative_variables)
        class(cellsystem               ), intent(inout) :: self
        class(time_increment_controller), intent(inout) :: controller
        class(configuration            ), intent(inout) :: config
        integer(int_kind               ), intent(in   ) :: num_conservative_variables

#ifdef _DEBUG
        call write_debuginfo("In initialize_time_increment_controller(), cellsystem.")
#endif
        call controller%initialize(config)
    end subroutine initialize_time_increment_controller

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

    module subroutine initialize_gradient_calculator(self, a_gradient_calculator, config, num_conservative_variables)
        class  (cellsystem         ), intent(inout) :: self
        class  (gradient_calculator), intent(inout) :: a_gradient_calculator
        class  (configuration      ), intent(inout) :: config
        integer(int_kind           ), intent(in   ) :: num_conservative_variables

#ifdef _DEBUG
        call write_debuginfo("In initialize_gradient_calculator(), cellsystem.")
#endif
        call a_gradient_calculator%initialize(config)
    end subroutine initialize_gradient_calculator

    module subroutine initialize_interpolator(self, a_interpolator, config, num_conservative_variables)
        class  (cellsystem   ), intent(inout) :: self
        class  (interpolator ), intent(inout) :: a_interpolator
        class  (configuration), intent(inout) :: config
        integer(int_kind     ), intent(in   ) :: num_conservative_variables

#ifdef _DEBUG
        call write_debuginfo("In initialize_interpolator(), cellsystem.")
#endif
        call a_interpolator%initialize(config)
    end subroutine initialize_interpolator

    module subroutine initialize_time_stepping(self, a_time_stepping, config, num_conservative_variables)
        class  (cellsystem   ), intent(inout) :: self
        class  (time_stepping), intent(inout) :: a_time_stepping
        class  (configuration), intent(inout) :: config
        integer(int_kind     ), intent(in   ) :: num_conservative_variables

#ifdef _DEBUG
        call write_debuginfo("In initialize_time_stepping(), cellsystem.")
#endif
        call a_time_stepping%initialize(config, self%num_cells, num_conservative_variables)
    end subroutine initialize_time_stepping

    module subroutine initialize_face_gradient_interpolator(self, a_face_gradient_interpolator, config, num_conservative_variables)
        class   (cellsystem                ), intent(inout) :: self
        class   (face_gradient_interpolator), intent(inout) :: a_face_gradient_interpolator
        class   (configuration             ), intent(inout) :: config
        integer (int_kind                  ), intent(in   ) :: num_conservative_variables

#ifdef _DEBUG
        call write_debuginfo("In initialize_face_gradient_interpolator(), cellsystem.")
#endif
        call a_face_gradient_interpolator%initialize(config)
    end subroutine initialize_face_gradient_interpolator

    module subroutine initialize_face_gradient_calculator(self, a_face_gradient_calculator, config, num_conservative_variables)
        class   (cellsystem              ), intent(inout) :: self
        class   (face_gradient_calculator), intent(inout) :: a_face_gradient_calculator
        class   (configuration           ), intent(inout) :: config
        integer (int_kind                ), intent(in   ) :: num_conservative_variables

#ifdef _DEBUG
        call write_debuginfo("In initialize_face_gradient_calculator(), cellsystem.")
#endif
        call a_face_gradient_calculator%initialize(config)
    end subroutine initialize_face_gradient_calculator

    module subroutine initialize_riemann_solver(self, a_riemann_solver, config, num_conservative_variables)
        class  (cellsystem    ), intent(inout) :: self
        class  (riemann_solver), intent(inout) :: a_riemann_solver
        class  (configuration ), intent(inout) :: config
        integer(int_kind      ), intent(in   ) :: num_conservative_variables

#ifdef _DEBUG
        call write_debuginfo("In initialize_riemann_solver(), cellsystem.")
#endif
        call a_riemann_solver%initialize(config)
    end subroutine initialize_riemann_solver

    module subroutine initialize_reconstructor(self, a_reconstructor, config, num_conservative_variables, a_reconstructor_generator)
        class  (cellsystem    ), intent(inout) :: self
        class  (reconstructor ), intent(inout) :: a_reconstructor
        class  (configuration ), intent(inout) :: config
        integer(int_kind      ), intent(in   ) :: num_conservative_variables

        class  (reconstructor_generator), optional, intent(inout) :: a_reconstructor_generator

#ifdef _DEBUG
        call write_debuginfo("In initialize_reconstructor(), cellsystem.")
#endif
        call a_reconstructor%initialize(config, a_reconstructor_generator)
    end subroutine initialize_reconstructor

    module subroutine initialize_result_writer(self, writer, config, num_conservative_variables)
        class  (cellsystem   ), intent(inout) :: self
        class  (result_writer), intent(inout) :: writer
        class  (configuration), intent(inout) :: config
        integer(int_kind     ), intent(in   ) :: num_conservative_variables

#ifdef _DEBUG
        call write_debuginfo("In initialize_result_writer(), cellsystem.")
#endif
        call writer%initialize(self%num_cells, self%num_points, self%is_real_cell, self%cell_geometries, self%cell_types, config)
    end subroutine initialize_result_writer
end submodule