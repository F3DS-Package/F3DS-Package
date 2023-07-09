submodule(class_cellsystem) boundary_condition_impl
    implicit none

    contains

    module subroutine apply_boundary_condition( &
        self, variables_set, num_variables ,    &
        compute_rotate_variables_function  ,    &
        compute_unrotate_variables_function,    &
        empty_condition_function           ,    &
        symmetric_condition_function       ,    &
        nonslip_wall_condition_function    ,    &
        slip_wall_condition_function       ,    &
        outflow_condition_function              &
    )
        class  (cellsystem  ), intent(inout) :: self
        real   (real_kind   ), intent(inout) :: variables_set(:,:)
        integer(int_kind    ), intent(in   ) :: num_variables

        procedure(compute_rotate_variables_function_interface  ) :: compute_rotate_variables_function
        procedure(compute_unrotate_variables_function_interface) :: compute_unrotate_variables_function

        procedure(boundary_condition_function_interface        ), optional :: empty_condition_function
        procedure(boundary_condition_function_interface        ), optional :: symmetric_condition_function
        procedure(boundary_condition_function_interface        ), optional :: nonslip_wall_condition_function
        procedure(boundary_condition_function_interface        ), optional :: slip_wall_condition_function
        procedure(boundary_condition_function_interface        ), optional :: outflow_condition_function

        integer :: i, face_index

#ifdef _DEBUG
        call write_debuginfo("In apply_boundary_condition(), cellsystem.")
#endif

        if (present(empty_condition_function)) call self%apply_empty_condition( &
            variables_set                                                     , &
            num_variables                                                     , &
            compute_rotate_variables_function                                 , &
            compute_unrotate_variables_function                               , &
            empty_condition_function                                            &
        )

        if (present(symmetric_condition_function)) call self%apply_symmetric_condition( &
            variables_set                                                             , &
            num_variables                                                             , &
            compute_rotate_variables_function                                         , &
            compute_unrotate_variables_function                                       , &
            symmetric_condition_function                                                &
        )

        if (present(nonslip_wall_condition_function)) call self%apply_nonslip_wall_condition( &
            variables_set                                                                   , &
            num_variables                                                                   , &
            compute_rotate_variables_function                                               , &
            compute_unrotate_variables_function                                             , &
            nonslip_wall_condition_function                                                   &
        )

        if (present(slip_wall_condition_function)) call self%apply_slip_wall_condition( &
            variables_set                                                             , &
            num_variables                                                             , &
            compute_rotate_variables_function                                         , &
            compute_unrotate_variables_function                                       , &
            slip_wall_condition_function                                                &
        )

        if (present(outflow_condition_function)) call self%apply_outflow_condition( &
            variables_set                                                         , &
            num_variables                                                         , &
            compute_rotate_variables_function                                     , &
            compute_unrotate_variables_function                                   , &
            outflow_condition_function                                              &
        )
    end subroutine apply_boundary_condition

    module subroutine apply_outflow_condition(self, variables_set, num_variables, &
        compute_rotate_variables_function, compute_unrotate_variables_function, boundary_condition_function)
        class  (cellsystem  ), intent(inout) :: self
        real   (real_kind   ), intent(inout) :: variables_set(:,:)
        integer(int_kind    ), intent(in   ) :: num_variables

        procedure(compute_rotate_variables_function_interface  ) :: compute_rotate_variables_function
        procedure(compute_unrotate_variables_function_interface) :: compute_unrotate_variables_function
        procedure(boundary_condition_function_interface        ) :: boundary_condition_function

        integer :: i, face_index

#ifdef _DEBUG
        call write_debuginfo("In apply_outflow_condition(), cellsystem.")
#endif

!$omp parallel do private(i, face_index)
        do i = 1, self%num_outflow_faces, 1
            face_index = self%outflow_face_indexes(i)
            call apply_boundary_condition_common_impl(          &
                variables_set                                 , &
                self%face_normal_vectors     (:,face_index)   , &
                self%face_tangential1_vectors(:,face_index)   , &
                self%face_tangential2_vectors(:,face_index)   , &
                self%face_to_cell_indexes    (:,face_index)   , &
                self%num_local_cells                          , &
                num_variables                                 , &
                compute_rotate_variables_function             , &
                compute_unrotate_variables_function           , &
                boundary_condition_function                     &
            )
        end do
    end subroutine apply_outflow_condition

    module subroutine apply_nonslip_wall_condition(self, variables_set, num_variables, &
        compute_rotate_variables_function, compute_unrotate_variables_function, boundary_condition_function)
        class  (cellsystem  ), intent(inout) :: self
        real   (real_kind   ), intent(inout) :: variables_set(:,:)
        integer(int_kind    ), intent(in   ) :: num_variables

        procedure(compute_rotate_variables_function_interface  ) :: compute_rotate_variables_function
        procedure(compute_unrotate_variables_function_interface) :: compute_unrotate_variables_function
        procedure(boundary_condition_function_interface        ) :: boundary_condition_function

        integer :: i, face_index

#ifdef _DEBUG
        call write_debuginfo("In apply_nonslip_wall_condition(), cellsystem.")
#endif

!$omp parallel do private(i, face_index)
        do i = 1, self%num_nonslip_wall_faces, 1
            face_index = self%nonslip_wall_face_indexes(i)
            call apply_boundary_condition_common_impl(         &
                variables_set                                , &
                self%face_normal_vectors     (:,face_index)  , &
                self%face_tangential1_vectors(:,face_index)  , &
                self%face_tangential2_vectors(:,face_index)  , &
                self%face_to_cell_indexes    (:,face_index)  , &
                self%num_local_cells                         , &
                num_variables                                , &
                compute_rotate_variables_function            , &
                compute_unrotate_variables_function          , &
                boundary_condition_function                    &
            )
        end do
    end subroutine apply_nonslip_wall_condition

    module subroutine apply_symmetric_condition(self, variables_set, num_variables, &
        compute_rotate_variables_function, compute_unrotate_variables_function, boundary_condition_function)
        class  (cellsystem  ), intent(inout) :: self
        real   (real_kind   ), intent(inout) :: variables_set(:,:)
        integer(int_kind    ), intent(in   ) :: num_variables

        procedure(compute_rotate_variables_function_interface  ) :: compute_rotate_variables_function
        procedure(compute_unrotate_variables_function_interface) :: compute_unrotate_variables_function
        procedure(boundary_condition_function_interface        ) :: boundary_condition_function

        integer :: i, face_index

#ifdef _DEBUG
        call write_debuginfo("In apply_slip_and_symmetric_condition(), cellsystem.")
#endif

!$omp parallel do private(i, face_index)
        do i = 1, self%num_symmetric_faces, 1
            face_index = self%symmetric_face_indexes(i)
            call apply_boundary_condition_common_impl(          &
                variables_set                                 , &
                self%face_normal_vectors     (:,face_index)   , &
                self%face_tangential1_vectors(:,face_index)   , &
                self%face_tangential2_vectors(:,face_index)   , &
                self%face_to_cell_indexes    (:,face_index)   , &
                self%num_local_cells                          , &
                num_variables                                 , &
                compute_rotate_variables_function             , &
                compute_unrotate_variables_function           , &
                boundary_condition_function                     &
            )
        end do
    end subroutine apply_symmetric_condition

    module subroutine apply_slip_wall_condition(self, variables_set, num_variables, &
        compute_rotate_variables_function, compute_unrotate_variables_function, boundary_condition_function)
        class  (cellsystem  ), intent(inout) :: self
        real   (real_kind   ), intent(inout) :: variables_set(:,:)
        integer(int_kind    ), intent(in   ) :: num_variables

        procedure(compute_rotate_variables_function_interface  ) :: compute_rotate_variables_function
        procedure(compute_unrotate_variables_function_interface) :: compute_unrotate_variables_function
        procedure(boundary_condition_function_interface        ) :: boundary_condition_function

        integer :: i, face_index

#ifdef _DEBUG
        call write_debuginfo("In apply_slip_wall_condition(), cellsystem.")
#endif

!$omp parallel do private(i, face_index)
        do i = 1, self%num_slip_wall_faces, 1
            face_index = self%slip_wall_face_indexes(i)
            call apply_boundary_condition_common_impl(          &
                variables_set                                 , &
                self%face_normal_vectors     (:,face_index)   , &
                self%face_tangential1_vectors(:,face_index)   , &
                self%face_tangential2_vectors(:,face_index)   , &
                self%face_to_cell_indexes    (:,face_index)   , &
                self%num_local_cells                          , &
                num_variables                                 , &
                compute_rotate_variables_function             , &
                compute_unrotate_variables_function           , &
                boundary_condition_function                     &
            )
        end do
    end subroutine apply_slip_wall_condition

    module subroutine apply_empty_condition(self, variables_set, num_variables, &
        compute_rotate_variables_function, compute_unrotate_variables_function, boundary_condition_function)
        class  (cellsystem  ), intent(inout) :: self
        real   (real_kind   ), intent(inout) :: variables_set(:,:)
        integer(int_kind    ), intent(in   ) :: num_variables

        procedure(compute_rotate_variables_function_interface  ) :: compute_rotate_variables_function
        procedure(compute_unrotate_variables_function_interface) :: compute_unrotate_variables_function
        procedure(boundary_condition_function_interface        ) :: boundary_condition_function

        integer :: i, face_index

#ifdef _DEBUG
        call write_debuginfo("In apply_empty_condition(), cellsystem.")
#endif

!$omp parallel do private(i, face_index)
        do i = 1, self%num_empty_faces, 1
            face_index = self%empty_face_indexes(i)
            call apply_boundary_condition_empty_impl(           &
                variables_set                                 , &
                self%face_normal_vectors     (:,face_index)   , &
                self%face_tangential1_vectors(:,face_index)   , &
                self%face_tangential2_vectors(:,face_index)   , &
                self%face_to_cell_indexes    (:,face_index)   , &
                self%num_local_cells                          , &
                num_variables                                 , &
                compute_rotate_variables_function             , &
                compute_unrotate_variables_function           , &
                boundary_condition_function                     &
            )
        end do
    end subroutine apply_empty_condition
end submodule