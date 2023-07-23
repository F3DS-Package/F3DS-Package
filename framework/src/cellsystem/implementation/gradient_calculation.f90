submodule(class_cellsystem) gradient_calculation_impl
    implicit none

    contains

    module subroutine compute_gradient_rank2(self, a_gradient_calculator, variables_set, gradient_variables_set, num_variables)
        class  (cellsystem         ), intent(inout) :: self
        class  (gradient_calculator), intent(inout) :: a_gradient_calculator
        real   (real_kind          ), intent(in   ) :: variables_set            (:,:)
        real   (real_kind          ), intent(inout) :: gradient_variables_set   (:,:)
        integer(int_kind           ), intent(in   ) :: num_variables

        integer(int_kind ) :: face_index, rhc_index, lhc_index, var_index

#ifdef _DEBUG
        call write_debuginfo("In compute_gradient_rank2(), cellsystem.")
#endif

!$omp parallel do private(face_index, rhc_index, lhc_index, var_index)
        do face_index = 1, self%num_faces, 1
            lhc_index = self%face_to_cell_indexes(self%num_local_cells - 0, face_index)
            rhc_index = self%face_to_cell_indexes(self%num_local_cells + 1, face_index)
            do var_index = 1, num_variables, 1
                associate(                                &
                    vec_start_index => 3*(var_index-1)+1, &
                    vec_end_index   => 3*(var_index-1)+3, &
                    residual        => a_gradient_calculator%compute_residual(                         &
                                           variables_set                      (var_index, lhc_index ), &
                                           variables_set                      (var_index, rhc_index ), &
                                           self%cell_centor_positions         (1:3      , lhc_index ), &
                                           self%cell_centor_positions         (1:3      , rhc_index ), &
                                           self%face_normal_vectors           (1:3      , face_index), &
                                           self%face_positions                (1:3      , face_index), &
                                           self%face_areas                               (face_index)  &
                                        )                                                              &
                )
                    gradient_variables_set(vec_start_index:vec_end_index, rhc_index) = gradient_variables_set(vec_start_index:vec_end_index, rhc_index)  &
                                                        - (1.d0 / self%cell_volumes(rhc_index)) * residual
                    gradient_variables_set(vec_start_index:vec_end_index, lhc_index) = gradient_variables_set(vec_start_index:vec_end_index, lhc_index)  &
                                                        + (1.d0 / self%cell_volumes(lhc_index)) * residual
                end associate
            end do
        end do
    end subroutine compute_gradient_rank2

    module subroutine compute_gradient_rank1(self, a_gradient_calculator, variable_set, gradient_variable_set)
        class(cellsystem         ), intent(inout) :: self
        class(gradient_calculator), intent(inout) :: a_gradient_calculator
        real (real_kind          ), intent(in   ) :: variable_set            (:)
        real (real_kind          ), intent(inout) :: gradient_variable_set   (:,:)

        integer(int_kind) :: face_index, rhc_index, lhc_index

#ifdef _DEBUG
        call write_debuginfo("In compute_gradient_rank1(), cellsystem.")
#endif

!$omp parallel do private(face_index, rhc_index, lhc_index)
        do face_index = 1, self%num_faces, 1
            lhc_index = self%face_to_cell_indexes(self%num_local_cells - 0, face_index)
            rhc_index = self%face_to_cell_indexes(self%num_local_cells + 1, face_index)

            associate(                                                                     &
                residual => a_gradient_calculator%compute_residual(                        &
                               variable_set                                  (lhc_index ), &
                               variable_set                                  (rhc_index ), &
                               self%cell_centor_positions         (1:3      , lhc_index ), &
                               self%cell_centor_positions         (1:3      , rhc_index ), &
                               self%face_normal_vectors           (1:3      , face_index), &
                               self%face_positions                (1:3      , face_index), &
                               self%face_areas                               (face_index)  &
                            )                                                              &
            )
                gradient_variable_set(1:3,rhc_index) = gradient_variable_set(1:3,rhc_index)                    &
                                                 - (1.d0 / self%cell_volumes(rhc_index)) * residual
                gradient_variable_set(1:3,lhc_index) = gradient_variable_set(1:3,lhc_index)                    &
                                                 + (1.d0 / self%cell_volumes(lhc_index)) * residual
            end associate
        end do
    end subroutine compute_gradient_rank1
end submodule