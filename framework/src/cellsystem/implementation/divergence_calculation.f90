submodule(class_cellsystem) divergence_calculation_impl
    implicit none

    contains

    module subroutine compute_divergence_rank2(self, a_interpolator, variables_set, divergence_variables_set, num_variables, gradient_variables_set, velosity_set)
        class  (cellsystem  ), intent(inout) :: self
        class  (interpolator), intent(inout) :: a_interpolator
        real   (real_kind   ), intent(in   ) :: variables_set           (:,:)
        real   (real_kind   ), intent(inout) :: divergence_variables_set(:,:)
        integer(int_kind    ), intent(in   ) :: num_variables
        real   (real_kind   ), intent(in   ), optional :: gradient_variables_set(:,:)
        real   (real_kind   ), intent(in   ), optional :: velosity_set          (:,:)

        integer(int_kind ) :: face_index, rhc_index, lhc_index, var_index, num_divergence_variables
        real   (real_kind) :: face_variables(num_variables)

#ifdef _DEBUG
        call write_debuginfo("In compute_divergence_rank2(), cellsystem.")
#endif

        num_divergence_variables = num_variables/3

!$omp parallel do private(face_index, rhc_index, lhc_index, var_index, face_variables)
        do face_index = 1, self%num_faces, 1
            lhc_index = self%face_to_cell_indexes(self%num_local_cells - 0, face_index)
            rhc_index = self%face_to_cell_indexes(self%num_local_cells + 1, face_index)

            face_variables(:) = a_interpolator%interpolate_face_variables(&
               variables_set                      (:   , :         ), &
               self%face_to_cell_indexes          (:   , face_index), &
               self%cell_centor_positions         (1:3 , :         ), &
               self%face_positions                (1:3 , face_index), &
               self%face_normal_vectors           (1:3 , face_index), &
               self%num_local_cells                                 , &
               num_variables                                        , &
               gradient_variables_set             (:   , :         ), &
               velosity_set                       (:   , :         )  &
            )

            do var_index = 1, num_divergence_variables, 1
                associate(face_variable=>face_variables(3*(var_index-1)+1:3*(var_index-1)+3))
                    divergence_variables_set(var_index, lhc_index) = divergence_variables_set(var_index, lhc_index) &
                                                        + (1.d0 / self%cell_volumes(lhc_index)) * vector_multiply(face_variable, self%face_normal_vectors(:, face_index) * self%face_areas(face_index))
                    divergence_variables_set(var_index, rhc_index) = divergence_variables_set(var_index, rhc_index) &
                                                        - (1.d0 / self%cell_volumes(rhc_index)) * vector_multiply(face_variable, self%face_normal_vectors(:, face_index) * self%face_areas(face_index))
                end associate
            end do
        end do
    end subroutine compute_divergence_rank2

    module subroutine compute_divergence_rank1(self, a_interpolator, variable_set, divergence_variable_set, gradient_variables_set, velosity_set)
        class(cellsystem   ), intent(inout) :: self
        class(interpolator ), intent(inout) :: a_interpolator
        real (real_kind    ), intent(in   ) :: variable_set              (:,:)
        real (real_kind    ), intent(inout) :: divergence_variable_set   (:)
        real   (real_kind   ), intent(in   ), optional :: gradient_variables_set(:,:)
        real   (real_kind   ), intent(in   ), optional :: velosity_set          (:,:)

        integer(int_kind ) :: face_index, rhc_index, lhc_index
        real   (real_kind) :: face_variables(3)

#ifdef _DEBUG
        call write_debuginfo("In compute_divergence_rank1(), cellsystem.")
#endif

!$omp parallel do private(face_index, rhc_index, lhc_index, face_variables)
        do face_index = 1, self%num_faces, 1
            lhc_index = self%face_to_cell_indexes(self%num_local_cells - 0, face_index)
            rhc_index = self%face_to_cell_indexes(self%num_local_cells + 1, face_index)

            face_variables(:) = a_interpolator%interpolate_face_variables( &
                variable_set(:,:), self%face_to_cell_indexes(:, face_index), self%cell_centor_positions(:,:), self%face_positions(:,face_index), self%face_normal_vectors(1:3, face_index), self%num_local_cells, 3, gradient_variables_set, velosity_set &
            )

            divergence_variable_set(lhc_index) = divergence_variable_set(lhc_index)               &
                                               + (1.d0 / self%cell_volumes(lhc_index)) * vector_multiply(face_variables(1:3), self%face_normal_vectors(1:3, face_index) * self%face_areas(face_index))
            divergence_variable_set(rhc_index) = divergence_variable_set(rhc_index)               &
                                               - (1.d0 / self%cell_volumes(rhc_index)) * vector_multiply(face_variables(1:3), self%face_normal_vectors(1:3, face_index) * self%face_areas(face_index))
        end do
    end subroutine compute_divergence_rank1

    module subroutine compute_divergence_godunov_rank2(self, a_reconstructor, a_riemann_solver, an_eos,                                         &
                                             primitive_variables_set, residual_set, num_conservative_variables, num_primitive_variables, &
                                             primitive_to_conservative_function, element_function)

        class  (cellsystem    ), intent(in   ) :: self
        class  (reconstructor ), intent(in   ) :: a_reconstructor
        class  (riemann_solver), intent(in   ) :: a_riemann_solver
        class  (eos           ), intent(in   ) :: an_eos
        real   (real_kind     ), intent(in   ) :: primitive_variables_set (:,:)
        real   (real_kind     ), intent(inout) :: residual_set            (:,:)
        integer(int_kind      ), intent(in   ) :: num_conservative_variables
        integer(int_kind      ), intent(in   ) :: num_primitive_variables

        procedure(primitive_to_conservative_function_interface   ) :: primitive_to_conservative_function
        procedure(element_provided_with_riemann_function_interface) :: element_function

        integer(int_kind ) :: i
        integer(int_kind ) :: lhc_index
        integer(int_kind ) :: rhc_index
        real   (real_kind) :: reconstructed_primitive_variables_lhc   (num_primitive_variables)
        real   (real_kind) :: reconstructed_primitive_variables_rhc   (num_primitive_variables)
        real   (real_kind) :: element                                 (num_conservative_variables, 1:2)

#ifdef _DEBUG
        call write_debuginfo("In compute_divergence_godunov_rank2 (), cellsystem.")
#endif

!$omp parallel do private(i,lhc_index,rhc_index,reconstructed_primitive_variables_lhc,reconstructed_primitive_variables_rhc,element)
        do i = 1, self%num_faces, 1
            lhc_index = self%face_to_cell_indexes(self%num_local_cells+0, i)
            rhc_index = self%face_to_cell_indexes(self%num_local_cells+1, i)

            reconstructed_primitive_variables_lhc(:) = a_reconstructor%reconstruct_lhc(                              &
                primitive_variables_set, self%face_to_cell_indexes, self%cell_centor_positions, self%face_positions, &
                i, self%num_local_cells, num_primitive_variables                                                     &
            )
            reconstructed_primitive_variables_rhc(:) = a_reconstructor%reconstruct_rhc(                              &
                primitive_variables_set, self%face_to_cell_indexes, self%cell_centor_positions, self%face_positions, &
                i, self%num_local_cells, num_primitive_variables                                                     &
            )

            element(:,:) = element_function(&
                an_eos                                      , &
                a_riemann_solver                            , &
                primitive_variables_set       (:,lhc_index) , &
                primitive_variables_set       (:,rhc_index) , &
                reconstructed_primitive_variables_lhc   (:) , &
                reconstructed_primitive_variables_rhc   (:) , &
                self%cell_volumes            (lhc_index)    , &
                self%cell_volumes            (rhc_index)    , &
                self%face_areas                  (i)        , &
                self%face_normal_vectors     (1:3,i)        , &
                self%face_tangential1_vectors(1:3,i)        , &
                self%face_tangential2_vectors(1:3,i)        , &
                num_conservative_variables                  , &
                num_primitive_variables                       &
            )

            residual_set(:,lhc_index) = residual_set(:,lhc_index) + element(:, 1)
            residual_set(:,rhc_index) = residual_set(:,rhc_index) + element(:, 2)
        end do
    end subroutine compute_divergence_godunov_rank2

    module subroutine compute_divergence_godunov_facegrad_rank2(self, a_reconstructor, a_riemann_solver, an_eos, a_face_gradient_interpolator, a_face_gradient_calculator,  &
                                        primitive_variables_set, gradient_primitive_variables_set, residual_set, num_conservative_variables, num_primitive_variables, &
                                        primitive_to_conservative_function, element_function)

        class  (cellsystem                ), intent(in   ) :: self
        class  (reconstructor             ), intent(in   ) :: a_reconstructor
        class  (riemann_solver            ), intent(in   ) :: a_riemann_solver
        class  (eos                       ), intent(in   ) :: an_eos
        class  (face_gradient_interpolator), intent(in   ) :: a_face_gradient_interpolator
        class  (face_gradient_calculator  ), intent(in   ) :: a_face_gradient_calculator
        real   (real_kind                 ), intent(in   ) :: primitive_variables_set          (:,:)
        real   (real_kind                 ), intent(in   ) :: gradient_primitive_variables_set (:,:)
        real   (real_kind                 ), intent(inout) :: residual_set                     (:,:)
        integer(int_kind                  ), intent(in   ) :: num_conservative_variables
        integer(int_kind                  ), intent(in   ) :: num_primitive_variables

        procedure(primitive_to_conservative_function_interface            ) :: primitive_to_conservative_function
        procedure(element_provided_with_riemann_facegrad_function_interface) :: element_function

        integer(int_kind ) :: i
        real   (real_kind) :: reconstructed_primitive_variables_lhc   (num_primitive_variables)
        real   (real_kind) :: reconstructed_primitive_variables_rhc   (num_primitive_variables)
        real   (real_kind) :: face_gradient_primitive_variables       (num_primitive_variables*3)
        real   (real_kind) :: element                                 (num_conservative_variables, 1:2)

#ifdef _DEBUG
        call write_debuginfo("In compute_divergence_godunov_facegrad_rank2(), cellsystem.")
#endif

!$omp parallel do private(i,reconstructed_primitive_variables_lhc,reconstructed_primitive_variables_rhc,face_gradient_primitive_variables,element)
        do i = 1, self%num_faces, 1
            associate(                                                             &
                lhc_index => self%face_to_cell_indexes(self%num_local_cells+0, i), &
                rhc_index => self%face_to_cell_indexes(self%num_local_cells+1, i)  &
            )
                reconstructed_primitive_variables_lhc(:) = a_reconstructor%reconstruct_lhc(                              &
                    primitive_variables_set, self%face_to_cell_indexes, self%cell_centor_positions, self%face_positions, &
                    i, self%num_local_cells, num_primitive_variables                                                     &
                )
                reconstructed_primitive_variables_rhc(:) = a_reconstructor%reconstruct_rhc(                              &
                    primitive_variables_set, self%face_to_cell_indexes, self%cell_centor_positions, self%face_positions, &
                    i, self%num_local_cells, num_primitive_variables                                                     &
                )

                if(.not. self%is_real_cell(rhc_index))then ! It's boundary face!
                    face_gradient_primitive_variables(:) = a_face_gradient_calculator%compute(            &
                        primitive_variables_set   (:,lhc_index), primitive_variables_set   (:,rhc_index), &
                        self%cell_centor_positions(:,lhc_index), self%cell_centor_positions(:,rhc_index), &
                        self%face_normal_vectors(:, i),                                                   &
                        num_primitive_variables                                                           &
                    )
                else
                    face_gradient_primitive_variables(:) = a_face_gradient_interpolator%interpolate(                  &
                        gradient_primitive_variables_set(:,lhc_index), gradient_primitive_variables_set(:,rhc_index), &
                        primitive_variables_set         (:,lhc_index), primitive_variables_set         (:,rhc_index), &
                        self%cell_centor_positions      (:,lhc_index), self%cell_centor_positions      (:,rhc_index), &
                        num_primitive_variables                                                                       &
                    )
                end if

                element(:,:) = element_function(        &
                    an_eos                                              , &
                    a_riemann_solver                                    , &
                    primitive_variables_set              (:,lhc_index)  , &
                    primitive_variables_set              (:,rhc_index)  , &
                    reconstructed_primitive_variables_lhc(:)            , &
                    reconstructed_primitive_variables_rhc(:)            , &
                    face_gradient_primitive_variables    (:)            , &
                    self%cell_volumes                      (lhc_index)  , &
                    self%cell_volumes                      (rhc_index)  , &
                    self%face_areas                  (i)                , &
                    self%face_normal_vectors     (1:3,i)                , &
                    self%face_tangential1_vectors(1:3,i)                , &
                    self%face_tangential2_vectors(1:3,i)                , &
                    num_conservative_variables                          , &
                    num_primitive_variables                               &
                )

                residual_set(:,lhc_index) = residual_set(:,lhc_index) + element(:, 1)
                residual_set(:,rhc_index) = residual_set(:,rhc_index) + element(:, 2)
            end associate
        end do
    end subroutine compute_divergence_godunov_facegrad_rank2
end submodule