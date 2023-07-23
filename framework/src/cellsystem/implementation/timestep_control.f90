submodule(class_cellsystem) timestep_control_impl
    implicit none

    contains

    module subroutine update_time_increment_eos_rank2(self, controller, an_eos, variables_set, spectral_radius_function)
        class  (cellsystem               ), intent(inout) :: self
        class  (time_increment_controller), intent(in   ) :: controller
        class  (eos                      ), intent(in   ) :: an_eos
        real   (real_kind                ), intent(in   ) :: variables_set(:,:)

        procedure(spectral_radius_function_eos_rank1_interface) :: spectral_radius_function

        integer(int_kind ) :: i

#ifdef _DEBUG
        call write_debuginfo("In update_time_increment_eos_rank2(), cellsystem.")
#endif

        if ( controller%returns_constant() ) then
            self%time_increment = controller%get_constant_dt()
            return
        end if

        self%time_increment = large_value
        do i = 1, self%num_faces, 1
            associate(                                                                                                             &
                lhc_q        => variables_set    (:, self%face_to_cell_indexes(self%num_local_cells+0, i)), &
                rhc_q        => variables_set    (:, self%face_to_cell_indexes(self%num_local_cells+1, i)), &
                lhc_v        => self%cell_volumes   (self%face_to_cell_indexes(self%num_local_cells+0, i)), &
                rhc_v        => self%cell_volumes   (self%face_to_cell_indexes(self%num_local_cells+1, i)), &
                s            => self%face_areas                                                       (i) , &
                lhc_is_real  => self%is_real_cell   (self%face_to_cell_indexes(self%num_local_cells+0, i)), &
                rhc_is_real  => self%is_real_cell   (self%face_to_cell_indexes(self%num_local_cells+1, i))  &
            )
                if(lhc_is_real)then
                    self%time_increment = min(controller%compute_local_dt(lhc_v, s, spectral_radius_function(an_eos, lhc_q, lhc_v / s)), self%time_increment)
                end if
                if(rhc_is_real)then
                    self%time_increment = min(controller%compute_local_dt(rhc_v, s, spectral_radius_function(an_eos, rhc_q, rhc_v / s)), self%time_increment)
                end if
            end associate
        end do
    end subroutine update_time_increment_eos_rank2

    module subroutine update_time_increment_rank1(self, controller, variables_set, spectral_radius_function)
        class  (cellsystem               ), intent(inout) :: self
        class  (time_increment_controller), intent(in   ) :: controller
        real   (real_kind                ), intent(in   ) :: variables_set(:)

        procedure(spectral_radius_function_rank0_interface) :: spectral_radius_function

        integer(int_kind ) :: i

#ifdef _DEBUG
        call write_debuginfo("In update_time_increment_eos_rank2(), cellsystem.")
#endif

        if ( controller%returns_constant() ) then
            self%time_increment = controller%get_constant_dt()
            return
        end if

        self%time_increment = large_value
        do i = 1, self%num_faces, 1
            associate(                                                                                                             &
                lhc_q        => variables_set    (self%face_to_cell_indexes(self%num_local_cells+0, i)), &
                rhc_q        => variables_set    (self%face_to_cell_indexes(self%num_local_cells+1, i)), &
                lhc_v        => self%cell_volumes(self%face_to_cell_indexes(self%num_local_cells+0, i)), &
                rhc_v        => self%cell_volumes(self%face_to_cell_indexes(self%num_local_cells+1, i)), &
                s            => self%face_areas                                                    (i) , &
                lhc_is_real  => self%is_real_cell(self%face_to_cell_indexes(self%num_local_cells+0, i)), &
                rhc_is_real  => self%is_real_cell(self%face_to_cell_indexes(self%num_local_cells+1, i))  &
            )
                if(lhc_is_real)then
                    self%time_increment = min(controller%compute_local_dt(lhc_v, s, spectral_radius_function(lhc_q, lhc_v / s)), self%time_increment)
                end if
                if(rhc_is_real)then
                    self%time_increment = min(controller%compute_local_dt(rhc_v, s, spectral_radius_function(rhc_q, rhc_v / s)), self%time_increment)
                end if
            end associate
        end do
    end subroutine update_time_increment_rank1
end submodule