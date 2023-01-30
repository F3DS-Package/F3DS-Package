module class_minmod_muscl3
    !**
    !* 3rd order TVD-MUSCL with minmod limmiter
    !*

    use typedef_module
    use abstract_reconstructor
    use abstract_configuration
    use stdio_module

    implicit none

    private

    type, public, extends(reconstructor) :: minmod_muscl3
        private

        real   (real_kind) :: kappa_ = 1.d0 / 3.d0

        contains

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: reconstruct_lhc
        procedure, public, pass(self) :: reconstruct_rhc

        procedure, pass(self) :: minmod
    end type minmod_muscl3

    contains

    pure function minmod(self, s) result(r)
        class(minmod_muscl3), intent(in) :: self
        real (real_kind    ), intent(in) :: s
        real (real_kind    )             :: r
        r = max(0.d0, min(1.d0, s))
    end function

    subroutine initialize(self, config)
        class(minmod_muscl3), intent(inout) :: self
        class(configuration), intent(inout) :: config
        logical :: found

        call config%get_real("Reconstructor.Minmod MUSCL3.Kappa", self%kappa_, found, 1.d0 / 3.d0)
        if(.not. found) call write_warring("'Reconstructor.Minmod MUSCL3.Kappa' is not found in configuration you set. To be set dafault value.")
    end subroutine initialize

    pure function reconstruct_lhc( &
        self                     , &
        primitive_variables_set  , &
        face_to_cell_index       , &
        cell_centor_positions    , &
        face_centor_positions    , &
        face_index               , &
        num_local_cells          , &
        num_primitive_variables      ) result(reconstructed_primitive_variables)

        class  (minmod_muscl3), intent(in) :: self
        integer(int_kind     ), intent(in) :: face_index
        real   (real_kind    ), intent(in) :: primitive_variables_set          (:, :)
        integer(int_kind     ), intent(in) :: face_to_cell_index               (:, :)
        real   (real_kind    ), intent(in) :: cell_centor_positions            (:, :)
        real   (real_kind    ), intent(in) :: face_centor_positions            (:, :)
        integer(int_kind     ), intent(in) :: num_local_cells
        integer(int_kind     ), intent(in) :: num_primitive_variables
        real   (real_kind    )             :: reconstructed_primitive_variables(num_primitive_variables)

        integer(int_kind ) :: i
        real   (real_kind) :: w(3), p(3), cell_pos_l(3), cell_pos_r(3), face_pos(3)
        real   (real_kind) :: cell_cell_distance, cell_face_distanse, s

        do i = 1, num_primitive_variables, 1
            associate(                                                                                       &
                v_m1 => primitive_variables_set(i, face_to_cell_index(num_local_cells-1, face_index)), &
                v    => primitive_variables_set(i, face_to_cell_index(num_local_cells+0, face_index)), &
                v_p1 => primitive_variables_set(i, face_to_cell_index(num_local_cells+1, face_index))  &
            )
                if((v - v_m1) == 0.d0 .or. (v_p1 - v) == 0.d0)then
                    reconstructed_primitive_variables(i) = v
                else
                    s = (v - v_m1) / (v_p1 - v)
                    reconstructed_primitive_variables(i) = v &
                                            + 0.25d0 * (1.d0 - self%kappa_) * self%minmod(1.d0 / s) * (v    - v_m1) &
                                            + 0.25d0 * (1.d0 - self%kappa_) * self%minmod(       s) * (v_p1 - v   )
                end if
            end associate
        end do
    end function reconstruct_lhc

    pure function reconstruct_rhc( &
        self                     , &
        primitive_variables_set  , &
        face_to_cell_index       , &
        cell_centor_positions    , &
        face_centor_positions    , &
        face_index               , &
        num_local_cells          , &
        num_primitive_variables      ) result(reconstructed_primitive_variables)

        class  (minmod_muscl3), intent(in) :: self
        integer(int_kind     ), intent(in) :: face_index
        real   (real_kind    ), intent(in) :: primitive_variables_set          (:, :)
        integer(int_kind     ), intent(in) :: face_to_cell_index               (:, :)
        real   (real_kind    ), intent(in) :: cell_centor_positions            (:, :)
        real   (real_kind    ), intent(in) :: face_centor_positions            (:, :)
        integer(int_kind     ), intent(in) :: num_local_cells
        integer(int_kind     ), intent(in) :: num_primitive_variables
        real   (real_kind    )             :: reconstructed_primitive_variables(num_primitive_variables)

        integer(int_kind ) :: i
        real   (real_kind) :: w(3), p(3), cell_pos_l(3), cell_pos_r(3), face_pos(3)
        real   (real_kind) :: cell_cell_distance, cell_face_distanse, s

        do i = 1, num_primitive_variables, 1
            associate(                                                                                 &
                v_m1 => primitive_variables_set(i, face_to_cell_index(num_local_cells+0, face_index)), &
                v    => primitive_variables_set(i, face_to_cell_index(num_local_cells+1, face_index)), &
                v_p1 => primitive_variables_set(i, face_to_cell_index(num_local_cells+2, face_index))  &
            )
                if((v - v_m1) == 0.d0 .or. (v_p1 - v) == 0.d0)then
                    reconstructed_primitive_variables(i) = v
                else
                    s = (v - v_m1) / (v_p1 - v)
                    reconstructed_primitive_variables(i) = v &
                                               - 0.25d0 * (1.d0 - self%kappa_) * self%minmod(1.d0 / s) * (v    - v_m1) &
                                               - 0.25d0 * (1.d0 - self%kappa_) * self%minmod(       s) * (v_p1 - v   )
                end if
            end associate
        end do
    end function reconstruct_rhc
end module class_minmod_muscl3