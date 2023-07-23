module class_mp
    use typedef_module
    use abstract_reconstructor
    use abstract_configuration
    use weno_utils
    use stdio_module

    implicit none

    private

    type, public, extends(reconstructor) :: mp
        private

        class(reconstructor), pointer :: primary_reconstructor_

        real(real_kind) :: alpha_
        real(real_kind) :: beta_

        contains

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: reconstruct_lhc
        procedure, public, pass(self) :: reconstruct_rhc

        procedure, pass(self) :: minmod
        procedure, pass(self) :: median
        procedure, pass(self) :: curvature_measure
    end type mp

    contains

    subroutine initialize(self, config, a_reconstructor_generator)
        class(mp                     ),           intent(inout) :: self
        class(configuration          ),           intent(inout) :: config
        class(reconstructor_generator), optional, intent(inout) :: a_reconstructor_generator
        logical :: found
        character(len=:), allocatable :: name

        call config%get_real("Reconstructor.MP.Alpha", self%alpha_, found, 10.d0)
        if(.not. found) call write_warring("'Reconstructor.MP.Alpha' is not found in configration you set. To be set dafault value.")

        call config%get_real("Reconstructor.MP.Beta", self%beta_, found, 4.d0)
        if(.not. found) call write_warring("'Reconstructor.MP.Beta' is not found in configration you set. To be set dafault value.")

        call config%get_char("Reconstructor.MP.Primary reconstructor", name, found, "WENO5-JS")
        if(.not. found) call write_warring("'Reconstructor.MP.Primary reconstructor' is not found in configuration you set. To be set dafault value.")

        if(.not. present(a_reconstructor_generator))then
            call call_error("MP requests a reconstructor-generator.")
        end if

        call a_reconstructor_generator%generate(self%primary_reconstructor_, name)

        call self%primary_reconstructor_%initialize(config, a_reconstructor_generator)
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

        class  (mp         ), intent(in) :: self
        integer(int_kind   ), intent(in) :: face_index
        real   (real_kind  ), intent(in) :: primitive_variables_set          (:, :)
        integer(int_kind   ), intent(in) :: face_to_cell_index               (:, :)
        real   (real_kind  ), intent(in) :: cell_centor_positions            (:, :)
        real   (real_kind  ), intent(in) :: face_centor_positions            (:, :)
        integer(int_kind   ), intent(in) :: num_local_cells
        integer(int_kind   ), intent(in) :: num_primitive_variables
        real   (real_kind  )             :: reconstructed_primitive_variables(num_primitive_variables)

        integer(int_kind   ) :: i

        real   (real_kind  ) :: d_p1, d, d_m1, dm4
        real   (real_kind  ) :: v_ul, v_md, v_lc, v_min, v_max

        reconstructed_primitive_variables(:) = self%primary_reconstructor_%reconstruct_lhc( &
            primitive_variables_set  , &
            face_to_cell_index       , &
            cell_centor_positions    , &
            face_centor_positions    , &
            face_index               , &
            num_local_cells          , &
            num_primitive_variables    &
        )

        do i = 1, num_primitive_variables, 1
            associate(                                                                                 &
                v_m2 => primitive_variables_set(i, face_to_cell_index(num_local_cells-2, face_index)), &
                v_m1 => primitive_variables_set(i, face_to_cell_index(num_local_cells-1, face_index)), &
                v    => primitive_variables_set(i, face_to_cell_index(num_local_cells+0, face_index)), &
                v_p1 => primitive_variables_set(i, face_to_cell_index(num_local_cells+1, face_index)), &
                v_p2 => primitive_variables_set(i, face_to_cell_index(num_local_cells+2, face_index))  &
            )
                d_m1 = self%curvature_measure(v_m2, v_m1, v   )
                d    = self%curvature_measure(v_m1, v   , v_p1)
                d_p1 = self%curvature_measure(v   , v_p1, v_p2)
                dm4  = self%minmod(self%minmod(4.d0 * d - d_p1, 4.d0 * d_p1 - d), self%minmod(d, d_p1))
                v_ul = v + self%alpha_ * (v - v_m1)
                v_md = 0.5d0 * (v + v_p1 - dm4)
                v_lc = v + 0.5d0 * (v - v_m1) + (self%beta_ / 3.d0) * dm4
                v_min = max(min(v, v_p1, v_md), min(v, v_ul, v_lc))
                v_max = min(max(v, v_p1, v_md), max(v, v_ul, v_lc))
                reconstructed_primitive_variables(i) = self%median(reconstructed_primitive_variables(i), v_min, v_max)
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

        class  (mp         ), intent(in) :: self
        integer(int_kind   ), intent(in) :: face_index
        real   (real_kind  ), intent(in) :: primitive_variables_set          (:, :)
        integer(int_kind   ), intent(in) :: face_to_cell_index               (:, :)
        real   (real_kind  ), intent(in) :: cell_centor_positions            (:, :)
        real   (real_kind  ), intent(in) :: face_centor_positions            (:, :)
        integer(int_kind   ), intent(in) :: num_local_cells
        integer(int_kind   ), intent(in) :: num_primitive_variables
        real   (real_kind  )             :: reconstructed_primitive_variables(num_primitive_variables)

        integer(int_kind   ) :: i

        real   (real_kind  ) :: d_p1, d, d_m1, dm4
        real   (real_kind  ) :: v_ul, v_md, v_lc, v_min, v_max

        reconstructed_primitive_variables(:) = self%primary_reconstructor_%reconstruct_rhc( &
            primitive_variables_set  , &
            face_to_cell_index       , &
            cell_centor_positions    , &
            face_centor_positions    , &
            face_index               , &
            num_local_cells          , &
            num_primitive_variables    &
        )

        do i = 1, num_primitive_variables, 1
            associate(                                                                                 &
                v_m2 => primitive_variables_set(i, face_to_cell_index(num_local_cells-1, face_index)), &
                v_m1 => primitive_variables_set(i, face_to_cell_index(num_local_cells+0, face_index)), &
                v    => primitive_variables_set(i, face_to_cell_index(num_local_cells+1, face_index)), &
                v_p1 => primitive_variables_set(i, face_to_cell_index(num_local_cells+2, face_index)), &
                v_p2 => primitive_variables_set(i, face_to_cell_index(num_local_cells+3, face_index))  &
            )
                d_m1 = self%curvature_measure(v_m2, v_m1, v   )
                d    = self%curvature_measure(v_m1, v   , v_p1)
                d_p1 = self%curvature_measure(v   , v_p1, v_p2)
                dm4  = self%minmod(self%minmod(4.d0 * d - d_m1, 4.d0 * d_m1 - d), self%minmod(d, d_m1))
                v_ul = v + self%alpha_ * (v - v_p1)
                v_md = 0.5d0 * (v + v_m1 - dm4)
                v_lc = v + 0.5d0 * (v - v_p1) + (self%beta_ / 3.d0) * dm4
                v_min = max(min(v, v_m1, v_md), min(v, v_ul, v_lc))
                v_max = min(max(v, v_m1, v_md), max(v, v_ul, v_lc))
                reconstructed_primitive_variables(i) = self%median(reconstructed_primitive_variables(i), v_min, v_max)
            end associate
        end do
    end function reconstruct_rhc

    pure function minmod(self, x, y) result(m)
        class(mp), intent(in) :: self
        real (real_kind  ), intent(in) :: x, y
        real (real_kind  )             :: m
        m = 0.5d0 * (sign(1.d0, x) + sign(1.d0, y)) * min(abs(x), abs(y))
    end function minmod

    pure function median(self, x, y, z) result(m)
        class(mp), intent(in) :: self
        real (real_kind  ), intent(in) :: x, y, z
        real (real_kind  )             :: m
        m = x + self%minmod(y - x, z - x)
    end function median

    pure function curvature_measure(self, v_m1, v, v_p1) result(c)
        class(mp), intent(in) :: self
        real (real_kind  ), intent(in) :: v_m1, v, v_p1
        real (real_kind  )             :: c
        c = v_p1 - 2.d0 * v + v_m1
    end function curvature_measure
end module class_mp