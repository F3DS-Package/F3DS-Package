module class_mp_weno5_js
    !**
    !* NOTE: If you use a grid that has curvature, variables are not interpolated by 5th order.
    !* Efficient Implementation of Weighted ENO Schemes (https://www.sciencedirect.com/science/article/pii/S0021999196901308)
    !*

    use typedef_module
    use abstract_reconstructor
    use abstract_configuration
    use weno_utils

    implicit none

    private

    type, public, extends(reconstructor) :: mp_weno5_js
        private

        ! parameter for WENO
        real(real_kind) :: epsilon_

        ! parameters for monotonicity-preserving (MP)
        real(real_kind) :: alpha_
        real(real_kind) :: beta_

        contains

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: reconstruct_lhc
        procedure, public, pass(self) :: reconstruct_rhc

        procedure, pass(self) :: minmod
        procedure, pass(self) :: median
        procedure, pass(self) :: curvature_measure
    end type mp_weno5_js

    contains

    subroutine initialize(self, config)
        class(mp_weno5_js  ), intent(inout) :: self
        class(configuration), intent(inout) :: config
        logical :: found

        call config%get_real("Reconstructor.MP-WENO5.Epsilon", self%epsilon_, found, 1.d-6)
        if(.not. found) call write_warring("'Reconstructor.MP-WENO5.Epsilon' is not found in configration you set. To be set dafault value.")

        call config%get_real("Reconstructor.MP-WENO5.Alpha", self%alpha_, found, 10.d0)
        if(.not. found) call write_warring("'Reconstructor.MP-WENO5.Alpha' is not found in configration you set. To be set dafault value.")

        call config%get_real("Reconstructor.MP-WENO5.Beta", self%beta_, found, 4.d0)
        if(.not. found) call write_warring("'Reconstructor.MP-WENO5.Beta' is not found in configration you set. To be set dafault value.")
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

        class  (mp_weno5_js), intent(in) :: self
        integer(int_kind   ), intent(in) :: face_index
        real   (real_kind  ), intent(in) :: primitive_variables_set          (:, :)
        integer(int_kind   ), intent(in) :: face_to_cell_index               (:, :)
        real   (real_kind  ), intent(in) :: cell_centor_positions            (:, :)
        real   (real_kind  ), intent(in) :: face_centor_positions            (:, :)
        integer(int_kind   ), intent(in) :: num_local_cells
        integer(int_kind   ), intent(in) :: num_primitive_variables
        real   (real_kind  )             :: reconstructed_primitive_variables(num_primitive_variables)

        integer(int_kind   ) :: i
        ! WENO
        real   (real_kind  ) :: w(3), p(3)
        real   (real_kind  ) :: ideal_w(3), b_coef(3,3), p_coef(3,2)
        ! MP
        real   (real_kind  ) :: d_p1, d, d_m1, dm4
        real   (real_kind  ) :: v_ul, v_md, v_lc, v_min, v_max

        associate(                                                                                 &
            p_m3 => cell_centor_positions(1:3, face_to_cell_index(num_local_cells-2, face_index)), &
            p_m2 => cell_centor_positions(1:3, face_to_cell_index(num_local_cells-1, face_index)), &
            p_m1 => cell_centor_positions(1:3, face_to_cell_index(num_local_cells+0, face_index)), &
            p_0  => cell_centor_positions(1:3, face_to_cell_index(num_local_cells+1, face_index)), &
            p_p1 => cell_centor_positions(1:3, face_to_cell_index(num_local_cells+2, face_index)), &
            p_p2 => cell_centor_positions(1:3, face_to_cell_index(num_local_cells+3, face_index))  &
        )
            ideal_w(:  ) = compute_ideal_weights_left_side            (p_m3, p_m2, p_m1, p_0, p_p1, p_p2)
            b_coef (:,:) = compute_js_indicator_coefficients_left_side(p_m3, p_m2, p_m1, p_0, p_p1, p_p2)
            p_coef (:,:) = compute_polynomials_coefficients_left_side (p_m3, p_m2, p_m1, p_0, p_p1, p_p2)
        end associate

        do i = 1, num_primitive_variables, 1
            associate(                                                                                 &
                v_m2 => primitive_variables_set(i, face_to_cell_index(num_local_cells-2, face_index)), &
                v_m1 => primitive_variables_set(i, face_to_cell_index(num_local_cells-1, face_index)), &
                v    => primitive_variables_set(i, face_to_cell_index(num_local_cells+0, face_index)), &
                v_p1 => primitive_variables_set(i, face_to_cell_index(num_local_cells+1, face_index)), &
                v_p2 => primitive_variables_set(i, face_to_cell_index(num_local_cells+2, face_index))  &
            )
                ! WENO
                w(1:3) = compute_weights_left_side    (v_m2, v_m1, v, v_p1, v_p2, ideal_w, b_coef, self%epsilon_)
                p(1:3) = compute_polynomials_left_side(v_m2, v_m1, v, v_p1, v_p2, p_coef)
                reconstructed_primitive_variables(i) = w(1) * p(1) + w(2) * p(2) + w(3) * p(3)
                ! MP
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

        class  (mp_weno5_js), intent(in) :: self
        integer(int_kind   ), intent(in) :: face_index
        real   (real_kind  ), intent(in) :: primitive_variables_set          (:, :)
        integer(int_kind   ), intent(in) :: face_to_cell_index               (:, :)
        real   (real_kind  ), intent(in) :: cell_centor_positions            (:, :)
        real   (real_kind  ), intent(in) :: face_centor_positions            (:, :)
        integer(int_kind   ), intent(in) :: num_local_cells
        integer(int_kind   ), intent(in) :: num_primitive_variables
        real   (real_kind  )             :: reconstructed_primitive_variables(num_primitive_variables)

        integer(int_kind   ) :: i
        ! WENO
        real   (real_kind  ) :: w(3), p(3)
        real   (real_kind  ) :: ideal_w(3), b_coef(3,3), p_coef(3,2)
        ! MP
        real   (real_kind  ) :: d_p1, d, d_m1, dm4
        real   (real_kind  ) :: v_ul, v_md, v_lc, v_min, v_max

        associate(                                                                                 &
            p_m3 => cell_centor_positions(1:3, face_to_cell_index(num_local_cells-2, face_index)), &
            p_m2 => cell_centor_positions(1:3, face_to_cell_index(num_local_cells-1, face_index)), &
            p_m1 => cell_centor_positions(1:3, face_to_cell_index(num_local_cells+0, face_index)), &
            p_0  => cell_centor_positions(1:3, face_to_cell_index(num_local_cells+1, face_index)), &
            p_p1 => cell_centor_positions(1:3, face_to_cell_index(num_local_cells+2, face_index)), &
            p_p2 => cell_centor_positions(1:3, face_to_cell_index(num_local_cells+3, face_index))  &
        )
            ideal_w(:  ) = compute_ideal_weights_right_side            (p_m3, p_m2, p_m1, p_0, p_p1, p_p2)
            b_coef (:,:) = compute_js_indicator_coefficients_right_side(p_m3, p_m2, p_m1, p_0, p_p1, p_p2)
            p_coef (:,:) = compute_polynomials_coefficients_right_side (p_m3, p_m2, p_m1, p_0, p_p1, p_p2)
        end associate

        do i = 1, num_primitive_variables, 1
            associate(                                                                                 &
                v_m2 => primitive_variables_set(i, face_to_cell_index(num_local_cells-1, face_index)), &
                v_m1 => primitive_variables_set(i, face_to_cell_index(num_local_cells+0, face_index)), &
                v    => primitive_variables_set(i, face_to_cell_index(num_local_cells+1, face_index)), &
                v_p1 => primitive_variables_set(i, face_to_cell_index(num_local_cells+2, face_index)), &
                v_p2 => primitive_variables_set(i, face_to_cell_index(num_local_cells+3, face_index))  &
            )
                ! WENO
                w(1:3) = compute_weights_right_side    (v_m2, v_m1, v, v_p1, v_p2, ideal_w, b_coef, self%epsilon_)
                p(1:3) = compute_polynomials_right_side(v_m2, v_m1, v, v_p1, v_p2, p_coef)
                reconstructed_primitive_variables(i) = w(1) * p(1) + w(2) * p(2) + w(3) * p(3)
                ! MP
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
        class(mp_weno5_js), intent(in) :: self
        real (real_kind  ), intent(in) :: x, y
        real (real_kind  )             :: m
        m = 0.5d0 * (sign(1.d0, x) + sign(1.d0, y)) * min(abs(x), abs(y))
    end function minmod

    pure function median(self, x, y, z) result(m)
        class(mp_weno5_js), intent(in) :: self
        real (real_kind  ), intent(in) :: x, y, z
        real (real_kind  )             :: m
        m = x + self%minmod(y - x, z - x)
    end function median

    pure function curvature_measure(self, v_m1, v, v_p1) result(c)
        class(mp_weno5_js), intent(in) :: self
        real (real_kind  ), intent(in) :: v_m1, v, v_p1
        real (real_kind  )             :: c
        c = v_p1 - 2.d0 * v + v_m1
    end function curvature_measure
end module class_mp_weno5_js