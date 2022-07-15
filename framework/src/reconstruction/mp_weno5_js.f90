module mp_weno5_js_module
    !**
    !* NOTE: If you use a grid that has curvature, variables are not interpolated by 5th order.
    !* Efficient Implementation of Weighted ENO Schemes (https://www.sciencedirect.com/science/article/pii/S0021999196901308)
    !*

    use typedef_module
    use vector_module
    use abstract_mixture_eos

    implicit none

    private

    integer(int_kind) :: num_ghost_cells_ = 3

    ! WENO
    real   (real_kind) :: epsilon_ = 1.d-6

    ! monotonicity-preserving (MP)
    real   (real_kind) :: alpha_ = 10.d0
    real   (real_kind) :: beta_  = 4.d0

    public :: reconstruct_mp_weno5_js
    public :: reconstruct_lhc_mp_weno5_js
    public :: reconstruct_rhc_mp_weno5_js

    contains

    pure function minmod(x, y) result(m)
        real(real_kind), intent(in) :: x, y
        real(real_kind)             :: m
        m = 0.5d0 * (sign(1.d0, x) + sign(1.d0, y)) * min(abs(x), abs(y))
    end function minmod

    pure function median(x, y, z) result(m)
        real(real_kind), intent(in) :: x, y, z
        real(real_kind)             :: m
        m = x + minmod(y - x, z - x)
    end function median

    pure function curvature_measure(v_m1, v, v_p1) result(c)
        real(real_kind), intent(in) :: v_m1, v, v_p1
        real(real_kind)             :: c
        c = v_p1 - 2.d0 * v + v_m1
    end function curvature_measure

    ! See Appendix C in [Coralic 2014, JCP].
    !  i-5/2     i-3/2     i-1/2     i+1/2     i+3/2     i+5/2
    !  p_m3      p_m2      p_m1       p        p_p1      p_p2
    !   *---------*---------*---------*---------*---------*
    !   |         |         |         |         |         |
    !   |    o    |    o    |    o    x    o    |    o    |
    !   |   v_m2  |   v_m1  |    v    |   v_p1  |   v_p2  |
    !   *---------*---------*---------*---------*---------*
    ! I will rewrite follows to pre-calculate coefficients and weights.
    pure function compute_ideal_weights_left_side(p_m3, p_m2, p_m1, p, p_p1, p_p2) result(g)
        real(real_kind), intent(in) :: p_m3(3), p_m2(3), p_m1(3), p(3), p_p1(3), p_p2(3) ! face position
        real(real_kind)             :: g(3)

        g(1) = vector_multiply((p    - p_m3), (p    - p_m2)) &
             / vector_multiply((p_p2 - p_m3), (p_p2 - p_m2))
        g(3) = vector_multiply((p_p1 - p   ), (p_p2 - p   )) &
             / vector_multiply((p_p1 - p_m3), (p_p2 - p_m3))
        g(2) = 1.0d0 - g(1) - g(3)
    end function compute_ideal_weights_left_side

    pure function compute_indicator_coefficients_left_side(p_m3, p_m2, p_m1, p, p_p1, p_p2) result(b_coef)
        real(real_kind), intent(in) :: p_m3(3), p_m2(3), p_m1(3), p(3), p_p1(3), p_p2(3) ! face position
        real(real_kind)             :: b_coef(3,3)

        b_coef(1,1) = (10.d0*vector_squared(p - p_m1) + vector_multiply((p - p_m1), (p_p1 - p_m1) + (p_p2 - p)) + vector_squared((p_p1 - p_m1) + (p_p2 - p))) &
                    / vector_multiply((p_p1 - p_m1), (p_p2 - p_m1))**2.d0
        b_coef(1,2) = -(19.d0*vector_squared(p - p_m1) - vector_multiply((p - p_m1), (p_p2 - p)) + 2.d0*vector_multiply((p_p1 - p_m1), ((p_p1 - p_m1) + (p_p2 - p)))) &
                    / (vector_multiply((p_p1 - p_m1), (p_p2 - p))*vector_squared(p_p2 - p_m1))
        b_coef(1,3) = (10.d0*vector_squared(p - p_m1) + vector_multiply((p - p_m1), (p_p1 - p)) + vector_squared(p_p1 - p)) &
                    / vector_multiply((p_p2 - p_m1), (p_p2 - p))**2.d0

        b_coef(2,1) = -(vector_multiply((p - p_m1), ((p_m1 - p_m2) + 20.d0*(p - p_m1))) - vector_multiply((p_p1 - p_m1), (2.d0*(p_m1 - p_m2) + (p - p_m1)))) &
                    / (vector_multiply((p - p_m2), (p_p1 - p_m1))*vector_squared(p_p1 - p_m2))
        b_coef(2,2) = (10.d0*vector_squared(p - p_m1) + vector_multiply((p - p_m1), (p_p1 - p)) + vector_squared(p_p1 - p)) &
                    / vector_multiply((p - p_m2), (p_p1 - p_m2))**2.d0
        b_coef(2,3) = (10.d0*vector_squared(p - p_m1) + vector_multiply((p_m1 - p_m2), (p - p_m1)) + vector_squared(p_m1 - p_m2)) &
                    / vector_multiply((p_p1 - p_m2), (p_p1 - p_m1))**2.d0

        b_coef(3,1) = (12.d0*vector_squared(p - p_m1) + 3.d0*vector_multiply((p - p_m1), ((p_m1 - p_m3) + (p_m1 - p_m2))) + vector_squared((p_m1 - p_m3) + (p_m1 - p_m2))) &
                    / vector_multiply((p - p_m3), (p - p_m2))**2.d0
        b_coef(3,2) = (-19.d0*vector_squared(p - p_m1) - vector_multiply((p_m1 - p_m3), (p - p_m1)) + 2.d0*vector_multiply((p - p_m2), (p_m1 - p_m3) + (p - p_m2))) &
                    / (vector_multiply((p_m1 - p_m3), (p - p_m2))*vector_squared(p - p_m3))
        b_coef(3,3) = (10.d0*vector_squared(p - p_m1) + vector_multiply((p_m1 - p_m2), (p - p_m1)) + vector_squared(p_m1 - p_m2)) &
                    / vector_multiply((p_m1 - p_m3), (p - p_m3))**2.d0

        b_coef(:,:) =  b_coef(:,:) * 4.d0 * vector_squared(p - p_m1)
    end function compute_indicator_coefficients_left_side

    pure function compute_polynomials_coefficients_left_side(p_m3, p_m2, p_m1, p, p_p1, p_p2) result(p_coef)
        real(real_kind), intent(in) :: p_m3(3), p_m2(3), p_m1(3), p(3), p_p1(3), p_p2(3) ! face position
        real(real_kind)             :: p_coef(3,2)

        p_coef(1,1) = + vector_multiply((p    - p_m1), ((p_p1 - p_m1) + (p_p2 - p))) &
                      / vector_multiply((p_p1 - p_m1), (p_p2 - p_m1))
        p_coef(1,2) = - vector_multiply((p    - p_m1), (p_p1 - p   )) &
                      / vector_multiply((p_p2 - p_m1), (p_p2 - p   ))
        p_coef(2,1) = + vector_multiply((p    - p_m1), (p_p1 - p   )) &
                      / vector_multiply((p    - p_m2), (p_p1 - p_m2))
        p_coef(2,2) = + vector_multiply((p    - p_m2), (p    - p_m1)) &
                      / vector_multiply((p_p1 - p_m2), (p_p1 - p_m1))
        p_coef(3,1) = - vector_multiply((p    - p_m2), (p    - p_m1)) &
                      / vector_multiply((p_m1 - p_m3), (p    - p_m3))
        p_coef(3,2) = + vector_multiply((p    - p_m1), ((p - p_m3) + (p - p_m2))) &
                      / vector_multiply((p    - p_m3), (p    - p_m2))
    end function compute_polynomials_coefficients_left_side

    pure function compute_weights_left_side(v_m2, v_m1, v, v_p1, v_p2, g, b_coef) result(w)
        real(real_kind), intent(in) :: v_m2, v_m1, v, v_p1, v_p2 ! value
        real(real_kind), intent(in) :: g(3), b_coef(3,3)
        real(real_kind)             :: w(3)
        real(real_kind)             :: b(3)
        real(real_kind)             :: total_w

        b(1) = b_coef(1,1) * (v_p1 - v)**2.d0 &
             + b_coef(1,2) * (v_p1 - v)*(v_p2 - v_p1) &
             + b_coef(1,3) * (v_p2 - v_p1)**2.d0
        b(2) = b_coef(2,1) * (v - v_m1)*(v_p1 - v) &
             + b_coef(2,2) * (v - v_m1)**2.d0 &
             + b_coef(2,3) * (v_p1 - v)**2.d0
        b(3) = b_coef(3,1) * (v - v_m1)**2.d0 &
             + b_coef(3,2) * (v_m1 - v_m2)*(v - v_m1) &
             + b_coef(3,3) * (v_m1 - v_m2)**2.d0

        w(1) = g(1) / (b(1) + epsilon_)**2.d0
        w(2) = g(2) / (b(2) + epsilon_)**2.d0
        w(3) = g(3) / (b(3) + epsilon_)**2.d0
        total_w = w(1) + w(2) + w(3)
        w(1) = w(1) / total_w
        w(2) = w(2) / total_w
        w(3) = w(3) / total_w
    end function compute_weights_left_side

    pure function compute_polynomials_left_side(v_m2, v_m1, v, v_p1, v_p2, p_coef) result(p)
        real(real_kind), intent(in) :: v_m2, v_m1, v, v_p1, v_p2       ! value
        real(real_kind), intent(in) :: p_coef(3,2)
        real(real_kind)             :: p(3)

        p(1) = v + p_coef(1,1) * (v_p1 - v)    + p_coef(1,2) * (v_p2 - v_p1)
        p(2) = v + p_coef(2,1) * (v - v_m1)    + p_coef(2,2) * (v_p1 - v)
        p(3) = v + p_coef(3,1) * (v_m1 - v_m2) + p_coef(3,2) * (v - v_m1)
    end function compute_polynomials_left_side

    ! See Appendix C in [Coralic 2014, JCP].
    !  i-5/2     i-3/2     i-1/2     i+1/2     i+3/2     i+5/2
    !  p_m3      p_m2      p_m1       p        p_p1      p_p2
    !   *---------*---------*---------*---------*---------*
    !   |         |         |         |         |         |
    !   |    o    |    o    x    o    |    o    |    o    |
    !   |   v_m2  |   v_m1  |    v    |   v_p1  |   v_p2  |
    !   *---------*---------*---------*---------*---------*
    ! I will rewrite follows to pre-calculate coefficients and weights.
    pure function compute_ideal_weights_right_side(p_m3, p_m2, p_m1, p, p_p1, p_p2) result(g)
        real(real_kind), intent(in) :: p_m3(3), p_m2(3), p_m1(3), p(3), p_p1(3), p_p2(3) ! face position
        real(real_kind)             :: g(3)

        g(1) = vector_multiply((p_m1 - p_m3), (p_m1 - p_m2)) &
             / vector_multiply((p_p2 - p_m3), (p_p2 - p_m2))
        g(3) = vector_multiply((p_p1 - p_m1), (p_p2 - p_m1)) &
             / vector_multiply((p_p1 - p_m3), (p_p2 - p_m3))
        g(2) = 1.0d0 - g(1) - g(3)
    end function compute_ideal_weights_right_side

    pure function compute_indicator_coefficients_right_side(p_m3, p_m2, p_m1, p, p_p1, p_p2) result(b_coef)
        real(real_kind), intent(in) :: p_m3(3), p_m2(3), p_m1(3), p(3), p_p1(3), p_p2(3) ! face position
        real(real_kind)             :: b_coef(3,3)

        b_coef(1,1) = (10.d0*vector_squared(p - p_m1) + vector_multiply((p - p_m1), (p_p1 - p_m1) + (p_p2 - p)) + vector_squared((p_p1 - p_m1) + (p_p2 - p))) &
                    / vector_multiply((p_p1 - p_m1), (p_p2 - p_m1))**2.d0
        b_coef(1,2) = -(19.d0*vector_squared(p - p_m1) - vector_multiply((p - p_m1), (p_p2 - p)) + 2.d0*vector_multiply((p_p1 - p_m1), ((p_p1 - p_m1) + (p_p2 - p)))) &
                    / (vector_multiply((p_p1 - p_m1), (p_p2 - p))*vector_squared(p_p2 - p_m1))
        b_coef(1,3) = (10.d0*vector_squared(p - p_m1) + vector_multiply((p - p_m1), (p_p1 - p)) + vector_squared(p_p1 - p)) &
                    / vector_multiply((p_p2 - p_m1), (p_p2 - p))**2.d0

        b_coef(2,1) = -(vector_multiply((p - p_m1), ((p_m1 - p_m2) + 20.d0*(p - p_m1))) - vector_multiply((p_p1 - p_m1), (2.d0*(p_m1 - p_m2) + (p - p_m1)))) &
                    / (vector_multiply((p - p_m2), (p_p1 - p_m1))*vector_squared(p_p1 - p_m2))
        b_coef(2,2) = (10.d0*vector_squared(p - p_m1) + vector_multiply((p - p_m1), (p_p1 - p)) + vector_squared(p_p1 - p)) &
                    / vector_multiply((p - p_m2), (p_p1 - p_m2))**2.d0
        b_coef(2,3) = (10.d0*vector_squared(p - p_m1) + vector_multiply((p_m1 - p_m2), (p - p_m1)) + vector_squared(p_m1 - p_m2)) &
                    / vector_multiply((p_p1 - p_m2), (p_p1 - p_m1))**2.d0

        b_coef(3,1) = (12.d0*vector_squared(p - p_m1) + 3.d0*vector_multiply((p - p_m1), ((p_m1 - p_m3) + (p_m1 - p_m2))) + vector_squared((p_m1 - p_m3) + (p_m1 - p_m2))) &
                    / vector_multiply((p - p_m3), (p - p_m2))**2.d0
        b_coef(3,2) = (-19.d0*vector_squared(p - p_m1) - vector_multiply((p_m1 - p_m3), (p - p_m1)) + 2.d0*vector_multiply((p - p_m2), (p_m1 - p_m3) + (p - p_m2))) &
                    / (vector_multiply((p_m1 - p_m3), (p - p_m2))*vector_squared(p - p_m3))
        b_coef(3,3) = (10.d0*vector_squared(p - p_m1) + vector_multiply((p_m1 - p_m2), (p - p_m1)) + vector_squared(p_m1 - p_m2)) &
                    / vector_multiply((p_m1 - p_m3), (p - p_m3))**2.d0

        b_coef(:,:) =  b_coef(:,:) * 4.d0 * vector_squared(p - p_m1)
    end function compute_indicator_coefficients_right_side

    pure function compute_polynomials_coefficients_right_side(p_m3, p_m2, p_m1, p, p_p1, p_p2) result(p_coef)
        real(real_kind), intent(in) :: p_m3(3), p_m2(3), p_m1(3), p(3), p_p1(3), p_p2(3) ! face position
        real(real_kind)             :: p_coef(3,2)

        p_coef(1,1) = - vector_multiply((p    - p_m1), ((p_p1 - p_m1) + (p_p2 - p))) &
                      / vector_multiply((p_p1 - p_m1), (p_p2 - p_m1))
        p_coef(1,2) = + vector_multiply((p    - p_m1), (p_p1 - p   )) &
                      / vector_multiply((p_p2 - p_m1), (p_p2 - p   ))
        p_coef(2,1) = - vector_multiply((p    - p_m1), (p_p1 - p   )) &
                      / vector_multiply((p    - p_m2), (p_p1 - p_m2))
        p_coef(2,2) = - vector_multiply((p    - p_m2), (p    - p_m1)) &
                      / vector_multiply((p_p1 - p_m2), (p_p1 - p_m1))
        p_coef(3,1) = + vector_multiply((p    - p_m2), (p    - p_m1)) &
                      / vector_multiply((p_m1 - p_m3), (p    - p_m3))
        p_coef(3,2) = - vector_multiply((p    - p_m1), ((p - p_m3) + (p - p_m2))) &
                      / vector_multiply((p    - p_m3), (p    - p_m2))
    end function compute_polynomials_coefficients_right_side

    pure function compute_weights_right_side(v_m2, v_m1, v, v_p1, v_p2, g, b_coef) result(w)
        real(real_kind), intent(in) :: v_m2, v_m1, v, v_p1, v_p2 ! value
        real(real_kind), intent(in) :: g(3), b_coef(3,3)
        real(real_kind)             :: w(3)
        real(real_kind)             :: b(3)
        real(real_kind)             :: total_w

        b(1) = b_coef(1,1) * (v_p1 - v)**2.d0 &
             + b_coef(1,2) * (v_p1 - v)*(v_p2 - v_p1) &
             + b_coef(1,3) * (v_p2 - v_p1)**2.d0
        b(2) = b_coef(2,1) * (v - v_m1)*(v_p1 - v) &
             + b_coef(2,2) * (v - v_m1)**2.d0 &
             + b_coef(2,3) * (v_p1 - v)**2.d0
        b(3) = b_coef(3,1) * (v - v_m1)**2.d0 &
             + b_coef(3,2) * (v_m1 - v_m2)*(v - v_m1) &
             + b_coef(3,3) * (v_m1 - v_m2)**2.d0

        w(1) = g(1) / (b(1) + epsilon_)**2.d0
        w(2) = g(2) / (b(2) + epsilon_)**2.d0
        w(3) = g(3) / (b(3) + epsilon_)**2.d0
        total_w = w(1) + w(2) + w(3)
        w(1) = w(1) / total_w
        w(2) = w(2) / total_w
        w(3) = w(3) / total_w
    end function compute_weights_right_side

    pure function compute_polynomials_right_side(v_m2, v_m1, v, v_p1, v_p2, p_coef) result(p)
        real(real_kind), intent(in) :: v_m2, v_m1, v, v_p1, v_p2       ! value
        real(real_kind), intent(in) :: p_coef(3,2)
        real(real_kind)             :: p(3)

        p(1) = v + p_coef(1,1) * (v_p1 - v)    + p_coef(1,2) * (v_p2 - v_p1)
        p(2) = v + p_coef(2,1) * (v - v_m1)    + p_coef(2,2) * (v_p1 - v)
        p(3) = v + p_coef(3,1) * (v_m1 - v_m2) + p_coef(3,2) * (v - v_m1)
    end function compute_polynomials_right_side

    pure function reconstruct_lhc_mp_weno5_js( &
        primitive_values_set             , &
        face_to_cell_index               , &
        cell_positions                   , &
        face_positions                   , &
        face_index                          ) result(reconstructed_primitive)

        real   (real_kind), intent(in) :: primitive_values_set      (:, :)
        integer(int_kind ), intent(in) :: face_to_cell_index        (:, :)
        real   (real_kind), intent(in) :: cell_positions            (:, :)
        real   (real_kind), intent(in) :: face_positions            (:, :)
        integer(int_kind ), intent(in) :: face_index
        real   (real_kind)             :: reconstructed_primitive   (size(primitive_values_set(1, :)))

        integer(int_kind ) :: n_primitives, i
        ! WENO
        real   (real_kind) :: w(3), p(3)
        real   (real_kind) :: ideal_w(3), b_coef(3,3), p_coef(3,2)
        real   (real_kind) :: p_m3(3), p_m2(3), p_m1(3), p_0(3), p_p1(3), p_p2(3)
        ! MP
        real   (real_kind) :: d_p1, d, d_m1, dm4
        real   (real_kind) :: v_ul, v_md, v_lc, v_min, v_max

        n_primitives = size(primitive_values_set(1, :))

        !p_m3 => 1.5d0 * cell_positions(face_to_cell_index(face_index, num_ghost_cells_-2), 1:3) + 0.5d0 * cell_positions(face_to_cell_index(face_index, num_ghost_cells_-1), 1:3), &
        !p_m2 => 0.5d0 * (cell_positions(face_to_cell_index(face_index, num_ghost_cells_-2), 1:3) + cell_positions(face_to_cell_index(face_index, num_ghost_cells_-1), 1:3)), &
        !p_m1 => 0.5d0 * (cell_positions(face_to_cell_index(face_index, num_ghost_cells_-1), 1:3) + cell_positions(face_to_cell_index(face_index, num_ghost_cells_+0), 1:3)), &
        !p    => 0.5d0 * (cell_positions(face_to_cell_index(face_index, num_ghost_cells_+0), 1:3) + cell_positions(face_to_cell_index(face_index, num_ghost_cells_+1), 1:3)), &
        !p_p1 => 0.5d0 * (cell_positions(face_to_cell_index(face_index, num_ghost_cells_+1), 1:3) + cell_positions(face_to_cell_index(face_index, num_ghost_cells_+2), 1:3)), &
        !p_p2 => 0.5d0 * (cell_positions(face_to_cell_index(face_index, num_ghost_cells_+2), 1:3) + cell_positions(face_to_cell_index(face_index, num_ghost_cells_+3), 1:3))  &
        p_m3(1:3) = cell_positions(face_to_cell_index(face_index, num_ghost_cells_-2), 1:3)
        p_m2(1:3) = cell_positions(face_to_cell_index(face_index, num_ghost_cells_-1), 1:3)
        p_m1(1:3) = cell_positions(face_to_cell_index(face_index, num_ghost_cells_+0), 1:3)
        p_0 (1:3) = cell_positions(face_to_cell_index(face_index, num_ghost_cells_+1), 1:3)
        p_p1(1:3) = cell_positions(face_to_cell_index(face_index, num_ghost_cells_+2), 1:3)
        p_p2(1:3) = cell_positions(face_to_cell_index(face_index, num_ghost_cells_+3), 1:3)
        ideal_w = compute_ideal_weights_left_side           (p_m3, p_m2, p_m1, p_0, p_p1, p_p2)
        b_coef  = compute_indicator_coefficients_left_side  (p_m3, p_m2, p_m1, p_0, p_p1, p_p2)
        p_coef  = compute_polynomials_coefficients_left_side(p_m3, p_m2, p_m1, p_0, p_p1, p_p2)


        do i = 1, n_primitives, 1
            associate(                                                                                              &
                v_m2 => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_-2), i), &
                v_m1 => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_-1), i), &
                v    => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_+0), i), &
                v_p1 => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_+1), i), &
                v_p2 => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_+2), i)  &
            )
                ! WENO-JS
                w(1:3) = compute_weights_left_side    (v_m2, v_m1, v, v_p1, v_p2, ideal_w, b_coef)
                p(1:3) = compute_polynomials_left_side(v_m2, v_m1, v, v_p1, v_p2, p_coef)
                reconstructed_primitive(i) = w(1) * p(1) + w(2) * p(2) + w(3) * p(3)
                ! MP
                d_m1 = curvature_measure(v_m2, v_m1, v   )
                d    = curvature_measure(v_m1, v   , v_p1)
                d_p1 = curvature_measure(v   , v_p1, v_p2)
                dm4  = minmod(minmod(4.d0 * d - d_p1, 4.d0 * d_p1 - d), minmod(d, d_p1))
                v_ul = v + alpha_ * (v - v_m1)
                v_md = 0.5d0 * (v + v_p1 - dm4)
                v_lc = v + 0.5d0 * (v - v_m1) + (beta_ / 3.d0) * dm4
                v_min = max(min(v, v_p1, v_md), min(v, v_ul, v_lc))
                v_max = min(max(v, v_p1, v_md), max(v, v_ul, v_lc))
                reconstructed_primitive(i) = median(reconstructed_primitive(i), v_min, v_max)
            end associate
        end do
    end function reconstruct_lhc_mp_weno5_js

    pure function reconstruct_rhc_mp_weno5_js( &
        primitive_values_set             , &
        face_to_cell_index               , &
        cell_positions                   , &
        face_positions                   , &
        face_index                          ) result(reconstructed_primitive)

        real   (real_kind), intent(in) :: primitive_values_set      (:, :)
        integer(int_kind ), intent(in) :: face_to_cell_index (:, :)
        real   (real_kind), intent(in) :: cell_positions            (:, :)
        real   (real_kind), intent(in) :: face_positions            (:, :)
        integer(int_kind ), intent(in) :: face_index
        real   (real_kind)             :: reconstructed_primitive   (size(primitive_values_set(1, :)))

        integer(int_kind ) :: n_primitives, i
        ! WENO
        real   (real_kind) :: w(3), p(3)
        real   (real_kind) :: ideal_w(3), b_coef(3,3), p_coef(3,2)
        real   (real_kind) :: p_m3(3), p_m2(3), p_m1(3), p_0(3), p_p1(3), p_p2(3)
        ! MP
        real   (real_kind) :: d_p1, d, d_m1, dm4
        real   (real_kind) :: v_ul, v_md, v_lc, v_min, v_max

        n_primitives = size(primitive_values_set(1, :))

        !p_m3 => 0.5d0 * (cell_positions(face_to_cell_index(face_index, num_ghost_cells_-2), 1:3) + cell_positions(face_to_cell_index(face_index, num_ghost_cells_-1), 1:3)), &
        !p_m2 => 0.5d0 * (cell_positions(face_to_cell_index(face_index, num_ghost_cells_-1), 1:3) + cell_positions(face_to_cell_index(face_index, num_ghost_cells_+0), 1:3)), &
        !p_m1 => 0.5d0 * (cell_positions(face_to_cell_index(face_index, num_ghost_cells_+0), 1:3) + cell_positions(face_to_cell_index(face_index, num_ghost_cells_+1), 1:3)), &
        !p    => 0.5d0 * (cell_positions(face_to_cell_index(face_index, num_ghost_cells_+1), 1:3) + cell_positions(face_to_cell_index(face_index, num_ghost_cells_+2), 1:3)), &
        !p_p1 => 0.5d0 * (cell_positions(face_to_cell_index(face_index, num_ghost_cells_+2), 1:3) + cell_positions(face_to_cell_index(face_index, num_ghost_cells_+3), 1:3)), &
        !p_p2 => 1.5d0 * cell_positions(face_to_cell_index(face_index, num_ghost_cells_+3), 1:3) + 0.5d0 * cell_positions(face_to_cell_index(face_index, num_ghost_cells_+2), 1:3) &

        p_m3(1:3) = cell_positions(face_to_cell_index(face_index, num_ghost_cells_-2), 1:3)
        p_m2(1:3) = cell_positions(face_to_cell_index(face_index, num_ghost_cells_-1), 1:3)
        p_m1(1:3) = cell_positions(face_to_cell_index(face_index, num_ghost_cells_+0), 1:3)
        p_0 (1:3) = cell_positions(face_to_cell_index(face_index, num_ghost_cells_+1), 1:3)
        p_p1(1:3) = cell_positions(face_to_cell_index(face_index, num_ghost_cells_+2), 1:3)
        p_p2(1:3) = cell_positions(face_to_cell_index(face_index, num_ghost_cells_+3), 1:3)

        ideal_w = compute_ideal_weights_right_side           (p_m3, p_m2, p_m1, p_0, p_p1, p_p2)
        b_coef  = compute_indicator_coefficients_right_side  (p_m3, p_m2, p_m1, p_0, p_p1, p_p2)
        p_coef  = compute_polynomials_coefficients_right_side(p_m3, p_m2, p_m1, p_0, p_p1, p_p2)

        do i = 1, n_primitives, 1
            associate(                                                                                              &
                v_m2       => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_-1), i), &
                v_m1       => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_+0), i), &
                v          => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_+1), i), &
                v_p1       => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_+2), i), &
                v_p2       => primitive_values_set(face_to_cell_index(face_index, num_ghost_cells_+3), i)  &
            )
                ! WENO5-JS
                w(1:3) = compute_weights_right_side    (v_m2, v_m1, v, v_p1, v_p2, ideal_w, b_coef)
                p(1:3) = compute_polynomials_right_side(v_m2, v_m1, v, v_p1, v_p2, p_coef)
                reconstructed_primitive(i) = w(1) * p(1) + w(2) * p(2) + w(3) * p(3)
                ! MP
                d_m1 = curvature_measure(v_m2, v_m1, v   )
                d    = curvature_measure(v_m1, v   , v_p1)
                d_p1 = curvature_measure(v   , v_p1, v_p2)
                dm4  = minmod(minmod(4.d0 * d - d_m1, 4.d0 * d_m1 - d), minmod(d, d_m1))
                v_ul = v + alpha_ * (v - v_p1)
                v_md = 0.5d0 * (v + v_m1 - dm4)
                v_lc = v + 0.5d0 * (v - v_p1) + (beta_ / 3.d0) * dm4
                v_min = max(min(v, v_m1, v_md), min(v, v_ul, v_lc))
                v_max = min(max(v, v_m1, v_md), max(v, v_ul, v_lc))
                reconstructed_primitive(i) = median(reconstructed_primitive(i), v_min, v_max)
            end associate
        end do
    end function reconstruct_rhc_mp_weno5_js

    pure function reconstruct_mp_weno5_js(  &
        primitive_values_set              , &
        cell_centor_positions             , &
        cell_volumes                      , &
        face_to_cell_index                , &
        face_centor_positions             , &
        face_normal_vectors               , &
        face_tangential1_vectors          , &
        face_tangential2_vectors          , &
        face_areas                        , &
        face_index                        , &
        n_conservative_values             , &
        n_derivative_values               , &
        eos                               , &
        flux_function                     , &
        primitive_to_conservative_function, &
        integrated_element_function             ) result(element)

        real   (real_kind  ), intent(in) :: primitive_values_set     (:, :)
        real   (real_kind  ), intent(in) :: cell_centor_positions    (:, :)
        real   (real_kind  ), intent(in) :: cell_volumes             (:)
        integer(int_kind   ), intent(in) :: face_to_cell_index       (:,:)
        real   (real_kind  ), intent(in) :: face_normal_vectors      (:,:)
        real   (real_kind  ), intent(in) :: face_tangential1_vectors (:,:)
        real   (real_kind  ), intent(in) :: face_tangential2_vectors (:,:)
        real   (real_kind  ), intent(in) :: face_centor_positions    (:,:)
        real   (real_kind  ), intent(in) :: face_areas               (:)
        integer(int_kind   ), intent(in) :: face_index
        integer(int_kind   ), intent(in) :: n_conservative_values
        integer(int_kind   ), intent(in) :: n_derivative_values
        class  (mixture_eos), intent(in) :: eos

        real   (real_kind)              :: element(2, n_conservative_values+n_derivative_values)

        interface
            pure function flux_function(       &
                left_conservative            , &
                left_main_velocity           , &
                left_density                 , &
                left_pressure                , &
                left_soundspeed              , &
                right_conservative           , &
                right_main_velocity          , &
                right_density                , &
                right_pressure               , &
                right_soundspeed             ) result(flux)

                use typedef_module
                real(real_kind), intent(in) :: left_conservative   (:)
                real(real_kind), intent(in) :: left_main_velocity
                real(real_kind), intent(in) :: left_density
                real(real_kind), intent(in) :: left_pressure
                real(real_kind), intent(in) :: left_soundspeed
                real(real_kind), intent(in) :: right_conservative   (:)
                real(real_kind), intent(in) :: right_main_velocity
                real(real_kind), intent(in) :: right_density
                real(real_kind), intent(in) :: right_pressure
                real(real_kind), intent(in) :: right_soundspeed
                real(real_kind)             :: flux(size(left_conservative))
            end function flux_function

            pure function primitive_to_conservative_function(primitive, eos) result(conservative)
                use typedef_module
                use abstract_mixture_eos
                real (real_kind  ), intent(in)  :: primitive   (:)
                class(mixture_eos), intent(in)  :: eos
                real (real_kind  ), allocatable :: conservative(:)
            end function primitive_to_conservative_function

            pure function integrated_element_function( &
                reconstructed_lhc_primitive       , &
                reconstructed_rhc_primitive       , &
                lhc_primitive                     , &
                rhc_primitive                     , &
                lhc_cell_volume                   , &
                rhc_cell_volume                   , &
                face_normal_vector                , &
                face_tangential1_vector           , &
                face_tangential2_vector           , &
                face_area                         , &
                n_conservative_values             , &
                n_derivative_values               , &
                eos                               , &
                flux_function                     , &
                primitive_to_conservative_function   ) result(element)

                use typedef_module
                use abstract_mixture_eos

                real   (real_kind  ), intent(in) :: reconstructed_lhc_primitive (:)
                real   (real_kind  ), intent(in) :: reconstructed_rhc_primitive (:)
                real   (real_kind  ), intent(in) :: lhc_primitive               (:)
                real   (real_kind  ), intent(in) :: rhc_primitive               (:)
                real   (real_kind  ), intent(in) :: lhc_cell_volume
                real   (real_kind  ), intent(in) :: rhc_cell_volume
                real   (real_kind  ), intent(in) :: face_normal_vector                (3)
                real   (real_kind  ), intent(in) :: face_tangential1_vector           (3)
                real   (real_kind  ), intent(in) :: face_tangential2_vector           (3)
                real   (real_kind  ), intent(in) :: face_area
                integer(int_kind   ), intent(in) :: n_conservative_values
                integer(int_kind   ), intent(in) :: n_derivative_values
                class  (mixture_eos), intent(in) :: eos
                real   (real_kind)               :: element        (2, n_conservative_values+n_derivative_values)

                interface
                    pure function flux_function(       &
                        left_conservative            , &
                        left_main_velocity           , &
                        left_density                 , &
                        left_pressure                , &
                        left_soundspeed              , &
                        right_conservative           , &
                        right_main_velocity          , &
                        right_density                , &
                        right_pressure               , &
                        right_soundspeed             ) result(flux)

                        use typedef_module
                        real(real_kind), intent(in) :: left_conservative   (:)
                        real(real_kind), intent(in) :: left_main_velocity
                        real(real_kind), intent(in) :: left_density
                        real(real_kind), intent(in) :: left_pressure
                        real(real_kind), intent(in) :: left_soundspeed
                        real(real_kind), intent(in) :: right_conservative   (:)
                        real(real_kind), intent(in) :: right_main_velocity
                        real(real_kind), intent(in) :: right_density
                        real(real_kind), intent(in) :: right_pressure
                        real(real_kind), intent(in) :: right_soundspeed
                        real(real_kind)             :: flux(size(left_conservative))
                    end function flux_function

                    pure function primitive_to_conservative_function(primitive, eos) result(conservative)
                        use typedef_module
                        use abstract_mixture_eos
                        real (real_kind  ), intent(in)  :: primitive   (:)
                        class(mixture_eos), intent(in)  :: eos
                        real (real_kind  ), allocatable :: conservative(:)
                    end function primitive_to_conservative_function
                end interface
            end function integrated_element_function
        end interface

        integer(int_kind )              :: lhc_index, rhc_index
        real   (real_kind), allocatable :: lhc_primitive   (:)
        real   (real_kind), allocatable :: rhc_primitive   (:)
        integer(int_kind )              :: n_primitives, i
        real   (real_kind)              :: delta

        lhc_index = face_to_cell_index(face_index, num_ghost_cells_+0)
        rhc_index = face_to_cell_index(face_index, num_ghost_cells_+1)

        lhc_primitive = reconstruct_lhc_mp_weno5_js(      &
            primitive_values_set                       , &
            face_to_cell_index                         , &
            cell_centor_positions                      , &
            face_centor_positions                      , &
            face_index                                   &
        )
        rhc_primitive = reconstruct_rhc_mp_weno5_js(     &
            primitive_values_set                       , &
            face_to_cell_index                         , &
            cell_centor_positions                      , &
            face_centor_positions                      , &
            face_index                                   &
        )

        element = integrated_element_function(  &
            lhc_primitive                            , &
            rhc_primitive                            , &
            primitive_values_set(lhc_index, :)       , &
            primitive_values_set(rhc_index, :)       , &
            cell_volumes(lhc_index)                  , &
            cell_volumes(rhc_index)                  , &
            face_normal_vectors     (face_index, 1:3), &
            face_tangential1_vectors(face_index, 1:3), &
            face_tangential2_vectors(face_index, 1:3), &
            face_areas              (face_index)     , &
            n_conservative_values                    , &
            n_derivative_values                      , &
            eos                                      , &
            flux_function                            , &
            primitive_to_conservative_function         &
        )
    end function reconstruct_mp_weno5_js
end module mp_weno5_js_module