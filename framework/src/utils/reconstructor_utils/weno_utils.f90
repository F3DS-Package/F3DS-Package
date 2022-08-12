module weno_utils
    use typedef_module
    use vector_module

    ! On left side (See Appendix C in [Coralic 2014, JCP]).
    !  i-5/2     i-3/2     i-1/2     i+1/2     i+3/2     i+5/2
    !  p_m3      p_m2      p_m1       p        p_p1      p_p2
    !   *---------*---------*---------*---------*---------*
    !   |         |         |         |         |         |
    !   |    o    |    o    |    o    x    o    |    o    |
    !   |   v_m2  |   v_m1  |    v    |   v_p1  |   v_p2  |
    !   *---------*---------*---------*---------*---------*

    ! On right side (See Appendix C in [Coralic 2014, JCP]).
    !  i-5/2     i-3/2     i-1/2     i+1/2     i+3/2     i+5/2
    !  p_m3      p_m2      p_m1       p        p_p1      p_p2
    !   *---------*---------*---------*---------*---------*
    !   |         |         |         |         |         |
    !   |    o    |    o    x    o    |    o    |    o    |
    !   |   v_m2  |   v_m1  |    v    |   v_p1  |   v_p2  |
    !   *---------*---------*---------*---------*---------*

    implicit none

    private

    public :: compute_ideal_weights_left_side
    public :: compute_js_indicator_coefficients_left_side
    public :: compute_polynomials_coefficients_left_side
    public :: compute_weights_left_side
    public :: compute_polynomials_left_side

    public :: compute_ideal_weights_right_side
    public :: compute_js_indicator_coefficients_right_side
    public :: compute_polynomials_coefficients_right_side
    public :: compute_weights_right_side
    public :: compute_polynomials_right_side

    contains

    pure function compute_ideal_weights_left_side(p_m3, p_m2, p_m1, p, p_p1, p_p2) result(g)
        real(real_kind), intent(in) :: p_m3(3), p_m2(3), p_m1(3), p(3), p_p1(3), p_p2(3) ! face position
        real(real_kind)             :: g(3)

        g(1) = vector_multiply((p    - p_m3), (p    - p_m2)) &
             / vector_multiply((p_p2 - p_m3), (p_p2 - p_m2))
        g(3) = vector_multiply((p_p1 - p   ), (p_p2 - p   )) &
             / vector_multiply((p_p1 - p_m3), (p_p2 - p_m3))
        g(2) = 1.0d0 - g(1) - g(3)
    end function compute_ideal_weights_left_side

    pure function compute_js_indicator_coefficients_left_side(p_m3, p_m2, p_m1, p, p_p1, p_p2) result(b_coef)
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
    end function compute_js_indicator_coefficients_left_side

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

    pure function compute_weights_left_side(v_m2, v_m1, v, v_p1, v_p2, g, b_coef, epsilon) result(w)
        real(real_kind), intent(in) :: v_m2, v_m1, v, v_p1, v_p2 ! value
        real(real_kind), intent(in) :: g(3), b_coef(3,3)
        real(real_kind), intent(in) :: epsilon
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

        w(1) = g(1) / (b(1) + epsilon)**2.d0
        w(2) = g(2) / (b(2) + epsilon)**2.d0
        w(3) = g(3) / (b(3) + epsilon)**2.d0
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

    pure function compute_ideal_weights_right_side(p_m3, p_m2, p_m1, p, p_p1, p_p2) result(g)
        real(real_kind), intent(in) :: p_m3(3), p_m2(3), p_m1(3), p(3), p_p1(3), p_p2(3) ! face position
        real(real_kind)             :: g(3)

        g(1) = vector_multiply((p_m1 - p_m3), (p_m1 - p_m2)) &
             / vector_multiply((p_p2 - p_m3), (p_p2 - p_m2))
        g(3) = vector_multiply((p_p1 - p_m1), (p_p2 - p_m1)) &
             / vector_multiply((p_p1 - p_m3), (p_p2 - p_m3))
        g(2) = 1.0d0 - g(1) - g(3)
    end function compute_ideal_weights_right_side

    pure function compute_js_indicator_coefficients_right_side(p_m3, p_m2, p_m1, p, p_p1, p_p2) result(b_coef)
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
    end function compute_js_indicator_coefficients_right_side

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

    pure function compute_weights_right_side(v_m2, v_m1, v, v_p1, v_p2, g, b_coef, epsilon) result(w)
        real(real_kind), intent(in) :: v_m2, v_m1, v, v_p1, v_p2 ! value
        real(real_kind), intent(in) :: g(3), b_coef(3,3)
        real(real_kind), intent(in) :: epsilon
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

        w(1) = g(1) / (b(1) + epsilon)**2.d0
        w(2) = g(2) / (b(2) + epsilon)**2.d0
        w(3) = g(3) / (b(3) + epsilon)**2.d0
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

end module weno_utils