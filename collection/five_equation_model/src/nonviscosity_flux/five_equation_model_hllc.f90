module five_equation_model_hllc_module
    use typedef_module
    implicit none

    private

    public :: compute_flux_five_equation_model_hllc

    contains

    subroutine initialize_five_equation_model_hllc()
    end subroutine initialize_five_equation_model_hllc

    pure function compute_flux_five_equation_model_hllc(   &
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

        ! Characters are follows:
        ! lhc = left hand cell
        ! rhc = right hand cell
        ! ave = average
        ! c   = soundspeed
        ! rho = density
        real(real_kind) :: lhc_rho, rhc_rho
        real(real_kind) :: ave_vel, ave_c
        ! hllc wave speeds
        real(real_kind) :: s_mid
        real(real_kind) :: s_muinus
        real(real_kind) :: s_puls
        real(real_kind) :: s_left
        real(real_kind) :: s_right
        ! intermediate state
        real(real_kind) :: q_star_left (size(left_conservative ))
        real(real_kind) :: q_star_right(size(right_conservative))
        ! intermediate state coef.
        real(real_kind) :: c_star_left
        real(real_kind) :: c_star_right
        ! numerical flux
        real(real_kind) :: f_left (size(left_conservative ))
        real(real_kind) :: f_right(size(right_conservative))

        associate(                                     &
            lhc_rho1_z1 => left_conservative (1), &
            lhc_rho2_z2 => left_conservative (2), &
            lhc_rho_u   => left_conservative (3), &
            lhc_rho_v   => left_conservative (4), &
            lhc_rho_w   => left_conservative (5), &
            lhc_e       => left_conservative (6), &
            lhc_z1      => left_conservative (7), &
            lhc_u       => left_main_velocity   , &
            lhc_rho     => left_density         , &
            lhc_p       => left_pressure        , &
            lhc_c       => left_soundspeed      , &
            rhc_rho1_z1 => right_conservative(1), &
            rhc_rho2_z2 => right_conservative(2), &
            rhc_rho_u   => right_conservative(3), &
            rhc_rho_v   => right_conservative(4), &
            rhc_rho_w   => right_conservative(5), &
            rhc_e       => right_conservative(6), &
            rhc_z1      => right_conservative(7), &
            rhc_u       => right_main_velocity  , &
            rhc_rho     => right_density        , &
            rhc_p       => right_pressure       , &
            rhc_c       => right_soundspeed       &
        )
            ave_vel = 0.5d0 * (lhc_u + rhc_u)
            ave_c   = 0.5d0 * (lhc_c + rhc_c)

            s_left   = min(ave_vel - ave_c, lhc_u - lhc_c)
            s_right  = max(ave_vel + ave_c, rhc_u + rhc_c)
            s_muinus = min(0.d0, s_left)
            s_puls   = max(0.d0, s_right)
            s_mid    = (rhc_p - lhc_p + lhc_rho * lhc_u * (s_left  - lhc_u)  &
                                      - rhc_rho * rhc_u * (s_right - rhc_u)) &
                     / (lhc_rho * (s_left  - lhc_u) - rhc_rho * (s_right - rhc_u))

            f_left(1) = lhc_rho1_z1     * lhc_u
            f_left(2) = lhc_rho2_z2     * lhc_u
            f_left(3) = lhc_rho_u       * lhc_u + lhc_p
            f_left(4) = lhc_rho_v       * lhc_u
            f_left(5) = lhc_rho_w       * lhc_u
            f_left(6) = (lhc_e + lhc_p) * lhc_u
            f_left(7) = lhc_z1          * lhc_u

            f_right(1) = rhc_rho1_z1     * rhc_u
            f_right(2) = rhc_rho2_z2     * rhc_u
            f_right(3) = rhc_rho_u       * rhc_u + rhc_p
            f_right(4) = rhc_rho_v       * rhc_u
            f_right(5) = rhc_rho_w       * rhc_u
            f_right(6) = (rhc_e + rhc_p) * rhc_u
            f_right(7) = rhc_z1          * rhc_u

            c_star_left = (s_left - lhc_u) / (s_left - s_mid)
            q_star_left(1) = c_star_left * lhc_rho1_z1
            q_star_left(2) = c_star_left * lhc_rho2_z2
            q_star_left(3) = c_star_left * lhc_rho * s_mid
            q_star_left(4) = c_star_left * lhc_rho_v
            q_star_left(5) = c_star_left * lhc_rho_w
            q_star_left(6) = c_star_left * (lhc_e + (s_mid - lhc_u) * (lhc_rho * s_mid + lhc_p / (s_left - lhc_u)))
            q_star_left(7) = c_star_left * lhc_z1

            c_star_right = (s_right - rhc_u) / (s_right - s_mid)
            q_star_right(1) = c_star_right * rhc_rho1_z1
            q_star_right(2) = c_star_right * rhc_rho2_z2
            q_star_right(3) = c_star_right * rhc_rho * s_mid
            q_star_right(4) = c_star_right * rhc_rho_v
            q_star_right(5) = c_star_right * rhc_rho_w
            q_star_right(6) = c_star_right * (rhc_e + (s_mid - rhc_u) * (rhc_rho * s_mid + rhc_p / (s_right - rhc_u)))
            q_star_right(7) = c_star_right * rhc_z1

            flux(:) = 0.5d0 * (1.d0 + sign(1.d0, s_mid))                      &
                    * (f_left(:) + s_muinus * (q_star_left(:) - left_conservative(:))) &
                    + 0.5d0 * (1.d0 - sign(1.d0, s_mid))                      &
                    * (f_right(:) + s_puls * (q_star_right(:) - right_conservative(:)))
        end associate
    end function compute_flux_five_equation_model_hllc
end module five_equation_model_hllc_module