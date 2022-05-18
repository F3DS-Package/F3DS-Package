module five_equation_space_model_module
    use vector_module
    use abstract_mixture_eos

    implicit none

    private

    public :: compute_space_element_five_equation_model

    contains

    pure function compute_space_element_five_equation_model(&
        reconstructed_lhc_primitive                 , &
        reconstructed_rhc_primitive                 , &
        lhc_primitive                               , &
        rhc_primitive                               , &
        lhc_cell_volume                             , &
        rhc_cell_volume                             , &
        face_normal_vector                                , &
        face_tangential1_vector                           , &
        face_tangential2_vector                           , &
        face_area                                         , &
        n_conservative_values                             , &
        n_derivative_values                               , &
        eos                                               , &
        flux_function                                     , &
        primitive_to_conservative_function                 ) result(element)

        use typedef_module

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
        real   (real_kind)               :: element        (2, n_conservative_values+n_derivative_values) ! 1 -> left, 2-> right

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
            end function
        end interface

        real(real_kind) :: local_coordinate_lhc_primitive(7)
        real(real_kind) :: local_coordinate_rhc_primitive(7)
        real(real_kind) :: lhc_conservative              (7)
        real(real_kind) :: rhc_conservative              (7)
        real(real_kind) :: nonviscosity_flux             (7)
        real(real_kind) :: lhc_soundspeed, lhc_pressure, lhc_density, lhc_main_velocity
        real(real_kind) :: rhc_soundspeed, rhc_pressure, rhc_density, rhc_main_velocity

        ! ## HLLC wave speeds
        real(real_kind) :: ave_vel, ave_c
        real(real_kind) :: s_mid
        real(real_kind) :: s_muinus
        real(real_kind) :: s_puls
        real(real_kind) :: s_left
        real(real_kind) :: s_right
        real(real_kind) :: numerical_velocity


        ! # compute primitive-variables of face local coordinate
        ! ## left-side
        local_coordinate_lhc_primitive(1:2) = reconstructed_lhc_primitive(1:2)
        local_coordinate_lhc_primitive(3:5) = rotate_vector( &
            reconstructed_lhc_primitive(3:5),        &
            face_normal_vector,             &
            face_tangential1_vector,        &
            face_tangential2_vector         &
        )
        local_coordinate_lhc_primitive(6:7) = reconstructed_lhc_primitive(6:7)
        ! ## right-side
        local_coordinate_rhc_primitive(1:2) = reconstructed_rhc_primitive(1:2)
        local_coordinate_rhc_primitive(3:5) = rotate_vector( &
            reconstructed_rhc_primitive(3:5),       &
            face_normal_vector,             &
            face_tangential1_vector,        &
            face_tangential2_vector         &
        )
        local_coordinate_rhc_primitive(6:7) = reconstructed_rhc_primitive(6:7)

        ! # compute conservative-variables
        lhc_conservative = primitive_to_conservative_function( &
            local_coordinate_lhc_primitive,                    &
            eos                                                &
        )
        rhc_conservative = primitive_to_conservative_function( &
            local_coordinate_rhc_primitive,                    &
            eos                                                &
        )

        ! # compute EoS, main velosity, and fluxs
        associate(                                         &
                rho1_z1 => local_coordinate_lhc_primitive(1), &
                rho2_z2 => local_coordinate_lhc_primitive(2), &
                u       => local_coordinate_lhc_primitive(3), &
                v       => local_coordinate_lhc_primitive(4), &
                w       => local_coordinate_lhc_primitive(5), &
                p       => local_coordinate_lhc_primitive(6), &
                z1      => local_coordinate_lhc_primitive(7)  &
            )
            lhc_density    = rho1_z1 + rho2_z2
            lhc_pressure   = p
            lhc_soundspeed = eos%compute_soundspeed(p, lhc_density, z1)
            lhc_main_velocity = u
        end associate
        associate(                        &
                rho1_z1 => local_coordinate_rhc_primitive(1), &
                rho2_z2 => local_coordinate_rhc_primitive(2), &
                u       => local_coordinate_rhc_primitive(3), &
                v       => local_coordinate_rhc_primitive(4), &
                w       => local_coordinate_rhc_primitive(5), &
                p       => local_coordinate_rhc_primitive(6), &
                z1      => local_coordinate_rhc_primitive(7)  &
            )
            rhc_density    = rho1_z1 + rho2_z2
            rhc_pressure   = p
            rhc_soundspeed = eos%compute_soundspeed(p, rhc_density, z1)
            rhc_main_velocity = u
        end associate

        ! # compute flux
        ! ## local coordinate flux
        nonviscosity_flux = flux_function( &
            lhc_conservative             , &
            lhc_main_velocity            , &
            lhc_density                  , &
            lhc_pressure                 , &
            lhc_soundspeed               , &
            rhc_conservative             , &
            rhc_main_velocity            , &
            rhc_density                  , &
            rhc_pressure                 , &
            rhc_soundspeed                 &
        )
        ! ## convert to global coordinate flux
        nonviscosity_flux(3:5) = reverse_vector(  &
            nonviscosity_flux(3:5)              , &
            face_normal_vector                  , &
            face_tangential1_vector             , &
            face_tangential2_vector               &
        )

        ! # summation nonviscosity-flux
        element(1, 1:7) = (-1.d0 / lhc_cell_volume) * nonviscosity_flux(:) * face_area
        element(2, 1:7) = (+1.d0 / rhc_cell_volume) * nonviscosity_flux(:) * face_area

        ! # Compute (- alpha1 - K) * div(u) according to [Schmidmayer 2020, JCP] using HLLC.
        ! # If you choice other Riemann solver, term (- alpha1 - K) * div(u) is computed by HLLC forcely. Sorry.
        associate(                                               &
            lhc_u       => local_coordinate_lhc_primitive(3)   , &
            lhc_z1      => local_coordinate_lhc_primitive(7)   , &
            lhc_rho     => lhc_density                         , &
            lhc_p       => lhc_pressure                        , &
            lhc_c       => lhc_soundspeed                      , &
            rhc_u       => local_coordinate_rhc_primitive(3)   , &
            rhc_z1      => local_coordinate_rhc_primitive(7)   , &
            rhc_rho     => rhc_density                         , &
            rhc_p       => rhc_pressure                        , &
            rhc_c       => rhc_soundspeed                        &
        )
        ! ## Compute HLLC valiables and numerical velocity
            ave_vel  = 0.5d0 * (lhc_u + rhc_u)
            ave_c    = 0.5d0 * (lhc_c + rhc_c)
            s_left   = min(ave_vel - ave_c, lhc_u - lhc_c)
            s_right  = max(ave_vel + ave_c, rhc_u + rhc_c)
            s_muinus = min(0.d0, s_left )
            s_puls   = max(0.d0, s_right)
            s_mid    = (rhc_p - lhc_p + lhc_rho * lhc_u * (s_left  - lhc_u)  &
                                      - rhc_rho * rhc_u * (s_right - rhc_u)) &
                     / (lhc_rho * (s_left  - lhc_u) - rhc_rho * (s_right - rhc_u))
            numerical_velocity &
                = 0.5d0 * (1.d0 + sign(1.d0, s_mid)) * (lhc_u + s_muinus * ((s_left  - lhc_u) / (s_left  - s_mid) - 1.d0)) &
                + 0.5d0 * (1.d0 - sign(1.d0, s_mid)) * (rhc_u + s_puls   * ((s_right - rhc_u) / (s_right - s_mid) - 1.d0))
        end associate

        ! # z1 * div(u)
        element(1, 7) = element(1, 7) &
                      - lhc_primitive(7) * (-1.d0 / lhc_cell_volume) * numerical_velocity * face_area
        element(2, 7) = element(2, 7) &
                      - rhc_primitive(7) * (+1.d0 / rhc_cell_volume) * numerical_velocity * face_area
    end function compute_space_element_five_equation_model
end module five_equation_space_model_module