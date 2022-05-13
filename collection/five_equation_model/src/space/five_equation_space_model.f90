module five_equation_space_model_module
    use vector_module

    implicit none

    private

    public :: compute_space_element_five_equation_model

    contains

    pure function compute_space_element_five_equation_model(&
        reconstructed_leftside_primitive                  , &
        reconstructed_rightside_primitive                 , &
        leftside_cell_volume                              , &
        rightside_cell_volume                             , &
        face_normal_vector                                , &
        face_tangential1_vector                           , &
        face_tangential2_vector                           , &
        face_area                                         , &
        n_conservative_values                             , &
        flux_function                                     , &
        eos_pressure_function                             , &
        eos_soundspeed_function                           , &
        primitive_to_conservative_function                 ) result(element_lef_and_right_side)

        use typedef_module

        real   (real_kind), intent(in ) :: reconstructed_leftside_primitive  (:)
        real   (real_kind), intent(in ) :: reconstructed_rightside_primitive (:)
        real   (real_kind), intent(in ) :: leftside_cell_volume
        real   (real_kind), intent(in ) :: rightside_cell_volume
        real   (real_kind), intent(in ) :: face_normal_vector                (3)
        real   (real_kind), intent(in ) :: face_tangential1_vector           (3)
        real   (real_kind), intent(in ) :: face_tangential2_vector           (3)
        real   (real_kind), intent(in ) :: face_area
        integer(int_kind ), intent(in ) :: n_conservative_values
        real   (real_kind)              :: element_lef_and_right_side        (2, n_conservative_values) ! 1 -> left, 2-> right

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

            pure function eos_pressure_function(specific_internal_energy, density, volume_fruction) result(pressure)
                use typedef_module
                real(real_kind), intent(in) :: specific_internal_energy
                real(real_kind), intent(in) :: density
                real(real_kind), intent(in) :: volume_fruction
                real(real_kind)             :: pressure
            end function eos_pressure_function

            pure function eos_soundspeed_function(specific_internal_energy, density, volume_fruction) result(soundspeed)
                use typedef_module
                real(real_kind), intent(in) :: specific_internal_energy
                real(real_kind), intent(in) :: density
                real(real_kind), intent(in) :: volume_fruction
                real(real_kind)             :: soundspeed
            end function eos_soundspeed_function

            pure function primitive_to_conservative_function(primitive) result(conservative)
                use typedef_module
                real(real_kind), intent(in)  :: primitive   (:)
                real(real_kind), allocatable :: conservative(:)
            end function
        end interface

        real(real_kind) :: local_coordinate_lhc_primitive(7)
        real(real_kind) :: local_coordinate_rhc_primitive(7)
        real(real_kind) :: lhc_conservative              (7)
        real(real_kind) :: rhc_conservative              (7)
        real(real_kind) :: nonviscosity_flux             (7)
        real(real_kind) :: lhc_soundspeed, lhc_pressure, lhc_density, lhc_main_velocity
        real(real_kind) :: rhc_soundspeed, rhc_pressure, rhc_density, rhc_main_velocity
        ! # for computing (- alpha1 - K) * div(u) using HLLC
        real(real_kind) :: lhc_rho1, rhc_rho1, lhc_rho2, rhc_rho2
        real(real_kind) :: ave_vel, ave_c
        ! ## HLLC wave speeds
        real(real_kind) :: s_mid
        real(real_kind) :: s_muinus
        real(real_kind) :: s_puls
        real(real_kind) :: s_left
        real(real_kind) :: s_right
        ! ## 
        real(real_kind) :: numerical_velocity(3)
        ! ## soundspeed (c) in each cells
        real(real_kind) :: lhc_c1, lhc_c2, lhc_wood_c
        real(real_kind) :: rhc_c1, rhc_c2, rhc_wood_c
        ! ## kapila`s factor
        real(real_kind) :: lhc_kapila, rhc_kapila

        ! # Allocate
        !allocate(local_coordinate_lhc_primitive(7))
        !allocate(local_coordinate_rhc_primitive(7))
        !allocate(lhc_conservative              (7))
        !allocate(rhc_conservative              (7))

        ! # compute primitive-variables of face local coordinate
        ! ## left-side
        local_coordinate_lhc_primitive(1:2) = reconstructed_leftside_primitive(1:2)
        local_coordinate_lhc_primitive(3:5) = rotate_vector( &
            reconstructed_leftside_primitive(3:5),        &
            face_normal_vector,             &
            face_tangential1_vector,        &
            face_tangential2_vector         &
        )
        local_coordinate_lhc_primitive(6:7) = reconstructed_leftside_primitive(6:7)
        ! ## right-side
        local_coordinate_rhc_primitive(1:2) = reconstructed_rightside_primitive(1:2)
        local_coordinate_rhc_primitive(3:5) = rotate_vector( &
            reconstructed_rightside_primitive(3:5),       &
            face_normal_vector,             &
            face_tangential1_vector,        &
            face_tangential2_vector         &
        )
        local_coordinate_rhc_primitive(6:7) = reconstructed_rightside_primitive(6:7)

        ! # sanity check
        ! ## left-hand cell
        if(local_coordinate_lhc_primitive(7) > 1.d0) local_coordinate_lhc_primitive(7) = 1.d0
        if(local_coordinate_lhc_primitive(7) < 0.d0) local_coordinate_lhc_primitive(7) = 0.d0
        ! ## right-hand cell
        if(local_coordinate_rhc_primitive(7) > 1.d0) local_coordinate_rhc_primitive(7) = 1.d0
        if(local_coordinate_rhc_primitive(7) < 0.d0) local_coordinate_rhc_primitive(7) = 0.d0

        ! # compute conservative-variables
        lhc_conservative = primitive_to_conservative_function( &
            local_coordinate_lhc_primitive                     &
        )
        rhc_conservative = primitive_to_conservative_function( &
            local_coordinate_rhc_primitive                     &
        )

        ! # compute EoS, main velosity, and fluxs
        associate(                                         &
                rho1_z1 => local_coordinate_lhc_primitive(1), &
                rho2_z2 => local_coordinate_lhc_primitive(2), &
                u       => local_coordinate_lhc_primitive(3), &
                v       => local_coordinate_lhc_primitive(4), &
                w       => local_coordinate_lhc_primitive(5), &
                ie      => local_coordinate_lhc_primitive(6), &
                z1      => local_coordinate_lhc_primitive(7)  &
            )
            lhc_density    = rho1_z1 + rho2_z2
            lhc_pressure   = eos_pressure_function  (ie, lhc_density, z1)
            lhc_soundspeed = eos_soundspeed_function(ie, lhc_density, z1)
            lhc_main_velocity = u

            lhc_c1 = eos_soundspeed_function(ie, lhc_density, 1.d0)
            lhc_c2 = eos_soundspeed_function(ie, lhc_density, 0.d0)
        end associate
        associate(                        &
                rho1_z1 => local_coordinate_rhc_primitive(1), &
                rho2_z2 => local_coordinate_rhc_primitive(2), &
                u       => local_coordinate_rhc_primitive(3), &
                v       => local_coordinate_rhc_primitive(4), &
                w       => local_coordinate_rhc_primitive(5), &
                ie      => local_coordinate_rhc_primitive(6), &
                z1      => local_coordinate_rhc_primitive(7)  &
            )
            rhc_density    = rho1_z1 + rho2_z2
            rhc_pressure   = eos_pressure_function  (ie, rhc_density, z1)
            rhc_soundspeed = eos_soundspeed_function(ie, rhc_density, z1)
            rhc_main_velocity = u

            rhc_c1 = eos_soundspeed_function(ie, rhc_density, 1.d0)
            rhc_c2 = eos_soundspeed_function(ie, rhc_density, 0.d0)
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

        ! # integrate nonviscosity-flux
        element_lef_and_right_side(1, :) = (-1.d0 / leftside_cell_volume ) * nonviscosity_flux(:) * face_area
        element_lef_and_right_side(2, :) = (+1.d0 / rightside_cell_volume) * nonviscosity_flux(:) * face_area

        ! # Compute (- alpha1 - K) * div(u) according to [Schmidmayer 2020, JCP] using HLLC.
        ! # If you choice other Riemann solver, term (- alpha1 - K) * div(u) is computed by HLLC forcely. Sorry.
        associate(                                               &
            lhc_rho1_z1 => local_coordinate_lhc_primitive(1)   , &
            lhc_rho2_z2 => local_coordinate_lhc_primitive(2)   , &
            lhc_u       => local_coordinate_lhc_primitive(3)   , &
            lhc_v       => local_coordinate_lhc_primitive(4)   , &
            lhc_w       => local_coordinate_lhc_primitive(5)   , &
            lhc_z1      => local_coordinate_lhc_primitive(7)   , &
            lhc_rho     => lhc_density                         , &
            lhc_p       => lhc_pressure                        , &
            lhc_c       => lhc_soundspeed                      , &
            rhc_rho1_z1 => local_coordinate_rhc_primitive(1)   , &
            rhc_rho2_z2 => local_coordinate_rhc_primitive(2)   , &
            rhc_u       => local_coordinate_rhc_primitive(3)   , &
            rhc_v       => local_coordinate_rhc_primitive(4)   , &
            rhc_w       => local_coordinate_rhc_primitive(5)   , &
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
            numerical_velocity(1) &
                = 0.5d0 * (1.d0 + sign(1.d0, s_mid)) * (lhc_u + s_muinus * ((s_left  - lhc_u) / (s_left  - s_mid) - 1.d0)) &
                + 0.5d0 * (1.d0 - sign(1.d0, s_mid)) * (rhc_u + s_puls   * ((s_right - rhc_u) / (s_right - s_mid) - 1.d0))
            numerical_velocity(2) &
                = 0.5d0 * (1.d0 + sign(1.d0, s_mid)) * lhc_v + 0.5d0 * (1.d0 - sign(1.d0, s_mid)) * rhc_v
            numerical_velocity(3) &
                = 0.5d0 * (1.d0 + sign(1.d0, s_mid)) * lhc_w + 0.5d0 * (1.d0 - sign(1.d0, s_mid)) * rhc_w
        ! ## convert to global coordinate numerical velosity
            numerical_velocity = reverse_vector(  &
                numerical_velocity              , &
                face_normal_vector              , &
                face_tangential1_vector         , &
                face_tangential2_vector           &
            )
        ! ## Compute kapira
            !if(lhc_z1 < 0.d0 + 1.d-8)then
            !    lhc_rho1 = 0.d0
            !    lhc_rho2 = lhc_rho2_z2
            !    lhc_wood_c = 1.d0 / (lhc_rho2 * lhc_c2**2)
            !    lhc_kapila = (lhc_rho2 * lhc_c2**2) * lhc_wood_c
            !else if (1.d0 - 1.d-8 < lhc_z1)then
            !    lhc_rho1 = lhc_rho1_z1
            !    lhc_rho2 = 0.d0
            !    lhc_wood_c = lhc_z1 / (lhc_rho1 * lhc_c1**2)
            !    lhc_kapila = (lhc_rho1 * lhc_c1**2) * lhc_wood_c
            !else
            !    lhc_rho1 = lhc_rho1_z1 / lhc_z1
            !    lhc_rho2 = lhc_rho2_z2 / (1.d0 - lhc_z1)
            !    lhc_wood_c = lhc_z1 / (lhc_rho1 * lhc_c1**2) + (1.d0 - lhc_z1) / (lhc_rho2 * lhc_c2**2)
            !    lhc_kapila = (lhc_rho2 * lhc_c2**2 - lhc_rho1 * lhc_c1**2) * lhc_wood_c
            !endif
            !if(rhc_z1 < 0.d0 + 1.d-8)then
            !    rhc_rho1 = 0.d0
            !    rhc_rho2 = rhc_rho2_z2
            !    rhc_wood_c = 1.d0 / (rhc_rho2 * rhc_c2**2)
            !    rhc_kapila = (rhc_rho2 * rhc_c2**2) * rhc_wood_c
            !else if (1.d0 - 1.d-8 < rhc_z1)then
            !    rhc_rho1 = rhc_rho1_z1
            !    rhc_rho2 = 0.d0
            !    rhc_wood_c = rhc_z1 / (rhc_rho1 * rhc_c1**2)
            !    rhc_kapila = (rhc_rho1 * rhc_c1**2) * rhc_wood_c
            !else
            !    rhc_rho1 = rhc_rho1_z1 / rhc_z1
            !    rhc_rho2 = rhc_rho2_z2 / (1.d0 - rhc_z1)
            !    rhc_wood_c = rhc_z1 / (rhc_rho1 * rhc_c1**2) + (1.d0 - rhc_z1) / (rhc_rho2 * rhc_c2**2)
            !    rhc_kapila = (rhc_rho2 * rhc_c2**2 - rhc_rho1 * rhc_c1**2) * rhc_wood_c
            !endif
        ! ## Integrate (- alpha1 - K) * div(u)
            element_lef_and_right_side(1, 7) = element_lef_and_right_side(1, 7) &
                                             + (1.d0 / leftside_cell_volume ) * (-lhc_z1) * multiply_vector(numerical_velocity, face_normal_vector) * face_area
            element_lef_and_right_side(2, 7) = element_lef_and_right_side(2, 7) &
                                             + (1.d0 / rightside_cell_volume) * (-rhc_z1) * multiply_vector(numerical_velocity, face_normal_vector) * face_area
        end associate
    end function compute_space_element_five_equation_model
end module five_equation_space_model_module