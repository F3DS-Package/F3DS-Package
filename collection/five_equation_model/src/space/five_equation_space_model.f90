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
        n_derivative_values                               , &
        flux_function                                     , &
        eos_pressure_function                             , &
        eos_soundspeed_function                           , &
        primitive_to_conservative_function                 ) result(element)

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
        integer(int_kind ), intent(in ) :: n_derivative_values
        real   (real_kind)              :: element        (2, n_conservative_values+n_derivative_values) ! 1 -> left, 2-> right

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
        element(1, 1:7) = (-1.d0 / leftside_cell_volume ) * nonviscosity_flux(:) * face_area
        element(2, 1:7) = (+1.d0 / rightside_cell_volume) * nonviscosity_flux(:) * face_area

        ! # z1 * div(u)
        element(1, 7) = element(1, 7) &
                      - reconstructed_leftside_primitive(7)  * (1.d0 / leftside_cell_volume ) * multiply_vector(reconstructed_leftside_primitive (3:5), +1.d0 * face_normal_vector(1:3) * face_area)
        element(2, 7) = element(2, 7) &
                      - reconstructed_rightside_primitive(7) * (1.d0 / rightside_cell_volume) * multiply_vector(reconstructed_rightside_primitive(3:5), -1.d0 * face_normal_vector(1:3) * face_area)
    end function compute_space_element_five_equation_model
end module five_equation_space_model_module