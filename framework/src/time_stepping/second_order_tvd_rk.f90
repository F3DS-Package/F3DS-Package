module second_order_tvd_rk_module
    use typedef_module

    implicit none

    private

    public :: compute_next_state_second_order_tvd_rk

    real(real_kind), allocatable :: residual_set (:,:)

    real(real_kind), allocatable :: stage1_conservative_variables_set (:,:)

    contains

    subroutine initialize_second_order_tvd_rk(conservative_variables_set)
        real(real_kind), intent(in), allocatable :: conservative_variables_set (:,:)

        allocate(residual_set                     , source=conservative_variables_set)
        allocate(stage1_conservative_variables_set, source=conservative_variables_set)
    end subroutine initialize_second_order_tvd_rk

    subroutine compute_next_state_second_order_tvd_rk( &
            conservative_variables_set               , &
            primitive_variables_set                  , &
            cell_centor_positions                    , &
            cell_volumes                             , &
            reference_cell_indexs_set                , &
            face_normal_vectors                      , &
            face_tangential1_vectors                 , &
            face_tangential2_vectors                 , &
            face_centor_positions                    , &
            face_areas                               , &
            n_cells                                  , &
            n_faces                                  , &
            time_increment                           , &
            reconstruction_function                  , &
            integrated_element_function              , &
            flux_function                            , &
            eos_pressure_function                    , &
            eos_soundspeed_function                  , &
            primitive_to_conservative_function       , &
            conservative_to_primitive_function         &
        )

        real   (real_kind), intent(inout) :: conservative_variables_set (:,:)
        real   (real_kind), intent(inout) :: primitive_variables_set    (:,:)
        real   (real_kind), intent(in   ) :: cell_centor_positions      (:,:)
        real   (real_kind), intent(in   ) :: cell_volumes               (:)
        integer(int_kind ), intent(in   ) :: reference_cell_indexs_set  (:,:)
        real   (real_kind), intent(in   ) :: face_normal_vectors        (:,:)
        real   (real_kind), intent(in   ) :: face_tangential1_vectors   (:,:)
        real   (real_kind), intent(in   ) :: face_tangential2_vectors   (:,:)
        real   (real_kind), intent(in   ) :: face_centor_positions      (:,:)
        real   (real_kind), intent(in   ) :: face_areas                 (:)
        integer(int_kind ), intent(in   ) :: n_cells
        integer(int_kind ), intent(in   ) :: n_faces
        real   (real_kind), intent(in   ) :: time_increment

        interface
            pure function reconstruction_function(  &
                primitive_values_set              , &
                cell_centor_positions             , &
                cell_volumes                      , &
                reference_cell_indexs_set         , &
                face_centor_positions             , &
                face_normal_vectors               , &
                face_tangential1_vectors          , &
                face_tangential2_vectors          , &
                face_areas                        , &
                face_index                        , &
                n_conservative_values             , &
                flux_function                     , &
                eos_pressure_function             , &
                eos_soundspeed_function           , &
                primitive_to_conservative_function, &
                integrated_element_function             ) result(element_lef_and_right_side)

                use typedef_module

                real   (real_kind), intent(in ) :: primitive_values_set     (:, :)
                real   (real_kind), intent(in ) :: cell_centor_positions    (:, :)
                real   (real_kind), intent(in ) :: cell_volumes             (:)
                integer(int_kind ), intent(in ) :: reference_cell_indexs_set(:,:)
                real   (real_kind), intent(in ) :: face_normal_vectors      (:,:)
                real   (real_kind), intent(in ) :: face_tangential1_vectors (:,:)
                real   (real_kind), intent(in ) :: face_tangential2_vectors (:,:)
                real   (real_kind), intent(in ) :: face_centor_positions    (:,:)
                real   (real_kind), intent(in ) :: face_areas               (:)
                integer(int_kind ), intent(in ) :: face_index
                integer(int_kind ), intent(in ) :: n_conservative_values

                real   (real_kind)                           :: element_lef_and_right_side(n_conservative_values, 2)

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

                    pure function integrated_element_function( &
                        reconstructed_leftside_primitive  , &
                        reconstructed_rightside_primitive , &
                        leftside_cell_volume              , &
                        rightside_cell_volume             , &
                        face_normal_vector                , &
                        face_tangential1_vector           , &
                        face_tangential2_vector           , &
                        face_area                         , &
                        n_conservative_values             , &
                        flux_function                     , &
                        eos_pressure_function             , &
                        eos_soundspeed_function           , &
                        primitive_to_conservative_function   ) result(element_lef_and_right_side)

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
                        real   (real_kind)              :: element_lef_and_right_side        (n_conservative_values, 2)

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
                            end function primitive_to_conservative_function
                        end interface
                    end function integrated_element_function
                end interface
            end function reconstruction_function

            pure function integrated_element_function( &
                reconstructed_leftside_primitive  , &
                reconstructed_rightside_primitive , &
                leftside_cell_volume              , &
                rightside_cell_volume             , &
                face_normal_vector                , &
                face_tangential1_vector           , &
                face_tangential2_vector           , &
                face_area                         , &
                n_conservative_values             , &
                flux_function                     , &
                eos_pressure_function             , &
                eos_soundspeed_function           , &
                primitive_to_conservative_function   ) result(element_lef_and_right_side)

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
                real   (real_kind)              :: element_lef_and_right_side        (n_conservative_values, 2)

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
                    end function primitive_to_conservative_function
                end interface
            end function integrated_element_function

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
            end function primitive_to_conservative_function

            pure function conservative_to_primitive_function(conservative) result(primitive)
                use typedef_module
                real(real_kind), intent(in)  :: conservative(:)
                real(real_kind), allocatable :: primitive   (:)
            end function conservative_to_primitive_function
        end interface

        ! code
        integer(int_kind ) :: i, j
        integer(int_kind ) :: lhc_index, rhc_index
        integer(int_kind ) :: n_conservative_values
        real   (real_kind) :: element_lef_and_right_side(size(conservative_variables_set(:,0)), 2)

        n_conservative_values = size(conservative_variables_set(:,0))

        do j = 0, n_faces, 1
            element_lef_and_right_side = reconstruction_function( &
                primitive_variables_set            , &
                cell_centor_positions              , &
                cell_volumes                       , &
                reference_cell_indexs_set          , &
                face_centor_positions              , &
                face_normal_vectors                , &
                face_tangential1_vectors           , &
                face_tangential2_vectors           , &
                face_areas                         , &
                j                                  , &
                n_conservative_values              , &
                flux_function                      , &
                eos_pressure_function              , &
                eos_soundspeed_function            , &
                primitive_to_conservative_function , &
                integrated_element_function          &
            )
            lhc_index = reference_cell_indexs_set(0, j)
            rhc_index = reference_cell_indexs_set(1, j)
            residual_set(:, lhc_index) = residual_set(:, lhc_index) + element_lef_and_right_side(:, 1)
            residual_set(:, rhc_index) = residual_set(:, rhc_index) + element_lef_and_right_side(:, 2)
        end do

        do i = 1, n_cells, 1
            stage1_conservative_variables_set(:, i) = conservative_variables_set(:, i) &
                + time_increment * residual_set(:, i)
            primitive_variables_set(:, i) = conservative_to_primitive_function(stage1_conservative_variables_set(:, i))
        end do

        do j = 0, n_faces, 1
            element_lef_and_right_side = reconstruction_function( &
                primitive_variables_set            , &
                cell_centor_positions              , &
                cell_volumes                       , &
                reference_cell_indexs_set          , &
                face_centor_positions              , &
                face_normal_vectors                , &
                face_tangential1_vectors           , &
                face_tangential2_vectors           , &
                face_areas                         , &
                j                                  , &
                n_conservative_values              , &
                flux_function                      , &
                eos_pressure_function              , &
                eos_soundspeed_function            , &
                primitive_to_conservative_function , &
                integrated_element_function          &
            )
            lhc_index = reference_cell_indexs_set(0, j)
            rhc_index = reference_cell_indexs_set(1, j)
            residual_set(:, lhc_index) = residual_set(:, lhc_index) + element_lef_and_right_side(:, 1)
            residual_set(:, rhc_index) = residual_set(:, rhc_index) + element_lef_and_right_side(:, 2)
        end do

        do i = 1, n_cells, 1
            conservative_variables_set(:, i) = conservative_variables_set(:, i) &
                + 0.5d0 * (stage1_conservative_variables_set(:, i) + time_increment * residual_set(:, i))
            primitive_variables_set(:, i) = conservative_to_primitive_function(conservative_variables_set(:, i))
        end do

    end subroutine compute_next_state_second_order_tvd_rk
end module second_order_tvd_rk_module