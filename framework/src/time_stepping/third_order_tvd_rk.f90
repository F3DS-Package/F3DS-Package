module third_order_tvd_rk_module
    use typedef_module
    use abstract_mixture_eos

    implicit none

    private

    public :: initialize_third_order_tvd_rk
    public :: compute_next_state_third_order_tvd_rk

    real(real_kind), allocatable :: residual_set (:,:)

    real(real_kind), allocatable :: stage1_conservative_variables_set (:,:)
    real(real_kind), allocatable :: stage2_conservative_variables_set (:,:)

    contains

    subroutine initialize_third_order_tvd_rk(conservative_variables_set)
        real(real_kind), intent(in), allocatable :: conservative_variables_set (:,:)

        allocate(residual_set                     , source=conservative_variables_set)
        allocate(stage1_conservative_variables_set, source=conservative_variables_set)
        allocate(stage2_conservative_variables_set, source=conservative_variables_set)

        residual_set(:,:) = 0.d0
    end subroutine initialize_third_order_tvd_rk

    subroutine compute_next_state_third_order_tvd_rk( &
            conservative_variables_set               , &
            primitive_variables_set                  , &
            derivative_variables_set                 , &
            cell_centor_positions                    , &
            cell_volumes                             , &
            face_to_cell_index                       , &
            face_normal_vectors                      , &
            face_tangential1_vectors                 , &
            face_tangential2_vectors                 , &
            face_centor_positions                    , &
            face_areas                               , &
            outflow_face_indexs                      , &
            slipwall_face_indexs                     , &
            symmetric_face_indexs                    , &
            n_cells                                  , &
            n_faces                                  , &
            n_ghost_cells                            , &
            n_outflow_faces                          , &
            n_slipwall_faces                         , &
            n_symmetric_faces                        , &
            time_increment                           , &
            eos                                      , &
            reconstruction_function                  , &
            integrated_element_function              , &
            flux_function                            , &
            primitive_to_conservative_function       , &
            conservative_to_primitive_function       , &
            set_boundary_condition_function            &
        )

        real   (real_kind)  , intent(inout) :: conservative_variables_set (:,:)
        real   (real_kind)  , intent(inout) :: primitive_variables_set    (:,:)
        real   (real_kind)  , intent(inout) :: derivative_variables_set   (:,:)
        real   (real_kind)  , intent(in   ) :: cell_centor_positions      (:,:)
        real   (real_kind)  , intent(in   ) :: cell_volumes               (:)
        integer(int_kind )  , intent(in   ) :: face_to_cell_index         (:,:)
        real   (real_kind)  , intent(in   ) :: face_normal_vectors        (:,:)
        real   (real_kind)  , intent(in   ) :: face_tangential1_vectors   (:,:)
        real   (real_kind)  , intent(in   ) :: face_tangential2_vectors   (:,:)
        real   (real_kind)  , intent(in   ) :: face_centor_positions      (:,:)
        real   (real_kind)  , intent(in   ) :: face_areas                 (:)
        integer(int_kind )  , intent(in   ) :: outflow_face_indexs        (:)
        integer(int_kind )  , intent(in   ) :: slipwall_face_indexs       (:)
        integer(int_kind )  , intent(in   ) :: symmetric_face_indexs      (:)
        integer(int_kind )  , intent(in   ) :: n_cells
        integer(int_kind )  , intent(in   ) :: n_faces
        integer(int_kind )  , intent(in   ) :: n_ghost_cells
        integer(int_kind )  , intent(in   ) :: n_outflow_faces
        integer(int_kind )  , intent(in   ) :: n_slipwall_faces
        integer(int_kind )  , intent(in   ) :: n_symmetric_faces
        real   (real_kind)  , intent(in   ) :: time_increment
        class  (mixture_eos), intent(in   ) :: eos

        interface
            pure function reconstruction_function(  &
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

                use typedef_module
                use abstract_mixture_eos

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

                real   (real_kind) :: element(2, n_conservative_values+n_derivative_values)

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

                        real   (real_kind  ), intent(in ) :: reconstructed_lhc_primitive (:)
                        real   (real_kind  ), intent(in ) :: reconstructed_rhc_primitive (:)
                        real   (real_kind  ), intent(in ) :: lhc_primitive               (:)
                        real   (real_kind  ), intent(in ) :: rhc_primitive               (:)
                        real   (real_kind  ), intent(in ) :: lhc_cell_volume
                        real   (real_kind  ), intent(in ) :: rhc_cell_volume
                        real   (real_kind  ), intent(in ) :: face_normal_vector                (3)
                        real   (real_kind  ), intent(in ) :: face_tangential1_vector           (3)
                        real   (real_kind  ), intent(in ) :: face_tangential2_vector           (3)
                        real   (real_kind  ), intent(in ) :: face_area
                        integer(int_kind   ), intent(in ) :: n_conservative_values
                        integer(int_kind   ), intent(in ) :: n_derivative_values
                        class  (mixture_eos), intent(in ) :: eos
                        real   (real_kind)                :: element        (2, n_conservative_values+n_derivative_values)

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
            end function reconstruction_function

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

                real   (real_kind  ), intent(in ) :: reconstructed_lhc_primitive (:)
                real   (real_kind  ), intent(in ) :: reconstructed_rhc_primitive (:)
                real   (real_kind  ), intent(in ) :: lhc_primitive               (:)
                real   (real_kind  ), intent(in ) :: rhc_primitive               (:)
                real   (real_kind  ), intent(in ) :: lhc_cell_volume
                real   (real_kind  ), intent(in ) :: rhc_cell_volume
                real   (real_kind  ), intent(in ) :: face_normal_vector                (3)
                real   (real_kind  ), intent(in ) :: face_tangential1_vector           (3)
                real   (real_kind  ), intent(in ) :: face_tangential2_vector           (3)
                real   (real_kind  ), intent(in ) :: face_area
                integer(int_kind   ), intent(in ) :: n_conservative_values
                integer(int_kind   ), intent(in ) :: n_derivative_values
                class  (mixture_eos), intent(in ) :: eos
                real   (real_kind)                :: element        (2, n_conservative_values+n_derivative_values)

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
                        real(real_kind), intent(in) :: left_conservative(:)
                        real(real_kind), intent(in) :: left_main_velocity
                        real(real_kind), intent(in) :: left_density
                        real(real_kind), intent(in) :: left_pressure
                        real(real_kind), intent(in) :: left_soundspeed
                        real(real_kind), intent(in) :: right_conservative(:)
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

            pure function conservative_to_primitive_function(conservative, eos) result(primitive)
                use typedef_module
                use abstract_mixture_eos
                real (real_kind  ), intent(in)  :: conservative(:)
                class(mixture_eos), intent(in)  :: eos
                real (real_kind  ), allocatable :: primitive   (:)
            end function conservative_to_primitive_function

            function set_boundary_condition_function( &
                    primitive_variables_set   , &
                    face_to_cell_index        , &
                    face_normal_vectors       , &
                    face_tangential1_vectors  , &
                    face_tangential2_vectors  , &
                    outflow_face_indexs       , &
                    slipwall_face_indexs      , &
                    symmetric_face_indexs     , &
                    num_outflow_faces         , &
                    num_slipwall_faces        , &
                    num_symmetric_faces       , &
                    n_ghost_cells               &
                ) result(error)
                use typedef_module
                real   (real_kind), intent(inout) :: primitive_variables_set    (:,:)
                integer(int_kind ), intent(in   ) :: face_to_cell_index         (:,:)
                real   (real_kind), intent(in   ) :: face_normal_vectors        (:,:)
                real   (real_kind), intent(in   ) :: face_tangential1_vectors   (:,:)
                real   (real_kind), intent(in   ) :: face_tangential2_vectors   (:,:)
                integer(int_kind ), intent(in   ) :: outflow_face_indexs        (:)
                integer(int_kind ), intent(in   ) :: slipwall_face_indexs       (:)
                integer(int_kind ), intent(in   ) :: symmetric_face_indexs      (:)
                integer(int_kind ), intent(in   ) :: num_outflow_faces
                integer(int_kind ), intent(in   ) :: num_slipwall_faces
                integer(int_kind ), intent(in   ) :: num_symmetric_faces
                integer(int_kind ), intent(in   ) :: n_ghost_cells
                logical                           :: error
            end function
        end interface

        ! code
        integer(int_kind ) :: i, j
        integer(int_kind ) :: lhc_index, rhc_index
        integer(int_kind ) :: n_conservative_values
        integer(int_kind ) :: n_derivative_values
        real   (real_kind) :: element(2, size(conservative_variables_set(1,:)) + size(derivative_variables_set(1,:)))
        logical :: err

        n_conservative_values = size(conservative_variables_set(1,:))
        n_derivative_values   = size(derivative_variables_set  (1,:))

        !derivative_variables_set(:, :) = 0.d0

        err = set_boundary_condition_function( &
            primitive_variables_set   , &
            face_to_cell_index , &
            face_normal_vectors       , &
            face_tangential1_vectors  , &
            face_tangential2_vectors  , &
            outflow_face_indexs       , &
            slipwall_face_indexs      , &
            symmetric_face_indexs     , &
            n_outflow_faces           , &
            n_slipwall_faces          , &
            n_symmetric_faces         , &
            n_ghost_cells               &
        )

!$omp parallel do private(j, lhc_index, rhc_index, element)
        do j = 1, n_faces, 1
            lhc_index = face_to_cell_index(j, n_ghost_cells+0)
            rhc_index = face_to_cell_index(j, n_ghost_cells+1)
            element = reconstruction_function( &
                primitive_variables_set            , &
                cell_centor_positions              , &
                cell_volumes                       , &
                face_to_cell_index                 , &
                face_centor_positions              , &
                face_normal_vectors                , &
                face_tangential1_vectors           , &
                face_tangential2_vectors           , &
                face_areas                         , &
                j                                  , &
                n_conservative_values              , &
                n_derivative_values                , &
                eos                                , &
                flux_function                      , &
                primitive_to_conservative_function , &
                integrated_element_function          &
            )
            residual_set(lhc_index, :) = residual_set(lhc_index, :) + element(1, 1:n_conservative_values)
            residual_set(rhc_index, :) = residual_set(rhc_index, :) + element(2, 1:n_conservative_values)
            derivative_variables_set(lhc_index, :) = derivative_variables_set(lhc_index, :) + element(1, n_conservative_values+1:n_conservative_values+n_derivative_values)
            derivative_variables_set(rhc_index, :) = derivative_variables_set(rhc_index, :) + element(2, n_conservative_values+1:n_conservative_values+n_derivative_values)
        end do

!$omp parallel do private(i)
        do i = 1, n_cells, 1
            stage1_conservative_variables_set(i, :) = conservative_variables_set(i, :) &
                + time_increment * residual_set(i, :)
            primitive_variables_set(i, :) = conservative_to_primitive_function(stage1_conservative_variables_set(i, :), eos)
            residual_set(i, :) = 0.d0
            derivative_variables_set(i, :) = 0.d0
        end do

        err = set_boundary_condition_function( &
            primitive_variables_set   , &
            face_to_cell_index , &
            face_normal_vectors       , &
            face_tangential1_vectors  , &
            face_tangential2_vectors  , &
            outflow_face_indexs       , &
            slipwall_face_indexs      , &
            symmetric_face_indexs     , &
            n_outflow_faces           , &
            n_slipwall_faces          , &
            n_symmetric_faces         , &
            n_ghost_cells               &
        )

!$omp parallel do private(j, lhc_index, rhc_index, element)
        do j = 1, n_faces, 1
            lhc_index = face_to_cell_index(j, n_ghost_cells+0)
            rhc_index = face_to_cell_index(j, n_ghost_cells+1)
            element = reconstruction_function( &
                primitive_variables_set            , &
                cell_centor_positions              , &
                cell_volumes                       , &
                face_to_cell_index          , &
                face_centor_positions              , &
                face_normal_vectors                , &
                face_tangential1_vectors           , &
                face_tangential2_vectors           , &
                face_areas                         , &
                j                                  , &
                n_conservative_values              , &
                n_derivative_values                , &
                eos                                , &
                flux_function                      , &
                primitive_to_conservative_function , &
                integrated_element_function          &
            )
            residual_set(lhc_index, :) = residual_set(lhc_index, :) + element(1, 1:n_conservative_values)
            residual_set(rhc_index, :) = residual_set(rhc_index, :) + element(2, 1:n_conservative_values)
            derivative_variables_set(lhc_index, :) = derivative_variables_set(lhc_index, :) + element(1, n_conservative_values+1:n_conservative_values+n_derivative_values)
            derivative_variables_set(rhc_index, :) = derivative_variables_set(rhc_index, :) + element(2, n_conservative_values+1:n_conservative_values+n_derivative_values)
        end do

!$omp parallel do private(i)
        do i = 1, n_cells, 1
            stage2_conservative_variables_set(i, :) = 0.75d0 * conservative_variables_set(i, :) + 0.25d0 * (stage1_conservative_variables_set(i, :) + time_increment * residual_set(i, :))
            primitive_variables_set(i, :) = conservative_to_primitive_function(stage2_conservative_variables_set(i, :), eos)
            residual_set(i, :) = 0.d0
            derivative_variables_set(i, :) = 0.d0
        end do

        err = set_boundary_condition_function( &
            primitive_variables_set   , &
            face_to_cell_index , &
            face_normal_vectors       , &
            face_tangential1_vectors  , &
            face_tangential2_vectors  , &
            outflow_face_indexs       , &
            slipwall_face_indexs      , &
            symmetric_face_indexs     , &
            n_outflow_faces           , &
            n_slipwall_faces          , &
            n_symmetric_faces         , &
            n_ghost_cells               &
        )

!$omp parallel do private(j, lhc_index, rhc_index, element)
        do j = 1, n_faces, 1
            lhc_index = face_to_cell_index(j, n_ghost_cells+0)
            rhc_index = face_to_cell_index(j, n_ghost_cells+1)
            element = reconstruction_function( &
                primitive_variables_set            , &
                cell_centor_positions              , &
                cell_volumes                       , &
                face_to_cell_index          , &
                face_centor_positions              , &
                face_normal_vectors                , &
                face_tangential1_vectors           , &
                face_tangential2_vectors           , &
                face_areas                         , &
                j                                  , &
                n_conservative_values              , &
                n_derivative_values                , &
                eos                                , &
                flux_function                      , &
                primitive_to_conservative_function , &
                integrated_element_function          &
            )
            residual_set(lhc_index, :) = residual_set(lhc_index, :) + element(1, 1:n_conservative_values)
            residual_set(rhc_index, :) = residual_set(rhc_index, :) + element(2, 1:n_conservative_values)
            derivative_variables_set(lhc_index, :) = derivative_variables_set(lhc_index, :) + element(1, n_conservative_values+1:n_conservative_values+n_derivative_values)
            derivative_variables_set(rhc_index, :) = derivative_variables_set(rhc_index, :) + element(2, n_conservative_values+1:n_conservative_values+n_derivative_values)
        end do

!$omp parallel do private(i)
        do i = 1, n_cells, 1
            conservative_variables_set(i, :) = (1.d0/3.d0) * conservative_variables_set(i, :) + (2.d0/3.d0) * (stage2_conservative_variables_set(i, :) + time_increment * residual_set(i, :))
            primitive_variables_set(i, :) = conservative_to_primitive_function(conservative_variables_set(i, :), eos)
            residual_set(i, :) = 0.d0
        end do

    end subroutine compute_next_state_third_order_tvd_rk
end module third_order_tvd_rk_module