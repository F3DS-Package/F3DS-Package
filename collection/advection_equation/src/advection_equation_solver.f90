module advection_equation_module
    use typedef_module
    use vector_module
    use stdio_module
    use class_cellsystem
    use abstract_configuration
    use abstract_result_writer

    implicit none
    private

    real(real_kind) :: speed(3)

    public :: initialize_model
    public :: write_result
    public :: spectral_radius
    public :: flux_calculation_operator
    public :: rotate_flux
    public :: unrotate_flux
    public :: empty_bc
    public :: symmetric_bc

contains

    subroutine initialize_model(a_configuration)
        class(configuration) :: a_configuration

        logical :: found

        call a_configuration%get_real  ( "Model.Speed.x", speed(1), found)
        if(.not. found) call call_error("'Model.Speed.x' is not found in the configuration file you set.")
        call a_configuration%get_real  ( "Model.Speed.y", speed(2), found)
        if(.not. found) call call_error("'Model.Speed.y' is not found in the configuration file you set.")
        call a_configuration%get_real  ( "Model.Speed.z", speed(3), found)
        if(.not. found) call call_error("'Model.Speed.z' is not found in the configuration file you set.")
    end subroutine initialize_model

    subroutine write_result(a_cellsystem, a_result_writer, conservative_variables_set)
        type (cellsystem       ), intent(inout) :: a_cellsystem
        class(result_writer    ), intent(inout) :: a_result_writer
        real (real_kind        ), intent(inout) :: conservative_variables_set   (:)
        call a_cellsystem%open_file   (a_result_writer)
        call a_cellsystem%write_scolar(a_result_writer, "Phi", conservative_variables_set)
        call a_cellsystem%close_file  (a_result_writer)
    end subroutine

    pure function spectral_radius(variable, length) result(r)
        real   (real_kind), intent(in) :: variable
        real   (real_kind), intent(in) :: length
        real   (real_kind) :: r
        r = vector_magnitude(speed(:))
    end function

    subroutine flux_calculation_operator(flux, phi)
        real(real_kind), intent(inout) :: flux(:)
        real(real_kind), intent(in   ) :: phi
        flux(:) = -1.0_real_kind * speed(:) * phi
    end subroutine

    pure function rotate_flux(      &
        global_coodinate_flux     , &
        face_normal_vector        , &
        face_tangential1_vector   , &
        face_tangential2_vector   , &
        num_flux_variables                ) result(face_coordinate_flux)

        real   (real_kind     ), intent(in)  :: global_coodinate_flux   (:)
        real   (real_kind     ), intent(in)  :: face_normal_vector      (3)
        real   (real_kind     ), intent(in)  :: face_tangential1_vector (3)
        real   (real_kind     ), intent(in)  :: face_tangential2_vector (3)
        integer(int_kind      ), intent(in)  :: num_flux_variables
        real   (real_kind     )              :: face_coordinate_flux    (num_flux_variables)

        face_coordinate_flux(1:3) = vector_rotate(global_coodinate_flux(1:3), face_normal_vector, face_tangential1_vector, face_tangential2_vector)
    end function rotate_flux

    pure function unrotate_flux(    &
        face_coordinate_flux      , &
        face_normal_vector        , &
        face_tangential1_vector   , &
        face_tangential2_vector   , &
        num_flux_variables                ) result(global_coordinate_flux)

        real   (real_kind     ), intent(in)  :: face_coordinate_flux    (:)
        real   (real_kind     ), intent(in)  :: face_normal_vector      (3)
        real   (real_kind     ), intent(in)  :: face_tangential1_vector (3)
        real   (real_kind     ), intent(in)  :: face_tangential2_vector (3)
        integer(int_kind      ), intent(in)  :: num_flux_variables
        real   (real_kind     )              :: global_coordinate_flux   (num_flux_variables)

        global_coordinate_flux(1:3) = vector_unrotate(face_coordinate_flux(1:3), face_normal_vector, face_tangential1_vector, face_tangential2_vector)
    end function unrotate_flux

    pure function empty_bc(inner_cell_flux, num_flux_variavles) result(ghost_cell_flux)
        real   (real_kind), intent(in) :: inner_cell_flux(:)
        integer(int_kind ), intent(in) :: num_flux_variavles
        real   (real_kind)             :: ghost_cell_flux(num_flux_variavles)
        ghost_cell_flux(1) =  1.0_real_kind * inner_cell_flux(1)
        ghost_cell_flux(2) =  1.0_real_kind * inner_cell_flux(2)
        ghost_cell_flux(3) =  1.0_real_kind * inner_cell_flux(3)
    end function

    pure function symmetric_bc(inner_cell_flux, num_flux_variavles) result(ghost_cell_flux)
        real   (real_kind), intent(in) :: inner_cell_flux(:)
        integer(int_kind ), intent(in) :: num_flux_variavles
        real   (real_kind)             :: ghost_cell_flux(num_flux_variavles)
        ghost_cell_flux(1) = -1.0_real_kind * inner_cell_flux(1)
        ghost_cell_flux(2) =  1.0_real_kind * inner_cell_flux(2)
        ghost_cell_flux(3) =  1.0_real_kind * inner_cell_flux(3)
    end function
end module

program advection_equation_solver
    use typedef_module
    use vector_module
    use stdio_module
    use class_cellsystem
    use class_second_order_tvd_rk
    use class_midpoint_interpolator
    use class_legacy_grid_parser
    use class_legacy_init_parser
    use class_vtk_result_writer
    use class_end_time_criterion
    use class_constant_time_increment_controller
    use class_json_configuration
    use advection_equation_module

    implicit none

    integer(int_kind), parameter :: num_conservative_variables = 1
    integer(int_kind), parameter :: flux_lenght = 3

    real(real_kind), allocatable :: conservative_variables_set   (:)
    real(real_kind), allocatable :: flux_set                  (:, :)
    real(real_kind), allocatable :: residual_set                 (:)

    type(cellsystem                        ) :: a_cellsystem
    type(legacy_grid_parser                     ) :: a_grid_parser
    type(legacy_init_parser                     ) :: a_init_parser
    type(json_configuration                ) :: a_configuration
    type(midpoint_interpolator             ) :: an_interpolator
    type(second_order_tvd_rk               ) :: a_time_stepper
    type(vtk_result_writer                 ) :: a_result_writer
    type(end_time_criterion                ) :: a_termination_criterion
    type(constant_time_increment_controller) :: a_time_increment_controller

    integer(int_kind) :: stage_num

    call a_configuration%parse("config.json")

    call initialize_model(a_configuration)

    call a_cellsystem%initialize(a_configuration)

    call a_cellsystem%read(a_grid_parser, a_configuration)

    call a_cellsystem%initialize(conservative_variables_set             )
    call a_cellsystem%initialize(flux_set                  , flux_lenght)
    call a_cellsystem%initialize(residual_set                           )

    call a_cellsystem%initialize(an_interpolator            , a_configuration, num_conservative_variables)
    call a_cellsystem%initialize(a_time_stepper             , a_configuration, num_conservative_variables)
    call a_cellsystem%initialize(a_result_writer            , a_configuration, num_conservative_variables)
    call a_cellsystem%initialize(a_termination_criterion    , a_configuration, num_conservative_variables)
    call a_cellsystem%initialize(a_time_increment_controller, a_configuration, num_conservative_variables)

    call a_cellsystem%read_initial_condition(a_init_parser, a_configuration, conservative_variables_set)

    do while ( .not. a_cellsystem%satisfy_termination_criterion(a_termination_criterion) )
        call a_cellsystem%update_time_increment(a_time_increment_controller, conservative_variables_set, spectral_radius)

        if ( a_cellsystem%is_writable(a_result_writer) ) then
            call write_result(a_cellsystem, a_result_writer, conservative_variables_set)
        end if

        call a_cellsystem%show_timestepping_infomation()

        call a_cellsystem%prepare_time_stepping(a_time_stepper, conservative_variables_set, residual_set)

        do stage_num = 1, a_cellsystem%get_number_of_stages(a_time_stepper), 1
            call a_cellsystem%operate_cellwise(flux_set, conservative_variables_set, flux_calculation_operator)

            call a_cellsystem%apply_boundary_condition(                  &
                flux_set                                               , &
                flux_lenght                                            , &
                rotate_flux                                            , &
                unrotate_flux                                          , &
                empty_condition_function        = empty_bc             , &
                symmetric_condition_function    = symmetric_bc           &
            )

            call a_cellsystem%compute_divergence(an_interpolator, flux_set, residual_set)

            call a_cellsystem%compute_next_stage(a_time_stepper, stage_num, conservative_variables_set, residual_set)
        end do

        call a_cellsystem%increment_time()
    end do

    if ( a_cellsystem%is_writable(a_result_writer) ) then
        call write_result(a_cellsystem, a_result_writer, conservative_variables_set)
    end if

    call a_result_writer%cleanup()
end program