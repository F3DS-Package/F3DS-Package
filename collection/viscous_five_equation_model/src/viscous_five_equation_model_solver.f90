program viscous_five_equation_model_solver
    ! Utils
    use typedef_module
    use stdio_module
    ! Cell system
    use class_cellsystem
    ! Configuration
    use class_json_configuration
    ! EoS
    use class_stiffened_gas_eos
    ! Gradient
    use class_green_gauss
    ! Face interpolation
    use class_midpoint_interpolator
    ! Face gradient
    use class_midpoint_face_gradient_interpolator
    ! Time stepping
    use abstract_time_stepping
    use time_stepping_generator_module
    ! Reconstructor
    use abstract_reconstructor
    use five_equation_model_reconstructor_generator_module
    ! Rieman solver
    use class_hllc
    ! Grid & initial condition reader
    use class_nlinit_parser
    use class_nlgrid_parser
    ! Result file output
    use class_vtk_result_writer
    ! Termination criteria
    use class_end_time_criterion
    ! Time increment control
    use abstract_time_increment_controller
    use time_increment_controller_generator_module
    ! Parallel computing support
    use class_openmp_parallelizer
    ! Measurements
    use class_line_plotter
    use class_control_volume_profiler
    ! Model
    use viscous_five_equation_model_utils_module
    use viscous_five_equation_model_boundary_condition_module
    use viscous_five_equation_model_module

    implicit none

    integer(int_kind), parameter :: num_conservative_variables       = 7
    integer(int_kind), parameter :: num_primitive_variables          = 8
    integer(int_kind), parameter :: num_gradient_primitive_variables = 24
    integer(int_kind), parameter :: num_surface_tension_variables    = 5

    ! Elm. 1) following variables are saved
    ! conservative_variables_set(1  , :)   = Z1*rho1   : Z1*density of fluid1
    ! conservative_variables_set(2  , :)   = Z2*rho2   : Z2*density of fluid2
    ! conservative_variables_set(3:5, :)   = rho*v     : momentum vector
    ! conservative_variables_set(6  , :)   = e         : energy density
    ! conservative_variables_set(7  , :)   = Z1        : volume fraction of fluid1
    ! Elm. 2) 1 : {@code num_cells}, cell index
    real(real_kind), allocatable :: conservative_variables_set(:,:)

    ! Elm. 1) following variables are saved
    ! primitive_variables_set(1  , :)   : density of fluid1
    ! primitive_variables_set(2  , :)   : density of fluid2
    ! primitive_variables_set(3:5, :)   : velocity vector (u,v,w)
    ! primitive_variables_set(6  , :)   : pressre
    ! primitive_variables_set(7  , :)   : volume fraction of fluid1
    ! primitive_variables_set(8  , :)   : interface curvature κ
    ! Elm. 2) 1 : {@code num_cells}, cell index
    real(real_kind), allocatable :: primitive_variables_set(:,:)

    ! Elm. 1) following variables are saved
    ! gradient_primitive_variables_set(1 :3 , :) : dρ1/dx, dρ1/dy, dρ1/dz
    ! gradient_primitive_variables_set(4 :6 , :) : dρ2/dx, dρ2/dy, dρ2/dz
    ! gradient_primitive_variables_set(7 :9 , :) : du/dx, du/dy, du/dz
    ! gradient_primitive_variables_set(10:12, :) : dv/dx, dv/dy, dv/dz
    ! gradient_primitive_variables_set(13:15, :) : dw/dx, dw/dy, dw/dz
    ! gradient_primitive_variables_set(16:18, :) : dp/dx, dp/dy, dp/dz
    ! gradient_primitive_variables_set(19:21, :) : dα1/dx, dα1/dy, dα1/dz
    ! gradient_primitive_variables_set(22:24, :) : dκ/dx, dκ/dy, dκ/dz
    ! Elm. 2) 1 : {@code num_cells}, cell index
    real(real_kind), allocatable :: gradient_primitive_variables_set(:,:)

    ! Elm. 1) following variables are saved
    ! conservative_variables_set(1  , :)   = Z1*rho1   : Z1*density of fluid1
    ! conservative_variables_set(2  , :)   = Z2*rho2   : Z2*density of fluid2
    ! conservative_variables_set(3:5, :)   = rho*v     : momentum vector
    ! conservative_variables_set(6  , :)   = e         : energy density
    ! conservative_variables_set(7  , :)   = Z1        : volume fraction of fluid1
    ! Elm. 2) 1 : {@code num_cells}, cell index
    real(real_kind), allocatable :: residual_set(:,:)

    ! Elm. 1) following variables are saved
    ! surface_tension_variables_set(1  , :) = smoothed volume fraction
    ! surface_tension_variables_set(2:4, :) = normalized gradient volume fraction
    ! surface_tension_variables_set(5  , :) = curvature
    ! Elm. 2) 1 : {@code num_cells}, cell index
    real(real_kind), allocatable :: surface_tension_variables_set(:,:)

    ! F3DS Flamework
    type(cellsystem                                   ) :: a_cellsystem
    type(nlgrid_parser                                ) :: a_grid_parser
    type(nlinit_parser                                ) :: an_initial_condition_parser
    type(json_configuration                           ) :: a_configuration
    type(stiffened_gas_eos                            ) :: an_eos
    type(hllc                                         ) :: a_riemann_solver
    type(green_gauss                                  ) :: a_gradient_calculator
    type(midpoint_interpolator                        ) :: a_interpolator
    type(midpoint_face_gradient_interpolator          ) :: a_face_gradient_interpolator
    type(vtk_result_writer                            ) :: a_result_writer
    type(end_time_criterion                           ) :: a_termination_criterion
    type(openmp_parallelizer                          ) :: a_parallelizer
    type(line_plotter                                 ) :: a_line_plotter
    type(control_volume_profiler                      ) :: a_control_volume_profiler

    ! These schemes and methods can be cahnge from configuration file you set.
    class(time_stepping            ), pointer :: a_time_stepping
    class(reconstructor            ), pointer :: a_reconstructor
    class(time_increment_controller), pointer :: a_time_increment_controller

    ! Loop index
    integer(int_kind) :: stage_num, cell_index, smooth_num

    ! Read config
    call a_configuration%parse("config.json")

    ! Allocate schemes
    call f3ds_time_stepping_generator               (a_configuration, a_time_stepping            )
    call five_equation_model_reconstructor_generator(a_configuration, a_reconstructor            )
    call f3ds_time_increment_controller_generator   (a_configuration, a_time_increment_controller)

    ! Read grid
    call a_cellsystem%read(a_grid_parser, a_configuration)

    ! Initialize variables
    call a_cellsystem%initialize(conservative_variables_set      , num_conservative_variables      )
    call a_cellsystem%initialize(residual_set                    , num_conservative_variables      )
    call a_cellsystem%initialize(primitive_variables_set         , num_primitive_variables         )
    call a_cellsystem%initialize(gradient_primitive_variables_set, num_gradient_primitive_variables)
    call a_cellsystem%initialize(surface_tension_variables_set   , num_surface_tension_variables   )
    ! Initialize schemes & utils
    call a_cellsystem%initialize(a_parallelizer              , a_configuration, num_conservative_variables, num_primitive_variables)
    call a_cellsystem%initialize(an_eos                      , a_configuration, num_conservative_variables, num_primitive_variables)
    call a_cellsystem%initialize(a_riemann_solver            , a_configuration, num_conservative_variables, num_primitive_variables)
    call a_cellsystem%initialize(a_time_stepping             , a_configuration, num_conservative_variables, num_primitive_variables)
    call a_cellsystem%initialize(a_reconstructor             , a_configuration, num_conservative_variables, num_primitive_variables)
    call a_cellsystem%initialize(a_gradient_calculator       , a_configuration, num_conservative_variables, num_primitive_variables)
    call a_cellsystem%initialize(a_interpolator              , a_configuration, num_conservative_variables, num_primitive_variables)
    call a_cellsystem%initialize(a_result_writer             , a_configuration, num_conservative_variables, num_primitive_variables)
    call a_cellsystem%initialize(a_termination_criterion     , a_configuration, num_conservative_variables, num_primitive_variables)
    call a_cellsystem%initialize(a_time_increment_controller , a_configuration, num_conservative_variables, num_primitive_variables)
    call a_cellsystem%initialize(a_line_plotter              , a_configuration, num_conservative_variables, num_primitive_variables)
    call a_cellsystem%initialize(a_control_volume_profiler   , a_configuration, num_conservative_variables, num_primitive_variables)
    call a_cellsystem%initialize(a_face_gradient_interpolator, a_configuration, num_conservative_variables, num_primitive_variables)

    ! Set initial condition
    call a_cellsystem%read_initial_condition(an_initial_condition_parser, a_configuration, conservative_variables_set)
    call a_cellsystem%conservative_to_primitive_variables_all(a_parallelizer, an_eos, conservative_variables_set, primitive_variables_set, num_primitive_variables, conservative_to_primitive)

    ! Initialize model
    call initialize_model(a_configuration, primitive_variables_set, a_cellsystem%get_number_of_cells())

    ! Timestepping loop
    do while ( .not. a_cellsystem%satisfies_termination_criterion(a_termination_criterion) )
        call a_cellsystem%update_time_increment(a_time_increment_controller, an_eos, primitive_variables_set, spectral_radius)

        if ( a_cellsystem%is_writable(a_result_writer) ) then
            call write_result(a_cellsystem, a_result_writer, primitive_variables_set, surface_tension_variables_set)
        end if

        if ( a_cellsystem%is_writable(a_line_plotter) ) then
            call a_cellsystem%write(a_line_plotter, primitive_variables_set)
        end if

        if ( a_cellsystem%is_writable(a_control_volume_profiler) ) then
            call a_cellsystem%write(a_control_volume_profiler, primitive_variables_set)
        end if

        call a_cellsystem%show_timestepping_infomation()

        call a_cellsystem%prepare_stepping(a_parallelizer, a_time_stepping, conservative_variables_set, primitive_variables_set, residual_set)

        do stage_num = 1, a_cellsystem%get_number_of_stages(a_time_stepping), 1
            call a_cellsystem%apply_empty_condition             (a_parallelizer, primitive_variables_set, num_primitive_variables, rotate_primitive, unrotate_primitive, empty_bc             )
            call a_cellsystem%apply_outflow_condition           (a_parallelizer, primitive_variables_set, num_primitive_variables, rotate_primitive, unrotate_primitive, outflow_bc           )
            call a_cellsystem%apply_nonslip_wall_condition      (a_parallelizer, primitive_variables_set, num_primitive_variables, rotate_primitive, unrotate_primitive, nonslip_wall_bc      )
            call a_cellsystem%apply_slip_and_symmetric_condition(a_parallelizer, primitive_variables_set, num_primitive_variables, rotate_primitive, unrotate_primitive, slip_and_symmetric_bc)

            ! Compute gradient primitive variables
            call a_cellsystem%compute_gradient(a_parallelizer, a_gradient_calculator, primitive_variables_set(:,:), gradient_primitive_variables_set(:,:), num_primitive_variables)

            if(apply_surface_tension())then
                ! Compute normarize gradient volume fraction
                call a_cellsystem%processes_variables_set(a_parallelizer, surface_tension_variables_set, primitive_variables_set, num_surface_tension_variables, compute_smoothed_volume_fraction)
                call a_cellsystem%compute_gradient(a_parallelizer, a_gradient_calculator, surface_tension_variables_set(1,:), surface_tension_variables_set(2:4,:))
                call a_cellsystem%processes_variables_set(a_parallelizer, surface_tension_variables_set, num_surface_tension_variables, normalize_gradient_volume_fraction)

                ! Apply BC for normalized volume flaction
                call a_cellsystem%apply_empty_condition             (a_parallelizer, surface_tension_variables_set(2:4, :), 3, rotate_gradient_value, unrotate_gradient_value, surface_normal_bc)
                call a_cellsystem%apply_outflow_condition           (a_parallelizer, surface_tension_variables_set(2:4, :), 3, rotate_gradient_value, unrotate_gradient_value, surface_normal_bc)
                call a_cellsystem%apply_nonslip_wall_condition      (a_parallelizer, surface_tension_variables_set(2:4, :), 3, rotate_gradient_value, unrotate_gradient_value, surface_normal_bc_nonslip_wall)
                call a_cellsystem%apply_slip_and_symmetric_condition(a_parallelizer, surface_tension_variables_set(2:4, :), 3, rotate_gradient_value, unrotate_gradient_value, surface_normal_bc)

                ! Compute a negative heavest fluid curvature
                call a_cellsystem%compute_divergence(a_parallelizer, a_interpolator, surface_tension_variables_set(2:4, :), surface_tension_variables_set(5, :))

                ! Smoothing curvature
                call a_cellsystem%smooth_variables(a_parallelizer, surface_tension_variables_set, curvature_smoothing_weight)

                ! Curvature is copied to {@code primitive_variables_set}.
                call a_cellsystem%processes_variables_set(a_parallelizer, primitive_variables_set, surface_tension_variables_set, num_primitive_variables, curvature_preprocessing)

                ! Re apply BC for curvature
                call a_cellsystem%apply_empty_condition             (a_parallelizer, primitive_variables_set, num_primitive_variables, rotate_primitive, unrotate_primitive, empty_bc             )
                call a_cellsystem%apply_outflow_condition           (a_parallelizer, primitive_variables_set, num_primitive_variables, rotate_primitive, unrotate_primitive, outflow_bc           )
                call a_cellsystem%apply_nonslip_wall_condition      (a_parallelizer, primitive_variables_set, num_primitive_variables, rotate_primitive, unrotate_primitive, nonslip_wall_bc      )
                call a_cellsystem%apply_slip_and_symmetric_condition(a_parallelizer, primitive_variables_set, num_primitive_variables, rotate_primitive, unrotate_primitive, slip_and_symmetric_bc)
            end if

            call a_cellsystem%compute_residual(   &
                a_parallelizer                  , &
                a_reconstructor                 , &
                a_riemann_solver                , &
                an_eos                          , &
                a_face_gradient_interpolator    , &
                primitive_variables_set         , &
                gradient_primitive_variables_set, &
                residual_set                    , &
                num_conservative_variables      , &
                num_primitive_variables         , &
                primitive_to_conservative       , &
                compute_residual_element          &
            )

            call a_cellsystem%compute_source_term(a_parallelizer, primitive_variables_set, residual_set, num_conservative_variables, compute_source_term)

            call a_cellsystem%compute_next_stage(a_parallelizer, a_time_stepping, an_eos, stage_num, conservative_variables_set, primitive_variables_set, residual_set, num_primitive_variables, conservative_to_primitive)
        end do

        call a_cellsystem%increment_time()
    end do

    if ( a_cellsystem%is_writable(a_result_writer) ) then
        call write_result(a_cellsystem, a_result_writer, primitive_variables_set, surface_tension_variables_set)
    end if

    if ( a_cellsystem%is_writable(a_line_plotter) ) then
        call a_cellsystem%write(a_line_plotter, primitive_variables_set)
    end if

    if ( a_cellsystem%is_writable(a_control_volume_profiler) ) then
        call a_cellsystem%write(a_control_volume_profiler, primitive_variables_set)
    end if

    call a_result_writer%cleanup()

    call write_message("fv5eq is successfully terminated. done...")
end program viscous_five_equation_model_solver