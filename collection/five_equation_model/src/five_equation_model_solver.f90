program five_eq_model_solver
    ! Utils
    use typedef_module
    ! Cell system
    use class_cellsystem
    use five_equation_model_variables_module
    use five_equation_model_boundary_condition_module
    use five_equation_model_space_discretization
    ! Configuration
    use class_json_configuration
    ! EoS
    use class_stiffened_gas_eos
    ! Gradient
    use class_green_gauss
    ! Divergence
    use class_gauss_divergence
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
    ! Time incriment control
    use class_constant_time_incriment_controller
    ! Parallel computing support
    use class_openmp_parallelizer
    ! Line plotter
    use class_line_plotter

    implicit none

    type(cellsystem                        ) :: a_cellsystem
    type(nlgrid_parser                     ) :: a_grid_parser
    type(nlinit_parser                     ) :: an_initial_condition_parser
    type(json_configuration                ) :: a_configuration
    type(stiffened_gas_eos                 ) :: an_eos
    type(hllc                              ) :: a_riemann_solver
    type(green_gauss                       ) :: a_gradient_calculator
    type(gauss_divergence                  ) :: a_divergence_calculator
    type(vtk_result_writer                 ) :: a_result_writer
    type(end_time_criterion                ) :: a_termination_criterion
    type(constant_time_incriment_controller) :: a_time_incriment_controller
    type(openmp_parallelizer               ) :: a_parallelizer
    type(line_plotter                      ) :: a_line_plotter

    ! These schemes can be cahnge from configuration file you set.
    class(time_stepping), pointer :: a_time_stepping
    class(reconstructor), pointer :: a_reconstructor

    ! Loop index
    integer(int_kind) :: state_num, cell_index

    ! Additional results
    real(real_kind), allocatable :: density(:)

    ! Read config
    call a_configuration%parse("config.json")

    ! Allocate schemes
    call default_time_stepping_generator            (a_configuration, a_time_stepping)
    call five_equation_model_reconstructor_generator(a_configuration, a_reconstructor)

    ! Support for parallel computing
    call a_parallelizer%initialize(a_configuration)

    ! Read grid
    call a_cellsystem%read(a_grid_parser, a_configuration)

    ! Initialize variables
    call a_cellsystem%initialize(conservative_variables_set   , num_conservative_variables   )
    call a_cellsystem%initialize(residual_set                 , num_conservative_variables   )
    call a_cellsystem%initialize(primitive_variables_set      , num_primitive_variables      )
    call a_cellsystem%initialize(surface_tension_variables_set, num_surface_tension_variables)
    ! Initialize schemes & utils
    call a_cellsystem%initialize(an_eos                     , a_configuration, num_conservative_variables, num_primitive_variables)
    call a_cellsystem%initialize(a_riemann_solver           , a_configuration, num_conservative_variables, num_primitive_variables)
    call a_cellsystem%initialize(a_time_stepping            , a_configuration, num_conservative_variables, num_primitive_variables)
    call a_cellsystem%initialize(a_reconstructor            , a_configuration, num_conservative_variables, num_primitive_variables)
    call a_cellsystem%initialize(a_gradient_calculator      , a_configuration, num_conservative_variables, num_primitive_variables)
    call a_cellsystem%initialize(a_divergence_calculator    , a_configuration, num_conservative_variables, num_primitive_variables)
    call a_cellsystem%initialize(a_result_writer            , a_configuration, num_conservative_variables, num_primitive_variables)
    call a_cellsystem%initialize(a_termination_criterion    , a_configuration, num_conservative_variables, num_primitive_variables)
    call a_cellsystem%initialize(a_time_incriment_controller, a_configuration, num_conservative_variables, num_primitive_variables)
    call a_cellsystem%initialize(a_line_plotter             , a_configuration, num_conservative_variables, num_primitive_variables)

    ! Set initial condition
    call a_cellsystem%read_initial_condition(an_initial_condition_parser, a_configuration, conservative_variables_set)
    call a_cellsystem%conservative_to_primitive_variables_all(an_eos, conservative_variables_set, primitive_variables_set, num_primitive_variables, conservative_to_primitive)

    ! Allocate
    allocate(density(a_cellsystem%get_number_of_cells()))

    ! Timestepping loop
    do while ( .not. a_cellsystem%satisfies_termination_criterion(a_termination_criterion) )
        call a_cellsystem%update_time_incriment(a_time_incriment_controller, an_eos, primitive_variables_set, spectral_radius)

        if ( a_cellsystem%is_writable(a_result_writer) ) then
            call a_cellsystem%open_file   (a_result_writer)

            print *, "Write a result "//a_cellsystem%get_filename(a_result_writer)//"..."

            do cell_index = 1, a_cellsystem%get_number_of_cells(), 1
                associate(                                          &
                    rho1 => primitive_variables_set(1, cell_index), &
                    rho2 => primitive_variables_set(2, cell_index), &
                    z    => primitive_variables_set(7, cell_index)  &
                )
                    density(cell_index) = z * rho1 + (1.d0 - z) * rho2
                end associate
            end do

            call a_cellsystem%write_scolar(a_result_writer, "Density"        , density                        )
            call a_cellsystem%write_scolar(a_result_writer, "Density 1"      , primitive_variables_set(1  , :))
            call a_cellsystem%write_scolar(a_result_writer, "Density 2"      , primitive_variables_set(2  , :))
            call a_cellsystem%write_vector(a_result_writer, "Velocity"       , primitive_variables_set(3:5, :))
            call a_cellsystem%write_scolar(a_result_writer, "Pressure"       , primitive_variables_set(6  , :))
            call a_cellsystem%write_scolar(a_result_writer, "Volume fraction", primitive_variables_set(7  , :))

            call a_cellsystem%write_vector(a_result_writer, "Gradient volume fraction", surface_tension_variables_set(1:3, :))
            call a_cellsystem%write_scolar(a_result_writer, "Curvature"               , surface_tension_variables_set(  4, :))

            call a_cellsystem%close_file  (a_result_writer)
        end if

        if ( a_cellsystem%is_writable(a_line_plotter) ) then
            call a_cellsystem%write(a_line_plotter, primitive_variables_set)
        end if

        print *, "Step "          , a_cellsystem%get_number_of_steps(), ", ", &
                 "Time incriment ", a_cellsystem%get_time_increment() , ", ", &
                 "Time "          , a_cellsystem%get_time()

        call a_cellsystem%prepare_stepping(a_time_stepping, conservative_variables_set, primitive_variables_set, residual_set)

        do state_num = 1, a_cellsystem%get_number_of_states(a_time_stepping), 1
            call a_cellsystem%apply_empty_condition    (primitive_variables_set, num_primitive_variables, rotate_primitive, unrotate_primitive, empty_bc    )
            call a_cellsystem%apply_outflow_condition  (primitive_variables_set, num_primitive_variables, rotate_primitive, unrotate_primitive, outflow_bc  )
            call a_cellsystem%apply_slipwall_condition (primitive_variables_set, num_primitive_variables, rotate_primitive, unrotate_primitive, slipwall_bc )
            call a_cellsystem%apply_symmetric_condition(primitive_variables_set, num_primitive_variables, rotate_primitive, unrotate_primitive, symmetric_bc)

            ! Compute normalized volume flaction
            call a_cellsystem%compute_gradient(a_gradient_calculator, primitive_variables_set(7,:), surface_tension_variables_set(1:3, :))
            call a_cellsystem%processes_variables_set(surface_tension_variables_set, primitive_variables_set, num_surface_tension_variables, normarize_gradient_volume_fraction)

            ! Apply BC for normalized volume flaction
            call a_cellsystem%apply_empty_condition    (surface_tension_variables_set(1:3, :), 3, rotate_gradient_value, unrotate_gradient_value, gradient_volume_fraction_bc)
            call a_cellsystem%apply_outflow_condition  (surface_tension_variables_set(1:3, :), 3, rotate_gradient_value, unrotate_gradient_value, gradient_volume_fraction_bc)
            call a_cellsystem%apply_slipwall_condition (surface_tension_variables_set(1:3, :), 3, rotate_gradient_value, unrotate_gradient_value, gradient_volume_fraction_bc)
            call a_cellsystem%apply_symmetric_condition(surface_tension_variables_set(1:3, :), 3, rotate_gradient_value, unrotate_gradient_value, gradient_volume_fraction_bc)

            ! Compute (negative) curvature
            call a_cellsystem%compute_divergence(a_divergence_calculator, surface_tension_variables_set(1:3, :), surface_tension_variables_set(4, :))

            ! Compute puressure jump induced surface tension effect
            call a_cellsystem%processes_variables_set(primitive_variables_set, surface_tension_variables_set, num_primitive_variables, compute_pressure_jump)

            call a_cellsystem%apply_empty_condition    (primitive_variables_set, num_primitive_variables, rotate_primitive, unrotate_primitive, empty_bc    )
            call a_cellsystem%apply_outflow_condition  (primitive_variables_set, num_primitive_variables, rotate_primitive, unrotate_primitive, outflow_bc  )
            call a_cellsystem%apply_slipwall_condition (primitive_variables_set, num_primitive_variables, rotate_primitive, unrotate_primitive, slipwall_bc )
            call a_cellsystem%apply_symmetric_condition(primitive_variables_set, num_primitive_variables, rotate_primitive, unrotate_primitive, symmetric_bc)

            call a_cellsystem%compute_residual(         &
                a_reconstructor                       , &
                a_riemann_solver                      , &
                an_eos                                , &
                primitive_variables_set               , &
                residual_set                          , &
                num_conservative_variables            , &
                num_primitive_variables               , &
                primitive_to_conservative             , &
                five_equation_model_residual_element    &
            )

            call a_cellsystem%compute_next_state(a_time_stepping, an_eos, state_num, conservative_variables_set, primitive_variables_set, residual_set, num_primitive_variables, conservative_to_primitive)
        end do

        call a_cellsystem%incriment_time()
    end do

    if ( a_cellsystem%is_writable(a_result_writer) ) then
        call a_cellsystem%open_file   (a_result_writer)

        print *, "Write a result "//a_cellsystem%get_filename(a_result_writer)//"..."

        do cell_index = 1, a_cellsystem%get_number_of_cells(), 1
            associate(                                          &
                rho1 => primitive_variables_set(1, cell_index), &
                rho2 => primitive_variables_set(2, cell_index), &
                z    => primitive_variables_set(7, cell_index)  &
            )
                density(cell_index) = z * rho1 + (1.d0 - z) * rho2
            end associate
        end do

        call a_cellsystem%write_scolar(a_result_writer, "Density"        , density                        )
        call a_cellsystem%write_scolar(a_result_writer, "Density 1"      , primitive_variables_set(1  , :))
        call a_cellsystem%write_scolar(a_result_writer, "Density 2"      , primitive_variables_set(2  , :))
        call a_cellsystem%write_vector(a_result_writer, "Velocity"       , primitive_variables_set(3:5, :))
        call a_cellsystem%write_scolar(a_result_writer, "Pressure"       , primitive_variables_set(6  , :))
        call a_cellsystem%write_scolar(a_result_writer, "Volume fraction", primitive_variables_set(7  , :))
        call a_cellsystem%write_vector(a_result_writer, "Gradient volume fraction", surface_tension_variables_set(1:3, :))
        call a_cellsystem%write_scolar(a_result_writer, "Curvature"               , surface_tension_variables_set(  4, :))
        call a_cellsystem%close_file  (a_result_writer)
    end if

    call a_result_writer%cleanup()

    print *, "done..."
end program five_eq_model_solver