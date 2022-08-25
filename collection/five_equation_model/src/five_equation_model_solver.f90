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
    use class_weighted_green_gasuu
    ! Divergence
    use class_weighted_gauss_divergence
    ! Time stepping
    use abstract_time_stepping
    use time_stepping_generator_module
    ! Reconstructor
    use abstract_reconstructor
    use reconstructor_generator_module
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

    implicit none

    type(cellsystem                        ) :: a_cellsystem
    type(nlgrid_parser                     ) :: a_grid_parser
    type(nlinit_parser                     ) :: an_initial_condition_parser
    type(json_configuration                ) :: a_configuration
    type(stiffened_gas_eos                 ) :: an_eos
    type(hllc                              ) :: a_riemann_solver
    type(weighted_green_gauss              ) :: a_gradient_calculator
    type(weighted_gauss_divergence         ) :: a_divergence_calculator
    type(vtk_result_writer                 ) :: a_result_writer
    type(end_time_criterion                ) :: a_termination_criterion
    type(constant_time_incriment_controller) :: a_time_incriment_controller
    type(openmp_parallelizer               ) :: a_parallelizer

    ! These schemes can be cahnge from configuration file you set.
    class(time_stepping), pointer :: a_time_stepping
    class(reconstructor), pointer :: a_reconstructor

    ! Loop index
    integer(int_kind) :: state_num

    ! Read config
    call a_configuration%parse("config.json")

    ! Allocate schemes
    call default_time_stepping_generator(a_configuration, a_time_stepping)
    call default_reconstructor_generator(a_configuration, a_reconstructor)

    ! Support for parallel computing
    call a_parallelizer%initialize(a_configuration)

    ! Read grid
    call a_cellsystem%read(a_grid_parser, a_configuration)

    ! Initialize variables
    call a_cellsystem%initialize(conservative_variables_set, num_conservative_variables)
    call a_cellsystem%initialize(residual_set              , num_conservative_variables)
    call a_cellsystem%initialize(primitive_variables_set   , num_primitive_variables)
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

    ! Set initial condition
    call a_cellsystem%read_initial_condition(an_initial_condition_parser, a_configuration, conservative_variables_set)
    call a_cellsystem%conservative_to_primitive_variables_all(an_eos, conservative_variables_set, primitive_variables_set, num_primitive_variables, conservative_to_primitive)

    ! Timestepping loop
    do while ( .not. a_cellsystem%satisfies_termination_criterion(a_termination_criterion) )
        call a_cellsystem%update_time_incriment(a_time_incriment_controller, an_eos, primitive_variables_set, spectral_radius)

        if ( a_cellsystem%is_writable(a_result_writer) ) then
            call a_cellsystem%open_file   (a_result_writer)

            print *, "Write a result "//a_cellsystem%get_filename(a_result_writer)//"..."

            call a_cellsystem%write_scolar(a_result_writer, "Density 1"      , primitive_variables_set(1  , :))
            call a_cellsystem%write_scolar(a_result_writer, "Density 2"      , primitive_variables_set(2  , :))
            call a_cellsystem%write_vector(a_result_writer, "Velocity"       , primitive_variables_set(3:5, :))
            call a_cellsystem%write_scolar(a_result_writer, "Pressre"        , primitive_variables_set(6  , :))
            call a_cellsystem%write_scolar(a_result_writer, "Volume flaction", primitive_variables_set(7  , :))
            call a_cellsystem%close_file  (a_result_writer)
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

        call a_cellsystem%write_scolar(a_result_writer, "Density 1"      , primitive_variables_set(1  , :))
        call a_cellsystem%write_scolar(a_result_writer, "Density 2"      , primitive_variables_set(2  , :))
        call a_cellsystem%write_vector(a_result_writer, "Velocity"       , primitive_variables_set(3:5, :))
        call a_cellsystem%write_scolar(a_result_writer, "Pressre"        , primitive_variables_set(6  , :))
        call a_cellsystem%write_scolar(a_result_writer, "Volume flaction", primitive_variables_set(7  , :))
        call a_cellsystem%close_file  (a_result_writer)
    end if

    call a_result_writer%cleanup()

    print *, "done..."
end program five_eq_model_solver