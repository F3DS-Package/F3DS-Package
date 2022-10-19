module class_cellsystem
    use typedef_module
    use stdio_module
    use math_constant_module
    use vector_module
    use matrix_module
    use face_type_module
    use boundary_condition_module
    use class_point_id_list
    use abstract_grid_parser
    use abstract_result_writer
    use abstract_configuration
    use abstract_time_stepping
    use abstract_eos
    use abstract_reconstructor
    use abstract_riemann_solver
    use abstract_gradient_calculator
    use abstract_divergence_calculator
    use abstract_termination_criterion
    use abstract_time_increment_controller
    use abstract_initial_condition_parser
    use abstract_face_gradient_interpolator
    use abstract_model
    use class_line_plotter
    use class_control_volume_profiler

    implicit none

    private

    ! Data structure and rule:
    ! (local cell index)          1         2         3         4         5         6
    !                             |         |         |         |         |         |
    !
    !                        *---------*---------*---------*---------*---------*---------*
    !                        |         |         |         | n       |         |         | "n" is a normal vector.
    !                        |    o    |    o    |    o    x--> o    |    o    |    o    | Normal vector must be oriented to the right-hand cell.
    ! (cell index, e.g.) ->  |   101   |   102   |   103   |   104   |   105   |   106   |
    !                        *---------*---------*---------*---------*---------*---------*
    !
    !         (inner cells) ______|_________|_________|    |    |_________|_________|________ (inner or ghost cells)
    !                                                      |             *If fase is boundary, ghost cells are stored right-side.
    !                                            (boundary/inner face)
    type, public :: cellsystem
        private

        ! Time
        integer(int_kind ) :: num_steps      = 0
        real   (real_kind) :: time           = 0.d0
        real   (real_kind) :: time_increment = 0.d0

        ! Real cells only (NOT include points of ghost/dummy cells)
        integer(int_kind) :: num_points

        ! Include ghost/dummy cells.
        integer(int_kind) :: num_cells

        ! Real cells only
        integer(int_kind) :: num_faces

        ! Also known as number of ghost cells
        integer(int_kind) :: num_local_cells

        ! Number of boundary faces
        integer(int_kind) :: num_outflow_faces
        integer(int_kind) :: num_wall_faces
        integer(int_kind) :: num_symmetric_faces
        integer(int_kind) :: num_empty_faces

        ! ### Faces ###

        ! Reference for cell index
        ! Elm. 1) -Number of local cells + 1: Number of local cells , Local cell index
        ! Elm. 2) 1 : Number of faces                               , Face number
        integer(int_kind), allocatable :: face_to_cell_indexes(:,:)

        ! Normal vectors that is directed right-hand cell (local index is one)
        ! Elm. 1) 1 : 3              , vector compornents
        ! Elm. 2) 1 : Number of faces, face number
        real(real_kind), allocatable :: face_normal_vectors(:,:)

        ! Tangental vectors (1st)
        ! Elm. 1) 1 : 3              , vector compornents
        ! Elm. 2) 1 : Number of faces, face number
        real(real_kind), allocatable :: face_tangential1_vectors(:,:)

        ! Tangental vectors (2nd)
        ! Elm. 1) 1 : 3              , vector compornents
        ! Elm. 2) 1 : Number of faces, face number
        real(real_kind), allocatable :: face_tangential2_vectors(:,:)

        ! Face centor position
        ! Elm. 1) 1 : 3              , vector compornents
        ! Elm. 2) 1 : Number of faces, face number
        real(real_kind), allocatable :: face_positions(:,:)

        ! Face area
        ! Elm. 1) 1 : Number of faces, face number
        real(real_kind), allocatable :: face_areas(:)

        ! ### Cells ###

        ! Cell centor positions
        ! Elm. 1) 1 : 3                , vector compornents
        ! Elm. 2) 1 : {@code num_cells}, cell index
        real(real_kind), allocatable :: cell_centor_positions(:,:)

        ! Cell volumes
        ! Elm. 1) 1 : {@code num_cells}, cell index
        real(real_kind), allocatable :: cell_volumes(:)

        ! Elm. 1) 1 : {@code num_cells}, cell index
        logical, allocatable :: is_real_cell(:)

        ! Elm. 1) 1:3 = x, y, z
        ! Elm. 2) 1:{@code num_points}
        real(real_kind), allocatable :: points(:, :)

        ! Elm. 1) 1 : {@code num_cells}, cell index
        class(point_id_list), allocatable :: cell_geometries(:)

        ! Followed VTK cell type
        ! Elm. 1) 1 : {@code num_cells}, cell index
        integer(type_kind), allocatable :: cell_types(:)

        ! ### Boundary Condition ###

        ! Elements are stored boundary face index.
        ! Elm. 1) 1 : num_{boundary condition}_faces
        integer(int_kind), allocatable :: outflow_face_indexes  (:)
        integer(int_kind), allocatable :: wall_face_indexes (:)
        integer(int_kind), allocatable :: symmetric_face_indexes(:)
        integer(int_kind), allocatable :: empty_face_indexes    (:)

        ! ### Status ###
        logical :: read_cellsystem = .false.

        contains
        ! ### Common ###
        procedure, public, pass(self) :: result_writer_initialize
        procedure, public, pass(self) :: time_stepping_initialize
        procedure, public, pass(self) :: reconstructor_initialize
        procedure, public, pass(self) :: riemann_solver_initialize
        procedure, public, pass(self) :: eos_initialize
        procedure, public, pass(self) :: gradient_calculator_initialize
        procedure, public, pass(self) :: variables_initialize
        procedure, public, pass(self) :: variables_1darray_initialize
        procedure, public, pass(self) :: termination_criterion_initialize
        procedure, public, pass(self) :: time_increment_controller_initialize
        procedure, public, pass(self) :: divergence_calculator_initialize
        procedure, public, pass(self) :: line_plotter_initialize
        procedure, public, pass(self) :: control_volume_profiler_initialize
        procedure, public, pass(self) :: model_initialize
        procedure, public, pass(self) :: face_gradient_interpolator_initialize
        generic  , public             :: initialize  => result_writer_initialize            , &
                                                        time_stepping_initialize            , &
                                                        reconstructor_initialize            , &
                                                        riemann_solver_initialize           , &
                                                        eos_initialize                      , &
                                                        gradient_calculator_initialize      , &
                                                        variables_initialize                , &
                                                        variables_1darray_initialize        , &
                                                        termination_criterion_initialize    , &
                                                        time_increment_controller_initialize, &
                                                        divergence_calculator_initialize    , &
                                                        line_plotter_initialize             , &
                                                        control_volume_profiler_initialize  , &
                                                        model_initialize                    , &
                                                        face_gradient_interpolator_initialize
        procedure, public, pass(self) :: result_writer_open_file
        generic  , public             :: open_file   => result_writer_open_file
        procedure, public, pass(self) :: result_writer_close_file
        generic  , public             :: close_file  => result_writer_close_file
        procedure, public, pass(self) :: result_writer_is_writable
        procedure, public, pass(self) :: line_plotter_is_writable
        procedure, public, pass(self) :: control_volume_profiler_is_writable
        generic  , public             :: is_writable => result_writer_is_writable , &
                                                        line_plotter_is_writable  , &
                                                        control_volume_profiler_is_writable
        procedure, public, pass(self) :: line_plotter_write
        procedure, public, pass(self) :: control_volume_profiler_write
        generic  , public             :: write => line_plotter_write  , &
                                                  control_volume_profiler_write

        ! ### Vatiables ###
        procedure, public, pass(self) :: read_initial_condition
        procedure, public, pass(self) :: conservative_to_primitive_variables_all
        procedure, public, pass(self) :: processes_variables_set_single_array
        procedure, public, pass(self) :: processes_variables_set_two_array
        generic  , public             :: processes_variables_set => processes_variables_set_single_array, &
                                                                    processes_variables_set_two_array

        ! ### Time increment control ###
        procedure, public, pass(self) :: update_time_increment_functional_api
        procedure, public, pass(self) :: update_time_increment_objective_api
        generic  , public             :: update_time_increment => update_time_increment_functional_api, &
                                                                  update_time_increment_objective_api

        ! ### Termination criterion ###
        procedure, public, pass(self) :: satisfies_termination_criterion

        ! ### Boundary Condition ###
        procedure, public, pass(self) :: apply_outflow_condition
        procedure, public, pass(self) :: apply_wall_condition
        procedure, public, pass(self) :: apply_symmetric_condition
        procedure, public, pass(self) :: apply_empty_condition

        ! ### Gradient Calculator ###
        procedure, public, pass(self) :: compute_gradient_1darray
        procedure, public, pass(self) :: compute_gradient_2darray
        generic  , public             :: compute_gradient => compute_gradient_1darray, &
                                                             compute_gradient_2darray


        ! ### Divergence Calculator ###
        procedure, public, pass(self) :: compute_divergence_1darray
        procedure, public, pass(self) :: compute_divergence_2darray
        generic  , public             :: compute_divergence => compute_divergence_1darray, &
                                                               compute_divergence_2darray

        ! ### Resudual ###
        procedure, public, pass(self) :: compute_residual_functional_api
        procedure, public, pass(self) :: compute_residual_objective_api
        generic  , public             :: compute_residual => compute_residual_functional_api, &
                                                             compute_residual_objective_api
        procedure, public, pass(self) :: compute_source_term_objective_api
        generic  , public             :: compute_source_term => compute_source_term_objective_api

        ! ### Time Stepping ###
        procedure, public, pass(self) :: compute_next_state
        procedure, public, pass(self) :: prepare_stepping
        procedure, public, pass(self) :: get_number_of_states

        ! ### Result Writer ###
        procedure, public, pass(self) :: write_scolar
        procedure, public, pass(self) :: write_vector
        procedure, public, pass(self) :: get_filename

        ! ### Cellsystem ###
        procedure, public, pass(self) :: read
        procedure, public, pass(self) :: increment_time

        ! ### Getter ###
        procedure, public, pass(self) :: get_number_of_faces
        procedure, public, pass(self) :: get_number_of_local_cells
        procedure, public, pass(self) :: get_number_of_points
        procedure, public, pass(self) :: get_number_of_cells
        procedure, public, pass(self) :: get_number_of_outflow_faces
        procedure, public, pass(self) :: get_number_of_wall_faces
        procedure, public, pass(self) :: get_number_of_symmetric_faces
        procedure, public, pass(self) :: get_number_of_empty_faces
        procedure, public, pass(self) :: get_number_of_steps
        procedure, public, pass(self) :: get_time
        procedure, public, pass(self) :: get_time_increment

        ! ### Finalizer ###
        procedure, public, pass(self) :: finalize

        ! ### Utils ###
        procedure, private, pass(self) :: initialise_faces
        procedure, private, pass(self) :: initialise_cells
        procedure, private, pass(self) :: initialise_boundary_references
        procedure, private, pass(self) :: finalize_faces
        procedure, private, pass(self) :: finalize_cells
        procedure, private, pass(self) :: finalize_boundary_references
        procedure, private, pass(self) :: assign_boundary
        procedure, private, pass(self) :: compute_boundary_gradient ! TODO: Move it! Make 'compute_face_gradient' class!
    end type

    contains

    ! ### Face gradient interpolator ###
    subroutine face_gradient_interpolator_initialize(self, a_face_gradient_interpolator, config, num_conservative_variables, num_primitive_variables)
        class   (cellsystem                ), intent(inout) :: self
        class   (face_gradient_interpolator), intent(inout) :: a_face_gradient_interpolator
        class   (configuration             ), intent(inout) :: config
        integer (int_kind                  ), intent(in   ) :: num_conservative_variables
        integer (int_kind                  ), intent(in   ) :: num_primitive_variables
#ifdef _DEBUG
        call write_debuginfo("In face_gradient_interpolator_initialize(), cellsystem.")
#endif
        call a_face_gradient_interpolator%initialize(config)
    end subroutine face_gradient_interpolator_initialize

    ! ### Model ###
    subroutine model_initialize(self, a_model, a_config, num_conservative_variables, num_primitive_variables)
        class   (cellsystem   ), intent(inout) :: self
        class   (model        ), intent(inout) :: a_model
        class   (configuration), intent(inout) :: a_config
        integer(int_kind      ), intent(in   ) :: num_conservative_variables
        integer(int_kind      ), intent(in   ) :: num_primitive_variables
#ifdef _DEBUG
        call write_debuginfo("In model_initialize(), cellsystem.")
#endif
        call a_model%initialize(a_config)
    end subroutine model_initialize

    ! ### Control volume profiler ###
    subroutine control_volume_profiler_initialize(self, plotter, config, num_conservative_variables, num_primitive_variables)
        class(cellsystem             ), intent(inout) :: self
        class(control_volume_profiler), intent(inout) :: plotter
        class(configuration          ), intent(inout) :: config
        integer(int_kind             ), intent(in   ) :: num_conservative_variables
        integer(int_kind             ), intent(in   ) :: num_primitive_variables
#ifdef _DEBUG
        call write_debuginfo("In control_volume_profiler_initialize(), cellsystem.")
#endif
        call plotter%initialize(config, self%cell_centor_positions, self%is_real_cell, self%num_cells)
    end subroutine control_volume_profiler_initialize

    subroutine control_volume_profiler_write(self, plotter, values_set)
        class(cellsystem             ), intent(inout) :: self
        class(control_volume_profiler), intent(inout) :: plotter
        real (real_kind              ), intent(in   ) :: values_set(:,:)
#ifdef _DEBUG
        call write_debuginfo("In control_volume_profiler_write(), cellsystem.")
#endif
        call plotter%write(self%time, values_set, self%cell_centor_positions, self%cell_volumes)
    end subroutine control_volume_profiler_write

    pure function control_volume_profiler_is_writable(self, plotter) result(juge)
        class(cellsystem             ), intent(in) :: self
        class(control_volume_profiler), intent(in) :: plotter
        logical                                    :: juge
        juge = plotter%is_writable(self%time)
    end function control_volume_profiler_is_writable

    ! ### Line plotter ###
    subroutine line_plotter_initialize(self, plotter, config, num_conservative_variables, num_primitive_variables)
        class(cellsystem   ), intent(inout) :: self
        class(line_plotter ), intent(inout) :: plotter
        class(configuration), intent(inout) :: config
        integer(int_kind   ), intent(in   ) :: num_conservative_variables
        integer(int_kind   ), intent(in   ) :: num_primitive_variables
#ifdef _DEBUG
        call write_debuginfo("In line_plotter_initialize(), cellsystem.")
#endif
        call plotter%initialize(config, self%cell_centor_positions, self%is_real_cell, self%num_cells)
    end subroutine line_plotter_initialize

    subroutine line_plotter_write(self, plotter, values_set)
        class(cellsystem  ), intent(inout) :: self
        class(line_plotter), intent(inout) :: plotter
        real (real_kind   ), intent(in   ) :: values_set(:,:)
#ifdef _DEBUG
        call write_debuginfo("In line_plotter_write(), cellsystem.")
#endif
        call plotter%write(self%time, values_set, self%cell_centor_positions)
    end subroutine line_plotter_write

    pure function line_plotter_is_writable(self, plotter) result(juge)
        class(cellsystem  ), intent(in) :: self
        class(line_plotter), intent(in) :: plotter
        logical                         :: juge
        juge = plotter%is_writable(self%time)
    end function line_plotter_is_writable

    ! ### Time increment control ###
    subroutine time_increment_controller_initialize(self, controller, config, num_conservative_variables, num_primitive_variables)
        class(cellsystem               ), intent(inout) :: self
        class(time_increment_controller), intent(inout) :: controller
        class(configuration            ), intent(inout) :: config
        integer(int_kind               ), intent(in   ) :: num_conservative_variables
        integer(int_kind               ), intent(in   ) :: num_primitive_variables
#ifdef _DEBUG
        call write_debuginfo("In time_increment_controller_initialize(), cellsystem.")
#endif
        call controller%initialize(config)
    end subroutine time_increment_controller_initialize

    subroutine update_time_increment_functional_api(self, controller, an_eos, variables_set, spectral_radius_function)
        class  (cellsystem               ), intent(inout) :: self
        class  (time_increment_controller), intent(in   ) :: controller
        class  (eos                      ), intent(in   ) :: an_eos
        real   (real_kind                ), intent(in   ) :: variables_set(:,:)

        interface
            pure function spectral_radius_function(an_eos, variables, length) result(r)
                use abstract_eos
                use typedef_module
                class  (eos      ), intent(in) :: an_eos
                real   (real_kind), intent(in) :: variables(:)
                real   (real_kind), intent(in) :: length
                real   (real_kind) :: r
            end function spectral_radius_function
        end interface

        integer(int_kind ) :: i

#ifdef _DEBUG
        call write_debuginfo("In update_time_increment_functional_api(), cellsystem.")
#endif

        if ( controller%returns_constant() ) then
            self%time_increment = controller%get_constant_dt()
            return
        end if

        self%time_increment = large_value
        do i = 1, self%num_faces, 1
            associate(                                                                                                             &
                lhc_q        => variables_set    (:, self%face_to_cell_indexes(self%num_local_cells+0, i)), &
                rhc_q        => variables_set    (:, self%face_to_cell_indexes(self%num_local_cells+1, i)), &
                lhc_v        => self%cell_volumes   (self%face_to_cell_indexes(self%num_local_cells+0, i)), &
                rhc_v        => self%cell_volumes   (self%face_to_cell_indexes(self%num_local_cells+1, i)), &
                s            => self%face_areas                                                       (i) , &
                lhc_is_real  => self%is_real_cell   (self%face_to_cell_indexes(self%num_local_cells+0, i)), &
                rhc_is_real  => self%is_real_cell   (self%face_to_cell_indexes(self%num_local_cells+1, i))  &
            )
                if(lhc_is_real)then
                    self%time_increment = min(controller%compute_local_dt(lhc_v, s, spectral_radius_function(an_eos, lhc_q, lhc_v / s)), self%time_increment)
                end if
                if(rhc_is_real)then
                    self%time_increment = min(controller%compute_local_dt(rhc_v, s, spectral_radius_function(an_eos, rhc_q, rhc_v / s)), self%time_increment)
                end if
            end associate
        end do
    end subroutine update_time_increment_functional_api

    subroutine update_time_increment_objective_api(self, controller, an_eos, variables_set, a_model)
        class  (cellsystem               ), intent(inout) :: self
        class  (time_increment_controller), intent(in   ) :: controller
        class  (eos                      ), intent(in   ) :: an_eos
        real   (real_kind                ), intent(in   ) :: variables_set(:,:)
        class  (model                    ), intent(in   ) :: a_model

        integer(int_kind ) :: i

#ifdef _DEBUG
        call write_debuginfo("In update_time_increment_objective_api(), cellsystem.")
#endif

        if ( controller%returns_constant() ) then
            self%time_increment = controller%get_constant_dt()
            return
        end if

        self%time_increment = large_value
        do i = 1, self%num_faces, 1
            associate(                                                                                                             &
                lhc_q        => variables_set    (:, self%face_to_cell_indexes(self%num_local_cells+0, i)), &
                rhc_q        => variables_set    (:, self%face_to_cell_indexes(self%num_local_cells+1, i)), &
                lhc_v        => self%cell_volumes   (self%face_to_cell_indexes(self%num_local_cells+0, i)), &
                rhc_v        => self%cell_volumes   (self%face_to_cell_indexes(self%num_local_cells+1, i)), &
                s            => self%face_areas                                                       (i) , &
                lhc_is_real  => self%is_real_cell   (self%face_to_cell_indexes(self%num_local_cells+0, i)), &
                rhc_is_real  => self%is_real_cell   (self%face_to_cell_indexes(self%num_local_cells+1, i))  &
            )
                if(lhc_is_real)then
                    self%time_increment = min(controller%compute_local_dt(lhc_v, s, a_model%spectral_radius(an_eos, lhc_q, lhc_v / s)), self%time_increment)
                end if
                if(rhc_is_real)then
                    self%time_increment = min(controller%compute_local_dt(rhc_v, s, a_model%spectral_radius(an_eos, rhc_q, rhc_v / s)), self%time_increment)
                end if
            end associate
        end do
    end subroutine update_time_increment_objective_api

    ! ### Termination criterion ###
    subroutine termination_criterion_initialize(self, criterion, config, num_conservative_variables, num_primitive_variables)
        class  (cellsystem           ), intent(inout) :: self
        class  (termination_criterion), intent(inout) :: criterion
        class  (configuration        ), intent(inout) :: config
        integer(int_kind             ), intent(in   ) :: num_conservative_variables
        integer(int_kind             ), intent(in   ) :: num_primitive_variables
#ifdef _DEBUG
        call write_debuginfo("In termination_criterion_initialize(), cellsystem.")
#endif
        call criterion%initialize(config)
    end subroutine termination_criterion_initialize

    pure function satisfies_termination_criterion(self, criterion) result(judge)
        class(cellsystem           ), intent(in) :: self
        class(termination_criterion), intent(in) :: criterion
        logical :: judge
        judge = criterion%is_satisfied(self%time, self%num_steps)
    end function satisfies_termination_criterion

    ! ### Boundary Condition ###
    subroutine apply_outflow_condition(self, primitive_variables_set, num_primitive_variables, &
        compute_rotate_primitive_variables_function, compute_unrotate_primitive_variables_function, boundary_condition_function)
        class  (cellsystem), intent(inout) :: self
        real   (real_kind ), intent(inout) :: primitive_variables_set(:,:)
        integer(int_kind  ), intent(in   ) :: num_primitive_variables

        interface
            pure function compute_rotate_primitive_variables_function( &
                primitive_variables                                  , &
                face_normal_vector                                   , &
                face_tangential1_vector                              , &
                face_tangential2_vector                              , &
                num_primitive_variables                                  ) result(rotate_primitive_variables)

                use typedef_module
                real   (real_kind     ), intent(in)  :: primitive_variables     (:)
                real   (real_kind     ), intent(in)  :: face_normal_vector      (3)
                real   (real_kind     ), intent(in)  :: face_tangential1_vector (3)
                real   (real_kind     ), intent(in)  :: face_tangential2_vector (3)
                integer(int_kind      ), intent(in)  :: num_primitive_variables
                real   (real_kind     )              :: rotate_primitive_variables(num_primitive_variables)
            end function compute_rotate_primitive_variables_function

            pure function compute_unrotate_primitive_variables_function( &
                primitive_variables                                    , &
                face_normal_vector                                     , &
                face_tangential1_vector                                , &
                face_tangential2_vector                                , &
                num_primitive_variables                                    ) result(unrotate_primitive_variables)

                use typedef_module
                real   (real_kind     ), intent(in)  :: primitive_variables     (:)
                real   (real_kind     ), intent(in)  :: face_normal_vector      (3)
                real   (real_kind     ), intent(in)  :: face_tangential1_vector (3)
                real   (real_kind     ), intent(in)  :: face_tangential2_vector (3)
                integer(int_kind      ), intent(in)  :: num_primitive_variables
                real   (real_kind     )              :: unrotate_primitive_variables(num_primitive_variables)
            end function compute_unrotate_primitive_variables_function

            pure function boundary_condition_function(inner_primitive_variables, num_primitive_variables) result(ghost_primitive_variables)
                use typedef_module
                real   (real_kind), intent(in) :: inner_primitive_variables(:)
                integer(int_kind ), intent(in) :: num_primitive_variables
                real   (real_kind)             :: ghost_primitive_variables(num_primitive_variables)
            end function boundary_condition_function
        end interface

        integer :: i, face_index

#ifdef _DEBUG
        call write_debuginfo("In apply_outflow_condition(), cellsystem.")
#endif

        !$omp parallel do private(face_index)
        do i = 1, self%num_outflow_faces, 1
            face_index = self%outflow_face_indexes(i)
            call apply_boundary_condition_common_impl(       &
                primitive_variables_set                    , &
                self%face_normal_vectors     (:,face_index), &
                self%face_tangential1_vectors(:,face_index), &
                self%face_tangential2_vectors(:,face_index), &
                self%face_to_cell_indexes    (:,face_index), &
                self%num_local_cells                       , &
                num_primitive_variables                    , &
                compute_rotate_primitive_variables_function   , &
                compute_unrotate_primitive_variables_function , &
                boundary_condition_function                &
            )
        end do
    end subroutine apply_outflow_condition

    subroutine apply_wall_condition(self, primitive_variables_set, num_primitive_variables, &
        compute_rotate_primitive_variables_function, compute_unrotate_primitive_variables_function, boundary_condition_function)
        class  (cellsystem), intent(inout) :: self
        real   (real_kind ), intent(inout) :: primitive_variables_set(:,:)
        integer(int_kind  ), intent(in   ) :: num_primitive_variables

        interface
            pure function compute_rotate_primitive_variables_function( &
                primitive_variables                                  , &
                face_normal_vector                                   , &
                face_tangential1_vector                              , &
                face_tangential2_vector                              , &
                num_primitive_variables                                  ) result(rotate_primitive_variables)

                use typedef_module
                real   (real_kind     ), intent(in)  :: primitive_variables     (:)
                real   (real_kind     ), intent(in)  :: face_normal_vector      (3)
                real   (real_kind     ), intent(in)  :: face_tangential1_vector (3)
                real   (real_kind     ), intent(in)  :: face_tangential2_vector (3)
                integer(int_kind      ), intent(in)  :: num_primitive_variables
                real   (real_kind     )              :: rotate_primitive_variables(num_primitive_variables)
            end function compute_rotate_primitive_variables_function

            pure function compute_unrotate_primitive_variables_function( &
                primitive_variables                                    , &
                face_normal_vector                                     , &
                face_tangential1_vector                                , &
                face_tangential2_vector                                , &
                num_primitive_variables                                    ) result(unrotate_primitive_variables)

                use typedef_module
                real   (real_kind     ), intent(in)  :: primitive_variables     (:)
                real   (real_kind     ), intent(in)  :: face_normal_vector      (3)
                real   (real_kind     ), intent(in)  :: face_tangential1_vector (3)
                real   (real_kind     ), intent(in)  :: face_tangential2_vector (3)
                integer(int_kind      ), intent(in)  :: num_primitive_variables
                real   (real_kind     )              :: unrotate_primitive_variables(num_primitive_variables)
            end function compute_unrotate_primitive_variables_function

            pure function boundary_condition_function(inner_primitive_variables, num_primitive_variables) result(ghost_primitive_variables)
                use typedef_module
                real   (real_kind), intent(in) :: inner_primitive_variables(:)
                integer(int_kind ), intent(in) :: num_primitive_variables
                real   (real_kind)             :: ghost_primitive_variables(num_primitive_variables)
            end function boundary_condition_function
        end interface

        integer :: i, face_index

#ifdef _DEBUG
        call write_debuginfo("In apply_wall_condition(), cellsystem.")
#endif

        !$omp parallel do private(face_index)
        do i = 1, self%num_wall_faces, 1
            face_index = self%wall_face_indexes(i)
            call apply_boundary_condition_common_impl(         &
                primitive_variables_set                      , &
                self%face_normal_vectors     (:,face_index)  , &
                self%face_tangential1_vectors(:,face_index)  , &
                self%face_tangential2_vectors(:,face_index)  , &
                self%face_to_cell_indexes    (:,face_index)  , &
                self%num_local_cells                         , &
                num_primitive_variables                      , &
                compute_rotate_primitive_variables_function  , &
                compute_unrotate_primitive_variables_function, &
                boundary_condition_function                    &
            )
        end do
    end subroutine apply_wall_condition

    subroutine apply_symmetric_condition(self, primitive_variables_set, num_primitive_variables, &
        compute_rotate_primitive_variables_function, compute_unrotate_primitive_variables_function, boundary_condition_function)
        class  (cellsystem), intent(inout) :: self
        real   (real_kind ), intent(inout) :: primitive_variables_set(:,:)
        integer(int_kind  ), intent(in   ) :: num_primitive_variables

        interface
            pure function compute_rotate_primitive_variables_function( &
                primitive_variables                                  , &
                face_normal_vector                                   , &
                face_tangential1_vector                              , &
                face_tangential2_vector                              , &
                num_primitive_variables                                  ) result(rotate_primitive_variables)

                use typedef_module
                real   (real_kind     ), intent(in)  :: primitive_variables     (:)
                real   (real_kind     ), intent(in)  :: face_normal_vector      (3)
                real   (real_kind     ), intent(in)  :: face_tangential1_vector (3)
                real   (real_kind     ), intent(in)  :: face_tangential2_vector (3)
                integer(int_kind      ), intent(in)  :: num_primitive_variables
                real   (real_kind     )              :: rotate_primitive_variables(num_primitive_variables)
            end function compute_rotate_primitive_variables_function

            pure function compute_unrotate_primitive_variables_function( &
                primitive_variables                                    , &
                face_normal_vector                                     , &
                face_tangential1_vector                                , &
                face_tangential2_vector                                , &
                num_primitive_variables                                    ) result(unrotate_primitive_variables)

                use typedef_module
                real   (real_kind     ), intent(in)  :: primitive_variables     (:)
                real   (real_kind     ), intent(in)  :: face_normal_vector      (3)
                real   (real_kind     ), intent(in)  :: face_tangential1_vector (3)
                real   (real_kind     ), intent(in)  :: face_tangential2_vector (3)
                integer(int_kind      ), intent(in)  :: num_primitive_variables
                real   (real_kind     )              :: unrotate_primitive_variables(num_primitive_variables)
            end function compute_unrotate_primitive_variables_function

            pure function boundary_condition_function(inner_primitive_variables, num_primitive_variables) result(ghost_primitive_variables)
                use typedef_module
                real   (real_kind), intent(in) :: inner_primitive_variables(:)
                integer(int_kind ), intent(in) :: num_primitive_variables
                real   (real_kind)             :: ghost_primitive_variables(num_primitive_variables)
            end function boundary_condition_function
        end interface

        integer :: i, face_index

#ifdef _DEBUG
        call write_debuginfo("In apply_symmetric_condition(), cellsystem.")
#endif

        !$omp parallel do private(face_index)
        do i = 1, self%num_symmetric_faces, 1
            face_index = self%symmetric_face_indexes(i)
            call apply_boundary_condition_common_impl(          &
                primitive_variables_set                       , &
                self%face_normal_vectors     (:,face_index)   , &
                self%face_tangential1_vectors(:,face_index)   , &
                self%face_tangential2_vectors(:,face_index)   , &
                self%face_to_cell_indexes    (:,face_index)   , &
                self%num_local_cells                          , &
                num_primitive_variables                       , &
                compute_rotate_primitive_variables_function   , &
                compute_unrotate_primitive_variables_function , &
                boundary_condition_function                     &
            )
        end do
    end subroutine apply_symmetric_condition

    subroutine apply_empty_condition(self, primitive_variables_set, num_primitive_variables, &
        compute_rotate_primitive_variables_function, compute_unrotate_primitive_variables_function, boundary_condition_function)
        class  (cellsystem), intent(inout) :: self
        real   (real_kind ), intent(inout) :: primitive_variables_set(:,:)
        integer(int_kind  ), intent(in   ) :: num_primitive_variables

        interface
            pure function compute_rotate_primitive_variables_function( &
                primitive_variables                                  , &
                face_normal_vector                                   , &
                face_tangential1_vector                              , &
                face_tangential2_vector                              , &
                num_primitive_variables                                  ) result(rotate_primitive_variables)

                use typedef_module
                real   (real_kind     ), intent(in)  :: primitive_variables     (:)
                real   (real_kind     ), intent(in)  :: face_normal_vector      (3)
                real   (real_kind     ), intent(in)  :: face_tangential1_vector (3)
                real   (real_kind     ), intent(in)  :: face_tangential2_vector (3)
                integer(int_kind      ), intent(in)  :: num_primitive_variables
                real   (real_kind     )              :: rotate_primitive_variables(num_primitive_variables)
            end function compute_rotate_primitive_variables_function

            pure function compute_unrotate_primitive_variables_function( &
                primitive_variables                                    , &
                face_normal_vector                                     , &
                face_tangential1_vector                                , &
                face_tangential2_vector                                , &
                num_primitive_variables                                    ) result(unrotate_primitive_variables)

                use typedef_module
                real   (real_kind     ), intent(in)  :: primitive_variables     (:)
                real   (real_kind     ), intent(in)  :: face_normal_vector      (3)
                real   (real_kind     ), intent(in)  :: face_tangential1_vector (3)
                real   (real_kind     ), intent(in)  :: face_tangential2_vector (3)
                integer(int_kind      ), intent(in)  :: num_primitive_variables
                real   (real_kind     )              :: unrotate_primitive_variables(num_primitive_variables)
            end function compute_unrotate_primitive_variables_function

            pure function boundary_condition_function(inner_primitive_variables, num_primitive_variables) result(ghost_primitive_variables)
                use typedef_module
                real   (real_kind), intent(in) :: inner_primitive_variables(:)
                integer(int_kind ), intent(in) :: num_primitive_variables
                real   (real_kind)             :: ghost_primitive_variables(num_primitive_variables)
            end function boundary_condition_function
        end interface

        integer :: i, face_index

#ifdef _DEBUG
        call write_debuginfo("In apply_empty_condition(), cellsystem.")
#endif

        !$omp parallel do private(face_index)
        do i = 1, self%num_empty_faces, 1
            face_index = self%empty_face_indexes(i)
            call apply_boundary_condition_empty_impl(           &
                primitive_variables_set                       , &
                self%face_normal_vectors     (:,face_index)   , &
                self%face_tangential1_vectors(:,face_index)   , &
                self%face_tangential2_vectors(:,face_index)   , &
                self%face_to_cell_indexes    (:,face_index)   , &
                self%num_local_cells                          , &
                num_primitive_variables                       , &
                compute_rotate_primitive_variables_function   , &
                compute_unrotate_primitive_variables_function , &
                boundary_condition_function                     &
            )
        end do
    end subroutine apply_empty_condition

    ! ### Variables ###
    subroutine variables_initialize(self, variables_set, num_variables)
        class  (cellsystem), intent(inout)              :: self
        real   (real_kind ), intent(inout), allocatable :: variables_set(:,:)
        integer(int_kind  ), intent(in   )              :: num_variables

#ifdef _DEBUG
        call write_debuginfo("In variables_initialize(), cellsystem.")
#endif

        if(.not. self%read_cellsystem) call call_error("'read' subroutine is not called. You should call with following steps: first you call 'read' subroutine, next you initialze variables with 'initialze' subroutine. Please check your cord.")

        if(allocated(variables_set))then
            call call_error("Array variables_set is allocated. But you call 'initialize' subroutine.")
        end if
        allocate(variables_set(1:num_variables, 1:self%num_cells))
        variables_set(:,:) = 0.d0
    end subroutine variables_initialize

    subroutine variables_1darray_initialize(self, variables_set)
        class  (cellsystem), intent(inout)              :: self
        real   (real_kind ), intent(inout), allocatable :: variables_set(:)

#ifdef _DEBUG
        call write_debuginfo("In variables_1darray_initialize(), cellsystem.")
#endif

        if(.not. self%read_cellsystem) call call_error("'read' subroutine is not called. You should call with following steps: first you call 'read' subroutine, next you initialze variables with 'initialze' subroutine. Please check your cord.")

        if(allocated(variables_set))then
            call call_error("Array variables_set is allocated. But you call 'initialize' subroutine.")
        end if
        allocate(variables_set(1:self%num_cells))
        variables_set(:) = 0.d0
    end subroutine variables_1darray_initialize

    subroutine read_initial_condition(self, an_initial_condition_parser, config, conservative_variables_set)
        class  (cellsystem              ), intent(inout) :: self
        class  (initial_condition_parser), intent(inout) :: an_initial_condition_parser
        class  (configuration           ), intent(inout) :: config
        real   (real_kind               ), intent(inout) :: conservative_variables_set(:,:)

        character(len=:), allocatable :: filepath
        logical          :: found
        integer          :: i

#ifdef _DEBUG
        call write_debuginfo("In read_initial_condition(), cellsystem.")
#endif

        ! TODO: Following lines move to {@code initial_condition_parser} class.
        call config%get_char("Initial condition.Filepath", filepath, found)
        if(.not. found) call call_error("'Initial condition.Filepath' is not found in configuration file you set.")

        call an_initial_condition_parser%parse(filepath)
        call an_initial_condition_parser%get_conservative_variables_set(conservative_variables_set)
        call an_initial_condition_parser%close()
    end subroutine read_initial_condition

    subroutine conservative_to_primitive_variables_all(self, an_eos, conservative_variables_set, primitive_variables_set, num_primitive_variables, conservative_to_primitive_function)
        class  (cellsystem              ), intent(inout) :: self
        class  (eos                     ), intent(in)    :: an_eos
        real   (real_kind               ), intent(in)    :: conservative_variables_set(:,:)
        real   (real_kind               ), intent(inout) :: primitive_variables_set(:,:)
        integer(int_kind                ), intent(in)    :: num_primitive_variables
        interface
            pure function conservative_to_primitive_function(an_eos, conservative_variables, num_primitive_variables) result(primitive_variables)
                use typedef_module
                use abstract_eos
                class  (eos      ), intent(in)  :: an_eos
                real   (real_kind), intent(in)  :: conservative_variables(:)
                integer(int_kind ), intent(in)  :: num_primitive_variables
                real   (real_kind)              :: primitive_variables(num_primitive_variables)
            end function conservative_to_primitive_function
        end interface

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In conservative_to_primitive_variables_all(), cellsystem.")
#endif

        !$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            primitive_variables_set(:,i) = conservative_to_primitive_function(an_eos, conservative_variables_set(:,i), num_primitive_variables)
        end do
    end subroutine conservative_to_primitive_variables_all

    subroutine processes_variables_set_single_array(self, variables_set, num_variables, processing_function)
        class  (cellsystem), intent(inout) :: self
        real   (real_kind ), intent(inout) :: variables_set(:,:)
        integer(int_kind  ), intent(in   ) :: num_variables
        interface
            pure function processing_function(variables, num_variables) result(destination_variables)
                use typedef_module
                real   (real_kind ), intent(in) :: variables(:)
                integer(int_kind  ), intent(in) :: num_variables
                real   (real_kind )             :: destination_variables(num_variables)
            end function processing_function
        end interface

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In processes_variables_set_single_array(), cellsystem.")
#endif

        !$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            variables_set(:,i) = processing_function(variables_set(:,i), num_variables)
        end do
    end subroutine processes_variables_set_single_array

    subroutine processes_variables_set_two_array(self, primary_variables_set, secondary_variables_set, num_variables, processing_function)
        class  (cellsystem), intent(inout) :: self
        real   (real_kind ), intent(inout) :: primary_variables_set(:,:), secondary_variables_set(:,:)
        integer(int_kind  ), intent(in   ) :: num_variables
        interface
            pure function processing_function(primary_variables, secondary_variables, num_variables) result(destination_variables)
                use typedef_module
                real   (real_kind ), intent(in) :: primary_variables(:), secondary_variables(:)
                integer(int_kind  ), intent(in) :: num_variables
                real   (real_kind )             :: destination_variables(num_variables)
            end function processing_function
        end interface

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In processes_variables_set_two_array(), cellsystem.")
#endif

        !$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            primary_variables_set(:,i) = processing_function(primary_variables_set(:,i), secondary_variables_set(:,i), num_variables)
        end do
    end subroutine processes_variables_set_two_array

    ! ### Gradient Calculator ###
    subroutine gradient_calculator_initialize(self, a_gradient_calculator, config, num_conservative_variables, num_primitive_variables)
        class  (cellsystem         ), intent(inout) :: self
        class  (gradient_calculator), intent(inout) :: a_gradient_calculator
        class  (configuration      ), intent(inout) :: config
        integer(int_kind           ), intent(in   ) :: num_conservative_variables
        integer(int_kind           ), intent(in   ) :: num_primitive_variables
#ifdef _DEBUG
        call write_debuginfo("In gradient_calculator_initialize(), cellsystem.")
#endif
        call a_gradient_calculator%initialize(config)
    end subroutine gradient_calculator_initialize

    subroutine compute_gradient_2darray(self, a_gradient_calculator, variables_set, gradient_variables_set, num_variables)
        class  (cellsystem         ), intent(inout) :: self
        class  (gradient_calculator), intent(inout) :: a_gradient_calculator
        real   (real_kind          ), intent(in   ) :: variables_set            (:,:)
        real   (real_kind          ), intent(inout) :: gradient_variables_set   (:,:)
        integer(int_kind           ), intent(in   ) :: num_variables

        integer(int_kind ) :: face_index, cell_index, rhc_index, lhc_index, var_index

#ifdef _DEBUG
        call write_debuginfo("In compute_gradient_2darray(), cellsystem.")
#endif

!$omp parallel do private(cell_index)
        do cell_index = 1, self%num_cells, 1
            gradient_variables_set(:,cell_index) = 0.d0
        end do

!$omp parallel do private(face_index, rhc_index, lhc_index, var_index)
        do face_index = 1, self%num_faces, 1
            lhc_index = self%face_to_cell_indexes(self%num_local_cells - 0, face_index)
            rhc_index = self%face_to_cell_indexes(self%num_local_cells + 1, face_index)
            do var_index = 1, num_variables, 1
                associate(                                &
                    vec_start_index => 3*(var_index-1)+1, &
                    vec_end_index   => 3*(var_index-1)+3, &
                    residual        => a_gradient_calculator%compute_residual(                         &
                                           variables_set                      (var_index, lhc_index ), &
                                           variables_set                      (var_index, rhc_index ), &
                                           self%cell_centor_positions         (1:3      , lhc_index ), &
                                           self%cell_centor_positions         (1:3      , rhc_index ), &
                                           self%face_normal_vectors           (1:3      , face_index), &
                                           self%face_positions                (1:3      , face_index), &
                                           self%face_areas                               (face_index)  &
                                        )                                                              &
                )
                    gradient_variables_set(vec_start_index:vec_end_index, rhc_index) = gradient_variables_set(vec_start_index:vec_end_index, rhc_index)  &
                                                        - (1.d0 / self%cell_volumes(rhc_index)) * residual
                    gradient_variables_set(vec_start_index:vec_end_index, lhc_index) = gradient_variables_set(vec_start_index:vec_end_index, lhc_index)  &
                                                        + (1.d0 / self%cell_volumes(lhc_index)) * residual
                end associate
            end do
        end do
    end subroutine compute_gradient_2darray

    subroutine compute_gradient_1darray(self, a_gradient_calculator, variable_set, gradient_variable_set)
        class(cellsystem         ), intent(inout) :: self
        class(gradient_calculator), intent(inout) :: a_gradient_calculator
        real (real_kind          ), intent(in   ) :: variable_set            (:)
        real (real_kind          ), intent(inout) :: gradient_variable_set   (:,:)

        integer(int_kind) :: face_index, cell_index, rhc_index, lhc_index

#ifdef _DEBUG
        call write_debuginfo("In compute_gradient_1darray(), cellsystem.")
#endif

!$omp parallel do private(cell_index)
        do cell_index = 1, self%num_cells, 1
            gradient_variable_set(:,cell_index) = 0.d0
        end do

!$omp parallel do private(face_index, rhc_index, lhc_index)
        do face_index = 1, self%num_faces, 1
            lhc_index = self%face_to_cell_indexes(self%num_local_cells - 0, face_index)
            rhc_index = self%face_to_cell_indexes(self%num_local_cells + 1, face_index)

            associate(                                                                     &
                residual => a_gradient_calculator%compute_residual(                        &
                               variable_set                                  (lhc_index ), &
                               variable_set                                  (rhc_index ), &
                               self%cell_centor_positions         (1:3      , lhc_index ), &
                               self%cell_centor_positions         (1:3      , rhc_index ), &
                               self%face_normal_vectors           (1:3      , face_index), &
                               self%face_positions                (1:3      , face_index), &
                               self%face_areas                               (face_index)  &
                            )                                                              &
            )
                gradient_variable_set(1:3,rhc_index) = gradient_variable_set(1:3,rhc_index)                    &
                                                 - (1.d0 / self%cell_volumes(rhc_index)) * residual
                gradient_variable_set(1:3,lhc_index) = gradient_variable_set(1:3,lhc_index)                    &
                                                 + (1.d0 / self%cell_volumes(lhc_index))  * residual
            end associate
        end do
    end subroutine compute_gradient_1darray

    ! ### Divergence Calculator ###
    subroutine divergence_calculator_initialize(self, a_divergence_calculator, config, num_conservative_variables, num_primitive_variables)
        class  (cellsystem           ), intent(inout) :: self
        class  (divergence_calculator), intent(inout) :: a_divergence_calculator
        class  (configuration        ), intent(inout) :: config
        integer(int_kind             ), intent(in   ) :: num_conservative_variables
        integer(int_kind             ), intent(in   ) :: num_primitive_variables
#ifdef _DEBUG
        call write_debuginfo("In divergence_calculator_initialize(), cellsystem.")
#endif
        call a_divergence_calculator%initialize(config)
    end subroutine divergence_calculator_initialize

    subroutine compute_divergence_2darray(self, a_divergence_calculator, variables_set, divergence_variables_set, num_variables)
        class  (cellsystem           ), intent(inout) :: self
        class  (divergence_calculator), intent(inout) :: a_divergence_calculator
        real   (real_kind            ), intent(in   ) :: variables_set           (:,:)
        real   (real_kind            ), intent(inout) :: divergence_variables_set(:,:)
        integer(int_kind             ), intent(in   ) :: num_variables

        integer(int_kind) :: face_index, cell_index, rhc_index, lhc_index, var_index, num_divergence_variables
        real  (real_kind) :: residual

#ifdef _DEBUG
        call write_debuginfo("In compute_divergence_2darray(), cellsystem.")
#endif

!$omp parallel do private(cell_index)
        do cell_index = 1, self%num_cells, 1
            divergence_variables_set(:,cell_index) = 0.d0
        end do

        num_divergence_variables = num_variables/3

!$omp parallel do private(face_index, rhc_index, lhc_index, var_index)
        do face_index = 1, self%num_faces, 1
            lhc_index = self%face_to_cell_indexes(self%num_local_cells - 0, face_index)
            rhc_index = self%face_to_cell_indexes(self%num_local_cells + 1, face_index)
            do var_index = 1, num_divergence_variables, 1
                associate(                                &
                    vec_start_index => 3*(var_index-1)+1, &
                    vec_end_index   => 3*(var_index-1)+3  &
                )
                    residual = a_divergence_calculator%compute_residual(                               &
                       variables_set                      (vec_start_index:vec_end_index, lhc_index ), &
                       variables_set                      (vec_start_index:vec_end_index, rhc_index ), &
                       self%cell_centor_positions         (1:3                          , lhc_index ), &
                       self%cell_centor_positions         (1:3                          , rhc_index ), &
                       self%face_normal_vectors           (1:3                          , face_index), &
                       self%face_positions                (1:3                          , face_index), &
                       self%face_areas                                                   (face_index)  &
                    )
                    divergence_variables_set(var_index, rhc_index) = divergence_variables_set(var_index, rhc_index) &
                                                        - (1.d0 / self%cell_volumes(rhc_index)) * residual
                    divergence_variables_set(var_index, lhc_index) = divergence_variables_set(var_index, lhc_index) &
                                                        + (1.d0 / self%cell_volumes(lhc_index)) * residual
                end associate
            end do
        end do
    end subroutine compute_divergence_2darray

    subroutine compute_divergence_1darray(self, a_divergence_calculator, variable_set, divergence_variable_set)
        class(cellsystem           ), intent(inout) :: self
        class(divergence_calculator), intent(inout) :: a_divergence_calculator
        real (real_kind            ), intent(in   ) :: variable_set              (:,:)
        real (real_kind            ), intent(inout) :: divergence_variable_set   (:)

        integer(int_kind) :: face_index, cell_index, rhc_index, lhc_index

#ifdef _DEBUG
        call write_debuginfo("In compute_divergence_1darray(), cellsystem.")
#endif

!$omp parallel do private(cell_index)
        do cell_index = 1, self%num_cells, 1
            divergence_variable_set(cell_index) = 0.d0
        end do

!$omp parallel do private(face_index, rhc_index, lhc_index)
        do face_index = 1, self%num_faces, 1
            lhc_index = self%face_to_cell_indexes(self%num_local_cells - 0, face_index)
            rhc_index = self%face_to_cell_indexes(self%num_local_cells + 1, face_index)
            associate(                                                   &
                residual => a_divergence_calculator%compute_residual(    &
                    variable_set                       (1:3,lhc_index ), &
                    variable_set                       (1:3,rhc_index ), &
                    self%cell_centor_positions         (1:3,lhc_index ), &
                    self%cell_centor_positions         (1:3,rhc_index ), &
                    self%face_normal_vectors           (1:3,face_index), &
                    self%face_positions                (1:3,face_index), &
                    self%face_areas                        (face_index)  &
                )                                                        &
            )
                divergence_variable_set(rhc_index)   = divergence_variable_set(rhc_index)               &
                                                     - (1.d0 / self%cell_volumes(rhc_index)) * residual
                divergence_variable_set(lhc_index)   = divergence_variable_set(lhc_index)               &
                                                     + (1.d0 / self%cell_volumes(lhc_index)) * residual
            end associate
        end do
    end subroutine compute_divergence_1darray

    ! ### EoS ###
    subroutine eos_initialize(self, an_eos, config, num_conservative_variables, num_primitive_variables)
        class  (cellsystem    ), intent(inout) :: self
        class  (eos           ), intent(inout) :: an_eos
        class  (configuration ), intent(inout) :: config
        integer(int_kind      ), intent(in   ) :: num_conservative_variables
        integer(int_kind      ), intent(in   ) :: num_primitive_variables
#ifdef _DEBUG
        call write_debuginfo("In eos_initialize(), cellsystem.")
#endif
        call an_eos%initialize(config)
    end subroutine eos_initialize

    ! ### Resudual ###
    subroutine compute_residual_functional_api(self, a_reconstructor, a_riemann_solver, an_eos,                                            &
                                           primitive_variables_set, residual_set, num_conservative_variables, num_primitive_variables, &
                                           primitive_to_conservative_function, residual_element_function                                )

        class  (cellsystem    ), intent(in   ) :: self
        class  (reconstructor ), intent(in   ) :: a_reconstructor
        class  (riemann_solver), intent(in   ) :: a_riemann_solver
        class  (eos           ), intent(in   ) :: an_eos
        real   (real_kind     ), intent(in   ) :: primitive_variables_set (:,:)
        real   (real_kind     ), intent(inout) :: residual_set            (:,:)
        integer(int_kind      ), intent(in   ) :: num_conservative_variables
        integer(int_kind      ), intent(in   ) :: num_primitive_variables

        interface
            pure function primitive_to_conservative_function(an_eos, primitive_variables, num_conservative_values) result(conservative_values)
                use typedef_module
                use abstract_eos
                class  (eos      ), intent(in) :: an_eos
                real   (real_kind), intent(in) :: primitive_variables(:)
                integer(int_kind ), intent(in) :: num_conservative_values
                real   (real_kind)             :: conservative_values(num_conservative_values)
            end function primitive_to_conservative_function

            pure function residual_element_function(      &
                an_eos                                  , &
                an_riemann_solver                       , &
                primitive_variables_lhc                 , &
                primitive_variables_rhc                 , &
                reconstructed_primitive_variables_lhc   , &
                reconstructed_primitive_variables_rhc   , &
                lhc_cell_volume                         , &
                rhc_cell_volume                         , &
                face_area                               , &
                face_normal_vector                      , &
                face_tangential1_vector                 , &
                face_tangential2_vector                 , &
                num_conservative_values                 , &
                num_primitive_values                        ) result(residual_element)

                use typedef_module
                use abstract_eos
                use abstract_riemann_solver

                class  (eos           ), intent(in) :: an_eos
                class  (riemann_solver), intent(in) :: an_riemann_solver
                real   (real_kind     ), intent(in) :: primitive_variables_lhc                  (:)
                real   (real_kind     ), intent(in) :: primitive_variables_rhc                  (:)
                real   (real_kind     ), intent(in) :: reconstructed_primitive_variables_lhc    (:)
                real   (real_kind     ), intent(in) :: reconstructed_primitive_variables_rhc    (:)
                real   (real_kind     ), intent(in) :: lhc_cell_volume
                real   (real_kind     ), intent(in) :: rhc_cell_volume
                real   (real_kind     ), intent(in) :: face_area
                real   (real_kind     ), intent(in) :: face_normal_vector     (3)
                real   (real_kind     ), intent(in) :: face_tangential1_vector(3)
                real   (real_kind     ), intent(in) :: face_tangential2_vector(3)
                integer(int_kind      ), intent(in) :: num_conservative_values
                integer(int_kind      ), intent(in) :: num_primitive_values

                real   (real_kind)                  :: residual_element(num_conservative_values, 1:2)
            end function residual_element_function
        end interface

        integer(int_kind ) :: i
        integer(int_kind ) :: lhc_index
        integer(int_kind ) :: rhc_index
        real   (real_kind) :: reconstructed_primitive_variables_lhc   (num_primitive_variables)
        real   (real_kind) :: reconstructed_primitive_variables_rhc   (num_primitive_variables)
        real   (real_kind) :: reconstructed_conservative_variables_lhc(num_conservative_variables)
        real   (real_kind) :: reconstructed_conservative_variables_rhc(num_conservative_variables)
        real   (real_kind) :: residual_element                        (num_conservative_variables, 1:2)

        logical :: lhc_is_dummy, rhc_is_dummy

#ifdef _DEBUG
        call write_debuginfo("In compute_residual_functional_api(), cellsystem.")
#endif

!$omp parallel do private(i,lhc_index,rhc_index,reconstructed_primitive_variables_lhc,reconstructed_primitive_variables_rhc,reconstructed_conservative_variables_lhc,reconstructed_conservative_variables_rhc,residual_element)
        do i = 1, self%num_faces, 1
            lhc_index = self%face_to_cell_indexes(self%num_local_cells+0, i)
            rhc_index = self%face_to_cell_indexes(self%num_local_cells+1, i)

            !lhc_is_dummy = (self%is_real_cell(lhc_index) .eqv. .false.)
            !rhc_is_dummy = (self%is_real_cell(lhc_index) .eqv. .false.)
            !if(lhc_is_dummy .and. rhc_is_dummy)cycle

            reconstructed_primitive_variables_lhc(:) = a_reconstructor%reconstruct_lhc(                              &
                primitive_variables_set, self%face_to_cell_indexes, self%cell_centor_positions, self%face_positions, &
                i, self%num_local_cells, num_primitive_variables                                                     &
            )
            reconstructed_primitive_variables_rhc(:) = a_reconstructor%reconstruct_rhc(                              &
                primitive_variables_set, self%face_to_cell_indexes, self%cell_centor_positions, self%face_positions, &
                i, self%num_local_cells, num_primitive_variables                                                     &
            )

            residual_element(:,:) = residual_element_function(&
                an_eos                                      , &
                a_riemann_solver                            , &
                primitive_variables_set       (:,lhc_index) , &
                primitive_variables_set       (:,rhc_index) , &
                reconstructed_primitive_variables_lhc   (:) , &
                reconstructed_primitive_variables_rhc   (:) , &
                self%cell_volumes            (lhc_index)    , &
                self%cell_volumes            (rhc_index)    , &
                self%face_areas                  (i)        , &
                self%face_normal_vectors     (1:3,i)        , &
                self%face_tangential1_vectors(1:3,i)        , &
                self%face_tangential2_vectors(1:3,i)        , &
                num_conservative_variables                  , &
                num_primitive_variables                       &
            )

            residual_set(:,lhc_index) = residual_set(:,lhc_index) + residual_element(:, 1)
            residual_set(:,rhc_index) = residual_set(:,rhc_index) + residual_element(:, 2)
        end do
    end subroutine compute_residual_functional_api

    subroutine compute_residual_objective_api(self, a_reconstructor, a_riemann_solver, an_eos, a_face_gradient_interpolator,                                                &
                                          primitive_variables_set, gradient_primitive_variables_set, residual_set, num_conservative_variables, num_primitive_variables, &
                                          primitive_to_conservative_function, a_model)

        class  (cellsystem                ), intent(in   ) :: self
        class  (reconstructor             ), intent(in   ) :: a_reconstructor
        class  (riemann_solver            ), intent(in   ) :: a_riemann_solver
        class  (eos                       ), intent(in   ) :: an_eos
        class  (face_gradient_interpolator), intent(in   ) :: a_face_gradient_interpolator
        real   (real_kind                 ), intent(in   ) :: primitive_variables_set          (:,:)
        real   (real_kind                 ), intent(in   ) :: gradient_primitive_variables_set (:,:)
        real   (real_kind                 ), intent(inout) :: residual_set                     (:,:)
        integer(int_kind                  ), intent(in   ) :: num_conservative_variables
        integer(int_kind                  ), intent(in   ) :: num_primitive_variables
        class  (model                     ), intent(in   ) :: a_model

        interface
            pure function primitive_to_conservative_function(an_eos, primitive_variables, num_conservative_values) result(conservative_values)
                use typedef_module
                use abstract_eos
                class  (eos      ), intent(in) :: an_eos
                real   (real_kind), intent(in) :: primitive_variables(:)
                integer(int_kind ), intent(in) :: num_conservative_values
                real   (real_kind)             :: conservative_values(num_conservative_values)
            end function primitive_to_conservative_function
        end interface

        integer(int_kind ) :: i
        real   (real_kind) :: reconstructed_primitive_variables_lhc   (num_primitive_variables)
        real   (real_kind) :: reconstructed_primitive_variables_rhc   (num_primitive_variables)
        real   (real_kind) :: face_gradient_primitive_variables       (num_primitive_variables*3)
        real   (real_kind) :: reconstructed_conservative_variables_lhc(num_conservative_variables)
        real   (real_kind) :: reconstructed_conservative_variables_rhc(num_conservative_variables)
        real   (real_kind) :: residual_element                        (num_conservative_variables, 1:2)

#ifdef _DEBUG
        call write_debuginfo("In compute_residual_objective_api(), cellsystem.")
#endif

!$omp parallel do private(i,reconstructed_primitive_variables_lhc,reconstructed_primitive_variables_rhc,face_gradient_primitive_variables,reconstructed_conservative_variables_lhc,reconstructed_conservative_variables_rhc,residual_element)
        do i = 1, self%num_faces, 1
            associate(                                                             &
                lhc_index => self%face_to_cell_indexes(self%num_local_cells+0, i), &
                rhc_index => self%face_to_cell_indexes(self%num_local_cells+1, i)  &
            )

                reconstructed_primitive_variables_lhc(:) = a_reconstructor%reconstruct_lhc(                              &
                    primitive_variables_set, self%face_to_cell_indexes, self%cell_centor_positions, self%face_positions, &
                    i, self%num_local_cells, num_primitive_variables                                                     &
                )
                reconstructed_primitive_variables_rhc(:) = a_reconstructor%reconstruct_rhc(                              &
                    primitive_variables_set, self%face_to_cell_indexes, self%cell_centor_positions, self%face_positions, &
                    i, self%num_local_cells, num_primitive_variables                                                     &
                )

                if(.not. self%is_real_cell(rhc_index))then
                    face_gradient_primitive_variables(:) = self%compute_boundary_gradient(                &
                        primitive_variables_set   (:,lhc_index), primitive_variables_set   (:,rhc_index), &
                        self%cell_centor_positions(:,lhc_index), self%cell_centor_positions(:,rhc_index), &
                        num_primitive_variables                                                           &
                    )
                else
                    face_gradient_primitive_variables(:) = a_face_gradient_interpolator%interpolate(                  &
                        gradient_primitive_variables_set(:,lhc_index), gradient_primitive_variables_set(:,rhc_index), &
                        primitive_variables_set         (:,lhc_index), primitive_variables_set         (:,rhc_index), &
                        self%cell_centor_positions      (:,lhc_index), self%cell_centor_positions      (:,rhc_index), &
                        num_primitive_variables                                                                       &
                    )
                end if

                residual_element(:,:) = a_model%compute_residual_element( &
                    an_eos                                              , &
                    a_riemann_solver                                    , &
                    primitive_variables_set              (:,lhc_index)  , &
                    primitive_variables_set              (:,rhc_index)  , &
                    reconstructed_primitive_variables_lhc(:)            , &
                    reconstructed_primitive_variables_rhc(:)            , &
                    face_gradient_primitive_variables    (:)            , &
                    self%cell_volumes                      (lhc_index)  , &
                    self%cell_volumes                      (rhc_index)  , &
                    self%face_areas                  (i)                , &
                    self%face_normal_vectors     (1:3,i)                , &
                    self%face_tangential1_vectors(1:3,i)                , &
                    self%face_tangential2_vectors(1:3,i)                , &
                    num_conservative_variables                          , &
                    num_primitive_variables                               &
                )

                residual_set(:,lhc_index) = residual_set(:,lhc_index) + residual_element(:, 1)
                residual_set(:,rhc_index) = residual_set(:,rhc_index) + residual_element(:, 2)
            end associate
        end do
    end subroutine compute_residual_objective_api

    subroutine compute_source_term_objective_api(self, variables_set, residual_set, num_conservative_variables, a_model)
        class  (cellsystem), intent(in   ) :: self
        real   (real_kind ), intent(in   ) :: variables_set(:,:)
        real   (real_kind ), intent(inout) :: residual_set (:,:)
        integer(int_kind  ), intent(in   ) :: num_conservative_variables
        class  (model     ), intent(in   ) :: a_model

        integer(int_kind ) :: i
        integer(int_kind ) :: element(num_conservative_variables)

#ifdef _DEBUG
        call write_debuginfo("In compute_source_term_objective_api(), cellsystem.")
#endif

!$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            element(:) = a_model%compute_source_term(variables_set(:,i), num_conservative_variables)
            residual_set(:,i) = residual_set(:,i) + element(:)
        end do
    end subroutine compute_source_term_objective_api

    ! ### Riemann solver ###
    subroutine riemann_solver_initialize(self, a_riemann_solver, config, num_conservative_variables, num_primitive_variables)
        class  (cellsystem    ), intent(inout) :: self
        class  (riemann_solver), intent(inout) :: a_riemann_solver
        class  (configuration ), intent(inout) :: config
        integer(int_kind      ), intent(in   ) :: num_conservative_variables
        integer(int_kind      ), intent(in   ) :: num_primitive_variables
#ifdef _DEBUG
        call write_debuginfo("In riemann_solver_initialize(), cellsystem.")
#endif
        call a_riemann_solver%initialize(config)
    end subroutine riemann_solver_initialize

    ! ### Reconstructor ###
    subroutine reconstructor_initialize(self, a_reconstructor, config, num_conservative_variables, num_primitive_variables)
        class  (cellsystem    ), intent(inout) :: self
        class  (reconstructor ), intent(inout) :: a_reconstructor
        class  (configuration ), intent(inout) :: config
        integer(int_kind      ), intent(in   ) :: num_conservative_variables
        integer(int_kind      ), intent(in   ) :: num_primitive_variables
#ifdef _DEBUG
        call write_debuginfo("In reconstructor_initialize(), cellsystem.")
#endif
        call a_reconstructor%initialize(config)
    end subroutine reconstructor_initialize

    ! ### Time stepping ###
    subroutine time_stepping_initialize(self, a_time_stepping, config, num_conservative_variables, num_primitive_variables)
        class  (cellsystem   ), intent(inout) :: self
        class  (time_stepping), intent(inout) :: a_time_stepping
        class  (configuration), intent(inout) :: config
        integer(int_kind     ), intent(in   ) :: num_conservative_variables
        integer(int_kind     ), intent(in   ) :: num_primitive_variables
#ifdef _DEBUG
        call write_debuginfo("In time_stepping_initialize(), cellsystem.")
#endif
        call a_time_stepping%initialize(config, self%num_cells, num_conservative_variables)
    end subroutine time_stepping_initialize

    subroutine compute_next_state(self, a_time_stepping, an_eos, state_num, conservative_variables_set, primitive_variables_set, residual_set, num_primitive_variables, &
        conservative_to_primitive_function)

        class  (cellsystem         ), intent(inout) :: self
        class  (time_stepping      ), intent(inout) :: a_time_stepping
        class  (eos                ), intent(in   ) :: an_eos
        integer(int_kind           ), intent(in   ) :: state_num
        real   (real_kind          ), intent(inout) :: conservative_variables_set(:,:)
        real   (real_kind          ), intent(inout) :: primitive_variables_set   (:,:)
        real   (real_kind          ), intent(inout) :: residual_set              (:,:)
        integer(int_kind           ), intent(in   ) :: num_primitive_variables

        interface
            pure function conservative_to_primitive_function(an_eos, conservative_variables, num_primitive_variables) result(primitive_variables)
                use typedef_module
                use abstract_eos
                class  (eos      ), intent(in)  :: an_eos
                real   (real_kind), intent(in)  :: conservative_variables(:)
                integer(int_kind ), intent(in)  :: num_primitive_variables
                real   (real_kind)              :: primitive_variables(num_primitive_variables)
            end function conservative_to_primitive_function
        end interface

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In time_stepping_initialize(), cellsystem.")
#endif

!$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            call a_time_stepping%compute_next_state(i, state_num, self%time_increment, conservative_variables_set(:,i), residual_set(:,i))
            primitive_variables_set(:,i) = conservative_to_primitive_function(an_eos, conservative_variables_set(:,i), num_primitive_variables)
        end do
    end subroutine compute_next_state

    subroutine prepare_stepping(    &
        self                      , &
        a_time_stepping           , &
        conservative_variables_set, &
        primitive_variables_set   , &
        residual_set                  )
        class(cellsystem   ), intent(inout) :: self
        class(time_stepping), intent(inout) :: a_time_stepping
        real (real_kind    ), intent(inout) :: conservative_variables_set(:,:)
        real (real_kind    ), intent(inout) :: primitive_variables_set   (:,:)
        real (real_kind    ), intent(inout) :: residual_set              (:,:)

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In prepare_stepping(), cellsystem.")
#endif

!$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            call a_time_stepping%prepare_stepping(i, conservative_variables_set(:,i), primitive_variables_set(:,i), residual_set(:,i))
        end do
    end subroutine

    pure function get_number_of_states(self, a_time_stepping) result(n)
        class  (cellsystem   ), intent(in) :: self
        class  (time_stepping), intent(in) :: a_time_stepping
        integer(int_kind     )                :: n
        n = a_time_stepping%get_number_of_states()
    end function get_number_of_states

    ! ###  Result writer ###
    subroutine result_writer_initialize(self, writer, config, num_conservative_variables, num_primitive_variables)
        class  (cellsystem   ), intent(inout) :: self
        class  (result_writer), intent(inout) :: writer
        class  (configuration), intent(inout) :: config
        integer(int_kind     ), intent(in   ) :: num_conservative_variables
        integer(int_kind     ), intent(in   ) :: num_primitive_variables
#ifdef _DEBUG
        call write_debuginfo("In result_writer_initialize(), cellsystem.")
#endif
        call writer%initialize(self%num_cells, self%num_points, self%is_real_cell, self%cell_geometries, self%cell_types, config)
    end subroutine result_writer_initialize

    subroutine result_writer_open_file(self, writer)
        class(cellsystem   ), intent(inout) :: self
        class(result_writer), intent(inout) :: writer
#ifdef _DEBUG
        call write_debuginfo("In result_writer_open_file(), cellsystem.")
#endif
        call writer%open_file(self%time, self%points)
    end subroutine result_writer_open_file

    subroutine result_writer_close_file(self, writer)
        class(cellsystem   ), intent(inout) :: self
        class(result_writer), intent(inout) :: writer
#ifdef _DEBUG
        call write_debuginfo("In result_writer_close_file(), cellsystem.")
#endif
        call writer%close_file()
    end subroutine result_writer_close_file

    pure function result_writer_is_writable(self, writer) result(yes)
        class(cellsystem   ), intent(in) :: self
        class(result_writer), intent(in) :: writer
        logical                          :: yes
        yes = writer%is_writable(self%time)
    end function result_writer_is_writable

    subroutine write_scolar(self, writer, name, scolar_variables)
        class    (cellsystem   ), intent(inout) :: self
        class    (result_writer), intent(inout) :: writer
        character(len=*        ), intent(in   ) :: name
        real     (real_kind    ), intent(in   ) :: scolar_variables(:)
#ifdef _DEBUG
        call write_debuginfo("In write_scolar(), cellsystem.")
#endif
        call writer%write_scolar(self%is_real_cell, name, scolar_variables)
    end subroutine write_scolar

    subroutine write_vector(self, writer, name, vector_variables)
        class    (cellsystem   ), intent(inout) :: self
        class    (result_writer), intent(inout) :: writer
        character(len=*        ), intent(in   ) :: name
        real     (real_kind    ), intent(in   ) :: vector_variables(:,:)
#ifdef _DEBUG
        call write_debuginfo("In write_vector(), cellsystem.")
#endif
        call writer%write_vector(self%is_real_cell, name, vector_variables)
    end subroutine write_vector

    function get_filename(self, writer) result(name)
        class    (cellsystem   ), intent(inout) :: self
        class    (result_writer), intent(inout) :: writer
        character(len=:        ), allocatable :: name
        name = writer%get_filename()
    end function get_filename

    ! ### Cellsystem ###
    subroutine read(self, parser, config)
        class(cellsystem         ), intent(inout) :: self
        class(grid_parser        ), intent(inout) :: parser
        class(configuration      ), intent(inout) :: config

        integer(kind(face_type)), allocatable :: face_types(:)

#ifdef _DEBUG
        call write_debuginfo("In read(), cellsystem.")
#endif

        ! parse grid file
        call parser%parse(config)

        ! initialize local valiables
        allocate(face_types(parser%get_number_of_faces()))

        ! allocate grid
        call self%initialise_faces              (parser%get_number_of_faces (), parser%get_number_of_ghost_cells())
        call self%initialise_cells              (parser%get_number_of_points(), parser%get_number_of_cells      ())
        call self%initialise_boundary_references(parser%get_number_of_boundary_faces(outflow_face_type  ), &
                                                 parser%get_number_of_boundary_faces(wall_face_type), &
                                                 parser%get_number_of_boundary_faces(symmetric_face_type), &
                                                 parser%get_number_of_boundary_faces(empty_face_type    )    )

        ! get grid data
        call parser%get_cells          (self%cell_centor_positions, self%cell_volumes, self%is_real_cell)
        call parser%get_faces          (self%face_to_cell_indexes, self%face_normal_vectors, self%face_tangential1_vectors, self%face_tangential2_vectors, self%face_positions, self%face_areas)
        call parser%get_cell_geometries(self%points, self%cell_geometries, self%cell_types)

        ! boundary condition
        call parser%get_boundaries(face_types)
        call self%assign_boundary (face_types)

        ! close file
        call parser%close()

#ifdef _DEBUG
        call write_debuginfo("Read following mesh.")
        print *, "Number of cells          : ", self%num_cells
        print *, "Number of faces          : ", self%num_faces
        print *, "Number of symmetric faces: ", self%num_symmetric_faces
        print *, "Number of empty faces    : ", self%num_empty_faces
#endif

        self%read_cellsystem = .true.
    end subroutine read

    subroutine increment_time(self)
        class(cellsystem), intent(inout) :: self
#ifdef _DEBUG
        call write_debuginfo("In increment_time(), cellsystem.")
#endif
        self%time = self%time + self%time_increment
        self%num_steps = self%num_steps + 1
    end subroutine increment_time

    ! ### Getter ###
    pure function get_number_of_faces(self) result(n)
        class(cellsystem), intent(in) :: self
        integer(int_kind) :: n
        n = self%num_faces
    end function get_number_of_faces

    pure function get_number_of_local_cells(self) result(n)
        class(cellsystem), intent(in) :: self
        integer(int_kind) :: n
        n = self%num_local_cells
    end function get_number_of_local_cells

    pure function get_number_of_cells(self) result(n)
        class(cellsystem), intent(in) :: self
        integer(int_kind) :: n
        n = self%num_cells
    end function get_number_of_cells

    pure function get_number_of_points(self) result(n)
        class(cellsystem), intent(in) :: self
        integer(int_kind) :: n
        n = self%num_points
    end function

    pure function get_number_of_outflow_faces(self) result(n)
        class(cellsystem), intent(in) :: self
        integer(int_kind) :: n
        n = self%num_outflow_faces
    end function get_number_of_outflow_faces

    pure function get_number_of_wall_faces(self) result(n)
        class(cellsystem), intent(in) :: self
        integer(int_kind) :: n
        n = self%num_wall_faces
    end function get_number_of_wall_faces

    pure function get_number_of_symmetric_faces(self) result(n)
        class(cellsystem), intent(in) :: self
        integer(int_kind) :: n
        n = self%num_symmetric_faces
    end function get_number_of_symmetric_faces

    pure function get_number_of_empty_faces(self) result(n)
        class(cellsystem), intent(in) :: self
        integer(int_kind) :: n
        n = self%num_empty_faces
    end function get_number_of_empty_faces

    pure function get_number_of_steps(self) result(n)
        class(cellsystem), intent(in) :: self
        integer(int_kind) :: n
        n = self%num_steps
    end function get_number_of_steps

    pure function get_time(self) result(t)
        class(cellsystem), intent(in) :: self
        real(real_kind) :: t
        t = self%time
    end function get_time

    pure function get_time_increment(self) result(dt)
        class(cellsystem), intent(in) :: self
        real(real_kind) :: dt
        dt = self%time_increment
    end function get_time_increment

    ! ### Finalizer ###
    subroutine finalize(self)
        class(cellsystem), intent(inout) :: self
#ifdef _DEBUG
        call write_debuginfo("In finalize(), cellsystem.")
#endif
        call self%finalize_cells()
        call self%finalize_faces()
        call self%finalize_boundary_references()
    end subroutine finalize

    ! ### Inner utils ###
    subroutine assign_boundary(self, face_types)
        class  (cellsystem     ), intent(inout) :: self
        integer(kind(face_type)), intent(in   ) :: face_types(:)

        integer(int_kind) :: index, outflow_index, wall_index, symmetric_index, empty_index

#ifdef _DEBUG
        call write_debuginfo("In assign_boundary(), cellsystem.")
#endif

        outflow_index   = 0
        wall_index       = 0
        symmetric_index = 0
        empty_index     = 0

        do index = 1, self%num_faces, 1
            if (face_types(index) == outflow_face_type) then
                outflow_index = outflow_index + 1
                self%outflow_face_indexes(outflow_index) = index
            else if (face_types(index) == wall_face_type) then
                wall_index = wall_index + 1
                self%wall_face_indexes(wall_index) = index
            else if (face_types(index) == symmetric_face_type) then
                symmetric_index = symmetric_index + 1
                self%symmetric_face_indexes(symmetric_index) = index
            else if (face_types(index) == empty_face_type) then
                empty_index = empty_index + 1
                self%empty_face_indexes(empty_index) = index
            else if (face_types(index) == face_type) then
                ! It is inner face.
            else
                call call_error("Unknown boundary type found.")
            end if
        end do
    end subroutine assign_boundary

    subroutine initialise_faces(self, num_faces, num_local_cells)
        class(cellsystem), intent(inout) :: self
        integer(int_kind), intent(in) :: num_faces
        integer(int_kind), intent(in) :: num_local_cells

#ifdef _DEBUG
        call write_debuginfo("In initialise_faces(), cellsystem.")
#endif

        if(num_faces < 1)then
            call call_error("Number of face must be set over zero.")
        end if
        self%num_faces = num_faces

        if(num_local_cells < 1)then
            call call_error("Number of ghost cells must be set over zero.")
        end if
        self%num_local_cells = num_local_cells

        if(allocated(self%face_to_cell_indexes))then
            call call_error("Array face_to_cell_indexes is already allocated. But you call the initialiser for faces.")
        end if
        allocate(self%face_to_cell_indexes(2*self%num_local_cells, self%num_faces))

        if(allocated(self%face_normal_vectors))then
            call call_error("Array face_normal_vectors is already allocated. But you call the initialiser for faces.")
        end if
        allocate(self%face_normal_vectors(3, self%num_faces))

        if(allocated(self%face_tangential1_vectors))then
            call call_error("Array face_tangential1_vectors is already allocated. But you call the initialiser for faces.")
        end if
        allocate(self%face_tangential1_vectors(3, self%num_faces))

        if(allocated(self%face_tangential2_vectors))then
            call call_error("Array face_tangential2_vectors is already allocated. But you call the initialiser for faces.")
        end if
        allocate(self%face_tangential2_vectors(3, self%num_faces))

        if(allocated(self%face_positions))then
            call call_error("Array face_positions is already allocated. But you call the initialiser for faces.")
        end if
        allocate(self%face_positions(3, self%num_faces))

        if(allocated(self%face_areas))then
            call call_error("Array face_areas is already allocated. But you call the initialiser for faces.")
        end if
        allocate(self%face_areas(self%num_faces))
    end subroutine initialise_faces

    subroutine initialise_cells(self, num_points, num_cells)
        class(cellsystem), intent(inout) :: self
        integer(int_kind), intent(in) :: num_points, num_cells

#ifdef _DEBUG
        call write_debuginfo("In initialise_cells(), cellsystem.")
#endif

        if(num_points < 1)then
            call call_error("Number of points must be set over zero.")
        end if
        self%num_points = num_points

        if(num_cells < 1)then
            call call_error("Number of cells must be set over zero.")
        end if
        self%num_cells = num_cells

        if(allocated(self%cell_centor_positions))then
            call call_error("Array cell_centor_positions is already allocated. But you call the initialiser for cells.")
        end if
        allocate(self%cell_centor_positions(3, self%num_cells))

        if(allocated(self%cell_volumes))then
            call call_error("Array cell_volumes is already allocated. But you call the initialiser for cells.")
        end if
        allocate(self%cell_volumes(self%num_cells))

        if(allocated(self%is_real_cell))then
            call call_error("Array is_real_cell is already allocated. But you call the initialiser for cells.")
        end if
        allocate(self%is_real_cell(self%num_cells))

        if(allocated(self%points))then
            call call_error("Array points is already allocated. But you call the initialiser for cells.")
        end if
        allocate(self%points(3, self%num_points))

        if(allocated(self%cell_geometries))then
            call call_error("Array cell_geometries is already allocated. But you call the initialiser for cells.")
        end if
        allocate(self%cell_geometries(self%num_cells))

        if(allocated(self%cell_types))then
            call call_error("Array cell_types is already allocated. But you call the initialiser for cells.")
        end if
        allocate(self%cell_types(self%num_cells))
    end subroutine initialise_cells

    subroutine initialise_boundary_references(self, num_outflow_faces, num_wall_faces, num_symetric_faces, num_empty_faces)
        class(cellsystem), intent(inout) :: self
        integer(int_kind), intent(in) :: num_outflow_faces
        integer(int_kind), intent(in) :: num_wall_faces
        integer(int_kind), intent(in) :: num_symetric_faces
        integer(int_kind), intent(in) :: num_empty_faces

#ifdef _DEBUG
        call write_debuginfo("In initialise_boundary_references(), cellsystem.")
#endif

        if(num_outflow_faces < 0)then
            call call_error("Number of outflow BC faces must be set over zero.")
        end if
        self%num_outflow_faces = num_outflow_faces

        if(num_wall_faces < 0)then
            call call_error("Number of wall BC faces must be set over zero.")
        end if
        self%num_wall_faces = num_wall_faces

        if(num_symetric_faces < 0)then
            call call_error("Number of symetric BC faces must be set over zero.")
        end if
        self%num_symmetric_faces = num_symetric_faces

        if(num_empty_faces < 0)then
            call call_error("Number of empty BC faces must be set over zero.")
        end if
        self%num_empty_faces = num_empty_faces

        if(allocated(self%outflow_face_indexes))then
            call call_error("Array outflow_face_indexes is already allocated. But you call the initialiser for boundary_reference.")
        end if
        allocate(self%outflow_face_indexes(self%num_outflow_faces))

        if(allocated(self%wall_face_indexes))then
            call call_error("Array wall_face_indexes is already allocated. But you call the initialiser for boundary_reference.")
        end if
        allocate(self%wall_face_indexes(self%num_wall_faces))

        if(allocated(self%symmetric_face_indexes))then
            call call_error("Array symmetric_face_indexes is already allocated. But you call the initialiser for boundary_reference.")
        end if
        allocate(self%symmetric_face_indexes(self%num_symmetric_faces))

        if(allocated(self%empty_face_indexes))then
            call call_error("Array empty_face_indexes is already allocated. But you call the initialiser for boundary_reference.")
        end if
        allocate(self%empty_face_indexes(self%num_empty_faces))
    end subroutine initialise_boundary_references

    subroutine finalize_faces(self)
        class(cellsystem), intent(inout) :: self

#ifdef _DEBUG
        call write_debuginfo("In finalize_faces(), cellsystem.")
#endif

        if(.not. allocated(self%face_to_cell_indexes))then
            call call_error("Array face_to_cell_indexes are not allocated. But you call the finalizer for faces.")
        end if
        deallocate(self%face_to_cell_indexes)

        if(.not. allocated(self%face_normal_vectors))then
            call call_error("Array face_normal_vectors is not allocated. But you call the finalizer for faces.")
        end if
        deallocate(self%face_normal_vectors)

        if(.not. allocated(self%face_tangential1_vectors))then
            call call_error("Array face_tangential1_vectors is not allocated. But you call the finalizer for faces.")
        end if
        deallocate(self%face_tangential1_vectors)

        if(.not. allocated(self%face_tangential2_vectors))then
            call call_error("Array face_tangential2_vectors is not allocated. But you call the finalizer for faces.")
        end if
        deallocate(self%face_tangential2_vectors)

        if(.not. allocated(self%face_positions))then
            call call_error("Array face_positions is not allocated. But you call the finalizer for faces.")
        end if
        deallocate(self%face_positions)

        if(.not. allocated(self%face_areas))then
            call call_error("Array face_areas is not allocated. But you call the finalizer for faces.")
        end if
        deallocate(self%face_areas)

        self%num_faces = 0
        self%num_local_cells = 0
    end subroutine finalize_faces

    subroutine finalize_cells(self)
        class(cellsystem), intent(inout) :: self

#ifdef _DEBUG
        call write_debuginfo("In finalize_cells(), cellsystem.")
#endif

        if(.not. allocated(self%cell_centor_positions))then
            call call_error("Array cell_centor_positions is not allocated. But you call the finalizer for cell_geometries_module.")
        end if
        deallocate(self%cell_centor_positions)

        if(.not. allocated(self%cell_volumes))then
            call call_error("Array cell_volumes is not allocated. But you call the finalizer for cell_geometries_module.")
        end if
        deallocate(self%cell_volumes)

        if(.not. allocated(self%is_real_cell))then
            call call_error("Array is_real_cell is not allocated. But you call the finalizer for cell_geometries_module.")
        end if
        deallocate(self%is_real_cell)

        if(.not. allocated(self%points))then
            call call_error("Array points is not allocated. But you call the finalizer for cell_geometries_module.")
        end if
        deallocate(self%points)

        if(.not. allocated(self%cell_geometries))then
            call call_error("Array cell_geometries is not allocated. But you call the finalizer for cell_geometries_module.")
        end if
        deallocate(self%cell_geometries)

        if(.not. allocated(self%cell_types))then
            call call_error("Array cell_types is not allocated. But you call the finalizer for cell_geometries_module.")
        end if
        deallocate(self%cell_types)

        self%num_points = 0
        self%num_cells  = 0
    end subroutine finalize_cells

    subroutine finalize_boundary_references(self)
        class(cellsystem), intent(inout) :: self

#ifdef _DEBUG
        call write_debuginfo("In finalize_boundary_references(), cellsystem.")
#endif

        if(.not. allocated(self%outflow_face_indexes))then
            call call_error("Array outflow_face_indexes is not allocated. But you call the finalizer for boundary_reference module.")
        end if
        deallocate(self%outflow_face_indexes)

        if(.not. allocated(self%wall_face_indexes))then
            call call_error("Array wall_face_indexes is not allocated. But you call the finalizer for boundary_reference module.")
        end if
        deallocate(self%wall_face_indexes)

        if(.not. allocated(self%symmetric_face_indexes))then
            call call_error("Array symmetric_face_indexes is not allocated. But you call the finalizer for boundary_reference module.")
        end if
        deallocate(self%symmetric_face_indexes)

        if(.not. allocated(self%empty_face_indexes))then
            call call_error("Array empty_face_indexes is not allocated. But you call the finalizer for boundary_reference module.")
        end if
        deallocate(self%empty_face_indexes)

        self%num_local_cells     = 0
        self%num_outflow_faces   = 0
        self%num_wall_faces  = 0
        self%num_symmetric_faces = 0
    end subroutine finalize_boundary_references

    pure function compute_boundary_gradient(self, lhc_variables, rhc_variables, lhc_position, rhc_position, num_variables) result(gradient_variables)
        class  (cellsystem), intent(in) :: self
        real   (real_kind ), intent(in) :: lhc_variables(:)
        real   (real_kind ), intent(in) :: rhc_variables(:)
        real   (real_kind ), intent(in) :: lhc_position(3)
        real   (real_kind ), intent(in) :: rhc_position(3)
        integer(int_kind  ), intent(in) :: num_variables
        real   (real_kind )             :: gradient_variables(num_variables*3)
        integer(int_kind  ) :: i
        do i = 1, num_variables, 1
            associate(                                           &
                grad => gradient_variables(3*(i-1)+1:3*(i-1)+3), &
                dphi => lhc_variables(i) - rhc_variables(i)    , &
                dx   => lhc_position - rhc_position              &
            )
                grad(1:3) = vector_normalize(dx) * (dphi / vector_magnitude(dx))
            end associate
        end do
    end function compute_boundary_gradient
end module class_cellsystem