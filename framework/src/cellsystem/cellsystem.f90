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
    use abstract_interpolator
    use abstract_termination_criterion
    use abstract_time_increment_controller
    use abstract_initial_condition_parser
    use abstract_face_gradient_interpolator
    use abstract_face_gradient_calculator
    use class_line_plotter
    use class_control_volume_profiler
    use class_surface_profiler
#ifdef _OPENMP
    use omp_lib
#endif

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

        ! Parallel computing
        integer, private :: num_threads = 1

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
        integer(int_kind) :: num_nonslip_wall_faces
        integer(int_kind) :: num_slip_wall_faces
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
        integer(int_kind), allocatable :: outflow_face_indexes      (:)
        integer(int_kind), allocatable :: nonslip_wall_face_indexes (:)
        integer(int_kind), allocatable :: slip_wall_face_indexes    (:)
        integer(int_kind), allocatable :: symmetric_face_indexes    (:)
        integer(int_kind), allocatable :: empty_face_indexes        (:)

        ! ### Status ###
        logical :: read_cellsystem = .false.

        contains
        ! ### Common ###
        procedure, public, pass(self) :: initilaize_cellsystem
        procedure, public, pass(self) :: initialize_result_writer
        procedure, public, pass(self) :: initialize_time_stepping
        procedure, public, pass(self) :: initialize_reconstructor
        procedure, public, pass(self) :: initialize_riemann_solver
        procedure, public, pass(self) :: initialize_eos
        procedure, public, pass(self) :: initialize_gradient_calculator
        procedure, public, pass(self) :: initialize_variables_rank2
        procedure, public, pass(self) :: initialize_variables_rank1
        procedure, public, pass(self) :: initialize_termination_ctiterion
        procedure, public, pass(self) :: initialize_time_increment_controller
        procedure, public, pass(self) :: initialize_interpolator
        procedure, public, pass(self) :: initialize_line_plotter
        procedure, public, pass(self) :: initialize_control_volume_profiler
        procedure, public, pass(self) :: initialize_surface_profiler
        procedure, public, pass(self) :: initialize_face_gradient_interpolator
        procedure, public, pass(self) :: initialize_face_gradient_calculator
        generic  , public             :: initialize  => initilaize_cellsystem                , &
                                                        initialize_result_writer             , &
                                                        initialize_time_stepping             , &
                                                        initialize_reconstructor             , &
                                                        initialize_riemann_solver            , &
                                                        initialize_eos                       , &
                                                        initialize_gradient_calculator       , &
                                                        initialize_variables_rank2           , &
                                                        initialize_variables_rank1           , &
                                                        initialize_termination_ctiterion     , &
                                                        initialize_time_increment_controller , &
                                                        initialize_interpolator              , &
                                                        initialize_line_plotter              , &
                                                        initialize_control_volume_profiler   , &
                                                        initialize_surface_profiler          , &
                                                        initialize_face_gradient_interpolator, &
                                                        initialize_face_gradient_calculator
        procedure, public, pass(self) :: result_writer_open_file
        generic  , public             :: open_file   => result_writer_open_file
        procedure, public, pass(self) :: result_writer_close_file
        generic  , public             :: close_file  => result_writer_close_file
        procedure, public, pass(self) :: result_writer_is_writable
        procedure, public, pass(self) :: line_plotter_is_writable
        procedure, public, pass(self) :: control_volume_profiler_is_writable
        procedure, public, pass(self) :: surface_profiler_is_writable
        generic  , public             :: is_writable => result_writer_is_writable          , &
                                                        line_plotter_is_writable           , &
                                                        control_volume_profiler_is_writable, &
                                                        surface_profiler_is_writable
        procedure, public, pass(self) :: line_plotter_write
        procedure, public, pass(self) :: control_volume_profiler_write
        procedure, public, pass(self) :: surface_profiler_write
        generic  , public             :: write => line_plotter_write           , &
                                                  control_volume_profiler_write, &
                                                  surface_profiler_write

        ! ### Show infomation ###
        procedure, public, pass(self) :: show_timestepping_infomation

        ! ### Vatiables ###
        procedure, public, pass(self) :: read_initial_condition
        procedure, public, pass(self) :: conservative_to_primitive_variables_all
        procedure, public, pass(self) :: operate_cellwise_rank2
        procedure, public, pass(self) :: operate_cellwise_subarray_rank2
        generic  , public             :: operate_cellwise => operate_cellwise_rank2, &
                                                             operate_cellwise_subarray_rank2
        procedure, public, pass(self) :: smooth_variables
        procedure, public, pass(self) :: substitute_rank1
        procedure, public, pass(self) :: substitute_rank2
        generic  , public             :: substitute => substitute_rank1, &
                                                       substitute_rank2
        procedure, public, pass(self) :: substitute_zeros_rank1
        procedure, public, pass(self) :: substitute_zeros_rank2
        generic  , public             :: substitute_zeros => substitute_zeros_rank1, &
                                                             substitute_zeros_rank2

        ! ### Time increment control ###
        procedure, public, pass(self) :: update_time_increment

        ! ### Termination criterion ###
        procedure, public, pass(self) :: satisfy_termination_criterion

        ! ### Boundary Condition ###
        procedure, public , pass(self) :: apply_boundary_condition
        procedure, private, pass(self) :: apply_outflow_condition
        procedure, private, pass(self) :: apply_nonslip_wall_condition
        procedure, private, pass(self) :: apply_symmetric_condition
        procedure, private, pass(self) :: apply_slip_wall_condition
        procedure, private, pass(self) :: apply_empty_condition

        ! ### Gradient Calculator ###
        procedure, public, pass(self) :: compute_gradient_rank1
        procedure, public, pass(self) :: compute_gradient_rank2
        generic  , public             :: compute_gradient => compute_gradient_rank1, &
                                                             compute_gradient_rank2

        ! ### Divergence Calculator ###
        ! godunov: Compute flux by Godunov scheme.
        ! facegrad: Provide a face gradient variable(s) to a element function.
        procedure, public, pass(self) :: compute_divergence_rank1
        procedure, public, pass(self) :: compute_divergence_rank2
        procedure, public, pass(self) :: compute_divergence_godunov_rank2
        procedure, public, pass(self) :: compute_divergence_godunov_facegrad_rank2
        generic  , public             :: compute_divergence => compute_divergence_rank1, &
                                                               compute_divergence_rank2, &
                                                               compute_divergence_godunov_rank2, &
                                                               compute_divergence_godunov_facegrad_rank2

        ! ### Sorce Term
        procedure, public, pass(self) :: compute_source_term

        ! ### Time Stepping ###
        procedure, public, pass(self) :: compute_next_stage_primitive_rank2
        procedure, public, pass(self) :: compute_next_stage_rank2
        procedure, public, pass(self) :: compute_next_stage_rank1
        generic  , public             :: compute_next_stage => compute_next_stage_primitive_rank2, &
                                                               compute_next_stage_rank2          , &
                                                               compute_next_stage_rank1
        procedure, public, pass(self) :: prepare_time_stepping
        procedure, public, pass(self) :: get_number_of_stages

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
        procedure, public, pass(self) :: get_number_of_nonslip_wall_faces
        procedure, public, pass(self) :: get_number_of_slip_wall_faces
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

    ! ### Procedure interfaces ###
    interface
        subroutine operator_subroutine_interface(operated_variables)
            use typedef_module
            real   (real_kind ), intent(inout) :: operated_variables(:)
        end subroutine operator_subroutine_interface

        subroutine operator_with_subarray_subroutine_interface(operated_variables, sub_variables)
            use typedef_module
            real   (real_kind ), intent(inout) :: operated_variables(:)
            real   (real_kind ), intent(in   ) :: sub_variables(:)
        end subroutine operator_with_subarray_subroutine_interface

        function weight_function_interface(own_cell_position, neighbor_cell_position, own_variables, neighbor_variables) result(weight)
            use typedef_module
            real   (real_kind ), intent(in) :: own_cell_position     (3)
            real   (real_kind ), intent(in) :: neighbor_cell_position(3)
            real   (real_kind ), intent(in) :: own_variables         (:)
            real   (real_kind ), intent(in) :: neighbor_variables    (:)
            real   (real_kind )             :: weight
        end function weight_function_interface

        pure function compute_rotate_variables_function_interface( &
            variables                                            , &
            face_normal_vector                                   , &
            face_tangential1_vector                              , &
            face_tangential2_vector                              , &
            num_variables                                            ) result(rotate_variables)

            use typedef_module
            real   (real_kind     ), intent(in)  :: variables               (:)
            real   (real_kind     ), intent(in)  :: face_normal_vector      (3)
            real   (real_kind     ), intent(in)  :: face_tangential1_vector (3)
            real   (real_kind     ), intent(in)  :: face_tangential2_vector (3)
            integer(int_kind      ), intent(in)  :: num_variables
            real   (real_kind     )              :: rotate_variables(num_variables)
        end function compute_rotate_variables_function_interface

        pure function compute_unrotate_variables_function_interface( &
            variables                                              , &
            face_normal_vector                                     , &
            face_tangential1_vector                                , &
            face_tangential2_vector                                , &
            num_variables                                              ) result(unrotate_variables)

            use typedef_module
            real   (real_kind     ), intent(in)  :: variables               (:)
            real   (real_kind     ), intent(in)  :: face_normal_vector      (3)
            real   (real_kind     ), intent(in)  :: face_tangential1_vector (3)
            real   (real_kind     ), intent(in)  :: face_tangential2_vector (3)
            integer(int_kind      ), intent(in)  :: num_variables
            real   (real_kind     )              :: unrotate_variables(num_variables)
        end function compute_unrotate_variables_function_interface

        function boundary_condition_function_interface(inner_variables, num_variables) result(ghost_variables)
            use typedef_module
            real   (real_kind), intent(in) :: inner_variables(:)
            integer(int_kind ), intent(in) :: num_variables
            real   (real_kind)             :: ghost_variables(num_variables)
        end function boundary_condition_function_interface

        pure function primitive_to_conservative_function_interface(an_eos, primitive_variables, num_conservative_values) result(conservative_values)
            use typedef_module
            use abstract_eos
            class  (eos      ), intent(in) :: an_eos
            real   (real_kind), intent(in) :: primitive_variables(:)
            integer(int_kind ), intent(in) :: num_conservative_values
            real   (real_kind)             :: conservative_values(num_conservative_values)
        end function primitive_to_conservative_function_interface

        pure function conservative_to_primitive_function_interface(an_eos, conservative_variables, num_primitive_variables) result(primitive_variables)
            use typedef_module
            use abstract_eos
            class  (eos      ), intent(in)  :: an_eos
            real   (real_kind), intent(in)  :: conservative_variables(:)
            integer(int_kind ), intent(in)  :: num_primitive_variables
            real   (real_kind)              :: primitive_variables(num_primitive_variables)
        end function conservative_to_primitive_function_interface

        function spectral_radius_function_interface(an_eos, variables, length) result(r)
            use abstract_eos
            use typedef_module
            class  (eos      ), intent(in) :: an_eos
            real   (real_kind), intent(in) :: variables(:)
            real   (real_kind), intent(in) :: length
            real   (real_kind) :: r
        end function spectral_radius_function_interface

        function compute_source_term_function_interface(variables, num_conservative_variables) result(source)
            use typedef_module
            use abstract_eos
            real   (real_kind), intent(in) :: variables(:)
            integer(int_kind ), intent(in) :: num_conservative_variables
            real   (real_kind)             :: source(num_conservative_variables)
        end function compute_source_term_function_interface

        function element_provided_with_riemann_function_interface( &
            an_eos                                               , &
            an_riemann_solver                                    , &
            primitive_variables_lhc                              , &
            primitive_variables_rhc                              , &
            reconstructed_primitive_variables_lhc                , &
            reconstructed_primitive_variables_rhc                , &
            lhc_cell_volume                                      , &
            rhc_cell_volume                                      , &
            face_area                                            , &
            face_normal_vector                                   , &
            face_tangential1_vector                              , &
            face_tangential2_vector                              , &
            num_conservative_values                              , &
            num_primitive_values                                     ) result(flux)
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
            real   (real_kind)                  :: flux(num_conservative_values, 1:2)
        end function element_provided_with_riemann_function_interface

        function element_provided_with_riemann_facegrad_function_interface(&
            an_eos                                                       , &
            an_riemann_solver                                            , &
            primitive_variables_lhc                                      , &
            primitive_variables_rhc                                      , &
            reconstructed_primitive_variables_lhc                        , &
            reconstructed_primitive_variables_rhc                        , &
            face_gradient_primitive_variables                            , &
            lhc_cell_volume                                              , &
            rhc_cell_volume                                              , &
            face_area                                                    , &
            face_normal_vector                                           , &
            face_tangential1_vector                                      , &
            face_tangential2_vector                                      , &
            num_conservative_variables                                   , &
            num_primitive_variables                                          ) result(element)
            use typedef_module
            use abstract_eos
            use abstract_riemann_solver
            class  (eos           ), intent(in) :: an_eos
            class  (riemann_solver), intent(in) :: an_riemann_solver
            real   (real_kind     ), intent(in) :: primitive_variables_lhc              (:)
            real   (real_kind     ), intent(in) :: primitive_variables_rhc              (:)
            real   (real_kind     ), intent(in) :: reconstructed_primitive_variables_lhc(:)
            real   (real_kind     ), intent(in) :: reconstructed_primitive_variables_rhc(:)
            real   (real_kind     ), intent(in) :: face_gradient_primitive_variables    (:)
            real   (real_kind     ), intent(in) :: lhc_cell_volume
            real   (real_kind     ), intent(in) :: rhc_cell_volume
            real   (real_kind     ), intent(in) :: face_area
            real   (real_kind     ), intent(in) :: face_normal_vector     (3)
            real   (real_kind     ), intent(in) :: face_tangential1_vector(3)
            real   (real_kind     ), intent(in) :: face_tangential2_vector(3)
            integer(int_kind      ), intent(in) :: num_conservative_variables
            integer(int_kind      ), intent(in) :: num_primitive_variables
            real   (real_kind     )             :: element(num_conservative_variables, 1:2)
        end function element_provided_with_riemann_facegrad_function_interface
    end interface

    contains

    ! ### Cellsystem ###
    subroutine initilaize_cellsystem(self, config)
        class(cellsystem   ), intent(inout) :: self
        class(configuration), intent(inout) :: config

        logical           :: found

#ifdef _OPENMP
        call config%get_int                ("Parallel computing.Number of threads", self%num_threads, found, 1)
        if(.not. found) call write_warring("'Parallel computing.Number of threads' is not found in configuration file you set. Solver is executed a single thread.")

        call omp_set_num_threads(self%num_threads)
#endif
    end subroutine initilaize_cellsystem

    subroutine show_timestepping_infomation(self)
        class(cellsystem), intent(inout) :: self
        print '(3(a, g0))', "Step: ", self%num_steps, ", Time increment: ", self%time_increment, ", Time: ", self%time
    end subroutine show_timestepping_infomation

    ! ### Face gradient interpolator ###
    subroutine initialize_face_gradient_interpolator(self, a_face_gradient_interpolator, config, num_conservative_variables)
        class   (cellsystem                ), intent(inout) :: self
        class   (face_gradient_interpolator), intent(inout) :: a_face_gradient_interpolator
        class   (configuration             ), intent(inout) :: config
        integer (int_kind                  ), intent(in   ) :: num_conservative_variables

#ifdef _DEBUG
        call write_debuginfo("In initialize_face_gradient_interpolator(), cellsystem.")
#endif
        call a_face_gradient_interpolator%initialize(config)
    end subroutine initialize_face_gradient_interpolator

    ! ### Face gradient calculator ###
    subroutine initialize_face_gradient_calculator(self, a_face_gradient_calculator, config, num_conservative_variables)
        class   (cellsystem              ), intent(inout) :: self
        class   (face_gradient_calculator), intent(inout) :: a_face_gradient_calculator
        class   (configuration           ), intent(inout) :: config
        integer (int_kind                ), intent(in   ) :: num_conservative_variables

#ifdef _DEBUG
        call write_debuginfo("In initialize_face_gradient_calculator(), cellsystem.")
#endif
        call a_face_gradient_calculator%initialize(config)
    end subroutine initialize_face_gradient_calculator

    ! ### Surface profiler ###
    subroutine initialize_surface_profiler(self, plotter, config, num_conservative_variables)
        class  (cellsystem      ), intent(inout) :: self
        class  (surface_profiler), intent(inout) :: plotter
        class  (configuration   ), intent(inout) :: config
        integer(int_kind        ), intent(in   ) :: num_conservative_variables

#ifdef _DEBUG
        call write_debuginfo("In initialize_surface_profiler(), cellsystem.")
#endif
        call plotter%initialize(config, self%face_to_cell_indexes, self%face_positions, self%face_normal_vectors, self%is_real_cell, self%num_faces, self%num_local_cells)
    end subroutine initialize_surface_profiler

    subroutine surface_profiler_write(self, plotter, values_set)
        class(cellsystem      ), intent(inout) :: self
        class(surface_profiler), intent(inout) :: plotter
        real (real_kind       ), intent(in   ) :: values_set(:,:)
#ifdef _DEBUG
        call write_debuginfo("In surface_profiler_write(), cellsystem.")
#endif
        call plotter%write(self%time, values_set, self%face_to_cell_indexes, self%face_areas, self%num_local_cells)
    end subroutine surface_profiler_write

    pure function surface_profiler_is_writable(self, plotter) result(juge)
        class(cellsystem      ), intent(in) :: self
        class(surface_profiler), intent(in) :: plotter
        logical                             :: juge
        juge = plotter%is_writable(self%time)
    end function surface_profiler_is_writable

    ! ### Control volume profiler ###
    subroutine initialize_control_volume_profiler(self, plotter, config, num_conservative_variables)
        class(cellsystem             ), intent(inout) :: self
        class(control_volume_profiler), intent(inout) :: plotter
        class(configuration          ), intent(inout) :: config
        integer(int_kind             ), intent(in   ) :: num_conservative_variables

#ifdef _DEBUG
        call write_debuginfo("In initialize_control_volume_profiler(), cellsystem.")
#endif
        call plotter%initialize(config, self%cell_centor_positions, self%is_real_cell, self%num_cells)
    end subroutine initialize_control_volume_profiler

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
    subroutine initialize_line_plotter(self, plotter, config, num_conservative_variables)
        class(cellsystem   ), intent(inout) :: self
        class(line_plotter ), intent(inout) :: plotter
        class(configuration), intent(inout) :: config
        integer(int_kind   ), intent(in   ) :: num_conservative_variables

#ifdef _DEBUG
        call write_debuginfo("In initialize_line_plotter(), cellsystem.")
#endif
        call plotter%initialize(config, self%cell_centor_positions, self%is_real_cell, self%num_cells)
    end subroutine initialize_line_plotter

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
    subroutine initialize_time_increment_controller(self, controller, config, num_conservative_variables)
        class(cellsystem               ), intent(inout) :: self
        class(time_increment_controller), intent(inout) :: controller
        class(configuration            ), intent(inout) :: config
        integer(int_kind               ), intent(in   ) :: num_conservative_variables

#ifdef _DEBUG
        call write_debuginfo("In initialize_time_increment_controller(), cellsystem.")
#endif
        call controller%initialize(config)
    end subroutine initialize_time_increment_controller

    subroutine update_time_increment(self, controller, an_eos, variables_set, spectral_radius_function)
        class  (cellsystem               ), intent(inout) :: self
        class  (time_increment_controller), intent(in   ) :: controller
        class  (eos                      ), intent(in   ) :: an_eos
        real   (real_kind                ), intent(in   ) :: variables_set(:,:)

        procedure(spectral_radius_function_interface) :: spectral_radius_function

        integer(int_kind ) :: i

#ifdef _DEBUG
        call write_debuginfo("In update_time_increment (), cellsystem.")
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
    end subroutine update_time_increment

    ! ### Termination criterion ###
    subroutine initialize_termination_ctiterion(self, criterion, config, num_conservative_variables)
        class  (cellsystem           ), intent(inout) :: self
        class  (termination_criterion), intent(inout) :: criterion
        class  (configuration        ), intent(inout) :: config
        integer(int_kind             ), intent(in   ) :: num_conservative_variables

#ifdef _DEBUG
        call write_debuginfo("In initialize_termination_ctiterion(), cellsystem.")
#endif
        call criterion%initialize(config)
    end subroutine initialize_termination_ctiterion

    pure function satisfy_termination_criterion(self, criterion) result(judge)
        class(cellsystem           ), intent(in) :: self
        class(termination_criterion), intent(in) :: criterion
        logical :: judge
        judge = criterion%is_satisfied(self%time, self%num_steps)
    end function satisfy_termination_criterion

    ! ### Boundary Condition ###
    subroutine apply_boundary_condition(self, variables_set, num_variables  , &
                                        compute_rotate_variables_function   , &
                                        compute_unrotate_variables_function , &
                                        empty_condition_function            , &
                                        symmetric_condition_function        , &
                                        nonslip_wall_condition_function     , &
                                        slip_wall_condition_function        , &
                                        outflow_condition_function               )
        class  (cellsystem  ), intent(inout) :: self
        real   (real_kind   ), intent(inout) :: variables_set(:,:)
        integer(int_kind    ), intent(in   ) :: num_variables

        procedure(compute_rotate_variables_function_interface  ) :: compute_rotate_variables_function
        procedure(compute_unrotate_variables_function_interface) :: compute_unrotate_variables_function

        procedure(boundary_condition_function_interface        ), optional :: empty_condition_function
        procedure(boundary_condition_function_interface        ), optional :: symmetric_condition_function
        procedure(boundary_condition_function_interface        ), optional :: nonslip_wall_condition_function
        procedure(boundary_condition_function_interface        ), optional :: slip_wall_condition_function
        procedure(boundary_condition_function_interface        ), optional :: outflow_condition_function

        integer :: i, face_index

#ifdef _DEBUG
        call write_debuginfo("In apply_boundary_condition(), cellsystem.")
#endif

        if (present(empty_condition_function)) call self%apply_empty_condition( &
            variables_set                                                     , &
            num_variables                                                     , &
            compute_rotate_variables_function                                 , &
            compute_unrotate_variables_function                               , &
            empty_condition_function                                            &
        )

        if (present(symmetric_condition_function)) call self%apply_symmetric_condition( &
            variables_set                                                             , &
            num_variables                                                             , &
            compute_rotate_variables_function                                         , &
            compute_unrotate_variables_function                                       , &
            symmetric_condition_function                                                &
        )

        if (present(nonslip_wall_condition_function)) call self%apply_nonslip_wall_condition( &
            variables_set                                                                   , &
            num_variables                                                                   , &
            compute_rotate_variables_function                                               , &
            compute_unrotate_variables_function                                             , &
            nonslip_wall_condition_function                                                   &
        )

        if (present(slip_wall_condition_function)) call self%apply_slip_wall_condition( &
            variables_set                                                             , &
            num_variables                                                             , &
            compute_rotate_variables_function                                         , &
            compute_unrotate_variables_function                                       , &
            slip_wall_condition_function                                                &
        )

        if (present(outflow_condition_function)) call self%apply_outflow_condition( &
            variables_set                                                         , &
            num_variables                                                         , &
            compute_rotate_variables_function                                     , &
            compute_unrotate_variables_function                                   , &
            outflow_condition_function                                              &
        )
    end subroutine apply_boundary_condition

    subroutine apply_outflow_condition(self, variables_set, num_variables, &
        compute_rotate_variables_function, compute_unrotate_variables_function, boundary_condition_function)
        class  (cellsystem  ), intent(inout) :: self
        real   (real_kind   ), intent(inout) :: variables_set(:,:)
        integer(int_kind    ), intent(in   ) :: num_variables

        procedure(compute_rotate_variables_function_interface  ) :: compute_rotate_variables_function
        procedure(compute_unrotate_variables_function_interface) :: compute_unrotate_variables_function
        procedure(boundary_condition_function_interface        ) :: boundary_condition_function

        integer :: i, face_index

#ifdef _DEBUG
        call write_debuginfo("In apply_outflow_condition(), cellsystem.")
#endif

!$omp parallel do private(i, face_index)
        do i = 1, self%num_outflow_faces, 1
            face_index = self%outflow_face_indexes(i)
            call apply_boundary_condition_common_impl(          &
                variables_set                                 , &
                self%face_normal_vectors     (:,face_index)   , &
                self%face_tangential1_vectors(:,face_index)   , &
                self%face_tangential2_vectors(:,face_index)   , &
                self%face_to_cell_indexes    (:,face_index)   , &
                self%num_local_cells                          , &
                num_variables                                 , &
                compute_rotate_variables_function             , &
                compute_unrotate_variables_function           , &
                boundary_condition_function                     &
            )
        end do
    end subroutine apply_outflow_condition

    subroutine apply_nonslip_wall_condition(self, variables_set, num_variables, &
        compute_rotate_variables_function, compute_unrotate_variables_function, boundary_condition_function)
        class  (cellsystem  ), intent(inout) :: self
        real   (real_kind   ), intent(inout) :: variables_set(:,:)
        integer(int_kind    ), intent(in   ) :: num_variables

        procedure(compute_rotate_variables_function_interface  ) :: compute_rotate_variables_function
        procedure(compute_unrotate_variables_function_interface) :: compute_unrotate_variables_function
        procedure(boundary_condition_function_interface        ) :: boundary_condition_function

        integer :: i, face_index

#ifdef _DEBUG
        call write_debuginfo("In apply_nonslip_wall_condition(), cellsystem.")
#endif

!$omp parallel do private(i, face_index)
        do i = 1, self%num_nonslip_wall_faces, 1
            face_index = self%nonslip_wall_face_indexes(i)
            call apply_boundary_condition_common_impl(         &
                variables_set                                , &
                self%face_normal_vectors     (:,face_index)  , &
                self%face_tangential1_vectors(:,face_index)  , &
                self%face_tangential2_vectors(:,face_index)  , &
                self%face_to_cell_indexes    (:,face_index)  , &
                self%num_local_cells                         , &
                num_variables                                , &
                compute_rotate_variables_function            , &
                compute_unrotate_variables_function          , &
                boundary_condition_function                    &
            )
        end do
    end subroutine apply_nonslip_wall_condition

    subroutine apply_symmetric_condition(self, variables_set, num_variables, &
        compute_rotate_variables_function, compute_unrotate_variables_function, boundary_condition_function)
        class  (cellsystem  ), intent(inout) :: self
        real   (real_kind   ), intent(inout) :: variables_set(:,:)
        integer(int_kind    ), intent(in   ) :: num_variables

        procedure(compute_rotate_variables_function_interface  ) :: compute_rotate_variables_function
        procedure(compute_unrotate_variables_function_interface) :: compute_unrotate_variables_function
        procedure(boundary_condition_function_interface        ) :: boundary_condition_function

        integer :: i, face_index

#ifdef _DEBUG
        call write_debuginfo("In apply_slip_and_symmetric_condition(), cellsystem.")
#endif

!$omp parallel do private(i, face_index)
        do i = 1, self%num_symmetric_faces, 1
            face_index = self%symmetric_face_indexes(i)
            call apply_boundary_condition_common_impl(          &
                variables_set                                 , &
                self%face_normal_vectors     (:,face_index)   , &
                self%face_tangential1_vectors(:,face_index)   , &
                self%face_tangential2_vectors(:,face_index)   , &
                self%face_to_cell_indexes    (:,face_index)   , &
                self%num_local_cells                          , &
                num_variables                                 , &
                compute_rotate_variables_function             , &
                compute_unrotate_variables_function           , &
                boundary_condition_function                     &
            )
        end do
    end subroutine apply_symmetric_condition

    subroutine apply_slip_wall_condition(self, variables_set, num_variables, &
        compute_rotate_variables_function, compute_unrotate_variables_function, boundary_condition_function)
        class  (cellsystem  ), intent(inout) :: self
        real   (real_kind   ), intent(inout) :: variables_set(:,:)
        integer(int_kind    ), intent(in   ) :: num_variables

        procedure(compute_rotate_variables_function_interface  ) :: compute_rotate_variables_function
        procedure(compute_unrotate_variables_function_interface) :: compute_unrotate_variables_function
        procedure(boundary_condition_function_interface        ) :: boundary_condition_function

        integer :: i, face_index

#ifdef _DEBUG
        call write_debuginfo("In apply_slip_wall_condition(), cellsystem.")
#endif

!$omp parallel do private(i, face_index)
        do i = 1, self%num_slip_wall_faces, 1
            face_index = self%slip_wall_face_indexes(i)
            call apply_boundary_condition_common_impl(          &
                variables_set                                 , &
                self%face_normal_vectors     (:,face_index)   , &
                self%face_tangential1_vectors(:,face_index)   , &
                self%face_tangential2_vectors(:,face_index)   , &
                self%face_to_cell_indexes    (:,face_index)   , &
                self%num_local_cells                          , &
                num_variables                                 , &
                compute_rotate_variables_function             , &
                compute_unrotate_variables_function           , &
                boundary_condition_function                     &
            )
        end do
    end subroutine apply_slip_wall_condition

    subroutine apply_empty_condition(self, variables_set, num_variables, &
        compute_rotate_variables_function, compute_unrotate_variables_function, boundary_condition_function)
        class  (cellsystem  ), intent(inout) :: self
        real   (real_kind   ), intent(inout) :: variables_set(:,:)
        integer(int_kind    ), intent(in   ) :: num_variables

        procedure(compute_rotate_variables_function_interface  ) :: compute_rotate_variables_function
        procedure(compute_unrotate_variables_function_interface) :: compute_unrotate_variables_function
        procedure(boundary_condition_function_interface        ) :: boundary_condition_function

        integer :: i, face_index

#ifdef _DEBUG
        call write_debuginfo("In apply_empty_condition(), cellsystem.")
#endif

!$omp parallel do private(i, face_index)
        do i = 1, self%num_empty_faces, 1
            face_index = self%empty_face_indexes(i)
            call apply_boundary_condition_empty_impl(           &
                variables_set                                 , &
                self%face_normal_vectors     (:,face_index)   , &
                self%face_tangential1_vectors(:,face_index)   , &
                self%face_tangential2_vectors(:,face_index)   , &
                self%face_to_cell_indexes    (:,face_index)   , &
                self%num_local_cells                          , &
                num_variables                                 , &
                compute_rotate_variables_function             , &
                compute_unrotate_variables_function           , &
                boundary_condition_function                     &
            )
        end do
    end subroutine apply_empty_condition

    ! ### Variables ###
    subroutine initialize_variables_rank2(self, variables_set, num_variables)
        class  (cellsystem), intent(inout)              :: self
        real   (real_kind ), intent(inout), allocatable :: variables_set(:,:)
        integer(int_kind  ), intent(in   )              :: num_variables

#ifdef _DEBUG
        call write_debuginfo("In initialize_variables_rank2(), cellsystem.")
#endif

        if(.not. self%read_cellsystem) call call_error("'read' subroutine is not called. You should call with following steps: first you call 'read' subroutine, next you initialze variables with 'initialze' subroutine. Please check your cord.")

        if(allocated(variables_set))then
            call call_error("Array variables_set is allocated. But you call 'initialize' subroutine.")
        end if
        allocate(variables_set(1:num_variables, 1:self%num_cells))
        variables_set(:,:) = 0.d0
    end subroutine initialize_variables_rank2

    subroutine initialize_variables_rank1(self, variables_set)
        class  (cellsystem), intent(inout)              :: self
        real   (real_kind ), intent(inout), allocatable :: variables_set(:)

#ifdef _DEBUG
        call write_debuginfo("In initialize_variables_rank1(), cellsystem.")
#endif

        if(.not. self%read_cellsystem) call call_error("'read' subroutine is not called. You should call with following steps: first you call 'read' subroutine, next you initialze variables with 'initialze' subroutine. Please check your cord.")

        if(allocated(variables_set))then
            call call_error("Array variables_set is allocated. But you call 'initialize' subroutine.")
        end if
        allocate(variables_set(1:self%num_cells))
        variables_set(:) = 0.d0
    end subroutine initialize_variables_rank1

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
        class  (eos                     ), intent(in   ) :: an_eos
        real   (real_kind               ), intent(in   ) :: conservative_variables_set(:,:)
        real   (real_kind               ), intent(inout) :: primitive_variables_set(:,:)
        integer(int_kind                ), intent(in   ) :: num_primitive_variables

        procedure(conservative_to_primitive_function_interface) :: conservative_to_primitive_function

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In conservative_to_primitive_variables_all(), cellsystem.")
#endif

!$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            primitive_variables_set(:,i) = conservative_to_primitive_function(an_eos, conservative_variables_set(:,i), num_primitive_variables)
        end do
    end subroutine conservative_to_primitive_variables_all

    subroutine operate_cellwise_rank2(self, variables_set, operator_subroutine)
        class  (cellsystem  ), intent(inout) :: self
        real   (real_kind   ), intent(inout) :: variables_set(:,:)

        procedure(operator_subroutine_interface) :: operator_subroutine

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In operate_cellwise_rank2(), cellsystem.")
#endif

!$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            call operator_subroutine(variables_set(:,i))
        end do
    end subroutine operate_cellwise_rank2

    subroutine operate_cellwise_subarray_rank2(self, primary_variables_set, secondary_variables_set, operator_with_subarray_subroutine)
        class  (cellsystem  ), intent(inout) :: self
        real   (real_kind   ), intent(inout) :: primary_variables_set  (:,:)
        real   (real_kind   ), intent(in   ) :: secondary_variables_set(:,:)

        procedure(operator_with_subarray_subroutine_interface) :: operator_with_subarray_subroutine

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In operate_cellwise_subarray_rank2(), cellsystem.")
#endif

!$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            call operator_with_subarray_subroutine(primary_variables_set(:,i), secondary_variables_set(:,i))
        end do
    end subroutine operate_cellwise_subarray_rank2

    subroutine smooth_variables(self, variables_set, weight_function)
        class  (cellsystem  ), intent(inout) :: self
        real   (real_kind   ), intent(inout) :: variables_set(:,:)
        procedure(weight_function_interface) :: weight_function

        integer(int_kind ) :: face_index, cell_index, rhc_index, lhc_index
        real   (real_kind) :: smoothed_variables_set(size(variables_set(:,1)), size(variables_set(1,:))), total_weight_set(size(variables_set(1,:)))
        real   (real_kind) :: lhc_w, rhc_w

#ifdef _DEBUG
        call write_debuginfo("In smooth_variables(), cellsystem.")
#endif

!$omp parallel do private(cell_index)
        do cell_index = 1, self%num_cells, 1
            smoothed_variables_set(:,cell_index) = 0.d0
            total_weight_set        (cell_index) = 0.d0
        end do

!$omp parallel do private(face_index, rhc_index, lhc_index, lhc_w, rhc_w)
        do face_index = 1, self%num_faces, 1
            lhc_index = self%face_to_cell_indexes(self%num_local_cells - 0, face_index)
            rhc_index = self%face_to_cell_indexes(self%num_local_cells + 1, face_index)

            associate(                                               &
                x_lhc => self%cell_centor_positions(1:3, lhc_index), &
                x_rhc => self%cell_centor_positions(1:3, rhc_index), &
                v_lhc => variables_set             ( : , lhc_index), &
                v_rhc => variables_set             ( : , rhc_index)  &
            )
                lhc_w = weight_function(x_lhc, x_rhc, v_lhc, v_rhc)
                if(self%is_real_cell(rhc_index))then
                    rhc_w = weight_function(x_lhc, x_rhc, v_lhc, v_rhc)
                else
                    rhc_w = 0.d0
                endif

                smoothed_variables_set(:, lhc_index) = smoothed_variables_set(:, lhc_index) + lhc_w * v_lhc
                smoothed_variables_set(:, rhc_index) = smoothed_variables_set(:, rhc_index) + rhc_w * v_rhc

                total_weight_set(lhc_index) = total_weight_set(lhc_index) + lhc_w
                total_weight_set(rhc_index) = total_weight_set(rhc_index) + rhc_w
            end associate
        end do

!$omp parallel do private(cell_index)
        do cell_index = 1, self%num_cells, 1
            if(total_weight_set(cell_index) > 0.d0)then
                variables_set(:,cell_index) = smoothed_variables_set(:,cell_index) / total_weight_set(cell_index)
            else
                variables_set(:,cell_index) = 0.d0
            end if
        end do
    end subroutine smooth_variables

    subroutine substitute_rank2(self, variables, val)
        class(cellsystem  ), intent(inout) :: self
        real (real_kind   ), intent(inout) :: variables(:,:)
        real (real_kind   ), intent(in   ) :: val

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In substitute_rank2(), cellsystem.")
#endif
!$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            variables(:,i) = val
        end do
    end subroutine substitute_rank2

    subroutine substitute_zeros_rank2(self, variables)
        class(cellsystem  ), intent(inout) :: self
        real (real_kind   ), intent(inout) :: variables(:,:)

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In substitute_zeros_rank2(), cellsystem.")
#endif
        call self%substitute_rank2(variables, 0._real_kind)
    end subroutine substitute_zeros_rank2

    subroutine substitute_rank1(self, variables, val)
        class(cellsystem  ), intent(inout) :: self
        real (real_kind   ), intent(inout) :: variables(:)
        real (real_kind   ), intent(in   ) :: val

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In substitute_rank1(), cellsystem.")
#endif
!$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            variables(i) = val
        end do
    end subroutine substitute_rank1

    subroutine substitute_zeros_rank1(self, variables)
        class(cellsystem  ), intent(inout) :: self
        real (real_kind   ), intent(inout) :: variables(:)

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In substitute_zeros_rank1(), cellsystem.")
#endif
        call self%substitute_rank1(variables, 0._real_kind)
    end subroutine substitute_zeros_rank1

    ! ### Gradient Calculator ###
    subroutine initialize_gradient_calculator(self, a_gradient_calculator, config, num_conservative_variables)
        class  (cellsystem         ), intent(inout) :: self
        class  (gradient_calculator), intent(inout) :: a_gradient_calculator
        class  (configuration      ), intent(inout) :: config
        integer(int_kind           ), intent(in   ) :: num_conservative_variables

#ifdef _DEBUG
        call write_debuginfo("In initialize_gradient_calculator(), cellsystem.")
#endif
        call a_gradient_calculator%initialize(config)
    end subroutine initialize_gradient_calculator

    subroutine compute_gradient_rank2(self, a_gradient_calculator, variables_set, gradient_variables_set, num_variables)
        class  (cellsystem         ), intent(inout) :: self
        class  (gradient_calculator), intent(inout) :: a_gradient_calculator
        real   (real_kind          ), intent(in   ) :: variables_set            (:,:)
        real   (real_kind          ), intent(inout) :: gradient_variables_set   (:,:)
        integer(int_kind           ), intent(in   ) :: num_variables

        integer(int_kind ) :: face_index, rhc_index, lhc_index, var_index

#ifdef _DEBUG
        call write_debuginfo("In compute_gradient_rank2(), cellsystem.")
#endif

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
    end subroutine compute_gradient_rank2

    subroutine compute_gradient_rank1(self, a_gradient_calculator, variable_set, gradient_variable_set)
        class(cellsystem         ), intent(inout) :: self
        class(gradient_calculator), intent(inout) :: a_gradient_calculator
        real (real_kind          ), intent(in   ) :: variable_set            (:)
        real (real_kind          ), intent(inout) :: gradient_variable_set   (:,:)

        integer(int_kind) :: face_index, rhc_index, lhc_index

#ifdef _DEBUG
        call write_debuginfo("In compute_gradient_rank1(), cellsystem.")
#endif

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
                                                 + (1.d0 / self%cell_volumes(lhc_index)) * residual
            end associate
        end do
    end subroutine compute_gradient_rank1

    ! ### Divergence Calculator ###
    subroutine initialize_interpolator(self, a_interpolator, config, num_conservative_variables)
        class  (cellsystem   ), intent(inout) :: self
        class  (interpolator ), intent(inout) :: a_interpolator
        class  (configuration), intent(inout) :: config
        integer(int_kind     ), intent(in   ) :: num_conservative_variables

#ifdef _DEBUG
        call write_debuginfo("In initialize_interpolator(), cellsystem.")
#endif
        call a_interpolator%initialize(config)
    end subroutine initialize_interpolator

    subroutine compute_divergence_rank2(self, a_interpolator, variables_set, divergence_variables_set, num_variables)
        class  (cellsystem  ), intent(inout) :: self
        class  (interpolator), intent(inout) :: a_interpolator
        real   (real_kind   ), intent(in   ) :: variables_set           (:,:)
        real   (real_kind   ), intent(inout) :: divergence_variables_set(:,:)
        integer(int_kind    ), intent(in   ) :: num_variables

        integer(int_kind ) :: face_index, rhc_index, lhc_index, var_index, num_divergence_variables
        real   (real_kind) :: face_variables(3)

#ifdef _DEBUG
        call write_debuginfo("In compute_divergence_rank2(), cellsystem.")
#endif

        num_divergence_variables = num_variables/3

!$omp parallel do private(face_index, rhc_index, lhc_index, var_index, face_variables)
        do face_index = 1, self%num_faces, 1
            lhc_index = self%face_to_cell_indexes(self%num_local_cells - 0, face_index)
            rhc_index = self%face_to_cell_indexes(self%num_local_cells + 1, face_index)

            do var_index = 1, num_divergence_variables, 1
                associate(                                &
                    vec_start_index => 3*(var_index-1)+1, &
                    vec_end_index   => 3*(var_index-1)+3  &
                )
                    face_variables(:) = a_interpolator%interpolate_face_variables(                     &
                       variables_set                      (vec_start_index:vec_end_index, :         ), &
                       self%face_to_cell_indexes          (:                            , face_index), &
                       self%cell_centor_positions         (1:3                          , :         ), &
                       self%face_positions                (1:3                          , face_index), &
                       self%num_local_cells                                                          , &
                       3                                                                               &
                    )
                    divergence_variables_set(var_index, lhc_index) = divergence_variables_set(var_index, lhc_index) &
                                                        + (1.d0 / self%cell_volumes(lhc_index)) * vector_multiply(face_variables, self%face_normal_vectors(:, face_index) * self%face_areas(face_index))
                    divergence_variables_set(var_index, rhc_index) = divergence_variables_set(var_index, rhc_index) &
                                                        - (1.d0 / self%cell_volumes(rhc_index)) * vector_multiply(face_variables, self%face_normal_vectors(:, face_index) * self%face_areas(face_index))
                end associate
            end do
        end do
    end subroutine compute_divergence_rank2

    subroutine compute_divergence_rank1(self, a_interpolator, variable_set, divergence_variable_set)
        class(cellsystem   ), intent(inout) :: self
        class(interpolator ), intent(inout) :: a_interpolator
        real (real_kind    ), intent(in   ) :: variable_set              (:,:)
        real (real_kind    ), intent(inout) :: divergence_variable_set   (:)

        integer(int_kind ) :: face_index, rhc_index, lhc_index
        real   (real_kind) :: face_variables(3)

#ifdef _DEBUG
        call write_debuginfo("In compute_divergence_rank1(), cellsystem.")
#endif

!$omp parallel do private(face_index, rhc_index, lhc_index, face_variables)
        do face_index = 1, self%num_faces, 1
            lhc_index = self%face_to_cell_indexes(self%num_local_cells - 0, face_index)
            rhc_index = self%face_to_cell_indexes(self%num_local_cells + 1, face_index)

            face_variables(:) = a_interpolator%interpolate_face_variables(variable_set(:,:), self%face_to_cell_indexes(:, face_index), self%cell_centor_positions(:,:), self%face_positions(:,face_index), self%num_local_cells, 3)

            divergence_variable_set(lhc_index) = divergence_variable_set(lhc_index)               &
                                               + (1.d0 / self%cell_volumes(lhc_index)) * vector_multiply(face_variables(1:3), self%face_normal_vectors(1:3, face_index) * self%face_areas(face_index))
            divergence_variable_set(rhc_index) = divergence_variable_set(rhc_index)               &
                                               - (1.d0 / self%cell_volumes(rhc_index)) * vector_multiply(face_variables(1:3), self%face_normal_vectors(1:3, face_index) * self%face_areas(face_index))
        end do
    end subroutine compute_divergence_rank1

    subroutine compute_divergence_godunov_rank2(self, a_reconstructor, a_riemann_solver, an_eos,                                         &
                                             primitive_variables_set, residual_set, num_conservative_variables, num_primitive_variables, &
                                             primitive_to_conservative_function, element_function)

        class  (cellsystem    ), intent(in   ) :: self
        class  (reconstructor ), intent(in   ) :: a_reconstructor
        class  (riemann_solver), intent(in   ) :: a_riemann_solver
        class  (eos           ), intent(in   ) :: an_eos
        real   (real_kind     ), intent(in   ) :: primitive_variables_set (:,:)
        real   (real_kind     ), intent(inout) :: residual_set            (:,:)
        integer(int_kind      ), intent(in   ) :: num_conservative_variables
        integer(int_kind      ), intent(in   ) :: num_primitive_variables

        procedure(primitive_to_conservative_function_interface   ) :: primitive_to_conservative_function
        procedure(element_provided_with_riemann_function_interface) :: element_function

        integer(int_kind ) :: i
        integer(int_kind ) :: lhc_index
        integer(int_kind ) :: rhc_index
        real   (real_kind) :: reconstructed_primitive_variables_lhc   (num_primitive_variables)
        real   (real_kind) :: reconstructed_primitive_variables_rhc   (num_primitive_variables)
        real   (real_kind) :: element                                 (num_conservative_variables, 1:2)

#ifdef _DEBUG
        call write_debuginfo("In compute_divergence_godunov_rank2 (), cellsystem.")
#endif

!$omp parallel do private(i,lhc_index,rhc_index,reconstructed_primitive_variables_lhc,reconstructed_primitive_variables_rhc,element)
        do i = 1, self%num_faces, 1
            lhc_index = self%face_to_cell_indexes(self%num_local_cells+0, i)
            rhc_index = self%face_to_cell_indexes(self%num_local_cells+1, i)

            reconstructed_primitive_variables_lhc(:) = a_reconstructor%reconstruct_lhc(                              &
                primitive_variables_set, self%face_to_cell_indexes, self%cell_centor_positions, self%face_positions, &
                i, self%num_local_cells, num_primitive_variables                                                     &
            )
            reconstructed_primitive_variables_rhc(:) = a_reconstructor%reconstruct_rhc(                              &
                primitive_variables_set, self%face_to_cell_indexes, self%cell_centor_positions, self%face_positions, &
                i, self%num_local_cells, num_primitive_variables                                                     &
            )

            element(:,:) = element_function(&
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

            residual_set(:,lhc_index) = residual_set(:,lhc_index) + element(:, 1)
            residual_set(:,rhc_index) = residual_set(:,rhc_index) + element(:, 2)
        end do
    end subroutine compute_divergence_godunov_rank2

    subroutine compute_divergence_godunov_facegrad_rank2(self, a_reconstructor, a_riemann_solver, an_eos, a_face_gradient_interpolator, a_face_gradient_calculator,  &
                                        primitive_variables_set, gradient_primitive_variables_set, residual_set, num_conservative_variables, num_primitive_variables, &
                                        primitive_to_conservative_function, element_function)

        class  (cellsystem                ), intent(in   ) :: self
        class  (reconstructor             ), intent(in   ) :: a_reconstructor
        class  (riemann_solver            ), intent(in   ) :: a_riemann_solver
        class  (eos                       ), intent(in   ) :: an_eos
        class  (face_gradient_interpolator), intent(in   ) :: a_face_gradient_interpolator
        class  (face_gradient_calculator  ), intent(in   ) :: a_face_gradient_calculator
        real   (real_kind                 ), intent(in   ) :: primitive_variables_set          (:,:)
        real   (real_kind                 ), intent(in   ) :: gradient_primitive_variables_set (:,:)
        real   (real_kind                 ), intent(inout) :: residual_set                     (:,:)
        integer(int_kind                  ), intent(in   ) :: num_conservative_variables
        integer(int_kind                  ), intent(in   ) :: num_primitive_variables

        procedure(primitive_to_conservative_function_interface            ) :: primitive_to_conservative_function
        procedure(element_provided_with_riemann_facegrad_function_interface) :: element_function

        integer(int_kind ) :: i
        real   (real_kind) :: reconstructed_primitive_variables_lhc   (num_primitive_variables)
        real   (real_kind) :: reconstructed_primitive_variables_rhc   (num_primitive_variables)
        real   (real_kind) :: face_gradient_primitive_variables       (num_primitive_variables*3)
        real   (real_kind) :: element                                 (num_conservative_variables, 1:2)

#ifdef _DEBUG
        call write_debuginfo("In compute_divergence_godunov_facegrad_rank2(), cellsystem.")
#endif

!$omp parallel do private(i,reconstructed_primitive_variables_lhc,reconstructed_primitive_variables_rhc,face_gradient_primitive_variables,element)
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

                if(.not. self%is_real_cell(rhc_index))then ! It's boundary face!
                    face_gradient_primitive_variables(:) = a_face_gradient_calculator%compute(            &
                        primitive_variables_set   (:,lhc_index), primitive_variables_set   (:,rhc_index), &
                        self%cell_centor_positions(:,lhc_index), self%cell_centor_positions(:,rhc_index), &
                        self%face_normal_vectors(:, i),                                                   &
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

                element(:,:) = element_function(        &
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

                residual_set(:,lhc_index) = residual_set(:,lhc_index) + element(:, 1)
                residual_set(:,rhc_index) = residual_set(:,rhc_index) + element(:, 2)
            end associate
        end do
    end subroutine compute_divergence_godunov_facegrad_rank2

    subroutine compute_source_term(self, variables_set, residual_set, num_conservative_variables, compute_source_term_function)
        class  (cellsystem  ), intent(in   ) :: self
        real   (real_kind   ), intent(in   ) :: variables_set(:,:)
        real   (real_kind   ), intent(inout) :: residual_set (:,:)
        integer(int_kind    ), intent(in   ) :: num_conservative_variables

        procedure(compute_source_term_function_interface) :: compute_source_term_function

        integer(int_kind ) :: i
        integer(int_kind ) :: element(num_conservative_variables)

#ifdef _DEBUG
        call write_debuginfo("In compute_source_term (), cellsystem.")
#endif

!$omp parallel do private(i, element)
        do i = 1, self%num_cells, 1
            element(:) = compute_source_term_function(variables_set(:,i), num_conservative_variables)
            residual_set(:,i) = residual_set(:,i) + element(:)
        end do
    end subroutine compute_source_term

    ! ### EoS ###
    subroutine initialize_eos(self, an_eos, config, num_conservative_variables)
        class  (cellsystem    ), intent(inout) :: self
        class  (eos           ), intent(inout) :: an_eos
        class  (configuration ), intent(inout) :: config
        integer(int_kind      ), intent(in   ) :: num_conservative_variables

#ifdef _DEBUG
        call write_debuginfo("In initialize_eos(), cellsystem.")
#endif
        call an_eos%initialize(config)
    end subroutine initialize_eos

    ! ### Riemann solver ###
    subroutine initialize_riemann_solver(self, a_riemann_solver, config, num_conservative_variables)
        class  (cellsystem    ), intent(inout) :: self
        class  (riemann_solver), intent(inout) :: a_riemann_solver
        class  (configuration ), intent(inout) :: config
        integer(int_kind      ), intent(in   ) :: num_conservative_variables

#ifdef _DEBUG
        call write_debuginfo("In initialize_riemann_solver(), cellsystem.")
#endif
        call a_riemann_solver%initialize(config)
    end subroutine initialize_riemann_solver

    ! ### Reconstructor ###
    subroutine initialize_reconstructor(self, a_reconstructor, config, num_conservative_variables, a_reconstructor_generator)
        class  (cellsystem    ), intent(inout) :: self
        class  (reconstructor ), intent(inout) :: a_reconstructor
        class  (configuration ), intent(inout) :: config
        integer(int_kind      ), intent(in   ) :: num_conservative_variables

        class  (reconstructor_generator), optional, intent(inout) :: a_reconstructor_generator

#ifdef _DEBUG
        call write_debuginfo("In initialize_reconstructor(), cellsystem.")
#endif
        call a_reconstructor%initialize(config, a_reconstructor_generator)
    end subroutine initialize_reconstructor

    ! ### Time stepping ###
    subroutine initialize_time_stepping(self, a_time_stepping, config, num_conservative_variables)
        class  (cellsystem   ), intent(inout) :: self
        class  (time_stepping), intent(inout) :: a_time_stepping
        class  (configuration), intent(inout) :: config
        integer(int_kind     ), intent(in   ) :: num_conservative_variables

#ifdef _DEBUG
        call write_debuginfo("In initialize_time_stepping(), cellsystem.")
#endif
        call a_time_stepping%initialize(config, self%num_cells, num_conservative_variables)
    end subroutine initialize_time_stepping

    subroutine compute_next_stage_primitive_rank2(self, a_time_stepping, an_eos, stage_num, conservative_variables_set, primitive_variables_set, residual_set, num_primitive_variables, &
        conservative_to_primitive_function)

        class  (cellsystem         ), intent(inout) :: self
        class  (time_stepping      ), intent(inout) :: a_time_stepping
        class  (eos                ), intent(in   ) :: an_eos
        integer(int_kind           ), intent(in   ) :: stage_num
        real   (real_kind          ), intent(inout) :: conservative_variables_set(:,:)
        real   (real_kind          ), intent(inout) :: primitive_variables_set   (:,:)
        real   (real_kind          ), intent(inout) :: residual_set              (:,:)
        integer(int_kind           ), intent(in   ) :: num_primitive_variables

        procedure(conservative_to_primitive_function_interface) :: conservative_to_primitive_function

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In compute_next_stage_primitive_rank2(), cellsystem.")
#endif

!$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            conservative_variables_set(:,i) = a_time_stepping%compute_next_stage(i, stage_num, self%time_increment, conservative_variables_set(:,i), residual_set(:,i))
            residual_set              (:,i) = 0.d0
            primitive_variables_set   (:,i) = conservative_to_primitive_function(an_eos, conservative_variables_set(:,i), num_primitive_variables)
        end do
    end subroutine compute_next_stage_primitive_rank2

    subroutine compute_next_stage_rank2(self, a_time_stepping, stage_num, conservative_variables_set, residual_set)

        class  (cellsystem         ), intent(inout) :: self
        class  (time_stepping      ), intent(inout) :: a_time_stepping
        integer(int_kind           ), intent(in   ) :: stage_num
        real   (real_kind          ), intent(inout) :: conservative_variables_set(:,:)
        real   (real_kind          ), intent(inout) :: residual_set              (:,:)

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In compute_next_stage_rank2(), cellsystem.")
#endif

!$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            conservative_variables_set(:,i) = a_time_stepping%compute_next_stage(i, stage_num, self%time_increment, conservative_variables_set(:,i), residual_set(:,i))
            residual_set              (:,i) = 0.d0
        end do
    end subroutine compute_next_stage_rank2

    subroutine compute_next_stage_rank1(self, a_time_stepping, stage_num, conservative_variable_set, residual_set)

        class  (cellsystem         ), intent(inout) :: self
        class  (time_stepping      ), intent(inout) :: a_time_stepping
        integer(int_kind           ), intent(in   ) :: stage_num
        real   (real_kind          ), intent(inout) :: conservative_variable_set(:)
        real   (real_kind          ), intent(inout) :: residual_set             (:)

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In compute_next_stage_rank1(), cellsystem.")
#endif

!$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            ! Fortran cannot implicitly convert a 1-dimensional array of 1 elements to a scalar.
            ! In the following, the built-in function sum is used to convert to a scalar.
            conservative_variable_set(i) = sum( a_time_stepping%compute_next_stage(i, stage_num, self%time_increment, [conservative_variable_set(i)], [residual_set(i)]) )
            residual_set             (i) = 0.d0
        end do
    end subroutine compute_next_stage_rank1

    subroutine prepare_time_stepping(    &
        self                      , &
        a_time_stepping           , &
        conservative_variables_set, &
        residual_set                  )
        class(cellsystem   ), intent(inout) :: self
        class(time_stepping), intent(inout) :: a_time_stepping
        real (real_kind    ), intent(inout) :: conservative_variables_set(:,:)
        real (real_kind    ), intent(inout) :: residual_set              (:,:)

        integer(int_kind) :: i

#ifdef _DEBUG
        call write_debuginfo("In prepare_time_stepping(), cellsystem.")
#endif

!$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            call a_time_stepping%prepare_time_stepping(i, conservative_variables_set(:,i), residual_set(:,i))
        end do
    end subroutine

    pure function get_number_of_stages(self, a_time_stepping) result(n)
        class  (cellsystem   ), intent(in) :: self
        class  (time_stepping), intent(in) :: a_time_stepping
        integer(int_kind     )                :: n
        n = a_time_stepping%get_number_of_stages()
    end function get_number_of_stages

    ! ###  Result writer ###
    subroutine initialize_result_writer(self, writer, config, num_conservative_variables)
        class  (cellsystem   ), intent(inout) :: self
        class  (result_writer), intent(inout) :: writer
        class  (configuration), intent(inout) :: config
        integer(int_kind     ), intent(in   ) :: num_conservative_variables

#ifdef _DEBUG
        call write_debuginfo("In initialize_result_writer(), cellsystem.")
#endif
        call writer%initialize(self%num_cells, self%num_points, self%is_real_cell, self%cell_geometries, self%cell_types, config)
    end subroutine initialize_result_writer

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

        integer(boundary_face_type_kind), allocatable :: face_types(:)

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
        call self%initialise_boundary_references(parser%get_number_of_boundary_faces(outflow_face_type     ), &
                                                 parser%get_number_of_boundary_faces(nonslip_wall_face_type), &
                                                 parser%get_number_of_boundary_faces(slip_wall_face_type   ), &
                                                 parser%get_number_of_boundary_faces(symmetric_face_type   ), &
                                                 parser%get_number_of_boundary_faces(empty_face_type       )    )

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

    pure function get_number_of_nonslip_wall_faces(self) result(n)
        class(cellsystem), intent(in) :: self
        integer(int_kind) :: n
        n = self%num_nonslip_wall_faces
    end function get_number_of_nonslip_wall_faces

    pure function get_number_of_slip_wall_faces(self) result(n)
        class(cellsystem), intent(in) :: self
        integer(int_kind) :: n
        n = self%num_slip_wall_faces
    end function get_number_of_slip_wall_faces

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
        integer(boundary_face_type_kind), intent(in   ) :: face_types(:)

        integer(int_kind) :: index, outflow_index, nonslip_wall_index, slip_wall_index, empty_index, symmetric_index

#ifdef _DEBUG
        call write_debuginfo("In assign_boundary(), cellsystem.")
#endif

        outflow_index            = 0
        nonslip_wall_index       = 0
        slip_wall_index          = 0
        symmetric_index          = 0
        empty_index              = 0

        do index = 1, self%num_faces, 1
            if (face_types(index) == outflow_face_type) then
                outflow_index = outflow_index + 1
                self%outflow_face_indexes(outflow_index) = index
            else if (face_types(index) == nonslip_wall_face_type) then
                nonslip_wall_index = nonslip_wall_index + 1
                self%nonslip_wall_face_indexes(nonslip_wall_index) = index
            else if (face_types(index) == slip_wall_face_type) then
                slip_wall_index = slip_wall_index + 1
                self%slip_wall_face_indexes(slip_wall_index) = index
            else if (face_types(index) == symmetric_face_type) then
                symmetric_index = symmetric_index + 1
                self%symmetric_face_indexes(symmetric_index) = index
            else if (face_types(index) == empty_face_type) then
                empty_index = empty_index + 1
                self%empty_face_indexes(empty_index) = index
            else if (face_types(index) == boundary_face_type) then
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

    subroutine initialise_boundary_references(self, num_outflow_faces, num_nonslip_wall_faces, num_slip_wall_faces, num_symetric_faces, num_empty_faces)
        class(cellsystem), intent(inout) :: self
        integer(int_kind), intent(in) :: num_outflow_faces
        integer(int_kind), intent(in) :: num_nonslip_wall_faces
        integer(int_kind), intent(in) :: num_slip_wall_faces
        integer(int_kind), intent(in) :: num_symetric_faces
        integer(int_kind), intent(in) :: num_empty_faces

#ifdef _DEBUG
        call write_debuginfo("In initialise_boundary_references(), cellsystem.")
#endif

        if(num_outflow_faces < 0)then
            call call_error("Number of outflow BC faces must be set over zero.")
        end if
        self%num_outflow_faces = num_outflow_faces

        if(num_nonslip_wall_faces < 0)then
            call call_error("Number of nonslip wall BC faces must be set over zero.")
        end if
        self%num_nonslip_wall_faces = num_nonslip_wall_faces

        if(num_slip_wall_faces < 0)then
            call call_error("Number of slip wall BC faces must be set over zero.")
        end if
        self%num_slip_wall_faces = num_slip_wall_faces

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

        if(allocated(self%nonslip_wall_face_indexes))then
            call call_error("Array nonslip_wall_face_indexes is already allocated. But you call the initialiser for boundary_reference.")
        end if
        allocate(self%nonslip_wall_face_indexes(self%num_nonslip_wall_faces))

        if(allocated(self%slip_wall_face_indexes))then
            call call_error("Array slip_wall_face_indexes is already allocated. But you call the initialiser for boundary_reference.")
        end if
        allocate(self%slip_wall_face_indexes(self%num_slip_wall_faces))

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

        if(.not. allocated(self%nonslip_wall_face_indexes))then
            call call_error("Array nonslip_wall_face_indexes is not allocated. But you call the finalizer for boundary_reference module.")
        end if
        deallocate(self%nonslip_wall_face_indexes)

        if(.not. allocated(self%slip_wall_face_indexes))then
            call call_error("Array slip_wall_face_indexes is not allocated. But you call the finalizer for boundary_reference module.")
        end if
        deallocate(self%slip_wall_face_indexes)

        if(.not. allocated(self%symmetric_face_indexes))then
            call call_error("Array symmetric_face_indexes is not allocated. But you call the finalizer for boundary_reference module.")
        end if
        deallocate(self%symmetric_face_indexes)

        if(.not. allocated(self%empty_face_indexes))then
            call call_error("Array empty_face_indexes is not allocated. But you call the finalizer for boundary_reference module.")
        end if
        deallocate(self%empty_face_indexes)

        self%num_local_cells          = 0
        self%num_outflow_faces        = 0
        self%num_nonslip_wall_faces   = 0
        self%num_slip_wall_faces      = 0
        self%num_symmetric_faces      = 0
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
