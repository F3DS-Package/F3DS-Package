module class_cellsystem
    use typedef_module
    use stdio_module
    use vector_module
    use matrix_module
    use face_type_module
    use class_point_id_list
    use abstract_grid_parser
    use abstract_result_writer
    use abstract_configuration
    use abstract_time_stepping
    use abstract_eos
    use abstract_reconstruction
    use abstract_riemann_solver
    use abstract_gradient_scheme

    implicit none

    private

    ! Data structure and rule:
    ! (local cell index)   -2        -1         0         1         2         3
    !                       |         |         |         |         |         |
    !
    !                  *---------*---------*---------*---------*---------*---------*
    !                  |         |         |         | n       |         |         | "n" is a normal vector.
    !                  |    o    |    o    |    o    x--> o    |    o    |    o    | Normal vector must be oriented to the right-hand cell.
    ! (cell index) ->  |    1    |    2    |    3    |    4    |    5    |    6    |
    !                  *---------*---------*---------*---------*---------*---------*
    !
    !   (inner cells) ______|_________|_________|    |    |_________|_________|________ (inner or ghost cells)
    !                                                |             *If fase is boundary, ghost cells are stored right-side.
    !                                      (boundary/inner face)
    type, public :: cellsystem
        private

        ! Termination condition
        real(real_kind) :: end_time

        ! Time
        real(real_kind) :: num_steps
        real(real_kind) :: time
        real(real_kind) :: time_increment

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
        integer(int_kind) :: num_slipwall_faces
        integer(int_kind) :: num_symmetric_faces
        integer(int_kind) :: num_empty_face_indexes

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

        ! ### Boundary References ###

        ! Elements are stored boundary face index.
        ! Elm. 1) 1 : num_{boundary condition}_faces
        integer(int_kind), public, allocatable :: outflow_face_indexes  (:)
        integer(int_kind), public, allocatable :: slipwall_face_indexes (:)
        integer(int_kind), public, allocatable :: symmetric_face_indexes(:)
        integer(int_kind), public, allocatable :: empty_face_indexes    (:)

        ! ### Status ###
        logical :: read_cellsystem = .false.
        logical :: initialized_variables = .false.

        contains
        ! ### Common ###
        procedure, public, pass(self) :: initialize  => result_writer_initialize  , &
                                                        time_stepping_initialize  , &
                                                        reconstruction_initialize , &
                                                        riemann_solver_initialize , &
                                                        eos_initialize            , &
                                                        gradient_scheme_initialize, &
                                                        variables_initialize
        procedure, public, pass(self) :: open_file   => result_writer_open_file
        procedure, public, pass(self) :: close_file  => result_writer_close_file
        procedure, public, pass(self) :: is_writable => result_writer_is_writable

        ! ### Gradient Scheme ###
        procedure, public, pass(self) :: compute_gradient => compute_gradient_2darray, &
                                                             compute_gradient_1darray

        ! ### Flux ###
        procedure, public, pass(self) :: integrate_flux

        ! ### Time Stepping ###
        procedure, public, pass(self) :: compute_next_state
        procedure, public, pass(self) :: prepare_stepping
        procedure, public, pass(self) :: get_number_of_stage

        ! ### Result Writer ###
        procedure, public, pass(self) :: write_scolar
        procedure, public, pass(self) :: write_vector
        procedure, public, pass(self) :: get_filename

        ! ### Cellsystem ###
        procedure, public, pass(self) :: read

        ! ### Getter ###
        procedure, public, pass(self) :: get_number_of_faces
        procedure, public, pass(self) :: get_number_of_local_cells
        procedure, public, pass(self) :: get_number_of_points
        procedure, public, pass(self) :: get_number_of_cells
        procedure, public, pass(self) :: get_number_of_outflow_faces
        procedure, public, pass(self) :: get_number_of_slipwall_faces
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
    end type

    contains

    ! ### Variables ###
    subroutine variables_initialize(self, conservative_variables_set, primitive_variables_set, residual_set, num_conservative_variables, num_primitive_variables)
        class  (cellsystem), intent(inout)              :: self
        real   (real_kind ), intent(inout), allocatable :: conservative_variables_set(:,:)
        real   (real_kind ), intent(inout), allocatable :: primitive_variables_set   (:,:)
        real   (real_kind ), intent(inout), allocatable :: residual_set              (:,:)
        integer(int_kind  ), intent(in   )              :: num_conservative_variables
        integer(int_kind  ), intent(in   )              :: num_primitive_variables

        if(.not. read_cellsystem) call call_error("'read' subroutine is not called. You should call with following steps: first you call 'read' subroutine, next you initialze variables with 'initialze' subroutine. Please check your cord.")

        if(allocated(conservative_variables_set))then
            call call_error("Array conservative_variables_set is allocated. But you call 'initialize' subroutine.")
        end if
        allocate(conservative_variables_set(1:num_conservative_variables, 1:self%num_cells))

        if(allocated(primitive_variables_set))then
            call call_error("Array primitive_variables_set is allocated. But you call 'initialize' subroutine.")
        end if
        allocate(primitive_variables_set(1:num_primitive_variables, 1:self%num_cells))

        if(allocated(residual_set))then
            call call_error("Array residual_set is allocated. But you call 'initialize' subroutine.")
        end if
        allocate(residual_set(1:num_conservative_variables, 1:self%num_cells))

        self%initialized_variables = .true.
    end subroutine variables_initialize

    ! ### Gradient Schemes ###
    subroutine gradient_scheme_initialize(self, a_gradient_scheme, config)
        class  (cellsystem     ), intent(inout) :: self
        class  (gradient_scheme), intent(inout) :: a_gradient_scheme
        class  (configuration  ), intent(inout) :: config

        call a_gradient_scheme%initialize(config)
    end subroutine gradient_scheme_initialize

    subroutine compute_gradient_2darray(self, a_gradient_scheme, variables_set, gradient_variables_set)
        class(cellsystem     ), intent(inout) :: self
        class(gradient_scheme), intent(inout) :: a_gradient_scheme
        real (real_kind      ), intent(in   ) :: variables_set            (:,:)
        real (real_kind      ), intent(in   ) :: gradient_variables_set   (:,:)

        call a_gradient_scheme%compute_gradient( &
            variables_set                      , &
            gradient_variables_set             , &
            self%cell_centor_positions         , &
            self%cell_volumes                  , &
            self%face_to_cell_indexes          , &
            self%face_normal_vectors           , &
            self%face_positions                , &
            self%face_areas                    , &
            self%num_faces                       &
        )
    end subroutine compute_gradient_2darray

    subroutine compute_gradient_1darray(self, a_gradient_scheme, variable_set, gradient_variable_set)
        class(cellsystem     ), intent(inout) :: self
        class(gradient_scheme), intent(inout) :: a_gradient_scheme
        real (real_kind      ), intent(in   ) :: variable_set            (:)
        real (real_kind      ), intent(in   ) :: gradient_variable_set   (:)

        call a_gradient_scheme%compute_gradient_1darray( &
            variable_set                               , &
            gradient_variable_set                      , &
            self%cell_centor_positions                 , &
            self%cell_volumes                          , &
            self%face_to_cell_indexes                  , &
            self%face_normal_vectors                   , &
            self%face_positions                        , &
            self%face_areas                            , &
            self%num_faces                               &
        )
    end subroutine compute_gradient_1darray

    ! ### EoS ###
    subroutine eos_initialize(self, an_eos, config)
        class  (cellsystem    ), intent(inout) :: self
        class  (eos           ), intent(inout) :: an_eos
        class  (configuration ), intent(inout) :: config

        call an_eos%initialize(config)
    end subroutine eos_initialize

    ! ### Flux ###
    subroutine integrate_flux(self, a_reconstruction, a_riemann_solver, an_eos,                                           &
                              primitive_variables_set, residual_set, num_conservative_variables, num_primitive_variables, &
                              primitive_to_conservative_function, compute_rotate_matrix_primitive_function,               &
                              compute_unrotate_matrix_conservative_function, model_function                                 )

        class  (cellsystem    ), intent(in   ) :: self
        class  (reconstruction), intent(in   ) :: a_reconstruction
        class  (riemann_solver), intent(in   ) :: a_riemann_solver
        class  (eos           ), intent(in   ) :: an_eos
        real   (real_kind     ), intent(in   ) :: primitive_variables_set (:,:)
        real   (real_kind     ), intent(inout) :: residual_set            (:,:)
        integer(int_kind      ), intent(in   ) :: num_conservative_variables
        integer(int_kind      ), intent(in   ) :: num_primitive_variables

        interface
            pure function primitive_to_conservative_function(primitive, an_eos) result(conservative)
                use typedef_module
                use abstract_eos
                real (real_kind), intent(in)  :: primitive   (:)
                class(eos      ), intent(in)  :: an_eos
                real (real_kind), allocatable :: conservative(:)
            end function primitive_to_conservative_function

            pure function compute_rotate_matrix_primitive_function( &
                face_normal_vector       , &
                face_tangential1_vector  , &
                face_tangential2_vector  ,   ) result(matrix)

                use typedef_module
                real   (real_kind     ), intent(in) :: face_normal_vector      (3)
                real   (real_kind     ), intent(in) :: face_tangential1_vector (3)
                real   (real_kind     ), intent(in) :: face_tangential2_vector (3)
            end function compute_rotate_matrix_primitive_function

            pure function compute_unrotate_matrix_conservative_function( &
                face_normal_vector       , &
                face_tangential1_vector  , &
                face_tangential2_vector  ,   ) result(matrix)

                use typedef_module
                real   (real_kind     ), intent(in) :: face_normal_vector      (3)
                real   (real_kind     ), intent(in) :: face_tangential1_vector (3)
                real   (real_kind     ), intent(in) :: face_tangential2_vector (3)
            end function compute_unrotate_matrix_conservative_function

            pure function model_function(                 &
                primitive_variables_lhc                 , &
                primitive_variables_rhc                 , &
                reconstructed_primitive_variables_lhc   , &
                reconstructed_primitive_variables_rhc   , &
                reconstructed_conservative_variables_lhc, &
                reconstructed_conservative_variables_rhc, &
                num_conservative_values                 , &
                an_eos                                  , &
                an_riemann_solver                        ) result(flux)

                use typedef_module
                use abstract_eos
                use abstract_riemann_solver

                real   (real_kind     ), intent(in) :: primitive_variables_lhc                  (:)
                real   (real_kind     ), intent(in) :: primitive_variables_rhc                  (:)
                real   (real_kind     ), intent(in) :: reconstructed_primitive_variables_lhc    (:)
                real   (real_kind     ), intent(in) :: reconstructed_primitive_variables_rhc    (:)
                real   (real_kind     ), intent(in) :: reconstructed_conservative_variables_lhc (:)
                real   (real_kind     ), intent(in) :: reconstructed_conservative_variables_rhc (:)
                integer(int_kind      ), intent(in) :: num_conservative_values
                class  (eos           ), intent(in) :: an_eos
                class  (riemann_solver), intent(in) :: an_riemann_solver

                real   (real_kind)                  :: flux(num_conservative_values)
            end function model_function
        end interface

        integer(int_kind ) :: i
        integer(int_kind ) :: lhc_index
        integer(int_kind ) :: rhc_index
        real   (real_kind) :: rotate_matrix                        (num_primitive_variables, num_primitive_variables)
        real   (real_kind) :: reconstructed_primitive_variables_lhc(num_primitive_variables)
        real   (real_kind) :: reconstructed_primitive_variables_rhc(num_primitive_variables)
        real   (real_kind) :: rotated_primitive_variables_lhc      (num_primitive_variables)
        real   (real_kind) :: rotated_primitive_variables_rhc      (num_primitive_variables)
        real   (real_kind) :: rotated_conservative_variables_lhc   (num_conservative_variables)
        real   (real_kind) :: rotated_conservative_variables_rhc   (num_conservative_variables)
        real   (real_kind) :: local_coordinate_flux                (num_conservative_variables)
        real   (real_kind) :: global_coordinate_flux               (num_conservative_variables)

!$omp parallel do private(i,lhc_index,rhc_index,rotate_matrix,reconstructed_primitive_variables_lhc,reconstructed_primitive_variables_rhc,rotated_primitive_variables_lhc,rotated_conservative_variables_rhc,local_coordinate_flux,global_coordinate_flux)
        do i = 1, self%num_faces, 1
            lhc_index = face_to_cell_index(self%num_local_cells+0, i)
            rhc_index = face_to_cell_index(self%num_local_cells+1, i)

            reconstructed_primitive_variables_lhc(:) = a_reconstruction%reconstruct_lhc(                                      &
                i, primitive_variables_set, self%face_to_cell_indexes, self%cell_centor_positions, self%face_centor_positions &
            )
            reconstructed_primitive_variables_rhc(:) = a_reconstruction%reconstruct_rhc(                                      &
                i, primitive_variables_set, self%face_to_cell_indexes, self%cell_centor_positions, self%face_centor_positions &
            )

            rotate_matrix(:,:) = compute_rotate_matrix_primitive_function(                                            &
                self%face_normal_vectors(:,i), self%face_tangential1_vectors(:,i), self%face_tangential2_vectors(:,i) &
            )

            rotated_primitive_variables_lhc(:) = matrix_multiply( &
                rotate_matrix,                                    &
                reconstructed_primitive_variables_lhc             &
            )
            rotated_primitive_variables_rhc(:) = matrix_multiply( &
                rotate_matrix,                                    &
                reconstructed_primitive_variables_rhc,            &
            )

            rotated_conservative_variables_lhc(:) = primitive_to_conservative_function(rotated_primitive_variables_lhc, an_eos)
            rotated_conservative_variables_rhc(:) = primitive_to_conservative_function(rotated_primitive_variables_rhc, an_eos)

            local_coordinate_flux(:) = model_function(            &
                primitive_variables_set           (:, lhc_index), &
                primitive_variables_set           (:, rhc_index), &
                rotated_primitive_variables_lhc   (:)           , &
                rotated_primitive_variables_rhc   (:)           , &
                rotated_conservative_variables_lhc(:)           , &
                rotated_conservative_variables_rhc(:)           , &
                num_conservative_variables                      , &
                an_eos                                          , &
                a_riemann_solver                                  &
            )

            global_coordinate_flux(:) = matrix_multiply(                                                                  &
                compute_unrotate_matrix_conservative_function(                                                            &
                    self%face_normal_vectors(:,i), self%face_tangential1_vectors(:,i), self%face_tangential2_vectors(:,i) &
                )                                                                                                         &
                local_coordinate_flux                                                                                     &
            )

            residual_set(:,lhc_index) = residual_set(:,lhc_index) - (1.d0 / self%cell_volumes(lhc_index)) * global_coordinate_flux(:) * self%face_areas(i)
            residual_set(:,rhc_index) = residual_set(:,rhc_index) + (1.d0 / self%cell_volumes(rhc_index)) * global_coordinate_flux(:) * self%face_areas(i)
        end do
    end subroutine

    subroutine riemann_solver_initialize(self, a_riemann_solver, config)
        class  (cellsystem    ), intent(inout) :: self
        class  (riemann_solver), intent(inout) :: a_riemann_solver
        class  (configuration ), intent(inout) :: config

        call a_riemann_solver%initialize(config)
    end subroutine time_stepping_initialize

    subroutine reconstruction_initialize(self, a_reconstruction, config)
        class  (cellsystem    ), intent(inout) :: self
        class  (reconstruction), intent(inout) :: a_reconstruction
        class  (configuration ), intent(inout) :: config

        call reconstruction%initialize(config)
    end subroutine time_stepping_initialize

    ! ### Time stepping ###
    subroutine time_stepping_initialize(self, a_time_stepping, config, num_conservative_variables)
        class  (cellsystem   ), intent(inout) :: self
        class  (time_stepping), intent(inout) :: a_time_stepping
        class  (configuration), intent(inout) :: config
        integer(int_kind     ), intent(in   ) :: num_conservative_variables

        call a_time_stepping%initialize(config, self%num_cells, num_conservative_variables)
    end subroutine time_stepping_initialize

    subroutine compute_next_state(self, a_time_stepping, an_eos, state_num, conservative_variables_set, primitive_variables_set, residual_set, &
        conservative_to_primitive_function)

        class  (cellsystem         ), intent(inout) :: self
        class  (time_stepping      ), intent(inout) :: a_time_stepping
        class  (eos                ), intent(in   ) :: an_eos
        integer(int_kind           ), intent(in   ) :: state_num
        real   (real_kind          ), intent(inout) :: conservative_variables_set(:,:)
        real   (real_kind          ), intent(inout) :: primitive_variables_set   (:,:)
        real   (real_kind          ), intent(inout) :: residual_set              (:,:)

        interface
            pure function conservative_to_primitive_function(conservative, an_eos) result(primitive)
                use typedef_module
                use abstract_eos
                real (real_kind), intent(in)  :: conservative(:)
                class(eos      ), intent(in)  :: an_eos
                real (real_kind), allocatable :: primitive   (:)
            end function conservative_to_primitive_function
        end interface

        integer(int_kind) :: i

!$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            call a_time_stepping%compute_next_state(i, state_num, self%time_increment, conservative_variables_set(:,i), residual_set(:,i))
            primitive_variables_set(:,i) = conservative_to_primitive_function(conservative_variables_set(:,i), an_eos)
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

!$omp parallel do private(i)
        do i = 1, self%num_cells, 1
            call a_time_stepping%prepare_stepping(i, conservative_variables_set(:,i), primitive_variables_set(:,i), residual_set(:,i))
        end do
    end subroutine

    pure function get_number_of_stage(self, a_time_stepping) result(n)
        class  (cellsystem   ), intent(inout) :: self
        class  (time_stepping), intent(inout) :: a_time_stepping
        integer(int_kind     )                :: n

        n = a_time_stepping%get_number_of_stage()
    end function get_number_of_stage

    ! ###  Result writer ###
    subroutine result_writer_initialize(self, writer, config)
        class  (cellsystem   ), intent(inout) :: self
        class  (result_writer), intent(inout) :: writer
        class  (configuration), intent(inout) :: config

        call writer%initialize(self%num_cells, self%num_points, self%is_real_cell, self%cell_geometries, self%cell_types, self%end_time, config)
    end subroutine result_writer_initialize

    subroutine result_writer_open_file(self, writer)
        class(cellsystem   ), intent(inout) :: self
        class(result_writer), intent(inout) :: writer

        call writer%open_file(self%time, self%points)
    end subroutine result_writer_open_file

    subroutine result_writer_close_file(self, writer)
        class(cellsystem   ), intent(inout) :: self
        class(result_writer), intent(inout) :: writer

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

        call writer%write_scolar(self%is_real_cell, name, scolar_variables)
    end subroutine write_scolar

    subroutine write_vector(self, writer, name, vector_variables)
        class    (cellsystem   ), intent(inout) :: self
        class    (result_writer), intent(inout) :: writer
        character(len=*        ), intent(in   ) :: name
        real     (real_kind    ), intent(in   ) :: vector_variables(:,:)

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

        ! parse grid file
        call parser%parse(config)

        ! initialize local valiables
        allocate(face_types(parser%get_number_of_faces()))

        ! allocate grid
        call self%initialise_faces              (parser%get_number_of_faces (), parser%get_number_of_ghost_cells())
        call self%initialise_cells              (parser%get_number_of_points(), parser%get_number_of_cells      ())
        call self%initialise_boundary_references(parser%get_number_of_boundary_faces(outflow_face_type  ), &
                                                 parser%get_number_of_boundary_faces(slipwall_face_type ), &
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

        self%read_cellsystem = .true.
    end subroutine read

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

    pure function get_number_of_slipwall_faces(self) result(n)
        class(cellsystem), intent(in) :: self
        integer(int_kind) :: n
        n = self%num_slipwall_faces
    end function get_number_of_slipwall_faces

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
        call self%finalize_cells()
        call self%finalize_faces()
        call self%finalize_boundary_references()
    end subroutine finalize

    ! ### Inner utils ###
    subroutine assign_boundary(self, face_types)
        class  (cellsystem     ), intent(inout) :: self
        integer(kind(face_type)), intent(in   ) :: face_types(:)

        integer(int_kind) :: index, outflow_index, slipwall_index, symmetric_index, empty_index

        outflow_index   = 0
        slipwall_index  = 0
        symmetric_index = 0
        empty_index     = 0

        do index = 1, self%num_faces, 1
            if (face_types(index) == outflow_face_type) then
                outflow_index = outflow_index + 1
                self%outflow_face_indexes(outflow_index) = index
            else if (face_types(index) == slip_wall_face_type) then
                slipwall_index = slipwall_index + 1
                self%slipwall_face_indexes(slipwall_index) = index
            else if (face_types(index) == symmetric_face_type) then
                symmetric_index = symmetric_index + 1
                self%symmetric_face_indexes(symmetric_index) = index
            else if (face_types(index) == empty_face_type) then
                empty_index = empty_index + 1
                self%empty_face_indexes(empty_index) = index
            else
                call call_error("Unknown boundary type found.")
            end if
        end do
    end subroutine assign_boundary

    subroutine initialise_faces(self, num_faces, num_local_cells)
        class(cellsystem), intent(inout) :: self
        integer(int_kind), intent(in) :: num_faces
        integer(int_kind), intent(in) :: num_local_cells

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
        allocate(self%face_to_cell_indexes(-self%num_local_cells + 1:self%num_local_cells, self%num_faces))

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

    subroutine initialise_boundary_references(self, num_outflow_faces, num_slipwall_faces, num_symetric_faces, num_empty_face_indexes)
        class(cellsystem), intent(inout) :: self
        integer(int_kind), intent(in) :: num_outflow_faces
        integer(int_kind), intent(in) :: num_slipwall_faces
        integer(int_kind), intent(in) :: num_symetric_faces
        integer(int_kind), intent(in) :: num_empty_face_indexes

        if(num_outflow_faces < 0)then
            call call_error("Number of outflow BC faces must be set over zero.")
        end if
        self%num_outflow_faces = num_outflow_faces

        if(num_slipwall_faces < 0)then
            call call_error("Number of slipwall BC faces must be set over zero.")
        end if
        self%num_slipwall_faces = num_slipwall_faces

        if(num_symetric_faces < 0)then
            call call_error("Number of symetric BC faces must be set over zero.")
        end if
        self%num_symmetric_faces = num_symetric_faces

        if(num_empty_face_indexes < 0)then
            call call_error("Number of empty BC faces must be set over zero.")
        end if
        self%num_empty_face_indexes = num_empty_face_indexes

        if(allocated(self%outflow_face_indexes))then
            call call_error("Array outflow_face_indexes is already allocated. But you call the initialiser for boundary_reference.")
        end if
        allocate(self%outflow_face_indexes(self%num_outflow_faces))

        if(allocated(self%slipwall_face_indexes))then
            call call_error("Array slipwall_face_indexes is already allocated. But you call the initialiser for boundary_reference.")
        end if
        allocate(self%slipwall_face_indexes(self%num_slipwall_faces))

        if(allocated(self%symmetric_face_indexes))then
            call call_error("Array symmetric_face_indexes is already allocated. But you call the initialiser for boundary_reference.")
        end if
        allocate(self%symmetric_face_indexes(self%num_symmetric_faces))

        if(allocated(self%empty_face_indexes))then
            call call_error("Array empty_face_indexes is already allocated. But you call the initialiser for boundary_reference.")
        end if
        allocate(self%empty_face_indexes(self%num_empty_face_indexes))
    end subroutine initialise_boundary_references

    subroutine finalize_faces(self)
        class(cellsystem), intent(inout) :: self

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

        if(.not. allocated(self%outflow_face_indexes))then
            call call_error("Array outflow_face_indexes is not allocated. But you call the finalizer for boundary_reference module.")
        end if
        deallocate(self%outflow_face_indexes)

        if(.not. allocated(self%slipwall_face_indexes))then
            call call_error("Array slipwall_face_indexes is not allocated. But you call the finalizer for boundary_reference module.")
        end if
        deallocate(self%slipwall_face_indexes)

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
        self%num_slipwall_faces  = 0
        self%num_symmetric_faces = 0
    end subroutine finalize_boundary_references
end module class_cellsystem