program five_eq_model_solver
    ! Utils
    use typedef_module
    ! Cell system
    use faces_module
    use boundary_reference_module
    use cell_geometries_module
    use five_equation_model_variables_module
    ! Scheme
    use second_order_tvd_rk_module
    use weno5_module
    use five_equation_space_model_module
    use five_equation_model_hllc_module
    ! Model
    use mixture_stiffened_eos_module
    ! Grid & initial condition reader
    use class_nlinit_parser
    use class_nlgrid_parser
    ! Result file output
    use penf
    use vtk_fortran, only : vtk_file

    implicit none

    real   (real_kind)  :: time_increment
    integer(int_kind )  :: max_timestep, timestep
    integer(int_kind )  :: index
    ! Grid & initial condition I/O
    type(nlgrid_parser) :: a_grid_parser
    type(nlinit_parser) :: a_init_parser
    ! Result file output
    type(vtk_file)              :: a_vtk_file
    integer  (I4P)              :: vtk_error
    integer  (int_kind )        :: file_output_counter, vtk_index, cell_point_index
    character(9)                :: vtk_filename
    integer  (I4P)              :: n_output_cells, n_output_points, n_cell_points
    real     (R4P), allocatable :: vtk_x(:), vtk_y(:), vtk_z(:)
    real     (R4P), allocatable :: vtk_scolar(:)
    integer  (I1P), allocatable :: vtk_cell_type(:)
    integer  (I4P), allocatable :: vtk_offset(:), vtk_connect(:)

    time_increment = 1.d0
    max_timestep   = 100

    ! parse grid file
    call a_grid_parser%parse("grid.nlgrid")
    ! allocate grid & variable data
    call initialise_faces             (a_grid_parser%get_number_of_faces()    , a_grid_parser%get_number_of_ghost_cells())
    call initialise_cell_geometries   (a_grid_parser%get_number_of_points()   , a_grid_parser%get_number_of_cells())
    call initialise_boundary_reference(a_grid_parser%get_number_of_ghost_cells(), &
                                       a_grid_parser%get_number_of_outflow_faces(), a_grid_parser%get_number_of_slipwall_faces(), a_grid_parser%get_number_of_symmetric_faces())
    call initialise_variables         (a_grid_parser%get_number_of_cells())
    ! get grid data
    call a_grid_parser%get_cells          (cells_centor_position, cells_volume, cells_is_real_cell)
    call a_grid_parser%get_faces          (faces_reference_cell_index, faces_normal_vector, faces_tangential1_vector, faces_tangential2_vector, faces_position, faces_area)
    call a_grid_parser%get_boundaries     (outflow_face_indexs, slipwall_face_indexs, symetric_face_indexs)
    call a_grid_parser%get_cell_geometries(points, cell_geometries)
    ! close file
    call a_grid_parser%close()

    ! parse init file
    call a_init_parser%parse("init.nlinit")
    call a_init_parser%get_conservative_variables_set(conservative_variables_set)
    call a_init_parser%close()
    do index = 1, get_number_of_cells(), 1
        primitive_variables_set(:, index) = conservative_to_primitive(conservative_variables_set(:, index))
    end do

    ! VTK
    file_output_counter = 0
    n_output_cells = 0
    do index = 1, get_number_of_cells(), 1
        if(cells_is_real_cell(index)) n_output_cells = n_output_cells + 1
    end do
    allocate(vtk_scolar   (n_output_cells    ))
    allocate(vtk_cell_type(n_output_cells    ))
    allocate(vtk_offset   (n_output_cells    ))
    allocate(vtk_connect  (n_output_cells * 8))
    vtk_index = 1
    do index = 1, get_number_of_cells(), 1
        if(cells_is_real_cell(index))then
            n_cell_points = cell_geometries(vtk_index)%get_number_of_points()
            do cell_point_index = 1, n_cell_points, 1
                vtk_connect((vtk_index - 1) * 8 + cell_point_index) = cell_geometries(vtk_index)%get_point_id(cell_point_index)
            end do
            vtk_index = vtk_index + 1
        end if
    end do
    n_output_points = get_number_of_points()
    allocate(vtk_x(n_output_points))
    allocate(vtk_y(n_output_points))
    allocate(vtk_z(n_output_points))
    do index = 1, get_number_of_points(), 1
        vtk_x(index) = points(1, index)
        vtk_y(index) = points(2, index)
        vtk_z(index) = points(3, index)
    end do

    call initialize_second_order_tvd_rk(conservative_variables_set)

    ! solver timestepping loop
    do timestep = 0, max_timestep, 1
        print *, "step ", timestep

        call compute_next_state_second_order_tvd_rk(   &
            conservative_variables_set               , &
            primitive_variables_set                  , &
            cells_centor_position                    , & ! <- TODO: We remove arguments in a future. We make the solver module and compute_next_state() routine.
            cells_volume                             , & ! <-
            faces_reference_cell_index               , & ! <-
            faces_normal_vector                      , & ! <-
            faces_tangential1_vector                 , & ! <-
            faces_tangential2_vector                 , & ! <-
            faces_position                           , & ! <-
            faces_area                               , & ! <-
            get_number_of_cells()                    , & ! <-
            get_number_of_faces()                    , & ! <-
            time_increment                           , &
            reconstruct_weno5                        , &
            compute_space_element_five_equation_model, &
            compute_flux_five_equation_model_hllc    , &
            compute_pressure_mixture_stiffened_eos   , &
            compute_soundspeed_mixture_stiffened_eos , &
            primitive_to_conservative                , &
            conservative_to_primitive                  &
        )

        if (.true.) then
            write(vtk_filename, "(i5, a)") file_output_counter, ".vtu"
            print *, "write vtk "//vtk_filename//"..."
            vtk_error = a_vtk_file%initialize                   (format="binary", filename=vtk_filename, mesh_topology="UnstructuredGrid")
            vtk_error = a_vtk_file%xml_writer%write_piece       (np=n_output_points, nc=n_output_cells)
            vtk_error = a_vtk_file%xml_writer%write_geo         (np=n_output_points, nc=n_output_cells, x=vtk_x, y=vtk_y, z=vtk_z)
            vtk_error = a_vtk_file%xml_writer%write_connectivity(nc=n_output_cells, connectivity=vtk_connect, offset=vtk_offset, cell_type=vtk_cell_type)
            vtk_error = a_vtk_file%xml_writer%write_dataarray   (location='node', action='open')
            vtk_index = 0
            do index = 1, get_number_of_cells(), 1
                if(.not. cells_is_real_cell(index))then
                    associate(                                     &
                        rho1 => primitive_variables_set(1, index), &
                        rho2 => primitive_variables_set(2, index), &
                        ie   => primitive_variables_set(6, index), &
                        z1   => primitive_variables_set(7, index))
                        vtk_scolar(vtk_index) = compute_pressure_mixture_stiffened_eos(ie, rho1 + rho2, z1)
                    end associate
                    vtk_index = vtk_index + 1
                end if
            end do
            vtk_error = a_vtk_file%xml_writer%write_dataarray   (data_name='pressure', x=vtk_scolar)
            vtk_error = a_vtk_file%xml_writer%write_dataarray   (location='node', action='close')
            vtk_error = a_vtk_file%xml_writer%write_piece()
            vtk_error = a_vtk_file%finalize()
            file_output_counter = file_output_counter + 1
        end if
    end do

    call finalize_boundary_reference()
    call finalize_faces()
    call finalize_cell_geometries()
    call finalize_variables()
end program five_eq_model_solver