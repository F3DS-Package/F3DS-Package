program five_eq_model_solver
    !$ use omp_lib
    ! Utils
    use typedef_module
    ! Cell system
    use faces_module
    use boundary_reference_module
    use cell_geometries_module
    use five_equation_model_variables_module
    ! Scheme
    use second_order_tvd_rk_module
    use third_order_tvd_rk_module
    use jiang_weno5_module
    use five_equation_model_rho_thinc_module
    use five_equation_space_model_module
    use five_equation_model_hllc_module
    ! Model
    use class_stiffened_gas_mixture_eos
    ! Grid & initial condition reader
    use class_nlinit_parser
    use class_nlgrid_parser
    ! Result file output
    use penf
    use vtk_fortran, only : vtk_file
    use sytem_call_module
    ! BC
    use five_equation_model_boundary_condition_module
    ! Measurement
    use sensor_class
    use measurement_surface_class
    use line_plot_class

    implicit none

    real   (real_kind)  :: time_increment, time
    integer(int_kind )  :: max_timestep, timestep
    integer(int_kind )  :: index
    ! Grid & initial condition I/O
    type(nlgrid_parser) :: a_grid_parser
    type(nlinit_parser) :: a_init_parser
    ! Result file output
    type(vtk_file)              :: a_vtk_file
    integer  (I4P)              :: vtk_error
    integer  (int_kind )        :: file_output_counter, vtk_index, cell_point_index
    character(22)               :: vtk_filename
    integer  (I4P)              :: n_output_cells, n_output_points, n_cell_points, offset_incriment, n_output_file
    real     (R4P), allocatable :: vtk_density(:), vtk_cell_id(:), vtk_soundspeed(:), vtk_wood_soundspeed(:)
    integer  (I1P), allocatable :: vtk_cell_type(:)
    integer  (I4P), allocatable :: vtk_offset(:), vtk_connect(:)

    type(stiffened_gas_mixture_eos) :: eos
    type(sensor             ) :: pressure_sensor
    type(measurement_surface) :: surface
    type(line_plot          ) :: line
    real(real_kind)           :: normal(3)
    integer(int_kind ), allocatable :: line_ids(:)
    integer(int_kind )              :: n_line_ids

    time_increment = 1.d-4
    max_timestep   = 3*10**4!5*10**5
    n_output_file  = 100
    time           = 0.d0

#ifdef _OPENMP
    call omp_set_num_threads(16)
#endif

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
    call a_grid_parser%get_faces          (faces_to_cell_index, faces_normal_vector, faces_tangential1_vector, faces_tangential2_vector, faces_position, faces_area)
    call a_grid_parser%get_boundaries     (outflow_face_indexs, slipwall_face_indexs, symmetric_face_indexs)
    call a_grid_parser%get_cell_geometries(points, cell_geometries)
    ! close file
    call a_grid_parser%close()

    ! parse init file
    call a_init_parser%parse("init.nlinit")
    call a_init_parser%get_conservative_variables_set(conservative_variables_set)
    call a_init_parser%close()

    ! sensor
    call pressure_sensor%initialize("sensor.dat", 5.d3, 13633)

    ! measurement surface
    normal(1) = 0.d0
    normal(2) = 1.d0
    normal(3) = 0.d0
    call surface%initialize("surface.dat", 5.d3, 10.d0+0.04d0, 10.d0-0.04d0, 1.d0, 0.d0, 7.d0, 0.d0, normal, faces_to_cell_index, faces_position, faces_normal_vector, cells_is_real_cell)

    ! line plot
    allocate(line_ids(get_number_of_cells()))
    n_line_ids = 0
    do index = 1, get_number_of_cells(), 1
        if  ( 0.0d0 < cells_centor_position(index, 1) .and. cells_centor_position(index, 1) < 10.d0  &
        .and. 0.d0  < cells_centor_position(index, 2) .and. cells_centor_position(index, 2) < 0.04d0 &
        .and. 0.d0  < cells_centor_position(index, 3) .and. cells_centor_position(index, 3) < 0.04d0 ) then
            n_line_ids = n_line_ids + 1
            line_ids(n_line_ids) = index
        end if
    end do
    call line%initialize("line", 0.05d3, line_ids(1:n_line_ids))
    deallocate(line_ids)

    ! VTK
    file_output_counter = 0
    n_output_cells = 0
    offset_incriment=8_I4P
    do index = 1, get_number_of_cells(), 1
        if(cells_is_real_cell(index)) n_output_cells = n_output_cells + 1
    end do
    allocate(vtk_density        (n_output_cells    ))
    allocate(vtk_soundspeed     (n_output_cells    ))
    allocate(vtk_wood_soundspeed(n_output_cells    ))
    allocate(vtk_cell_id        (n_output_cells    ))
    allocate(vtk_cell_type      (n_output_cells    ))
    allocate(vtk_offset         (n_output_cells    ))
    allocate(vtk_connect        (n_output_cells * 8))
    vtk_index = 1
    do index = 1, get_number_of_cells(), 1
        if(cells_is_real_cell(index))then
            n_cell_points = cell_geometries(index)%get_number_of_points()
            do cell_point_index = 1, n_cell_points, 1
                vtk_connect((vtk_index - 1) * 8 + cell_point_index) = cell_geometries(index)%get_point_id(cell_point_index)
            end do
            vtk_cell_type(vtk_index) = 12_I1P
            vtk_offset   (vtk_index) = offset_incriment
            offset_incriment = offset_incriment + 8_I4P
            vtk_index = vtk_index + 1
        end if
    end do
    n_output_points = get_number_of_points()
    call make_dir("result/field")

    ! EoS and primitive-valiables
    !call eos%initialize(1.4d0, 2.35d0, 0.d0, 7.142d3) ![Pelanti 2019]
    call eos%initialize(1.4d0, 4.4d0, 0.d0, 4.286d3) ![Garick 2019]
    !call eos%initialize(1.4d0, 4.4d0, 0.d0, 4.3d-2) ![Garick 2019]
    !call eos%initialize(1.4d0, 7.0d0, 0.d0, 2.142d3) ![Lou 2004]
    !call eos%initialize(1.4d0, 6.12d0, 0.d0, 2.450d3) ! [Nishida]
    do index = 1, get_number_of_cells(), 1
        primitive_variables_set(index, :) = conservative_to_primitive(conservative_variables_set(index, :), eos)
    end do

    ! time-stepping
    call initialize_second_order_tvd_rk(conservative_variables_set)

    ! solver timestepping loop
    do timestep = 0, max_timestep, 1
        print *, "step ", timestep

        if (mod(timestep, max_timestep / n_output_file) == 0) then
            write(vtk_filename, "(a, i5.5, a)") "result/field/", file_output_counter, ".vtu"
            print *, "write vtk "//vtk_filename//"..."
            vtk_error = a_vtk_file%initialize                   (format="binary", filename=vtk_filename, mesh_topology="UnstructuredGrid")
            vtk_error = a_vtk_file%xml_writer%write_piece       (np=n_output_points, nc=n_output_cells)
            vtk_error = a_vtk_file%xml_writer%write_geo         (np=n_output_points, nc=n_output_cells, x=points(:, 1), y=points(:, 2), z=points(:, 3))
            vtk_error = a_vtk_file%xml_writer%write_connectivity(nc=n_output_cells, connectivity=vtk_connect, offset=vtk_offset, cell_type=vtk_cell_type)
            vtk_error = a_vtk_file%xml_writer%write_dataarray   (location='cell', action='open')
            vtk_index = 1
            do index = 1, get_number_of_cells(), 1
                if(cells_is_real_cell(index))then
                    associate(                                     &
                        rho1   => primitive_variables_set(index, 1), &
                        rho2   => primitive_variables_set(index, 2), &
                        p      => primitive_variables_set(index, 6), &
                        z1     => primitive_variables_set(index, 7))
                        vtk_density   (vtk_index) = rho1 * z1 + rho2 * (1.d0 - z1)
                        vtk_soundspeed(vtk_index) = eos%compute_soundspeed(p, rho1 * z1 + rho2 * (1.d0 - z1), z1)
                        vtk_wood_soundspeed(vtk_index) = eos%compute_wood_soundspeed(p, rho1, rho2, z1)
                        vtk_cell_id   (vtk_index) = index
                    end associate
                    vtk_index = vtk_index + 1
                end if
            end do
            vtk_error = a_vtk_file%xml_writer%write_dataarray(data_name='density', x=vtk_density)
            vtk_error = a_vtk_file%xml_writer%write_dataarray(data_name='density1', x=pack(primitive_variables_set(:, 1), mask=cells_is_real_cell))
            vtk_error = a_vtk_file%xml_writer%write_dataarray(data_name='density2', x=pack(primitive_variables_set(:, 2), mask=cells_is_real_cell))
            vtk_error = a_vtk_file%xml_writer%write_dataarray(data_name='velocity', x=pack(primitive_variables_set(:, 3), mask=cells_is_real_cell), y=pack(primitive_variables_set(:, 4), mask=cells_is_real_cell), z=pack(primitive_variables_set(:, 5), mask=cells_is_real_cell))
            vtk_error = a_vtk_file%xml_writer%write_dataarray(data_name='pressure', x=pack(primitive_variables_set(:, 6), mask=cells_is_real_cell))
            vtk_error = a_vtk_file%xml_writer%write_dataarray(data_name='volume fruction', x=pack(primitive_variables_set(:, 7), mask=cells_is_real_cell))
            vtk_error = a_vtk_file%xml_writer%write_dataarray(data_name='sound speed', x=vtk_soundspeed)
            vtk_error = a_vtk_file%xml_writer%write_dataarray(data_name='wood sound speed', x=vtk_wood_soundspeed)
            vtk_error = a_vtk_file%xml_writer%write_dataarray(data_name='cell id', x=vtk_cell_id)
            vtk_error = a_vtk_file%xml_writer%write_dataarray(data_name='cell position', x=pack(cells_centor_position(:, 1), mask=cells_is_real_cell), y=pack(cells_centor_position(:, 2), mask=cells_is_real_cell), z=pack(cells_centor_position(:, 3), mask=cells_is_real_cell))
            vtk_error = a_vtk_file%xml_writer%write_dataarray(location='cell', action='close')
            vtk_error = a_vtk_file%xml_writer%write_piece()
            vtk_error = a_vtk_file%finalize()
            file_output_counter = file_output_counter + 1
        end if

        call pressure_sensor%write(time, primitive_variables_set)
        call surface%write(time, primitive_variables_set, faces_to_cell_index, faces_area)
        call line%write(time, primitive_variables_set, cells_centor_position)

        call compute_next_state_second_order_tvd_rk(   &
            conservative_variables_set               , &
            primitive_variables_set                  , &
            derivative_variables_set                 , &
            cells_centor_position                    , & ! <-
            cells_volume                             , & ! <-
            faces_to_cell_index                      , & ! <-
            faces_normal_vector                      , & ! <-
            faces_tangential1_vector                 , & ! <-
            faces_tangential2_vector                 , & ! <-
            faces_position                           , & ! <-
            faces_area                               , & ! <-
            outflow_face_indexs                      , & ! <-
            slipwall_face_indexs                     , & ! <-
            symmetric_face_indexs                    , & ! <-
            get_number_of_cells()                    , & ! <-
            get_number_of_faces()                    , & ! <-
            get_number_of_ghost_cells()              , & ! <-
            get_number_of_outflow_faces()            , &
            get_number_of_slipwall_faces()           , &
            get_number_of_symmetric_faces()          , &
            time_increment                           , &
            eos                                      , &
            reconstruct_rho_thinc                    , &
            compute_space_element_five_equation_model, &
            compute_flux_five_equation_model_hllc    , &
            primitive_to_conservative                , &
            conservative_to_primitive                , &
            five_equation_model_set_boundary_condition &
        )

        time = time + time_increment
    end do

    call finalize_boundary_reference()
    call finalize_faces()
    call finalize_cell_geometries()
    call finalize_variables()
end program five_eq_model_solver