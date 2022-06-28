module class_vtk_result_writer
    use vtk_fortran
    use typedef_module
    use stdio_module
    use system_call_module
    use class_point_id_list
    use abstract_result_writer
    use abstract_configuration

    implicit none

    private

    type, public, extends(result_writer) :: vtk_result_writer
        private

        type(vtk_file) :: a_vtk_file

        real   (real_kind), allocatable :: cell_id  (:)
        integer(type_kind), allocatable :: cell_type(:)
        integer(int_kind ), allocatable :: offset   (:)
        integer(int_kind ), allocatable :: connect  (:)

        integer(int_kind) :: n_output_cells
        integer(int_kind) :: n_output_points

        integer(int_kind) :: n_output_file

        real(real_kind) :: output_timespan
        real(real_kind) :: next_output_time

        logical :: file_is_opened
        character(len=:), allocatable :: vtk_filename

        integer(int_kind) :: vtk_error

        contains

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: open_file
        procedure, public, pass(self) :: close_file
        procedure, public, pass(self) :: write_scolar
        procedure, public, pass(self) :: write_vector
        procedure, public, pass(self) :: is_writable
        procedure, public, pass(self) :: get_filename
    end type vtk_result_writer

    contains

    subroutine initialize(self, num_cells, num_points, is_real_cell, cell_geometries, cell_types, config)
        class  (vtk_result_writer), intent(inout) :: self
        integer(int_kind         ), intent(in   ) :: num_cells
        integer(int_kind         ), intent(in   ) :: num_points
        logical                   , intent(in   ) :: is_real_cell   (:)
        class  (point_id_list    ), intent(in   ) :: cell_geometries(:)
        integer(type_kind        ), intent(in   ) :: cell_types     (:)
        class  (configuration    ), intent(inout) :: config

        integer(int_kind ) :: index, vtk_index, cell_point_index
        integer(int_kind ) :: num_cell_points, offset_incriment
        real   (real_kind) :: frequency
        logical :: found

        ! count number of cells
        self%n_output_cells = 0
        do index = 1, num_cells, 1
            if(is_real_cell(index)) self%n_output_cells = self%n_output_cells + 1
        end do

        ! allocate
        allocate(self%cell_id  (self%n_output_cells    ))
        allocate(self%cell_type(self%n_output_cells    ))
        allocate(self%offset   (self%n_output_cells    ))
        allocate(self%connect  (self%n_output_cells * 8))

        ! make vtk variables
        vtk_index = 1
        offset_incriment = 0
        do index = 1, num_cells, 1
            if(is_real_cell(index))then
                num_cell_points = cell_geometries(index)%get_number_of_points()
                do cell_point_index = 1, num_cell_points, 1
                    self%connect((vtk_index - 1) * 8 + cell_point_index) = cell_geometries(index)%get_point_id(cell_point_index)
                end do
                self%cell_type(vtk_index) = cell_types(index)
                offset_incriment = offset_incriment + num_cell_points
                self%offset   (vtk_index) = offset_incriment
                self%cell_id  (vtk_index) = index
                vtk_index = vtk_index + 1
            end if
        end do
        self%n_output_points = num_points
        self%n_output_file   = 0
        call config%get_real("Result writer.Output frequency", frequency, found, 1.d3)
        if(.not. found) call write_warring("'Result writer.Output frequency' is not found in configration you set. To be set dafault value.")
        self%output_timespan = 1.d0 / frequency
        self%next_output_time = 0.d0

        ! make directory
        call make_dir("result/field")

        ! initialize
        self%file_is_opened = .false.
    end subroutine initialize

    subroutine open_file(self, time, points)
        class(vtk_result_writer), intent(inout) :: self
        real (real_kind        ), intent(in   ) :: time
        real (real_kind        ), intent(in   ) :: points(:, :)

        if(self%file_is_opened) call call_error("VTK result file '"//self%vtk_filename//"' is already opened. But you call open_file method.")

        if(time > self%next_output_time)then
            write(self%vtk_filename, "(a, i5.5, a)") "result/field/", self%n_output_file, ".vtu"

            self%vtk_error = self%a_vtk_file%initialize                   (format="binary", filename=self%vtk_filename, mesh_topology="UnstructuredGrid")
            self%vtk_error = self%a_vtk_file%xml_writer%write_piece       (np=self%n_output_points, nc=self%n_output_cells)
            self%vtk_error = self%a_vtk_file%xml_writer%write_geo         (np=self%n_output_points, nc=self%n_output_cells, x=points(:, 1), y=points(:, 2), z=points(:, 3))
            self%vtk_error = self%a_vtk_file%xml_writer%write_connectivity(nc=self%n_output_cells, connectivity=self%connect, offset=self%offset, cell_type=self%cell_type)
            self%vtk_error = self%a_vtk_file%xml_writer%write_dataarray   (location="cell", action="open")
            self%vtk_error = self%a_vtk_file%xml_writer%write_dataarray(data_name="cell id", x=self%cell_id)

            self%file_is_opened = .true.
        end if
    end subroutine open_file

    subroutine close_file(self)
        class(vtk_result_writer), intent(inout) :: self

        if(self%file_is_opened)then
            self%vtk_error = self%a_vtk_file%xml_writer%write_dataarray(location='cell', action='close')
            self%vtk_error = self%a_vtk_file%xml_writer%write_piece()
            self%vtk_error = self%a_vtk_file%finalize()

            self%n_output_file = self%n_output_file + 1
        end if
    end subroutine

    subroutine write_scolar(self, is_real_cell, name, scolar_variables)
        class    (vtk_result_writer), intent(inout) :: self
        logical                     , intent(in   ) :: is_real_cell(:)
        character(len=*            ), intent(in   ) :: name
        real     (real_kind        ), intent(in   ) :: scolar_variables(:)

        if(self%file_is_opened)then
            self%vtk_error = self%a_vtk_file%xml_writer%write_dataarray(data_name=name, x=pack(scolar_variables(:), mask=is_real_cell))
        end if
    end subroutine write_scolar

    subroutine write_vector(self, is_real_cell, name, vector_variables)
        class    (vtk_result_writer), intent(inout) :: self
        logical                     , intent(in   ) :: is_real_cell(:)
        character(len=*            ), intent(in   ) :: name
        real     (real_kind        ), intent(in   ) :: vector_variables(:,:)

        if(self%file_is_opened)then
            self%vtk_error = self%a_vtk_file%xml_writer%write_dataarray(data_name=name, x=pack(vector_variables(:, 1), mask=is_real_cell), y=pack(vector_variables(:, 2), mask=is_real_cell), z=pack(vector_variables(:, 3), mask=is_real_cell))
        end if
    end subroutine write_vector

    pure function is_writable(self, time) result(yes)
        class(vtk_result_writer), intent(in) :: self
        real (real_kind        ), intent(in) :: time
        logical                              :: yes

        if(time > self%next_output_time)then
            yes = .true.
        else
            yes = .false.
        end if
    end function is_writable

    function get_filename(self) result(name)
        class    (vtk_result_writer), intent(in)  :: self
        character(len=:            ), allocatable :: name
        if(self%file_is_opened)then
            name = self%vtk_filename
        else
            call call_error("Result file is not opened. But you call get_filename method.")
        end if
    end function get_filename
end module class_vtk_result_writer