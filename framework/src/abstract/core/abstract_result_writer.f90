module abstract_result_writer
    ! """
    ! This module defines an abstract interface for result writers.
    ! The `result_writer` type is an abstract type that defines the interface for result writers. It contains deferred procedures for
    ! initializing, cleaning up, opening and closing files, writing scalar and vector variables, checking if a writer is writable at a
    ! given time, and getting the filename associated with the writer.
    ! The `initialize_interface` subroutine initializes the result writer with the number of cells, number of points, information
    ! about whether cells are real or not, cell geometries, cell types, and a configuration object.
    ! The `cleanup_interface` subroutine cleans up the result writer.
    ! The `open_file_interface` subroutine opens a file for writing results, given a time and a matrix of points.
    ! The `close_file_interface` subroutine closes the file associated with the result writer.
    ! The `write_scalar_interface` subroutine writes scalar variables to the file, given information about whether cells are real or
    ! not, a variable name, and an array of scalar variables.
    ! The `write_vector_interface` subroutine writes vector variables to the file, given information about whether cells are real or
    ! not, a variable name, and a matrix of vector variables.
    ! The `is_writable_interface` function checks if the result writer is writable at a given time.
    ! The `get_filename_interface` function returns the filename associated with the result writer.
    ! Note: This module is written in Fortran.
    ! """
    implicit none

    private

    type, public, abstract :: result_writer
        ! """
        ! Abstract base type for result writers.
        ! Methods:
        ! - `initialize`: Initialize the result writer.
        ! - `cleanup`: Clean up any resources used by the result writer.
        ! - `open_file`: Open the file for writing.
        ! - `close_file`: Close the file.
        ! - `write_scolar`: Write a scalar value to the file.
        ! - `write_vector`: Write a vector of values to the file.
        ! - `is_writable`: Check if the result writer is writable.
        ! - `get_filename`: Get the filename associated with the result writer.
        ! """
        contains
        procedure(initialize_interface  ), pass(self), deferred :: initialize
        procedure(cleanup_interface     ), pass(self), deferred :: cleanup
        procedure(open_file_interface   ), pass(self), deferred :: open_file
        procedure(close_file_interface  ), pass(self), deferred :: close_file
        procedure(write_scolar_interface), pass(self), deferred :: write_scolar
        procedure(write_vector_interface), pass(self), deferred :: write_vector
        procedure(is_writable_interface ), pass(self), deferred :: is_writable
        procedure(get_filename_interface), pass(self), deferred :: get_filename
    end type result_writer

    abstract interface
        subroutine initialize_interface(self, num_cells, num_points, is_real_cell, cell_geometries, cell_types, config)
            ! This subroutine initializes the interface of the result_writer class.
            !
            ! Parameters:
            ! - self: The instance of the result_writer class.
            ! - num_cells: The number of cells.
            ! - num_points: The number of points.
            ! - is_real_cell: An array of logical values indicating whether each cell is real or not.
            ! - cell_geometries: An array of point_id_list objects representing the geometries of each cell.
            ! - cell_types: An array of integer values representing the types of each cell.
            ! - config: The configuration object.
            use typedef_module
            use class_point_id_list
            use abstract_configuration
            import result_writer
            class  (result_writer), intent(inout) :: self
            integer(int_kind     ), intent(in   ) :: num_cells
            integer(int_kind     ), intent(in   ) :: num_points
            logical               , intent(in   ) :: is_real_cell   (:)
            class  (point_id_list), intent(in   ) :: cell_geometries(:)
            integer(type_kind    ), intent(in   ) :: cell_types     (:)
            class  (configuration), intent(inout) :: config
        end subroutine initialize_interface

        subroutine cleanup_interface(self)
            ! -----------------------------------------------------------------------
            ! Subroutine cleanup_interface
            ! -----------------------------------------------------------------------
            ! Cleans up the interface of the result_writer object.
            !
            ! Parameters:
            ! self: inout, class(result_writer)
            ! The result_writer object to be cleaned up.
            ! -----------------------------------------------------------------------
            import result_writer
            class  (result_writer), intent(inout) :: self
        end subroutine cleanup_interface

        subroutine open_file_interface(self, time, points)
            ! -----------------------------------------------------------------------
            ! This subroutine is used to open a file interface.
            !
            ! Parameters:
            ! self: class(result_writer), intent(inout)
            ! - The instance of the result_writer class.
            ! time: real (real_kind), intent(in)
            ! - The time value.
            ! points: real (real_kind), intent(in)
            ! - The array of points.
            !
            ! -----------------------------------------------------------------------
            use typedef_module
            import result_writer
            class(result_writer), intent(inout) :: self
            real (real_kind    ), intent(in)    :: time
            real (real_kind    ), intent(in)    :: points(:, :)
        end subroutine open_file_interface

        subroutine close_file_interface(self)
            ! """
            ! Closes the file associated with the result_writer object.
            ! Parameters:
            !     self (result_writer): The result_writer object whose file needs to be closed.
            ! """
            import result_writer
            class(result_writer), intent(inout) :: self
        end subroutine close_file_interface

        subroutine write_scolar_interface(self, is_real_cell, name, scolar_variables)
            ! -------------------------------------------------------------------------
            ! Writes the scalar interface for a given result writer object.
            !
            ! Parameters:
            ! self: result_writer
            ! The result writer object to write the scalar interface for.
            ! is_real_cell: logical, dimension(:)
            ! An array indicating whether each cell is real or not.
            ! name: character(len=*)
            ! The name of the scalar interface.
            ! scolar_variables: real(real_kind), dimension(:)
            ! An array of scalar variables to be written.
            !
            ! -------------------------------------------------------------------------
            use typedef_module
            import result_writer
            class    (result_writer), intent(inout) :: self
            logical                 , intent(in   ) :: is_real_cell(:)
            character(len=*        ), intent(in   ) :: name
            real     (real_kind    ), intent(in   ) :: scolar_variables(:)
        end subroutine write_scolar_interface

        subroutine write_vector_interface(self, is_real_cell, name, vector_variables)
            ! ---------------------------------------------------------------------------
            ! Writes vector variables to the result file.
            !
            ! Parameters:
            ! self: result_writer
            ! The result writer object.
            ! is_real_cell: logical, dimension(:)
            ! Array indicating whether each cell is real or not.
            ! name: character(len=*)
            ! The name of the vector variables.
            ! vector_variables: real(real_kind), dimension(:,:)
            ! The vector variables to be written.
            ! ---------------------------------------------------------------------------
            use typedef_module
            import result_writer
            class    (result_writer), intent(inout) :: self
            logical                 , intent(in   ) :: is_real_cell(:)
            character(len=*        ), intent(in   ) :: name
            real     (real_kind    ), intent(in   ) :: vector_variables(:,:)
        end subroutine write_vector_interface

        pure function is_writable_interface(self, time) result(yes)
            ! -----------------------------------------------------------------------
            ! Function: is_writable_interface
            ! -----------------------------------------------------------------------
            ! Checks if the given result_writer object is writable at the specified time.
            !
            ! Parameters:
            ! - self: The result_writer object to be checked.
            ! - time: The time at which the writability is being checked.
            !
            ! Returns:
            ! - yes: A logical value indicating whether the object is writable or not.
            ! -----------------------------------------------------------------------
            use typedef_module
            import result_writer
            class(result_writer), intent(in) :: self
            real (real_kind    ), intent(in) :: time
            logical                          :: yes
        end function is_writable_interface

        function get_filename_interface(self) result(name)
            ! > Returns the filename associated with the given `result_writer` object.
            ! Parameters:
            ! - self: The `result_writer` object for which the filename is to be retrieved.
            ! Returns:
            ! - name: The filename associated with the `result_writer` object.
            use typedef_module
            import result_writer
            class    (result_writer), intent(in)  :: self
            character(len=:        ), allocatable :: name
        end function get_filename_interface
    end interface
end module abstract_result_writer