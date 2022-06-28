module abstract_result_writer
    implicit none

    private

    type, public, abstract :: result_writer
        contains
        procedure(initialize_interface  ), pass(self), deferred :: initialize
        procedure(open_file_interface   ), pass(self), deferred :: open_file
        procedure(close_file_interface  ), pass(self), deferred :: close_file
        procedure(write_scolar_interface), pass(self), deferred :: write_scolar
        procedure(write_vector_interface), pass(self), deferred :: write_vector
        procedure(is_writable_interface ), pass(self), deferred :: is_writable
        procedure(get_filename_interface), pass(self), deferred :: get_filename
    end type result_writer

    abstract interface
        subroutine initialize_interface(self, num_cells, num_points, is_real_cell, cell_geometries, cell_types, config)
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

        subroutine open_file_interface(self, time, points)
            use typedef_module
            import result_writer
            class(result_writer), intent(inout) :: self
            real (real_kind    ), intent(in)    :: time
            real (real_kind    ), intent(in)    :: points(:, :)
        end subroutine open_file_interface

        subroutine close_file_interface(self)
            import result_writer
            class(result_writer), intent(inout) :: self
        end subroutine close_file_interface

        subroutine write_scolar_interface(self, is_real_cell, name, scolar_variables)
            use typedef_module
            import result_writer
            class    (result_writer), intent(inout) :: self
            logical                 , intent(in   ) :: is_real_cell(:)
            character(len=*        ), intent(in   ) :: name
            real     (real_kind    ), intent(in   ) :: scolar_variables(:)
        end subroutine write_scolar_interface

        subroutine write_vector_interface(self, is_real_cell, name, vector_variables)
            use typedef_module
            import result_writer
            class    (result_writer), intent(inout) :: self
            logical                 , intent(in   ) :: is_real_cell(:)
            character(len=*        ), intent(in   ) :: name
            real     (real_kind    ), intent(in   ) :: vector_variables(:,:)
        end subroutine write_vector_interface

        pure function is_writable_interface(self, time) result(yes)
            use typedef_module
            import result_writer
            class(result_writer), intent(in) :: self
            real (real_kind    ), intent(in) :: time
            logical                          :: yes
        end function is_writable_interface

        function get_filename_interface(self) result(name)
            use typedef_module
            import result_writer
            class    (result_writer), intent(in)  :: self
            character(len=:        ), allocatable :: name
        end function get_filename_interface
    end interface
end module abstract_result_writer