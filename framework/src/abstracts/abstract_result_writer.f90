module abstract_result_writer
    implicit none

    private

    type, public, abstract :: result_writer
        contains
        procedure(initialize_interface       ), pass(self), deferred :: initialize
        procedure(write_interface            ), pass(self), deferred :: write
        procedure(next_writing_time_interface), pass(self), deferred :: next_writing_time
    end type result_writer

    abstract interface
        subroutine initialize_interface(self, is_real_cell, cell_geometries, cell_types, json)
            use typedef_module
            use class_point_id_list
            use json_module
            import result_writer
            class  (result_writer), intent(inout) :: self
            logical               , intent(in   ) :: is_real_cell   (:)
            class  (point_id_list), intent(in   ) :: cell_geometries(:)
            integer(type_kind    ), intent(in   ) :: cell_types     (:)
            class  (json_file    ), intent(in   ) :: json
        end subroutine initialize_interface

        subroutine write_interface(self, scolar_variables, vector_variables, time)
            use typedef_module
            import result_writer
            class(result_writer), intent(in) :: self
            real (real_kind    ), intent(in) :: scolar_variables(:,:)
            real (real_kind    ), intent(in) :: vector_variables(:,:)
            real (real_kind    ), intent(in) :: time
        end subroutine write_interface

        pure function next_writing_time_interface(self) result(t)
            use typedef_module
            import result_writer
            class(result_writer), intent(in) :: self
            real (real_kind    )             :: t
        end function
    end interface
end module abstract_result_writer