submodule(class_cellsystem) write_result_impl
    implicit none

    contains

    module subroutine result_writer_open_file(self, writer)
        class(cellsystem   ), intent(inout) :: self
        class(result_writer), intent(inout) :: writer
#ifdef _DEBUG
        call write_debuginfo("In result_writer_open_file(), cellsystem.")
#endif
        call writer%open_file(self%time, self%points)
    end subroutine result_writer_open_file

    module subroutine result_writer_close_file(self, writer)
        class(cellsystem   ), intent(inout) :: self
        class(result_writer), intent(inout) :: writer
#ifdef _DEBUG
        call write_debuginfo("In result_writer_close_file(), cellsystem.")
#endif
        call writer%close_file()
    end subroutine result_writer_close_file

    module pure function result_writer_is_writable(self, writer) result(yes)
        class(cellsystem   ), intent(in) :: self
        class(result_writer), intent(in) :: writer
        logical                          :: yes
        yes = writer%is_writable(self%time)
    end function result_writer_is_writable

    module subroutine write_scolar(self, writer, name, scolar_variables)
        class    (cellsystem   ), intent(inout) :: self
        class    (result_writer), intent(inout) :: writer
        character(len=*        ), intent(in   ) :: name
        real     (real_kind    ), intent(in   ) :: scolar_variables(:)
#ifdef _DEBUG
        call write_debuginfo("In write_scolar(), cellsystem.")
#endif
        call writer%write_scolar(self%is_real_cell, name, scolar_variables)
    end subroutine write_scolar

    module subroutine write_vector(self, writer, name, vector_variables)
        class    (cellsystem   ), intent(inout) :: self
        class    (result_writer), intent(inout) :: writer
        character(len=*        ), intent(in   ) :: name
        real     (real_kind    ), intent(in   ) :: vector_variables(:,:)
#ifdef _DEBUG
        call write_debuginfo("In write_vector(), cellsystem.")
#endif
        call writer%write_vector(self%is_real_cell, name, vector_variables)
    end subroutine write_vector

    module function get_filename(self, writer) result(name)
        class    (cellsystem   ), intent(inout) :: self
        class    (result_writer), intent(inout) :: writer
        character(len=:        ), allocatable :: name
        name = writer%get_filename()
    end function get_filename
end submodule