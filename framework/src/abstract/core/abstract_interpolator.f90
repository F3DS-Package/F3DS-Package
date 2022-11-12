module abstract_interpolator
    implicit none

    private

    type, public, abstract :: interpolator
        contains

        procedure(initialize_interface                ), pass(self), deferred :: initialize
        procedure(interpolate_face_variables_interface), pass(self), deferred :: interpolate_face_variables
    end type interpolator

    abstract interface
        subroutine initialize_interface(self, config)
            use abstract_configuration
            import interpolator

            class(interpolator ), intent(inout) :: self
            class(configuration), intent(in   ) :: config
        end subroutine initialize_interface

        pure function interpolate_face_variables_interface(self, variables, face_to_cell_index, cell_positions, face_position, num_local_cells, num_variables) result(face_variables)
            use typedef_module
            import interpolator

            class  (interpolator), intent(in) :: self
            real   (real_kind   ), intent(in) :: variables         (:,:)
            integer(int_kind    ), intent(in) :: face_to_cell_index(:)
            real   (real_kind   ), intent(in) :: cell_positions    (:,:)
            real   (real_kind   ), intent(in) :: face_position     (3)
            integer(int_kind    ), intent(in) :: num_local_cells
            integer(int_kind    ), intent(in) :: num_variables
            real   (real_kind   )             :: face_variables(num_variables)
        end function interpolate_face_variables_interface
    end interface
end module abstract_interpolator