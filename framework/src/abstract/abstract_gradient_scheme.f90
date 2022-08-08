module abstract_gradient_scheme
    implicit none

    private

    type, public, abstract :: gradient_scheme
        contains

        procedure(initialize_interface              ), pass(self), deferred :: initialize
        procedure(compute_gradient_interface        ), pass(self), deferred :: compute_gradient
        procedure(compure_gradient_1darray_interface), pass(self), deferred :: compure_gradient_1darray
    end type

    abstract interface

        subroutine initialize_interface(self, a_configuration)
            use abstract_configuration
            import gradient_scheme
            class(gradient_scheme), intent(inout) :: self
            class(configuration  ), intent(inout) :: a_configuration
        end subroutine initialize_interface

        subroutine compute_gradient_interface(self, variables_set, gradient_variables_set, cell_centor_positions, cell_volumes, face_to_cell_indexes, face_normal_vectors, face_positions, face_areas, num_faces)
            use typedef_module
            import gradient_scheme
            class  (gradient_scheme), intent(inout) :: self
            real   (real_kind      ), intent(in   ) :: variables_set            (:,:)
            real   (real_kind      ), intent(in   ) :: gradient_variables_set   (:,:)
            real   (real_kind      ), intent(in   ) :: cell_centor_positions    (:,:)
            real   (real_kind      ), intent(in   ) :: cell_volumes             (:)
            integer(int_kind       ), intent(in   ) :: face_to_cell_indexes     (:,:)
            real   (real_kind      ), intent(in   ) :: face_normal_vectors      (:,:)
            real   (real_kind      ), intent(in   ) :: face_positions           (:,:)
            real   (real_kind      ), intent(in   ) :: face_areas               (:)
            integer(int_kind       ), intent(in   ) :: num_faces
        end subroutine compute_gradient_interface

        subroutine compure_gradient_1darray_interface(self, variable_set, gradient_variable_set, cell_centor_positions, cell_volumes, face_to_cell_indexes, face_normal_vectors, face_positions, face_areas, num_faces)
            use typedef_module
            import gradient_scheme
            class  (gradient_scheme), intent(inout) :: self
            real   (real_kind      ), intent(in   ) :: variable_set             (:)
            real   (real_kind      ), intent(in   ) :: gradient_variable_set    (:)
            real   (real_kind      ), intent(in   ) :: cell_centor_positions    (:,:)
            real   (real_kind      ), intent(in   ) :: cell_volumes             (:)
            integer(int_kind       ), intent(in   ) :: face_to_cell_indexes     (:,:)
            real   (real_kind      ), intent(in   ) :: face_normal_vectors      (:,:)
            real   (real_kind      ), intent(in   ) :: face_positions           (:,:)
            real   (real_kind      ), intent(in   ) :: face_areas               (:)
            integer(int_kind       ), intent(in   ) :: num_faces
        end subroutine compure_gradient_1darray_interface

    end interface
end module abstract_gradient_scheme