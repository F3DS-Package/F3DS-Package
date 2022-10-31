module abstract_parallelizer
    implicit none

    private

    type, public, abstract :: parallelizer
        contains
        procedure(initialize_interface                         ), pass(self), deferred :: initialize
        procedure(get_node_number_interface                    ), pass(self), deferred :: get_node_number
        procedure(get_number_of_threads_interface              ), pass(self), deferred :: get_number_of_threads
        procedure(get_number_of_face_indexes_interface         ), pass(self), deferred :: get_number_of_face_indexes
        procedure(get_number_of_boundary_face_indexes_interface), pass(self), deferred :: get_number_of_boundary_face_indexes
        procedure(get_number_of_cell_indexes_interface         ), pass(self), deferred :: get_number_of_cell_indexes
        procedure(get_face_index_interface                     ), pass(self), deferred :: get_face_index
        procedure(get_boundary_face_index_interface            ), pass(self), deferred :: get_boundary_face_index
        procedure(get_cell_index_interface                     ), pass(self), deferred :: get_cell_index
    end type parallelizer

    abstract interface
        subroutine initialize_interface(self, config, num_cells, num_faces, &
                                        num_outflow_faces, num_wall_faces, num_symmetric_faces, num_empty_faces)
            use abstract_configuration
            use typedef_module
            import parallelizer

            class  (parallelizer ), intent(inout) :: self
            class  (configuration), intent(inout) :: config
            integer(int_kind     ), intent(in   ) :: num_cells
            integer(int_kind     ), intent(in   ) :: num_faces
            integer(int_kind     ), intent(in   ) :: num_outflow_faces
            integer(int_kind     ), intent(in   ) :: num_wall_faces
            integer(int_kind     ), intent(in   ) :: num_symmetric_faces
            integer(int_kind     ), intent(in   ) :: num_empty_faces
        end subroutine initialize_interface

        pure function get_node_number_interface(self) result(n)
            use    typedef_module
            import parallelizer

            class  (parallelizer ), intent(in) :: self
            integer(int_kind     )             :: n
        end function get_node_number_interface

        pure function get_number_of_threads_interface(self, node_number) result(n)
            use    typedef_module
            import parallelizer

            class  (parallelizer ), intent(in) :: self
            integer(int_kind     ), intent(in) :: node_number
            integer(int_kind     )             :: n
        end function get_number_of_threads_interface

        pure function get_number_of_face_indexes_interface(self, node_number, thread_number) result(n)
            use    typedef_module
            import parallelizer

            class  (parallelizer   ), intent(in) :: self
            integer(int_kind       ), intent(in) :: node_number
            integer(int_kind       ), intent(in) :: thread_number
            integer(int_kind       )             :: n
        end function get_number_of_face_indexes_interface

        pure function get_number_of_boundary_face_indexes_interface(self, node_number, thread_number, a_face_type) result(n)
            use    typedef_module
            use    face_type_module
            import parallelizer

            class  (parallelizer           ), intent(in) :: self
            integer(int_kind               ), intent(in) :: node_number
            integer(int_kind               ), intent(in) :: thread_number
            integer(boundary_face_type_kind), intent(in) :: a_face_type
            integer(int_kind               )             :: n
        end function get_number_of_boundary_face_indexes_interface

        pure function get_number_of_cell_indexes_interface(self, node_number, thread_number) result(n)
            use    typedef_module
            import parallelizer

            class  (parallelizer ), intent(in) :: self
            integer(int_kind     ), intent(in) :: node_number
            integer(int_kind     ), intent(in) :: thread_number
            integer(int_kind     )             :: n
        end function get_number_of_cell_indexes_interface

        pure function get_cell_index_interface(self, node_number, thread_number, local_index) result(index)
            use    typedef_module
            import parallelizer

            class  (parallelizer ), intent(in) :: self
            integer(int_kind     ), intent(in) :: node_number
            integer(int_kind     ), intent(in) :: thread_number
            integer(int_kind     ), intent(in) :: local_index
            integer(int_kind     )             :: index
        end function get_cell_index_interface

         pure function get_face_index_interface(self, node_number, thread_number, local_index) result(index)
            use    typedef_module
            import parallelizer

            class  (parallelizer   ), intent(in) :: self
            integer(int_kind       ), intent(in) :: node_number
            integer(int_kind       ), intent(in) :: thread_number
            integer(int_kind       ), intent(in) :: local_index
            integer(int_kind       )             :: index
        end function get_face_index_interface

        pure function get_boundary_face_index_interface(self, node_number, thread_number, local_index, a_boundary_face_type) result(index)
            use    typedef_module
            use    face_type_module
            import parallelizer

            class  (parallelizer           ), intent(in) :: self
            integer(int_kind               ), intent(in) :: node_number
            integer(int_kind               ), intent(in) :: thread_number
            integer(int_kind               ), intent(in) :: local_index
            integer(boundary_face_type_kind), intent(in) :: a_boundary_face_type
            integer(int_kind               )             :: index
        end function get_boundary_face_index_interface
    end interface
end module abstract_parallelizer