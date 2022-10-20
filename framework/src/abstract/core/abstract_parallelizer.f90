module abstract_parallelizer
    implicit none

    private

    type, public, abstract :: parallelizer
        contains
        procedure(initialize_interface                ), pass(self), deferred :: initialize
        procedure(get_node_number_interface           ), pass(self), deferred :: get_node_number
        procedure(get_number_of_threads_interface     ), pass(self), deferred :: get_number_of_threads
        procedure(get_number_of_face_indexes_interface), pass(self), deferred :: get_number_of_face_indexes
        procedure(get_number_of_cell_indexes_interface), pass(self), deferred :: get_number_of_cell_indexes
        procedure(get_face_index_interface            ), pass(self), deferred :: get_face_index
        procedure(get_cell_index_interface            ), pass(self), deferred :: get_cell_index
    end type parallelizer

    abstract interface
        subroutine initialize_interface(self, config, number_of_cells, number_of_faces)
            use abstract_configuration
            use typedef_module
            import parallelizer

            class  (parallelizer ), intent(inout) :: self
            class  (configuration), intent(inout) :: config
            integer(int_kind     ), intent(in   ) :: number_of_cells
            integer(int_kind     ), intent(in   ) :: number_of_faces
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

            class  (parallelizer ), intent(in) :: self
            integer(int_kind     ), intent(in) :: node_number
            integer(int_kind     ), intent(in) :: thread_number
            integer(int_kind     )             :: n
        end function get_number_of_face_indexes_interface

        pure function get_number_of_cell_indexes_interface(self, node_number, thread_number) result(n)
            use    typedef_module
            import parallelizer

            class  (parallelizer ), intent(in) :: self
            integer(int_kind     ), intent(in) :: node_number
            integer(int_kind     ), intent(in) :: thread_number
            integer(int_kind     )             :: n
        end function get_number_of_cell_indexes_interface

        pure function get_face_index_interface(self, node_number, thread_number, local_index) result(index)
            use    typedef_module
            import parallelizer

            class  (parallelizer ), intent(in) :: self
            integer(int_kind     ), intent(in) :: node_number
            integer(int_kind     ), intent(in) :: thread_number
            integer(int_kind     ), intent(in) :: local_index
            integer(int_kind     )             :: index
        end function get_face_index_interface

        pure function get_cell_index_interface(self, node_number, thread_number, local_index) result(index)
            use    typedef_module
            import parallelizer

            class  (parallelizer ), intent(in) :: self
            integer(int_kind     ), intent(in) :: node_number
            integer(int_kind     ), intent(in) :: thread_number
            integer(int_kind     ), intent(in) :: local_index
            integer(int_kind     )             :: index
        end function get_cell_index_interface
    end interface
end module abstract_parallelizer