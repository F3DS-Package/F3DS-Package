module class_openmp_parallelizer
    use typedef_module
    use abstract_parallelizer
    use abstract_configuration
    use stdio_module
#ifdef _OPENMP
    use omp_lib
#endif

    implicit none

    private

    type, public, extends(parallelizer) :: openmp_parallelizer
        integer(int_kind) :: num_nodes   = 1
        integer(int_kind) :: num_threads = 1

        ! store number of cell index in each theads; the range is i ∈ [1, {@code num_threads}].
        integer(int_kind), allocatable :: num_cell_indexes(:)

        ! store number of face index in each theads; the range is i ∈ [1, {@code num_threads}].
        integer(int_kind), allocatable :: num_face_indexes(:)

        ! store cell index  at first.
        ! range of elm. 1 is i ∈ [1, {@code num_threads}].
        integer(int_kind), allocatable :: start_cell_index(:)

        ! store face index at first.
        ! range of elm. 1 is i ∈ [1, {@code num_threads}].
        integer(int_kind), allocatable :: start_face_index(:)

        contains

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: get_node_number
        procedure, public, pass(self) :: get_number_of_threads
        procedure, public, pass(self) :: get_number_of_face_indexes
        procedure, public, pass(self) :: get_number_of_cell_indexes
        procedure, public, pass(self) :: get_face_index
        procedure, public, pass(self) :: get_cell_index
    end type openmp_parallelizer

    contains

    subroutine initialize(self, config, number_of_cells, number_of_faces)
        class  (openmp_parallelizer), intent(inout) :: self
        class  (configuration      ), intent(inout) :: config
        integer(int_kind           ), intent(in   ) :: number_of_cells
        integer(int_kind           ), intent(in   ) :: number_of_faces

        logical           :: found
        integer(int_kind) :: remainder, thread, local_index, start_index

#ifdef _OPENMP
        call config%get_int                ("Parallel computing.Number of threads", self%num_threads, found, 1)
        if(.not. found) call write_warring("'Parallel computing.Number of threads' is not found in configuration file you set. Solver is executed a single thread.")

        call omp_set_num_threads(self%num_threads)
#endif

        allocate(self%num_cell_indexes(self%num_threads))
        allocate(self%num_face_indexes(self%num_threads))
        allocate(self%start_cell_index(self%num_threads))
        allocate(self%start_face_index(self%num_threads))

        ! decompose cell indexes
        self%num_cell_indexes(:) = number_of_cells / self%num_threads
        remainder = mod(number_of_cells, self%num_threads)
        if (remainder > 0) then
            do thread = 1, remainder, 1
                self%num_cell_indexes(thread) = self%num_cell_indexes(thread) + 1
            end do
        end if
        start_index = 1
        do thread = 1, self%num_threads, 1
            self%start_cell_index(thread) = start_index
            start_index = start_index + self%num_cell_indexes(thread)
        end do

        ! decompose face indexes
        self%num_face_indexes(:) = number_of_faces / self%num_threads
        remainder = mod(number_of_faces, self%num_threads)
        if (remainder > 0) then
            do thread = 1, remainder, 1
                self%num_face_indexes(thread) = self%num_face_indexes(thread) + 1
            end do
        end if
        start_index = 1
        do thread = 1, self%num_threads, 1
            self%start_face_index(thread) = start_index
            start_index = start_index + self%num_face_indexes(thread)
        end do
    end subroutine initialize

    pure function get_node_number(self) result(n)
        class  (openmp_parallelizer), intent(in) :: self
        integer(int_kind           )             :: n

        n = self%num_nodes
    end function get_node_number

    pure function get_number_of_threads(self, node_number) result(n)
        class  (openmp_parallelizer), intent(in) :: self
        integer(int_kind           ), intent(in) :: node_number
        integer(int_kind           )             :: n

        n = self%num_threads
    end function get_number_of_threads

    pure function get_number_of_face_indexes(self, node_number, thread_number) result(n)
        class  (openmp_parallelizer), intent(in) :: self
        integer(int_kind           ), intent(in) :: node_number
        integer(int_kind           ), intent(in) :: thread_number
        integer(int_kind           )             :: n

        n = self%num_face_indexes(thread_number)
    end function get_number_of_face_indexes

    pure function get_number_of_cell_indexes(self, node_number, thread_number) result(n)
        class  (openmp_parallelizer), intent(in) :: self
        integer(int_kind           ), intent(in) :: node_number
        integer(int_kind           ), intent(in) :: thread_number
        integer(int_kind           )             :: n

        n = self%num_cell_indexes(thread_number)
    end function get_number_of_cell_indexes

    pure function get_face_index(self, node_number, thread_number, local_index) result(index)
        class  (openmp_parallelizer), intent(in) :: self
        integer(int_kind           ), intent(in) :: node_number
        integer(int_kind           ), intent(in) :: thread_number
        integer(int_kind           ), intent(in) :: local_index
        integer(int_kind           )             :: index

        index = self%start_face_index(thread_number) + (local_index - 1)
    end function get_face_index

    pure function get_cell_index(self, node_number, thread_number, local_index) result(index)
        class  (openmp_parallelizer), intent(in) :: self
        integer(int_kind           ), intent(in) :: node_number
        integer(int_kind           ), intent(in) :: thread_number
        integer(int_kind           ), intent(in) :: local_index
        integer(int_kind           )             :: index

        index = self%start_cell_index(thread_number) + (local_index - 1)
    end function get_cell_index
end module class_openmp_parallelizer