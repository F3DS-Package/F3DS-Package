module class_openmp_parallelizer
    use typedef_module
    use abstract_parallelizer
    use abstract_configuration
    use stdio_module
    use face_type_module
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

        ! store number of face index in each theads
        ! elm. 1: face type (see utils/cellsystem/face_type.f90)
        ! elm. 2: thiread number i ∈ [1, {@code num_threads}].
        integer(int_kind), allocatable :: num_boundary_face_indexes(:,:)

        ! store cell index at first.
        ! range of elm. 1 is i ∈ [1, {@code num_threads}].
        integer(int_kind), allocatable :: start_cell_index(:)

        ! store face index at first.
        ! elm. 1: thiread number i ∈ [1, {@code num_threads}].
        integer(int_kind), allocatable :: start_face_index(:)

        ! store boundary face index at first.
        ! elm. 1: face type (see utils/cellsystem/face_type.f90)
        ! elm. 2: thiread number i ∈ [1, {@code num_threads}].
        integer(int_kind), allocatable :: start_boundary_face_index(:,:)

        contains

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: get_node_number
        procedure, public, pass(self) :: get_number_of_threads
        procedure, public, pass(self) :: get_number_of_face_indexes
        procedure, public, pass(self) :: get_number_of_cell_indexes
        procedure, public, pass(self) :: get_number_of_boundary_face_indexes
        procedure, public, pass(self) :: get_face_index
        procedure, public, pass(self) :: get_cell_index
        procedure, public, pass(self) :: get_boundary_face_index

        procedure, pass(self) :: make_index_list
    end type openmp_parallelizer

    contains

    subroutine make_index_list(self, num_indexes, start_indexes, num_objs)
        class  (openmp_parallelizer), intent(in   ) :: self
        integer(int_kind           ), intent(inout) :: num_indexes  (:)
        integer(int_kind           ), intent(inout) :: start_indexes(:)
        integer(int_kind           ), intent(in   ) :: num_objs

        integer(int_kind) :: remainder, thread, local_index, start_index

        num_indexes(:) = num_objs / self%num_threads
        remainder = mod(num_objs, self%num_threads)
        if (remainder > 0) then
            do thread = 1, remainder, 1
                num_indexes(thread) = num_indexes(thread) + 1
            end do
        end if
        start_index = 1
        do thread = 1, self%num_threads, 1
            start_indexes(thread) = start_index
            start_index = start_index + num_indexes(thread)
        end do
    end subroutine make_index_list

    subroutine initialize(self, config, num_cells, num_faces, &
                          num_outflow_faces, num_wall_faces, num_symmetric_faces, num_empty_faces)
        class  (openmp_parallelizer), intent(inout) :: self
        class  (configuration      ), intent(inout) :: config
        integer(int_kind           ), intent(in   ) :: num_cells
        integer(int_kind           ), intent(in   ) :: num_faces
        integer(int_kind           ), intent(in   ) :: num_outflow_faces
        integer(int_kind           ), intent(in   ) :: num_wall_faces
        integer(int_kind           ), intent(in   ) :: num_symmetric_faces
        integer(int_kind           ), intent(in   ) :: num_empty_faces

        logical           :: found

#ifdef _OPENMP
        call config%get_int                ("Parallel computing.Number of threads", self%num_threads, found, 1)
        if(.not. found) call write_warring("'Parallel computing.Number of threads' is not found in configuration file you set. Solver is executed a single thread.")

        call omp_set_num_threads(self%num_threads)
#endif

        allocate(self%num_cell_indexes(self%num_threads))
        allocate(self%start_cell_index(self%num_threads))
        call self%make_index_list(self%num_cell_indexes, self%start_cell_index, num_cells)

        allocate(self%num_face_indexes(self%num_threads))
        allocate(self%start_face_index(self%num_threads))
        call self%make_index_list(self%num_face_indexes, self%start_face_index, num_faces)

        allocate(self%num_boundary_face_indexes(number_of_boundary_face_type, self%num_threads))
        allocate(self%start_boundary_face_index(number_of_boundary_face_type, self%num_threads))
        call self%make_index_list(self%num_boundary_face_indexes(outflow_face_type  , :), self%start_boundary_face_index(outflow_face_type  , :), num_outflow_faces  )
        call self%make_index_list(self%num_boundary_face_indexes(symmetric_face_type, :), self%start_boundary_face_index(symmetric_face_type, :), num_wall_faces     )
        call self%make_index_list(self%num_boundary_face_indexes(wall_face_type     , :), self%start_boundary_face_index(wall_face_type     , :), num_symmetric_faces)
        call self%make_index_list(self%num_boundary_face_indexes(empty_face_type    , :), self%start_boundary_face_index(empty_face_type    , :), num_empty_faces    )

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

    pure function get_number_of_cell_indexes(self, node_number, thread_number) result(n)
        class  (openmp_parallelizer), intent(in) :: self
        integer(int_kind           ), intent(in) :: node_number
        integer(int_kind           ), intent(in) :: thread_number
        integer(int_kind           )             :: n

        n = self%num_cell_indexes(thread_number)
    end function get_number_of_cell_indexes

    pure function get_number_of_face_indexes(self, node_number, thread_number) result(n)
        class  (openmp_parallelizer), intent(in) :: self
        integer(int_kind           ), intent(in) :: node_number
        integer(int_kind           ), intent(in) :: thread_number
        integer(int_kind           )             :: n

        n = self%num_face_indexes(thread_number)
    end function get_number_of_face_indexes

    pure function get_number_of_boundary_face_indexes(self, node_number, thread_number, a_face_type) result(n)
        class  (openmp_parallelizer    ), intent(in) :: self
        integer(int_kind               ), intent(in) :: node_number
        integer(int_kind               ), intent(in) :: thread_number
        integer(boundary_face_type_kind), intent(in) :: a_face_type
        integer(int_kind               )             :: n

        n = self%num_boundary_face_indexes(a_face_type, thread_number)
    end function get_number_of_boundary_face_indexes

    pure function get_cell_index(self, node_number, thread_number, local_index) result(index)
        class  (openmp_parallelizer), intent(in) :: self
        integer(int_kind           ), intent(in) :: node_number
        integer(int_kind           ), intent(in) :: thread_number
        integer(int_kind           ), intent(in) :: local_index
        integer(int_kind           )             :: index

        index = self%start_cell_index(thread_number) + (local_index - 1)
    end function get_cell_index

    pure function get_face_index(self, node_number, thread_number, local_index) result(index)
        class  (openmp_parallelizer), intent(in) :: self
        integer(int_kind           ), intent(in) :: node_number
        integer(int_kind           ), intent(in) :: thread_number
        integer(int_kind           ), intent(in) :: local_index
        integer(int_kind           )             :: index

        index = self%start_face_index(thread_number) + (local_index - 1)
    end function get_face_index

    pure function get_boundary_face_index(self, node_number, thread_number, local_index, a_boundary_face_type) result(index)
        class  (openmp_parallelizer    ), intent(in) :: self
        integer(int_kind               ), intent(in) :: node_number
        integer(int_kind               ), intent(in) :: thread_number
        integer(int_kind               ), intent(in) :: local_index
        integer(boundary_face_type_kind), intent(in) :: a_boundary_face_type
        integer(int_kind               )             :: index

        index = self%start_boundary_face_index(a_boundary_face_type, thread_number) + (local_index - 1)
    end function get_boundary_face_index
end module class_openmp_parallelizer