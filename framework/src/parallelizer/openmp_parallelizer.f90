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
        integer, private :: num_threads = 1

        contains

        procedure, public, pass(self) :: initialize
    end type openmp_parallelizer

    contains

    subroutine initialize(self, config)
        class  (openmp_parallelizer), intent(inout) :: self
        class  (configuration      ), intent(inout) :: config

        logical           :: found

#ifdef _OPENMP
        call config%get_int                ("Parallel computing.Number of threads", self%num_threads, found, 1)
        if(.not. found) call write_warring("'Parallel computing.Number of threads' is not found in configuration file you set. Solver is executed a single thread.")

        call omp_set_num_threads(self%num_threads)
#endif
    end subroutine initialize
end module class_openmp_parallelizer