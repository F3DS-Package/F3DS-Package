module abstract_parallelizer
    implicit none

    private

    type, public, abstract :: parallelizer
        contains
        procedure(initialize_interface), pass(self), deferred :: initialize
    end type parallelizer

    abstract interface
        subroutine initialize_interface(self, config)
            use abstract_configuration
            import parallelizer

            class(parallelizer ), intent(inout) :: self
            class(configuration), intent(inout) :: config
        end subroutine initialize_interface
    end interface
end module abstract_parallelizer