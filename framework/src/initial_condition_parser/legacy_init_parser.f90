module class_legacy_init_parser
    use abstract_initial_condition_parser
    use typedef_module
    use stdio_module
    use grid_utils_module
    use abstract_configuration

    implicit none

    private

    type, public, extends(initial_condition_parser) :: legacy_init_parser
        private

        integer(int_kind) :: num_ghost_cells_ = 3

        real(real_kind), allocatable :: conservative_variables_set_(:,:,:,:)
        real(real_kind), allocatable :: primitive_variables_set_(:,:,:,:)

        integer(int_kind) :: imax, jmax, kmax, imin, jmin, kmin
        integer(int_kind) :: n_conservative_variables, n_primitive_variables
        logical :: parsed = .false.
        logical :: nlinit_extension = .false.

        contains

        procedure, public, pass(self) :: parse
        procedure, public, pass(self) :: close
        procedure, public, pass(self) :: get_conservative_variables_set
    end type legacy_init_parser

    contains

    subroutine parse(self, config)
        class(legacy_init_parser), intent(inout) :: self
        class(configuration     ), intent(inout) :: config

        integer  (int_kind) :: i,j,k,n
        integer  (int_kind) :: unit_number
        logical             :: found
        character(len=:), allocatable :: filepath

        if (self%parsed) then
            call call_error("'parse' method of legacy_init_parser is already called. But you call 'parse' method.")
        end if

        call config%get_bool("Grid.NL Extension", self%nlinit_extension, found, .false.)

        call config%get_char("Initial condition.Filepath", filepath, found)
        if(.not. found) call call_error("'Initial condition.Filepath' is not found in configuration file you set.")

        open(newunit=unit_number, file=filepath, access = 'stream', form = 'unformatted', status = 'old')
        read(unit_number) self%imin, self%jmin, self%kmin
        read(unit_number) self%imax, self%jmax, self%kmax
        if (self%nlinit_extension) then
            read(unit_number) self%n_primitive_variables, self%n_conservative_variables
        else
            read(unit_number) self%n_conservative_variables
        endif

        if (self%nlinit_extension) ALLOCATE (self%primitive_variables_set_   (self%imin-2:self%imax+2, self%jmin-2:self%jmax+2, self%kmin-2:self%kmax+2, 1:self%n_primitive_variables   ) )
        if (self%nlinit_extension) read(unit_number) ((((self%primitive_variables_set_   (i,j,k,n), i = self%imin-2, self%imax+2), j = self%jmin-2, self%jmax+2), k = self%kmin-2, self%kmax+2), n = 1, self%n_primitive_variables   )

        ALLOCATE (self%conservative_variables_set_(self%imin-2:self%imax+2, self%jmin-2:self%jmax+2, self%kmin-2:self%kmax+2, 1:self%n_conservative_variables) )
        read(unit_number) ((((self%conservative_variables_set_(i,j,k,n), i = self%imin-2, self%imax+2), j = self%jmin-2, self%jmax+2), k = self%kmin-2, self%kmax+2), n = 1, self%n_conservative_variables)

        close(unit_number)

        self%parsed = .true.
    end subroutine parse

    subroutine close(self)
        class    (legacy_init_parser), intent(inout) :: self

        if (.not. self%parsed) then
            call call_error("'parse' method of legacy_init_parser is not called. But you call 'closs' method.")
        end if

        deallocate(self%conservative_variables_set_)
        deallocate(self%primitive_variables_set_)
        self%parsed = .false.
    end subroutine close

    subroutine get_conservative_variables_set(self, conservative_variables_set)
        class(legacy_init_parser), intent(inout) :: self
        real (real_kind    ), intent(inout) :: conservative_variables_set(:,:)

        integer(int_kind) :: i,j,k,n
        integer(int_kind) :: variables_index

        if (.not. self%parsed) then
            call call_error("'parse' method of legacy_init_parser is not called. But you call 'get_conservative_variables_set' method.")
        end if

        if(.not. (size(conservative_variables_set(:, 1)) == self%n_conservative_variables))then
            call call_error("Conservative-variables size you parse is not matched.")
        end if

        do k = self%kmin, self%kmax, 1
            do j = self%jmin, self%jmax, 1
                do i = self%imin, self%imax, 1
                    n = convert_structure_index_to_unstructure_index(i, j, k,                            &
                        self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
                        self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
                    do variables_index = 1, self%n_conservative_variables, 1
                        conservative_variables_set(variables_index, n) = self%conservative_variables_set_(i, j, k, variables_index)
                    end do
                end do
            end do
        end do
    end subroutine get_conservative_variables_set
end module class_legacy_init_parser