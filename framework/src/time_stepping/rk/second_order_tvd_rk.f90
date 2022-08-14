module class_second_order_tvd_rk
    use typedef_module
    use abstract_configuration
    use abstract_eos

    implicit none

    private

    type, public, extends(time_stepping) :: second_order_tvd_rk
        private

        integer(int_kind ), parameter :: nmu_stage_ = 2
        real   (real_kind), parameter :: alpha_(1:3) = [real(real_kind) :: 1, 1/2]
        real   (real_kind), parameter :: beta_ (1:3) = [real(real_kind) :: 1, 1/2]

        ! Stored init conservative-variables of each in cells
        ! Elm. 1) 1:number of conservative variables
        ! Elm. 2) 1:number of cells
        real(real_kind), allocatable :: init_conservative_variables_set(:,:)

        contains

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: compute_next_state
        procedure, public, pass(self) :: prepare_stepping
        procedure, public, pass(self) :: get_number_of_states
    end type second_order_tvd_rk

    contains

    subroutine initialize(self, config, num_cells, num_conservative_variables)
        class  (second_order_tvd_rk), intent(inout) :: self
        class  (configuration     ), intent(in   ) :: config
        integer(int_kind          ), intent(in   ) :: num_cells
        integer(int_kind          ), intent(in   ) :: num_conservative_variables

        allocate(init_conservative_variables_set(1:n_conservative_variables, 1:n_cells))
    end subroutine initialize

    subroutine compute_next_state(   &
        self                       , &
        state_num                  , &
        time_increment             , &
        conservative_variables     , &
        residuals                      )

        class  (second_order_tvd_rk), intent(inout) :: self
        integer(int_kind           ), intent(in   ) :: cell_index
        integer(int_kind           ), intent(in   ) :: state_num
        real   (real_kind          ), intent(in   ) :: time_increment
        real   (real_kind          ), intent(inout) :: conservative_variables(:)
        real   (real_kind          ), intent(inout) :: residuals             (:)


        conservative_variables_set(:) = self%alpha_(state_num) * self%init_conservative_variables_set(:, cell_index) &
                                      + self%beta_ (state_num) * (conservative_variables_set(:) + time_increment * residual_set(:))
        residual_set              (:) = 0.d0
    end subroutine compute_next_state

    subroutine prepare_stepping(   &
        self                     , &
        cell_index               , &
        conservative_variables   , &
        primitive_variables      , &
        residuals                    )

        class  (second_order_tvd_rk), intent(inout) :: self
        integer(int_kind           ), intent(in   ) :: cell_index
        real   (real_kind          ), intent(inout) :: conservative_variables(:)
        real   (real_kind          ), intent(inout) :: primitive_variables   (:)
        real   (real_kind          ), intent(inout) :: residuals             (:)

        self%init_conservative_variables_set(:, cell_index) = conservative_variables_set(:)
    end subroutine prepare_stepping

    pure function get_number_of_states(self) result(n)
        class  (second_order_tvd_rk), intent(inout) :: self
        integer(int_kind          )                :: n
        n = self%nmu_stage_
    end function get_number_of_states
end module class_second_order_tvd_rk