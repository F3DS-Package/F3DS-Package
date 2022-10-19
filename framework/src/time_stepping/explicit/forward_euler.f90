module class_forward_euler
    use typedef_module
    use abstract_time_stepping
    use abstract_configuration
    use abstract_eos

    implicit none

    private

    integer(int_kind ), parameter :: nmu_stage_ = 1

    type, public, extends(time_stepping) :: forward_euler
        private

        contains

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: compute_next_state
        procedure, public, pass(self) :: prepare_stepping
        procedure, public, pass(self) :: get_number_of_states
    end type forward_euler

    contains

    subroutine initialize(self, config, num_cells, num_conservative_variables)
        class  (forward_euler), intent(inout) :: self
        class  (configuration     ), intent(in   ) :: config
        integer(int_kind          ), intent(in   ) :: num_cells
        integer(int_kind          ), intent(in   ) :: num_conservative_variables

        return
    end subroutine initialize

    subroutine compute_next_state(   &
        self                       , &
        cell_index                 , &
        state_num                  , &
        time_increment             , &
        conservative_variables     , &
        residuals                      )

        class  (forward_euler), intent(inout) :: self
        integer(int_kind           ), intent(in   ) :: cell_index
        integer(int_kind           ), intent(in   ) :: state_num
        real   (real_kind          ), intent(in   ) :: time_increment
        real   (real_kind          ), intent(inout) :: conservative_variables(:)
        real   (real_kind          ), intent(inout) :: residuals             (:)

        conservative_variables(:) = conservative_variables(:) + time_increment * residuals(:)
        residuals             (:) = 0.d0
    end subroutine compute_next_state

    subroutine prepare_stepping(   &
        self                     , &
        cell_index               , &
        conservative_variables   , &
        primitive_variables      , &
        residuals                    )

        class  (forward_euler), intent(inout) :: self
        integer(int_kind           ), intent(in   ) :: cell_index
        real   (real_kind          ), intent(inout) :: conservative_variables(:)
        real   (real_kind          ), intent(inout) :: primitive_variables   (:)
        real   (real_kind          ), intent(inout) :: residuals             (:)

        return
    end subroutine prepare_stepping

    pure function get_number_of_states(self) result(n)
        class  (forward_euler), intent(in) :: self
        integer(int_kind          )              :: n
        n = nmu_stage_
    end function get_number_of_states
end module class_forward_euler