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
        procedure, public, pass(self) :: compute_next_stage
        procedure, public, pass(self) :: prepare_time_stepping
        procedure, public, pass(self) :: get_number_of_stages
    end type forward_euler

    contains

    subroutine initialize(self, config, num_cells, num_conservative_variables)
        class  (forward_euler), intent(inout) :: self
        class  (configuration     ), intent(in   ) :: config
        integer(int_kind          ), intent(in   ) :: num_cells
        integer(int_kind          ), intent(in   ) :: num_conservative_variables

        return
    end subroutine initialize

    pure function compute_next_stage(   &
        self                          , &
        cell_index                    , &
        stage_num                     , &
        time_increment                , &
        conservative_variables        , &
        residuals                         ) result(updated_conservative_variables)

        class  (forward_euler      ), intent(in) :: self
        integer(int_kind           ), intent(in) :: cell_index
        integer(int_kind           ), intent(in) :: stage_num
        real   (real_kind          ), intent(in) :: time_increment
        real   (real_kind          ), intent(in) :: conservative_variables        (:)
        real   (real_kind          ), intent(in) :: residuals                     (:)
        real   (real_kind          )             :: updated_conservative_variables(1:size(conservative_variables))

        updated_conservative_variables(:) = conservative_variables(:) + time_increment * residuals(:)
    end function compute_next_stage

    subroutine prepare_time_stepping(   &
        self                          , &
        cell_index                    , &
        conservative_variables        , &
        residuals                         )

        class  (forward_euler), intent(inout) :: self
        integer(int_kind     ), intent(in   ) :: cell_index
        real   (real_kind    ), intent(in   ) :: conservative_variables(:)
        real   (real_kind    ), intent(in   ) :: residuals             (:)

        return
    end subroutine prepare_time_stepping

    pure function get_number_of_stages(self) result(n)
        class  (forward_euler), intent(in) :: self
        integer(int_kind          )              :: n
        n = nmu_stage_
    end function get_number_of_stages
end module class_forward_euler