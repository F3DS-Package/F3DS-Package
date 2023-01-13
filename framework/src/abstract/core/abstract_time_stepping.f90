module abstract_time_stepping
    implicit none

    private

    type, public, abstract :: time_stepping
        contains

        procedure(initialize_interface          ), pass(self), deferred :: initialize
        procedure(compute_next_stage_interface  ), pass(self), deferred :: compute_next_stage
        procedure(prepare_stepping_interface    ), pass(self), deferred :: prepare_stepping
        !procedure(cleanup_stepping_interface    ), pass(self), deferred :: cleanup_stepping
        procedure(get_number_of_stages_interface), pass(self), deferred :: get_number_of_stages
    end type

    abstract interface
        subroutine initialize_interface(self, config, num_cells, num_conservative_variables)
            use abstract_configuration
            use typedef_module
            import time_stepping
            class  (time_stepping), intent(inout) :: self
            class  (configuration), intent(in   ) :: config
            integer(int_kind     ), intent(in   ) :: num_cells
            integer(int_kind     ), intent(in   ) :: num_conservative_variables
        end subroutine initialize_interface

        subroutine compute_next_stage_interface( &
            self                               , &
            cell_index                         , &
            stage_num                          , &
            time_increment                     , &
            conservative_variables             , &
            residuals                              )

            use typedef_module
            use abstract_eos
            import time_stepping

            class  (time_stepping      ), intent(inout) :: self
            integer(int_kind           ), intent(in   ) :: cell_index
            integer(int_kind           ), intent(in   ) :: stage_num
            real   (real_kind          ), intent(in   ) :: time_increment
            real   (real_kind          ), intent(inout) :: conservative_variables(:)
            real   (real_kind          ), intent(inout) :: residuals             (:)
        end subroutine compute_next_stage_interface

        subroutine prepare_stepping_interface(   &
            self                               , &
            cell_index                         , &
            conservative_variables             , &
            primitive_variables                , &
            residuals                              )
            use typedef_module
            import time_stepping
            class  (time_stepping), intent(inout) :: self
            integer(int_kind     ), intent(in   ) :: cell_index
            real   (real_kind    ), intent(inout) :: conservative_variables(:)
            real   (real_kind    ), intent(inout) :: primitive_variables   (:)
            real   (real_kind    ), intent(inout) :: residuals             (:)
        end subroutine prepare_stepping_interface

        pure function get_number_of_stages_interface(self) result(n)
            use typedef_module
            import time_stepping

            class  (time_stepping), intent(in   ) :: self
            integer(int_kind     )                :: n
        end function get_number_of_stages_interface
    end interface
end module abstract_time_stepping