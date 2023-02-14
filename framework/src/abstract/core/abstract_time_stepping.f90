module abstract_time_stepping
    implicit none

    private

    type, public, abstract :: time_stepping_generator
        contains
        generic, public :: generate => generate_from_configuration, generate_from_name
        procedure(generate_from_configuration_interface), pass(self), deferred :: generate_from_configuration
        procedure(generate_from_name_interface         ), pass(self), deferred :: generate_from_name
    end type time_stepping_generator

    type, public, abstract :: time_stepping
        contains

        procedure(initialize_interface           ), pass(self), deferred :: initialize
        procedure(compute_next_stage_interface   ), pass(self), deferred :: compute_next_stage
        procedure(prepare_time_stepping_interface), pass(self), deferred :: prepare_time_stepping
        !procedure(cleanup_stepping_interface    ), pass(self), deferred :: cleanup_stepping
        procedure(get_number_of_stages_interface ), pass(self), deferred :: get_number_of_stages
    end type time_stepping

    abstract interface
        subroutine generate_from_configuration_interface(self, a_time_stepping, a_config)
            use abstract_configuration
            import time_stepping_generator
            import time_stepping
            class(time_stepping_generator),          intent(inout) :: self
            class(time_stepping          ), pointer, intent(inout) :: a_time_stepping
            class(configuration          ),          intent(inout) :: a_config
        end subroutine generate_from_configuration_interface

        subroutine generate_from_name_interface(self, a_time_stepping, name)
            use abstract_configuration
            import time_stepping_generator
            import time_stepping
            class    (time_stepping_generator),              intent(inout) :: self
            class    (time_stepping          ), pointer    , intent(inout) :: a_time_stepping
            character(len=:                  ), allocatable, intent(in   ) :: name
        end subroutine generate_from_name_interface
    end interface

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

        pure function compute_next_stage_interface( &
            self                                  , &
            cell_index                            , &
            stage_num                             , &
            time_increment                        , &
            conservative_variables                , &
            residuals                                 ) result(updated_conservative_variables)

            use typedef_module
            import time_stepping

            class  (time_stepping      ), intent(in) :: self
            integer(int_kind           ), intent(in) :: cell_index
            integer(int_kind           ), intent(in) :: stage_num
            real   (real_kind          ), intent(in) :: time_increment
            real   (real_kind          ), intent(in) :: conservative_variables        (:)
            real   (real_kind          ), intent(in) :: residuals                     (:)
            real   (real_kind          )             :: updated_conservative_variables(1:size(conservative_variables))
        end function compute_next_stage_interface

        subroutine prepare_time_stepping_interface(   &
            self                               , &
            cell_index                         , &
            conservative_variables             , &
            residuals                              )
            use typedef_module
            import time_stepping
            class  (time_stepping), intent(inout) :: self
            integer(int_kind     ), intent(in   ) :: cell_index
            real   (real_kind    ), intent(inout) :: conservative_variables(:)
            real   (real_kind    ), intent(inout) :: residuals             (:)
        end subroutine prepare_time_stepping_interface

        pure function get_number_of_stages_interface(self) result(n)
            use typedef_module
            import time_stepping

            class  (time_stepping), intent(in   ) :: self
            integer(int_kind     )                :: n
        end function get_number_of_stages_interface
    end interface
end module abstract_time_stepping