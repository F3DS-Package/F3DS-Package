module abstract_time_increment_controller
    implicit none

    private

    type, public, abstract :: time_increment_controller
        contains

        procedure(initialize_interface   ), pass(self), deferred :: initialize
        procedure(update_global_interface), pass(self), deferred :: update_global
    end type time_increment_controller

    abstract interface
        subroutine initialize_interface(self, config)
            use abstract_configuration
            import time_increment_controller

            class(time_increment_controller), intent(inout) :: self
            class(configuration            ), intent(inout) :: config
        end subroutine initialize_interface

        pure function update_global_interface(self, an_eos, variables_set, cell_volumes, num_cells) result(dt)
            use typedef_module
            use abstract_eos
            import time_increment_controller

            class  (time_increment_controller), intent(in) :: self
            class  (eos                      ), intent(in) :: an_eos
            real   (real_kind                ), intent(in) :: variables_set(:,:)
            real   (real_kind                ), intent(in) :: cell_volumes (:)
            integer(int_kind                 ), intent(in) :: num_cells
            real   (real_kind                )             :: dt

            !interface
            !    pure function spectral_radius_function(an_eos, variables) result(radius)
            !        use typedef_module
            !        use abstract_eos
            !        class  (eos      ), intent(in) :: an_eos
            !        real   (real_kind), intent(in) :: variables(:)
            !        real   (real_kind)             :: radius
            !    end function spectral_radius_function
            !end interface
        end function update_global_interface
    end interface
end module abstract_time_increment_controller