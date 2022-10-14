module abstract_time_increment_controller
    implicit none

    private

    type, public, abstract :: time_increment_controller
        contains

        procedure(initialize_interface      ), pass(self), deferred :: initialize
        procedure(get_constant_dt_interface ), pass(self), deferred :: get_constant_dt
        procedure(compute_local_dt_interface), pass(self), deferred :: compute_local_dt
        procedure(returns_constant_interface), pass(self), deferred :: returns_constant
    end type time_increment_controller

    abstract interface
        subroutine initialize_interface(self, config)
            use abstract_configuration
            import time_increment_controller

            class(time_increment_controller), intent(inout) :: self
            class(configuration            ), intent(inout) :: config
        end subroutine initialize_interface

        pure function get_constant_dt_interface(self) result(dt)
            use typedef_module
            import time_increment_controller

            class  (time_increment_controller), intent(in) :: self
            real   (real_kind                )             :: dt
        end function get_constant_dt_interface

        pure function compute_local_dt_interface(self, cell_volume, spectral_radius) result(dt)
            use typedef_module
            import time_increment_controller

            class  (time_increment_controller), intent(in) :: self
            real   (real_kind                ), intent(in) :: cell_volume
            real   (real_kind                ), intent(in) :: spectral_radius
            real   (real_kind                )             :: dt
        end function compute_local_dt_interface

        pure function returns_constant_interface(self) result(ret)
            import time_increment_controller

            class  (time_increment_controller), intent(in) :: self
            logical                                        :: ret
        end function returns_constant_interface
    end interface
end module abstract_time_increment_controller