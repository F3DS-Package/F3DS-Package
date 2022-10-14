module class_constant_time_increment_controller
    use typedef_module
    use abstract_time_increment_controller
    use abstract_configuration
    use abstract_eos
    use stdio_module

    implicit none

    private

    type, public, extends(time_increment_controller) :: constant_time_increment_controller
        private

        real(real_kind) :: constant_dt

        contains

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: get_constant_dt
        procedure, public, pass(self) :: compute_local_dt
        procedure, public, pass(self) :: returns_constant
    end type constant_time_increment_controller

    contains

    subroutine initialize(self, config)
        class(constant_time_increment_controller), intent(inout) :: self
        class(configuration                     ), intent(inout) :: config

        logical :: found

        call config%get_real("Time increment control.dt", self%constant_dt, found)
        if(.not. found) call call_error("'Time increment control.dt' is not found in configuration file you set.")
    end subroutine initialize

    pure function get_constant_dt(self) result(dt)
        class  (constant_time_increment_controller), intent(in) :: self
        real   (real_kind                         )             :: dt

        dt = self%constant_dt
    end function get_constant_dt

    pure function compute_local_dt(self, cell_volume, spectral_radius) result(dt)
        class  (constant_time_increment_controller), intent(in) :: self
        real   (real_kind                         ), intent(in) :: cell_volume
        real   (real_kind                         ), intent(in) :: spectral_radius
        real   (real_kind                         )             :: dt

        dt = self%constant_dt
    end function compute_local_dt

    pure function returns_constant(self) result(ret)
        class  (constant_time_increment_controller), intent(in) :: self
        logical                                                 :: ret

        ret = .true.
    end function returns_constant
end module class_constant_time_increment_controller