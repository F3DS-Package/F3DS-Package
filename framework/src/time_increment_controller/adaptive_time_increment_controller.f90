module class_adaptive_time_increment_controller
    use typedef_module
    use abstract_time_increment_controller
    use abstract_configuration
    use abstract_eos
    use stdio_module

    implicit none

    private

    type, public, extends(time_increment_controller) :: adaptive_time_increment_controller
        private

        real(real_kind) :: courant_number

        contains

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: get_constant_dt
        procedure, public, pass(self) :: compute_local_dt
        procedure, public, pass(self) :: returns_constant
    end type adaptive_time_increment_controller

    contains

    subroutine initialize(self, config)
        class(adaptive_time_increment_controller), intent(inout) :: self
        class(configuration                     ), intent(inout) :: config

        logical :: found

        call config%get_real("Time increment control.Adaptive controller.Courant number", self%courant_number, found)
        if(.not. found) call call_error("'Time increment control.Adaptive controller.Courant number' is not found in configuration file you set.")
    end subroutine initialize

    function get_constant_dt(self) result(dt)
        class  (adaptive_time_increment_controller), intent(in) :: self
        real   (real_kind                         )             :: dt

        call call_error("Adaptive time increment controller does not return a constant dt. Please check your code.")

        dt = 0.d0
    end function get_constant_dt

    pure function compute_local_dt(self, cell_volume, spectral_radius) result(dt)
        class  (adaptive_time_increment_controller), intent(in) :: self
        real   (real_kind                         ), intent(in) :: cell_volume
        real   (real_kind                         ), intent(in) :: spectral_radius
        real   (real_kind                         )             :: dt

        dt = self%courant_number * cell_volume / spectral_radius
    end function compute_local_dt

    pure function returns_constant(self) result(ret)
        class  (adaptive_time_increment_controller), intent(in) :: self
        logical                                                 :: ret

        ret = .false.
    end function returns_constant
end module class_adaptive_time_increment_controller