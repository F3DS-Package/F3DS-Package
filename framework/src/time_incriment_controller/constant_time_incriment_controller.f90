module class_constant_time_incriment_controller
    use typedef_module
    use abstract_time_incriment_controller
    use abstract_configuration
    use abstract_eos
    use stdio_module

    implicit none

    private

    type, public, extends(time_incriment_controller) :: constant_time_incriment_controller
        private

        real(real_kind) :: global_dt

        contains

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: update_global
    end type constant_time_incriment_controller

    contains

    subroutine initialize(self, config)
        class(constant_time_incriment_controller), intent(inout) :: self
        class(configuration                     ), intent(inout) :: config

        logical :: found

        call config%get_real("Time incriment control.dt", self%global_dt, found)
        if(.not. found) call call_error("'Time incriment control.dt' is not found in configuration file you set.")
    end subroutine initialize

    pure function update_global(self, an_eos, variables_set, cell_volumes, num_cells) result(dt)
        class  (constant_time_incriment_controller), intent(in) :: self
        class  (eos                               ), intent(in) :: an_eos
        real   (real_kind                         ), intent(in) :: variables_set(:,:)
        real   (real_kind                         ), intent(in) :: cell_volumes (:)
        integer(int_kind                          ), intent(in) :: num_cells
        real   (real_kind                         )             :: dt

        !interface
        !    pure function spectral_radius_function(an_eos, variables) result(radius)
        !        use abstract_eos
        !        use typedef_module
        !        class  (eos      ), intent(in) :: an_eos
        !        real   (real_kind), intent(in) :: variables(:)
        !        real   (real_kind)             :: radius
        !    end function spectral_radius_function
        !end interface

        dt = self%global_dt
    end function update_global
end module class_constant_time_incriment_controller