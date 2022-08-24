module abstract_time_incriment_controller
    implicit none

    private

    type, public, abstract :: time_incriment_controller
        contains

        procedure(initialize_interface), pass(self), deferred :: initialize
        procedure(update_global_interface), pass(self), deferred :: update_global
    end type time_incriment_controller

    abstract interface
        subroutine initialize_interface(self, config)
            use abstract_configuration
            import time_incriment_controller

            class(time_incriment_controller), intent(inout) :: self
            class(configuration            ), intent(inout) :: config
        end subroutine initialize_interface

        pure function update_global_interface(self, an_eos, variables_set, cell_volumes, num_cells, spectral_radius_function) result(dt)
            use abstract_eos
            use typedef_module
            import time_incriment_controller

            class  (time_incriment_controller), intent(in) :: self
            class  (eos                      ), intent(in) :: an_eos
            real   (real_kind                ), intent(in) :: variables_set(:,:)
            real   (real_kind                ), intent(in) :: cell_volumes (:)
            integer(int_kind                 ), intent(in) :: num_cells

            real   (real_kind                )             :: dt

            interface
                pure function spectral_radius_function(an_eos, variables)
                    use abstract_eos
                    use typedef_module
                    class  (eos      ) :: an_eos
                    real   (real_kind) :: variables(:)
                end function spectral_radius_function
            end interface
        end function update_global_interface
    end interface
end module abstract_time_incriment_controller