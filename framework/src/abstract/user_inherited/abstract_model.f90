module abstract_model
    implicit none

    private

    type, public, abstract :: model
        contains

        procedure(initialize_interface              ), pass(self), deferred :: initialize
        procedure(compute_residual_element_interface), pass(self), deferred :: compute_residual_element
        procedure(compute_source_term_interface     ), pass(self), deferred :: compute_source_term
        procedure(spectral_radius_interface         ), pass(self), deferred :: spectral_radius
    end type model

    abstract interface
        subroutine initialize_interface(self, a_configuration)
            use abstract_configuration
            import model

            class(model        ), intent(inout) :: self
            class(configuration), intent(inout) :: a_configuration
        end subroutine initialize_interface

        pure function compute_residual_element_interface(   &
            self                                          , &
            an_eos                                        , &
            an_riemann_solver                             , &
            primitive_variables_lhc                       , &
            primitive_variables_rhc                       , &
            reconstructed_primitive_variables_lhc         , &
            reconstructed_primitive_variables_rhc         , &
            face_gradient_primitive_variables             , &
            lhc_cell_volume                               , &
            rhc_cell_volume                               , &
            face_area                                     , &
            face_normal_vector                            , &
            face_tangential1_vector                       , &
            face_tangential2_vector                       , &
            num_conservative_values                       , &
            num_primitive_values                              ) result(residual_element)

            use typedef_module
            use abstract_eos
            use abstract_riemann_solver
            import model

            class  (model         ), intent(in) :: self
            class  (eos           ), intent(in) :: an_eos
            class  (riemann_solver), intent(in) :: an_riemann_solver
            real   (real_kind     ), intent(in) :: primitive_variables_lhc              (:)
            real   (real_kind     ), intent(in) :: primitive_variables_rhc              (:)
            real   (real_kind     ), intent(in) :: reconstructed_primitive_variables_lhc(:)
            real   (real_kind     ), intent(in) :: reconstructed_primitive_variables_rhc(:)
            real   (real_kind     ), intent(in) :: face_gradient_primitive_variables    (:)
            real   (real_kind     ), intent(in) :: lhc_cell_volume
            real   (real_kind     ), intent(in) :: rhc_cell_volume
            real   (real_kind     ), intent(in) :: face_area
            real   (real_kind     ), intent(in) :: face_normal_vector     (3)
            real   (real_kind     ), intent(in) :: face_tangential1_vector(3)
            real   (real_kind     ), intent(in) :: face_tangential2_vector(3)
            integer(int_kind      ), intent(in) :: num_conservative_values
            integer(int_kind      ), intent(in) :: num_primitive_values

            real   (real_kind)                  :: residual_element(num_conservative_values, 1:2)
        end function compute_residual_element_interface

        pure function compute_source_term_interface(self, variables, num_conservative_values) result(source)
            use typedef_module
            use abstract_eos
            import model

            class  (model    ), intent(in) :: self
            real   (real_kind), intent(in) :: variables(:)
            integer(int_kind ), intent(in) :: num_conservative_values
            real   (real_kind)             :: source(num_conservative_values)
        end function compute_source_term_interface

        pure function spectral_radius_interface(self, an_eos, primitive_variables, length) result(r)
            use typedef_module
            use abstract_eos
            import model

            class(model    ), intent(in) :: self
            class(eos      ), intent(in) :: an_eos
            real (real_kind), intent(in) :: primitive_variables(:)
            real (real_kind), intent(in) :: length
            real (real_kind) :: r
        end function spectral_radius_interface
    end interface
end module abstract_model