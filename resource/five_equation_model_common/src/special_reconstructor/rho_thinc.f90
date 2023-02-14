module class_rho_thinc
    use typedef_module
    use abstract_reconstructor
    use abstract_configuration
    use stdio_module

    implicit none

    private
    type, public, extends(reconstructor) :: rho_thinc
        private

        real (real_kind    ) :: specified_slope_parameter_
        real (real_kind    ) :: epsilon_
        class(reconstructor), pointer :: primary_reconstructor_

        contains

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: reconstruct_lhc
        procedure, public, pass(self) :: reconstruct_rhc
    end type

    contains

    subroutine initialize(self, config, a_reconstructor_generator)
        class(rho_thinc              ),           intent(inout) :: self
        class(configuration          ),           intent(inout) :: config
        class(reconstructor_generator), optional, intent(inout) :: a_reconstructor_generator

        logical :: found
        character(len=:), allocatable :: name

        call config%get_real("Reconstructor.rho-THINC.Specified slope parameter", self%specified_slope_parameter_, found, 0.5d0)
        if(.not. found) call write_warring("'Reconstructor.rho-THINC.Specified slope parameter' is not found in configuration you set. To be set dafault value.")

        call config%get_real("Reconstructor.rho-THINC.Epsilon", self%epsilon_, found, 5.d-5)
        if(.not. found) call write_warring("'Reconstructor.rho-THINC.Epsilon' is not found in configuration you set. To be set dafault value.")

        call config%get_char("Reconstructor.rho-THINC.Primary reconstructor", name, found, "Minmod MUSCL")
        if(.not. found) call write_warring("'Reconstructor.rho-THINC.Primary reconstructor' is not found in configuration you set. To be set dafault value.")

        if(.not. present(a_reconstructor_generator))then
            call call_error("rho-THINC requests a reconstructor-generator.")
        end if

        call a_reconstructor_generator%generate(self%primary_reconstructor_, name)
    end subroutine initialize

    pure function reconstruct_lhc( &
        self                     , &
        primitive_variables_set  , &
        face_to_cell_index       , &
        cell_centor_positions    , &
        face_centor_positions    , &
        face_index               , &
        num_local_cells          , &
        num_primitive_variables      ) result(reconstructed_primitive_variables)

        class  (rho_thinc    ), intent(in) :: self
        integer(int_kind     ), intent(in) :: face_index
        real   (real_kind    ), intent(in) :: primitive_variables_set          (:, :)
        integer(int_kind     ), intent(in) :: face_to_cell_index               (:, :)
        real   (real_kind    ), intent(in) :: cell_centor_positions            (:, :)
        real   (real_kind    ), intent(in) :: face_centor_positions            (:, :)
        integer(int_kind     ), intent(in) :: num_local_cells
        integer(int_kind     ), intent(in) :: num_primitive_variables
        real   (real_kind    )             :: reconstructed_primitive_variables(num_primitive_variables)

        integer(int_kind )              :: lhc_index, rhc_index, llhc_index, rrhc_index

        real   (real_kind)              :: lhc_interface_location, lhc_sign, lhc_a, lhc_b, lhc_d, lhc_e
        logical                         :: lhc_is_monotonicity

        llhc_index = face_to_cell_index(num_local_cells-1, face_index)
        lhc_index  = face_to_cell_index(num_local_cells+0, face_index)
        rhc_index  = face_to_cell_index(num_local_cells+1, face_index)
        rrhc_index = face_to_cell_index(num_local_cells+2, face_index)

        reconstructed_primitive_variables(:) = self%primary_reconstructor_%reconstruct_lhc( &
            primitive_variables_set  , &
            face_to_cell_index       , &
            cell_centor_positions    , &
            face_centor_positions    , &
            face_index               , &
            num_local_cells          , &
            num_primitive_variables    &
        )

        lhc_sign  = sign(1.d0, primitive_variables_set(7, rhc_index) - primitive_variables_set(7, llhc_index))
        lhc_is_monotonicity = (primitive_variables_set(7, rhc_index) - primitive_variables_set(7, lhc_index)) * (primitive_variables_set(7, lhc_index) - primitive_variables_set(7,llhc_index)) > 0.d0
        associate(                                          &
            lhc_rho1 => primitive_variables_set(1, lhc_index), &
            lhc_rho2 => primitive_variables_set(2, lhc_index), &
            lhc_z1   => primitive_variables_set(7, lhc_index)  &
        )
            if (lhc_z1 < self%epsilon_) then
                reconstructed_primitive_variables(7) = 0.d0
            else if (lhc_z1 > 1.d0 - self%epsilon_) then
                reconstructed_primitive_variables(7) = 1.d0
            else if (lhc_is_monotonicity) then
                lhc_a    = exp (2.d0 * lhc_sign * self%specified_slope_parameter_)
                lhc_b    = exp (2.d0 * lhc_sign * self%specified_slope_parameter_ * lhc_z1)
                lhc_interface_location = 1.d0 / (2.d0 * self%specified_slope_parameter_) * log((lhc_b - 1.d0) / (lhc_a - lhc_b))
                lhc_d = exp(2.d0 * self%specified_slope_parameter_ * lhc_interface_location)
                lhc_e = log((lhc_a * lhc_d + 1.d0)**2.d0 / (lhc_a * (lhc_d + 1.d0)**2.d0))
                reconstructed_primitive_variables(1) =  (4.d0 * lhc_sign * self%specified_slope_parameter_ * lhc_rho1 * lhc_z1) &
                                                     /  (lhc_e + 2.d0 * lhc_sign * self%specified_slope_parameter_)
                reconstructed_primitive_variables(2) = -(4.d0 * lhc_sign * self%specified_slope_parameter_ * lhc_rho2 * (1.d0 - lhc_z1)) &
                                                     /  (lhc_e - 2.d0 * lhc_sign * self%specified_slope_parameter_)
                reconstructed_primitive_variables(7) = 0.5d0 * (1.d0 + tanh(self%specified_slope_parameter_ * (lhc_sign + lhc_interface_location)))
            end if
        end associate
    end function reconstruct_lhc

    pure function reconstruct_rhc( &
        self                     , &
        primitive_variables_set  , &
        face_to_cell_index       , &
        cell_centor_positions    , &
        face_centor_positions    , &
        face_index               , &
        num_local_cells          , &
        num_primitive_variables      ) result(reconstructed_primitive_variables)

        class  (rho_thinc    ), intent(in) :: self
        integer(int_kind     ), intent(in) :: face_index
        real   (real_kind    ), intent(in) :: primitive_variables_set          (:, :)
        integer(int_kind     ), intent(in) :: face_to_cell_index               (:, :)
        real   (real_kind    ), intent(in) :: cell_centor_positions            (:, :)
        real   (real_kind    ), intent(in) :: face_centor_positions            (:, :)
        integer(int_kind     ), intent(in) :: num_local_cells
        integer(int_kind     ), intent(in) :: num_primitive_variables
        real   (real_kind    )             :: reconstructed_primitive_variables(num_primitive_variables)

        integer(int_kind )              :: lhc_index, rhc_index, llhc_index, rrhc_index
        real   (real_kind)              :: rhc_interface_location, rhc_sign, rhc_a, rhc_b, rhc_d, rhc_e
        logical                         :: rhc_is_monotonicity

        llhc_index = face_to_cell_index(num_local_cells-1, face_index)
        lhc_index  = face_to_cell_index(num_local_cells+0, face_index)
        rhc_index  = face_to_cell_index(num_local_cells+1, face_index)
        rrhc_index = face_to_cell_index(num_local_cells+2, face_index)

        reconstructed_primitive_variables(:) = self%primary_reconstructor_%reconstruct_rhc( &
            primitive_variables_set  , &
            face_to_cell_index       , &
            cell_centor_positions    , &
            face_centor_positions    , &
            face_index               , &
            num_local_cells          , &
            num_primitive_variables    &
        )

        rhc_sign  = sign(1.d0, primitive_variables_set(7, rrhc_index) - primitive_variables_set(7, lhc_index))
        rhc_is_monotonicity = (primitive_variables_set(7, rrhc_index) - primitive_variables_set(7, rhc_index)) * (primitive_variables_set(7, rhc_index) - primitive_variables_set(7, lhc_index)) > 0.d0
        associate(                                          &
            rhc_rho1 => primitive_variables_set(1, rhc_index), &
            rhc_rho2 => primitive_variables_set(2, rhc_index), &
            rhc_z1   => primitive_variables_set(7, rhc_index)  &
        )
            if (rhc_z1 < self%epsilon_) then
                reconstructed_primitive_variables(7) = 0.d0
            else if (rhc_z1 > 1.d0 - self%epsilon_) then
                reconstructed_primitive_variables(7) = 1.d0
            else if (rhc_is_monotonicity) then
                rhc_a    = exp (2.d0 * rhc_sign * self%specified_slope_parameter_)
                rhc_b    = exp (2.d0 * rhc_sign * self%specified_slope_parameter_ * rhc_z1)
                rhc_interface_location = 1.d0 / (2.d0 * self%specified_slope_parameter_) * log((rhc_b - 1.d0) / (rhc_a - rhc_b))
                rhc_d = exp(2.d0 * self%specified_slope_parameter_ * rhc_interface_location)
                rhc_e = log((rhc_a * rhc_d + 1.d0)**2.d0 / (rhc_a * (rhc_d + 1.d0)**2.d0))
                reconstructed_primitive_variables(1) =  (4.d0 * rhc_sign * self%specified_slope_parameter_ * rhc_rho1 * rhc_z1) &
                                                     /  (rhc_e + 2.d0 * rhc_sign * self%specified_slope_parameter_)
                reconstructed_primitive_variables(2) = -(4.d0 * rhc_sign * self%specified_slope_parameter_ * rhc_rho2 * (1.d0 - rhc_z1)) &
                                                     /  (rhc_e - 2.d0 * rhc_sign * self%specified_slope_parameter_)
                reconstructed_primitive_variables(7) = 0.5d0 * (1.d0 + tanh(self%specified_slope_parameter_ * rhc_interface_location))
            end if
        end associate
    end function reconstruct_rhc
end module class_rho_thinc