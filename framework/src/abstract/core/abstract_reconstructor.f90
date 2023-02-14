module abstract_reconstructor
    implicit none

    private

    type, public, abstract :: reconstructor_generator
        contains
        generic, public :: generate => generate_from_configuration, generate_from_name
        procedure(generate_from_configuration_interface), pass(self), deferred :: generate_from_configuration
        procedure(generate_from_name_interface         ), pass(self), deferred :: generate_from_name
    end type reconstructor_generator

    type, public, abstract :: reconstructor
        contains

        procedure(initialize_interface     ), pass(self), deferred :: initialize
        procedure(reconstruct_lhc_interface), pass(self), deferred :: reconstruct_lhc
        procedure(reconstruct_rhc_interface), pass(self), deferred :: reconstruct_rhc
    end type

    abstract interface
        subroutine generate_from_configuration_interface(self, a_reconstructor, a_config)
            use abstract_configuration
            import reconstructor_generator
            import reconstructor
            class(reconstructor_generator),          intent(inout) :: self
            class(reconstructor          ), pointer, intent(inout) :: a_reconstructor
            class(configuration          ),          intent(inout) :: a_config
        end subroutine generate_from_configuration_interface

        subroutine generate_from_name_interface(self, a_reconstructor, name)
            use abstract_configuration
            import reconstructor_generator
            import reconstructor
            class    (reconstructor_generator),              intent(inout) :: self
            class    (reconstructor          ), pointer    , intent(inout) :: a_reconstructor
            character(len=:                  ), allocatable, intent(in   ) :: name
        end subroutine generate_from_name_interface
    end interface

    abstract interface
        subroutine initialize_interface(self, config, a_reconstructor_generator)
            use abstract_configuration
            import reconstructor
            import reconstructor_generator
            class(reconstructor          ),           intent(inout) :: self
            class(configuration          ),           intent(inout) :: config
            class(reconstructor_generator), optional, intent(inout) :: a_reconstructor_generator
        end subroutine initialize_interface

        pure function reconstruct_lhc_interface( &
            self                               , &
            primitive_variables_set            , &
            face_to_cell_index                 , &
            cell_centor_positions              , &
            face_centor_positions              , &
            face_index                         , &
            num_local_cells                    , &
            num_primitive_variables                ) result(reconstructed_primitive_variables)

            use typedef_module
            import reconstructor

            class  (reconstructor), intent(in) :: self
            integer(int_kind     ), intent(in) :: face_index
            real   (real_kind    ), intent(in) :: primitive_variables_set          (:, :)
            integer(int_kind     ), intent(in) :: face_to_cell_index               (:, :)
            real   (real_kind    ), intent(in) :: cell_centor_positions            (:, :)
            real   (real_kind    ), intent(in) :: face_centor_positions            (:, :)
            integer(int_kind     ), intent(in) :: num_local_cells
            integer(int_kind     ), intent(in) :: num_primitive_variables
            real   (real_kind    )             :: reconstructed_primitive_variables(num_primitive_variables)
        end function reconstruct_lhc_interface

        pure function reconstruct_rhc_interface( &
            self                               , &
            primitive_variables_set            , &
            face_to_cell_index                 , &
            cell_centor_positions              , &
            face_centor_positions              , &
            face_index                         , &
            num_local_cells                    , &
            num_primitive_variables                ) result(reconstructed_primitive_variables)

            use typedef_module
            import reconstructor

            class  (reconstructor), intent(in) :: self
            integer(int_kind     ), intent(in) :: face_index
            real   (real_kind    ), intent(in) :: primitive_variables_set          (:, :)
            integer(int_kind     ), intent(in) :: face_to_cell_index               (:, :)
            real   (real_kind    ), intent(in) :: cell_centor_positions            (:, :)
            real   (real_kind    ), intent(in) :: face_centor_positions            (:, :)
            integer(int_kind     ), intent(in) :: num_local_cells
            integer(int_kind     ), intent(in) :: num_primitive_variables
            real   (real_kind    )             :: reconstructed_primitive_variables(num_primitive_variables)
        end function reconstruct_rhc_interface
    end interface
end module abstract_reconstructor