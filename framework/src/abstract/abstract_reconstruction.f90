module abstract_reconstruction
    implicit none

    private

    type, public, abstract :: reconstruction
        contains

        procedure(initialize_interface     ), pass(self), deferred :: initialize
        procedure(reconstruct_lhc_interface), pass(self), deferred :: reconstruct_lhc
        procedure(reconstruct_rhc_interface), pass(self), deferred :: reconstruct_rhc
    end type

    abstract interface
        subroutine initialize_interface(self, config)
            use abstract_configuration
            import reconstruction
            class(reconstruction), intent(inout) :: self
            class(configuration ), intent(inout) :: config
        end subroutine initialize_interface

        pure function reconstruct_lhc_interface( &
            self                               , &
            face_index                         , &
            primitive_variables_set            , &
            face_to_cell_index                 , &
            cell_centor_positions              , &
            face_centor_positions                  ) result(reconstructed_primitive_variables)

            use typedef_module
            import reconstruction

            class  (reconstruction), intent(in) :: self
            integer(int_kind      ), intent(in) :: face_index
            real   (real_kind     ), intent(in) :: primitive_variables_set          (:, :)
            integer(int_kind      ), intent(in) :: face_to_cell_index               (:, :)
            real   (real_kind     ), intent(in) :: cell_centor_positions            (:, :)
            real   (real_kind     ), intent(in) :: face_centor_positions            (:, :)
            real   (real_kind     )             :: reconstructed_primitive_variables(size(primitive_variables_set(:,1)))
        end function reconstruct_lhc_interface

        pure function reconstruct_rhc_interface( &
            self                               , &
            face_index                         , &
            primitive_variables_set            , &
            face_to_cell_index                 , &
            cell_centor_positions              , &
            face_centor_positions                  ) result(reconstructed_primitive_variables)

            use typedef_module
            import reconstruction

            class  (reconstruction), intent(in) :: self
            integer(int_kind      ), intent(in) :: face_index
            real   (real_kind     ), intent(in) :: primitive_variables_set          (:, :)
            integer(int_kind      ), intent(in) :: face_to_cell_index               (:, :)
            real   (real_kind     ), intent(in) :: cell_centor_positions            (:, :)
            real   (real_kind     ), intent(in) :: face_centor_positions            (:, :)
            real   (real_kind     )             :: reconstructed_primitive_variables(size(primitive_variables_set(:,1)))
        end function reconstruct_rhc_interface
    end interface
end module abstract_reconstruction