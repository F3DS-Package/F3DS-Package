module abstract_measurement_tool
    implicit none

    private

    type, public, abstract :: measurement_tool_generator
        contains
        generic, public :: generate => generate_from_configuration, generate_from_name
        procedure(generate_from_configuration_interface), pass(self), deferred :: generate_from_configuration
        procedure(generate_from_name_interface         ), pass(self), deferred :: generate_from_name
    end type measurement_tool_generator

    type, public, abstract :: measurement_tool
        contains
        procedure(initialize_interface ), pass(self), deferred :: initialize
        procedure(write_interface      ), pass(self), deferred :: write
        procedure(is_writable_interface), pass(self), deferred :: is_writable
    end type measurement_tool

    ! ### Generator interfaces ###
    abstract interface
        subroutine generate_from_configuration_interface(self, an_measurement_tool, a_config)
            use abstract_configuration
            import measurement_tool_generator
            import measurement_tool
            class(measurement_tool_generator),          intent(inout) :: self
            class(measurement_tool          ), pointer, intent(inout) :: an_measurement_tool
            class(configuration             ),          intent(inout) :: a_config
        end subroutine generate_from_configuration_interface

        subroutine generate_from_name_interface(self, an_measurement_tool, name)
            use abstract_configuration
            import measurement_tool_generator
            import measurement_tool
            class    (measurement_tool_generator),              intent(inout) :: self
            class    (measurement_tool          ), pointer    , intent(inout) :: an_measurement_tool
            character(len=:                     ), allocatable, intent(in   ) :: name
        end subroutine generate_from_name_interface
    end interface

    ! ### Measurement tool interfaces ###
    abstract interface
        subroutine initialize_interface( &
            self,                        &
            a_config,                    &
            cell_positions,              &
            cell_volumes,                &
            is_real_cell,                &
            face_to_cell_index,          &
            face_centor_positions,       &
            face_normal_vectors,         &
            face_areas,                  &
            num_cells,                   &
            num_faces,                   &
            num_local_cells              &
        )
            use abstract_configuration
            use typedef_module
            import measurement_tool
            class  (measurement_tool), intent(inout) :: self
            class  (configuration   ), intent(inout) :: a_config
            real   (real_kind       ), intent(in   ) :: cell_positions       (:, :)
            real   (real_kind       ), intent(in   ) :: cell_volumes         (:)
            logical                  , intent(in   ) :: is_real_cell         (:)
            integer(int_kind        ), intent(in   ) :: face_to_cell_index   (:, :)
            real   (real_kind       ), intent(in   ) :: face_centor_positions(:, :)
            real   (real_kind       ), intent(in   ) :: face_normal_vectors  (:, :)
            real   (real_kind       ), intent(in   ) :: face_areas           (:)
            integer(int_kind        ), intent(in   ) :: num_cells
            integer(int_kind        ), intent(in   ) :: num_faces
            integer(int_kind        ), intent(in   ) :: num_local_cells
        end subroutine initialize_interface

        subroutine write_interface( &
            self,                   &
            time,                   &
            values_set,             &
            cell_positions,         &
            cell_volumes,           &
            face_to_cell_index,     &
            face_areas,             &
            num_local_cells         &
        )
            use typedef_module
            import measurement_tool
            class  (measurement_tool), intent(inout) :: self
            real   (real_kind       ), intent(in   ) :: time
            real   (real_kind       ), intent(in   ) :: values_set        (:,:)
            real   (real_kind       ), intent(in   ) :: cell_positions    (:,:)
            real   (real_kind       ), intent(in   ) :: cell_volumes      (:)
            integer(int_kind        ), intent(in   ) :: face_to_cell_index(:,:)
            real   (real_kind       ), intent(in   ) :: face_areas        (:)
            integer(int_kind        ), intent(in   ) :: num_local_cells
        end subroutine write_interface

        pure function is_writable_interface(self, time) result(juge)
            use typedef_module
            import measurement_tool
            class  (measurement_tool), intent(in) :: self
            real   (real_kind       ), intent(in) :: time
            logical                               :: juge
        end function is_writable_interface
    end interface
end module abstract_measurement_tool