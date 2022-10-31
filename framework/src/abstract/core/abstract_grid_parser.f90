module abstract_grid_parser
    implicit none

    private

    type, public, abstract :: grid_parser
        contains
        procedure(parse_interface), pass(self), deferred :: parse
        procedure(close_interface), pass(self), deferred :: close

        procedure(get_number_of_cells_interface)      , pass(self), deferred :: get_number_of_cells
        procedure(get_number_of_faces_interface)      , pass(self), deferred :: get_number_of_faces
        procedure(get_number_of_ghost_cells_interface), pass(self), deferred :: get_number_of_ghost_cells
        procedure(get_number_of_points_interface)     , pass(self), deferred :: get_number_of_points

        procedure(get_number_of_boundary_faces_interface)  , pass(self), deferred :: get_number_of_boundary_faces

        procedure(get_cells_interface          ), pass(self), deferred :: get_cells
        procedure(get_faces_interface          ), pass(self), deferred :: get_faces
        procedure(get_boundaries_interface     ), pass(self), deferred :: get_boundaries
        procedure(get_cell_geometries_interface), pass(self), deferred :: get_cell_geometries
    end type grid_parser

    abstract interface
        subroutine parse_interface(self, config)
            use abstract_configuration
            import grid_parser
            class(grid_parser  ), intent(inout) :: self
            class(configuration), intent(inout) :: config
        end subroutine parse_interface

        subroutine close_interface(self)
            import grid_parser
            class(grid_parser), intent(inout) :: self
        end subroutine close_interface

        function get_number_of_cells_interface(self) result(n)
            use typedef_module
            import grid_parser
            class  (grid_parser), intent(in) :: self
            integer(int_kind)                :: n
        end function get_number_of_cells_interface

        function get_number_of_faces_interface(self) result(n)
            use typedef_module
            import grid_parser
            class  (grid_parser), intent(in) :: self
            integer(int_kind)                :: n
        end function get_number_of_faces_interface

        function get_number_of_ghost_cells_interface(self) result(n)
            use typedef_module
            import grid_parser
            class  (grid_parser), intent(in) :: self
            integer(int_kind)                :: n
        end function get_number_of_ghost_cells_interface

        function get_number_of_boundary_faces_interface(self, type) result(n)
            use typedef_module
            use face_type_module
            import grid_parser
            class  (grid_parser        ), intent(in) :: self
            integer(boundary_face_type_kind), intent(in) :: type
            integer(int_kind)                :: n
        end function get_number_of_boundary_faces_interface

        function get_number_of_points_interface(self) result(n)
            use typedef_module
            import grid_parser
            class  (grid_parser), intent(in) :: self
            integer(int_kind)                :: n
        end function get_number_of_points_interface

        subroutine get_cells_interface(self, centor_positions, volumes, is_real_cell)
            use typedef_module
            import grid_parser
            class(grid_parser), intent(in   ) :: self
            real (real_kind  ), intent(inout) :: centor_positions(:,:)
            real (real_kind  ), intent(inout) :: volumes         (:)
            logical           , intent(inout) :: is_real_cell  (:)
        end subroutine get_cells_interface

        subroutine get_faces_interface(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas)
            use typedef_module
            import grid_parser
            class  (grid_parser), intent(in   ) :: self
            integer(int_kind   ), intent(inout) :: reference_cell_indexs(:,:)
            real   (real_kind  ), intent(inout) :: normal_vectors       (:,:)
            real   (real_kind  ), intent(inout) :: tangential1_vectors  (:,:)
            real   (real_kind  ), intent(inout) :: tangential2_vectors  (:,:)
            real   (real_kind  ), intent(inout) :: positions            (:,:)
            real   (real_kind  ), intent(inout) :: areas                (:)
        end subroutine get_faces_interface

        subroutine get_boundaries_interface(self, face_types)
            use typedef_module
            use face_type_module
            import grid_parser
            class  (grid_parser        ), intent(in   ) :: self
            integer(boundary_face_type_kind), intent(inout) :: face_types(:)
        end subroutine get_boundaries_interface

        subroutine get_cell_geometries_interface(self, points, cell_geometries, cell_types)
            use typedef_module
            use class_point_id_list
            import grid_parser
            class  (grid_parser  ), intent(in   ) :: self
            real   (real_kind    ), intent(inout) :: points         (:, :)
            class  (point_id_list), intent(inout) :: cell_geometries(:)
            integer(type_kind    ), intent(inout) :: cell_types     (:)
        end subroutine get_cell_geometries_interface
    end interface
end module abstract_grid_parser