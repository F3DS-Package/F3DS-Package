module abstract_grid_parser
    ! """
    ! This module defines an abstract grid parser for parsing grid data.
    ! The `grid_parser` type is an abstract type that contains deferred procedures for parsing and closing the grid data. It also
    ! contains deferred procedures for retrieving various information about the grid, such as the number of cells, faces, ghost cells,
    ! points, and boundary faces.
    ! The module also defines an abstract interface that specifies the signatures of the deferred procedures. These procedures include
    ! `parse_interface` for parsing the grid data, `close_interface` for closing the grid data, `get_number_of_cells_interface` for
    ! retrieving the number of cells, `get_number_of_faces_interface` for retrieving the number of faces,
    ! `get_number_of_ghost_cells_interface` for retrieving the number of ghost cells, `get_number_of_boundary_faces_interface` for
    ! retrieving the number of boundary faces, `get_number_of_points_interface` for retrieving the number of points,
    ! `get_cells_interface` for retrieving cell data, `get_faces_interface` for retrieving face data, `get_boundaries_interface` for
    ! retrieving boundary face data, and `get_cell_geometries_interface` for retrieving cell geometry data.
    ! The module is intended to be used as a base for implementing specific grid parsers for different grid formats.
    ! """
    implicit none

    private

    type, public, abstract :: grid_parser
        ! This type represents a grid parser object.
        ! > Parse the grid file.
        ! > Close the grid file.
        ! > Get the number of cells in the grid.
        ! > Get the number of faces in the grid.
        ! > Get the number of ghost cells in the grid.
        ! > Get the number of points in the grid.
        ! > Get the number of boundary faces in the grid.
        ! > Get the cells of the grid.
        ! > Get the faces of the grid.
        ! > Get the boundaries of the grid.
        ! > Get the geometries of the cells in the grid.
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
            ! -----------------------------------------------------------------------
            ! Subroutine: parse_interface
            ! -----------------------------------------------------------------------
            ! This subroutine is used to parse the interface of a grid parser object
            ! and a configuration object. It takes in a grid parser object and a
            ! configuration object as arguments and updates them accordingly.
            !
            ! Parameters:
            ! - self: A grid parser object. It is passed by reference and modified
            ! within the subroutine.
            ! - config: A configuration object. It is passed by reference and
            ! modified within the subroutine.
            ! -----------------------------------------------------------------------
            use abstract_configuration
            import grid_parser
            class(grid_parser  ), intent(inout) :: self
            class(configuration), intent(inout) :: config
        end subroutine parse_interface

        subroutine close_interface(self)
            ! -----------------------------------------------------------------------
            ! close_interface subroutine
            ! -----------------------------------------------------------------------
            ! Closes the interface of the grid_parser object.
            !
            ! Parameters:
            ! self: grid_parser object
            ! The grid_parser object whose interface is to be closed.
            ! -----------------------------------------------------------------------
            import grid_parser
            class(grid_parser), intent(inout) :: self
        end subroutine close_interface

        function get_number_of_cells_interface(self) result(n)
            ! > Returns the number of cells in the grid.
            ! This function is an interface for the `get_number_of_cells` method of the `grid_parser` class.
            !
            ! @param self: The `grid_parser` object.
            ! @return n: The number of cells in the grid.
            use typedef_module
            import grid_parser
            class  (grid_parser), intent(in) :: self
            integer(int_kind)                :: n
        end function get_number_of_cells_interface

        function get_number_of_faces_interface(self) result(n)
            ! > Returns the number of faces in the grid_parser object.
            ! This function is an interface to the get_number_of_faces_interface function.
            !
            ! @param self: The grid_parser object.
            ! @return n: The number of faces in the grid_parser object.
            use typedef_module
            import grid_parser
            class  (grid_parser), intent(in) :: self
            integer(int_kind)                :: n
        end function get_number_of_faces_interface

        function get_number_of_ghost_cells_interface(self) result(n)
            ! -----------------------------------------------------------------------
            ! Function: get_number_of_ghost_cells_interface
            ! -----------------------------------------------------------------------
            ! Description:
            ! This function returns the number of ghost cells for the given grid_parser object.
            !
            ! Parameters:
            ! self: grid_parser
            ! An object of type grid_parser representing the grid.
            !
            ! Returns:
            ! n: integer(int_kind)
            ! The number of ghost cells for the grid.
            ! -----------------------------------------------------------------------
            use typedef_module
            import grid_parser
            class  (grid_parser), intent(in) :: self
            integer(int_kind)                :: n
        end function get_number_of_ghost_cells_interface

        function get_number_of_boundary_faces_interface(self, type) result(n)
            ! -----------------------------------------------------------------------
            ! Function: get_number_of_boundary_faces_interface
            ! -----------------------------------------------------------------------
            ! Description:
            ! This function returns the number of boundary faces of a specific type
            ! in the given grid.
            !
            ! Parameters:
            ! self: grid_parser
            ! An object of the grid_parser class representing the grid.
            ! type: boundary_face_type_kind
            ! An integer representing the type of boundary face.
            !
            ! Returns:
            ! n: int_kind
            ! An integer representing the number of boundary faces of the
            ! specified type in the grid.
            ! -----------------------------------------------------------------------
            use typedef_module
            use face_type_module
            import grid_parser
            class  (grid_parser        ), intent(in) :: self
            integer(boundary_face_type_kind), intent(in) :: type
            integer(int_kind)                :: n
        end function get_number_of_boundary_faces_interface

        function get_number_of_points_interface(self) result(n)
            ! -----------------------------------------------------------------------
            ! Function: get_number_of_points_interface
            ! -----------------------------------------------------------------------
            ! Description:
            ! This function returns the number of points in the grid_parser object.
            !
            ! Parameters:
            ! self: grid_parser
            ! An object of type grid_parser.
            !
            ! Returns:
            ! n: integer(int_kind)
            ! The number of points in the grid_parser object.
            ! -----------------------------------------------------------------------
            use typedef_module
            import grid_parser
            class  (grid_parser), intent(in) :: self
            integer(int_kind)                :: n
        end function get_number_of_points_interface

        subroutine get_cells_interface(self, centor_positions, volumes, is_real_cell)
            ! ---------------------------------------------------------------------------
            ! This subroutine retrieves the cell interface information from the grid parser object.
            !
            ! Parameters:
            ! self: grid_parser
            ! The grid parser object containing the grid information.
            ! centor_positions: real (real_kind), intent(inout)
            ! The array of cell center positions.
            ! volumes: real (real_kind), intent(inout)
            ! The array of cell volumes.
            ! is_real_cell: logical, intent(inout)
            ! The array indicating whether each cell is real or not.
            !
            ! ---------------------------------------------------------------------------
            use typedef_module
            import grid_parser
            class(grid_parser), intent(in   ) :: self
            real (real_kind  ), intent(inout) :: centor_positions(:,:)
            real (real_kind  ), intent(inout) :: volumes         (:)
            logical           , intent(inout) :: is_real_cell  (:)
        end subroutine get_cells_interface

        subroutine get_faces_interface(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas)
            ! This subroutine calculates the interface properties of the faces in the grid.
            !
            ! Parameters:
            ! self: grid_parser object
            ! The grid_parser object containing the grid information.
            ! reference_cell_indexs: 2D integer array
            ! The array to store the reference cell indices of the faces.
            ! normal_vectors: 2D real array
            ! The array to store the normal vectors of the faces.
            ! tangential1_vectors: 2D real array
            ! The array to store the first tangential vectors of the faces.
            ! tangential2_vectors: 2D real array
            ! The array to store the second tangential vectors of the faces.
            ! positions: 2D real array
            ! The array to store the positions of the faces.
            ! areas: 1D real array
            ! The array to store the areas of the faces.
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
            ! -----------------------------------------------------------------------
            ! Subroutine: get_boundaries_interface
            ! -----------------------------------------------------------------------
            ! Description:
            ! This subroutine is used to get the boundary face types for a given grid.
            !
            ! Parameters:
            ! self: grid_parser (input)
            ! - An object of type grid_parser representing the grid.
            ! face_types: integer(boundary_face_type_kind) (input/output)
            ! - An array of boundary face types.
            !
            ! -----------------------------------------------------------------------
            use typedef_module
            use face_type_module
            import grid_parser
            class  (grid_parser        ), intent(in   ) :: self
            integer(boundary_face_type_kind), intent(inout) :: face_types(:)
        end subroutine get_boundaries_interface

        subroutine get_cell_geometries_interface(self, points, cell_geometries, cell_types)
            ! ---------------------------------------------------------------------------------------------
            ! This subroutine retrieves the cell geometries and types from the grid_parser object.
            !
            ! Parameters:
            ! self            - The grid_parser object containing the grid information.
            ! points          - A 2D array of real numbers representing the points.
            ! cell_geometries - An array of point_id_list objects representing the cell geometries.
            ! cell_types      - An array of integers representing the cell types.
            ! ---------------------------------------------------------------------------------------------
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