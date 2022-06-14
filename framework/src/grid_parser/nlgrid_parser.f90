module class_nlgrid_parser
    use, intrinsic :: iso_fortran_env
    use json_module
    use typedef_module
    use stdio_module
    use boundary_type_module
    use grid_utils_module
    use abstract_grid_parser
    use class_point_id_list

    implicit none

    private

    type, public, extends(grid_parser) :: nlgrid_parser
        integer(int_kind), private :: num_ghost_cells_ = 3

        integer(int_kind), private :: imax
        integer(int_kind), private :: jmax
        integer(int_kind), private :: kmax
        integer(int_kind), private :: imin
        integer(int_kind), private :: jmin
        integer(int_kind), private :: kmin

        real(real_kind), private, allocatable :: x_poss_cell_cent(:,:,:)
        real(real_kind), private, allocatable :: y_poss_cell_cent(:,:,:)
        real(real_kind), private, allocatable :: z_poss_cell_cent(:,:,:)
        real(real_kind), private, allocatable :: cell_vols(:,:,:)
        real(real_kind), private, allocatable :: areas_i_drct_face(:,:,:)
        real(real_kind), private, allocatable :: areas_j_drct_face(:,:,:)
        real(real_kind), private, allocatable :: areas_k_drct_face(:,:,:)
        real(real_kind), private, allocatable :: normal_vecs_i_drct_face(:,:,:,:)
        real(real_kind), private, allocatable :: normal_vecs_j_drct_face(:,:,:,:)
        real(real_kind), private, allocatable :: normal_vecs_k_drct_face(:,:,:,:)
        real(real_kind), private, allocatable :: tan1_vecs_i_drct_face(:,:,:,:)
        real(real_kind), private, allocatable :: tan1_vecs_j_drct_face(:,:,:,:)
        real(real_kind), private, allocatable :: tan1_vecs_k_drct_face(:,:,:,:)
        real(real_kind), private, allocatable :: tan2_vecs_i_drct_face(:,:,:,:)
        real(real_kind), private, allocatable :: tan2_vecs_j_drct_face(:,:,:,:)
        real(real_kind), private, allocatable :: tan2_vecs_k_drct_face(:,:,:,:)
        real(real_kind), private, allocatable :: drvts_i_drct(:,:,:,:)
        real(real_kind), private, allocatable :: drvts_j_drct(:,:,:,:)
        real(real_kind), private, allocatable :: drvts_k_drct(:,:,:,:)
        real(real_kind), private, allocatable :: x_poss_cell_edge(:,:,:)
        real(real_kind), private, allocatable :: y_poss_cell_edge(:,:,:)
        real(real_kind), private, allocatable :: z_poss_cell_edge(:,:,:)

        integer(kind(boundary_type)), private :: x_plus_direction, x_minus_direction
        integer(kind(boundary_type)), private :: y_plus_direction, y_minus_direction
        integer(kind(boundary_type)), private :: z_plus_direction, z_minus_direction

        logical, private :: parsed = .false.

        contains

        procedure, pass(self) :: parse
        procedure, pass(self) :: close

        procedure, pass(self) :: get_number_of_cells
        procedure, pass(self) :: get_number_of_faces
        procedure, pass(self) :: get_number_of_ghost_cells

        procedure, pass(self) :: get_number_of_outflow_faces
        procedure, pass(self) :: get_number_of_slipwall_faces
        procedure, pass(self) :: get_number_of_symmetric_faces
        procedure, pass(self) :: get_number_of_points

        procedure, pass(self) :: get_cells
        procedure, pass(self) :: get_faces
        procedure, pass(self) :: get_boundaries
        procedure, pass(self) :: get_cell_geometries

        procedure, private, pass(self) :: assign_upper_x_face
        procedure, private, pass(self) :: assign_upper_y_face
        procedure, private, pass(self) :: assign_upper_z_face
        procedure, private, pass(self) :: assign_lower_x_face
        procedure, private, pass(self) :: assign_lower_y_face
        procedure, private, pass(self) :: assign_lower_z_face
        procedure, private, pass(self) :: assign_boundary
    end type nlgrid_parser

    contains

    subroutine parse(self, filepath)
        class    (nlgrid_parser), intent(inout) :: self
        character(len=*)        , intent(in   ) :: filepath

        integer(int_kind)  :: i,j,k,n
        integer(int_kind)  :: unit_number
        real   (real_kind) :: dtheta

        type(json_file) :: json
        character(len=:), allocatable :: error_msg
        logical :: status_ok, found
        character(kind=json_CK,len=:), allocatable :: string_value

        if(self%parsed)then
            call call_error("'parse' method of nlgrid_parser is already called. But you call 'parse' method.")
        end if

        ! Real a grid file
        open(newunit=unit_number, file=filepath, access = 'stream', form = 'unformatted', status = 'old')
        read(unit_number) self%imin, self%jmin, self%kmin
        read(unit_number) self%imax, self%jmax, self%kmax
        read(unit_number) dtheta

#ifdef _DEBUG
        print *, "DEBUG: nlgrid parser: grid size:"
        print *, "i = [", self%imin, ", ", self%imax, "]"
        print *, "j = [", self%jmin, ", ", self%jmax, "]"
        print *, "k = [", self%kmin, ", ", self%kmax, "]"
        print *, "---------------------------------"
#endif

        allocate(self%x_poss_cell_cent       (self%imin-2:self%imax+2, self%jmin-2:self%jmax+2, self%kmin-2:self%kmax+2     ))
        allocate(self%y_poss_cell_cent       (self%imin-2:self%imax+2, self%jmin-2:self%jmax+2, self%kmin-2:self%kmax+2     ))
        allocate(self%z_poss_cell_cent       (self%imin-2:self%imax+2, self%jmin-2:self%jmax+2, self%kmin-2:self%kmax+2     ))
        allocate(self%cell_vols              (self%imin-1:self%imax+1, self%jmin-1:self%jmax+1, self%kmin-1:self%kmax+1     ))
        allocate(self%areas_i_drct_face      (self%imin-2:self%imax+1, self%jmin-1:self%jmax+1, self%kmin-1:self%kmax+1     ))
        allocate(self%areas_j_drct_face      (self%imin-1:self%imax+1, self%jmin-2:self%jmax+1, self%kmin-1:self%kmax+1     ))
        allocate(self%areas_k_drct_face      (self%imin-1:self%imax+1, self%jmin-1:self%jmax+1, self%kmin-2:self%kmax+1     ))
        allocate(self%normal_vecs_i_drct_face(self%imin-2:self%imax+1, self%jmin-1:self%jmax+1, self%kmin-1:self%kmax+1, 1:3))
        allocate(self%normal_vecs_j_drct_face(self%imin-1:self%imax+1, self%jmin-2:self%jmax+1, self%kmin-1:self%kmax+1, 1:3))
        allocate(self%normal_vecs_k_drct_face(self%imin-1:self%imax+1, self%jmin-1:self%jmax+1, self%kmin-2:self%kmax+1, 1:3))
        allocate(self%tan1_vecs_i_drct_face  (self%imin-1:self%imax,   self%jmin  :self%jmax,   self%kmin  :self%kmax,   1:3))
        allocate(self%tan1_vecs_j_drct_face  (self%imin  :self%imax,   self%jmin-1:self%jmax,   self%kmin  :self%kmax,   1:3))
        allocate(self%tan1_vecs_k_drct_face  (self%imin  :self%imax,   self%jmin  :self%jmax,   self%kmin-1:self%kmax,   1:3))
        allocate(self%tan2_vecs_i_drct_face  (self%imin-1:self%imax,   self%jmin  :self%jmax,   self%kmin  :self%kmax,   1:3))
        allocate(self%tan2_vecs_j_drct_face  (self%imin  :self%imax,   self%jmin-1:self%jmax,   self%kmin  :self%kmax,   1:3))
        allocate(self%tan2_vecs_k_drct_face  (self%imin  :self%imax,   self%jmin  :self%jmax,   self%kmin-1:self%kmax,   1:3))
        allocate(self%drvts_i_drct           (self%imin-1:self%imax+1, self%jmin-1:self%jmax+1, self%kmin-1:self%kmax+1, 1:3))
        allocate(self%drvts_j_drct           (self%imin-1:self%imax+1, self%jmin-1:self%jmax+1, self%kmin-1:self%kmax+1, 1:3))
        allocate(self%drvts_k_drct           (self%imin-1:self%imax+1, self%jmin-1:self%jmax+1, self%kmin-1:self%kmax+1, 1:3))
        allocate(self%x_poss_cell_edge       (self%imin-3:self%imax+2, self%jmin-3:self%jmax+2, self%kmin-3:self%kmax+2     ))
        allocate(self%y_poss_cell_edge       (self%imin-3:self%imax+2, self%jmin-3:self%jmax+2, self%kmin-3:self%kmax+2     ))
        allocate(self%z_poss_cell_edge       (self%imin-3:self%imax+2, self%jmin-3:self%jmax+2, self%kmin-3:self%kmax+2     ))

        read(unit_number) ((( self%x_poss_cell_cent       (i,j,k)  , i = self%imin-2, self%imax+2), j = self%jmin-2, self%jmax+2), k = self%kmin-2, self%kmax+2)
        read(unit_number) ((( self%y_poss_cell_cent       (i,j,k)  , i = self%imin-2, self%imax+2), j = self%jmin-2, self%jmax+2), k = self%kmin-2, self%kmax+2)
        read(unit_number) ((( self%z_poss_cell_cent       (i,j,k)  , i = self%imin-2, self%imax+2), j = self%jmin-2, self%jmax+2), k = self%kmin-2, self%kmax+2)
        read(unit_number) ((( self%cell_vols              (i,j,k)  , i = self%imin-1, self%imax+1), j = self%jmin-1, self%jmax+1), k = self%kmin-1, self%kmax+1)
        read(unit_number) ((( self%areas_i_drct_face      (i,j,k)  , i = self%imin-2, self%imax+1), j = self%jmin-1, self%jmax+1), k = self%kmin-1, self%kmax+1)
        read(unit_number) ((( self%areas_j_drct_face      (i,j,k)  , i = self%imin-1, self%imax+1), j = self%jmin-2, self%jmax+1), k = self%kmin-1, self%kmax+1)
        read(unit_number) ((( self%areas_k_drct_face      (i,j,k)  , i = self%imin-1, self%imax+1), j = self%jmin-1, self%jmax+1), k = self%kmin-2, self%kmax+1)
        read(unit_number) ((((self%normal_vecs_i_drct_face(i,j,k,n), i = self%imin-2, self%imax+1), j = self%jmin-1, self%jmax+1), k = self%kmin-1, self%kmax+1), n = 1, 3)
        read(unit_number) ((((self%normal_vecs_j_drct_face(i,j,k,n), i = self%imin-1, self%imax+1), j = self%jmin-2, self%jmax+1), k = self%kmin-1, self%kmax+1), n = 1, 3)
        read(unit_number) ((((self%normal_vecs_k_drct_face(i,j,k,n), i = self%imin-1, self%imax+1), j = self%jmin-1, self%jmax+1), k = self%kmin-2, self%kmax+1), n = 1, 3)
        read(unit_number) ((((self%tan1_vecs_i_drct_face  (i,j,k,n), i = self%imin-1, self%imax  ), j = self%jmin  , self%jmax  ), k = self%kmin  , self%kmax  ), n = 1, 3)
        read(unit_number) ((((self%tan1_vecs_j_drct_face  (i,j,k,n), i = self%imin  , self%imax  ), j = self%jmin-1, self%jmax  ), k = self%kmin  , self%kmax  ), n = 1, 3)
        read(unit_number) ((((self%tan1_vecs_k_drct_face  (i,j,k,n), i = self%imin  , self%imax  ), j = self%jmin  , self%jmax  ), k = self%kmin-1, self%kmax  ), n = 1, 3)
        read(unit_number) ((((self%tan2_vecs_i_drct_face  (i,j,k,n), i = self%imin-1, self%imax  ), j = self%jmin  , self%jmax  ), k = self%kmin  , self%kmax  ), n = 1, 3)
        read(unit_number) ((((self%tan2_vecs_j_drct_face  (i,j,k,n), i = self%imin  , self%imax  ), j = self%jmin-1, self%jmax  ), k = self%kmin  , self%kmax  ), n = 1, 3)
        read(unit_number) ((((self%tan2_vecs_k_drct_face  (i,j,k,n), i = self%imin  , self%imax  ), j = self%jmin  , self%jmax  ), k = self%kmin-1, self%kmax  ), n = 1, 3)
        read(unit_number) ((((self%drvts_i_drct           (i,j,k,n), i = self%imin-1, self%imax+1), j = self%jmin-1, self%jmax+1), k = self%kmin-1, self%kmax+1), n = 1, 3)
        read(unit_number) ((((self%drvts_j_drct           (i,j,k,n), i = self%imin-1, self%imax+1), j = self%jmin-1, self%jmax+1), k = self%kmin-1, self%kmax+1), n = 1, 3)
        read(unit_number) ((((self%drvts_k_drct           (i,j,k,n), i = self%imin-1, self%imax+1), j = self%jmin-1, self%jmax+1), k = self%kmin-1, self%kmax+1), n = 1, 3)
        read(unit_number) ((( self%x_poss_cell_edge       (i,j,k)  , i = self%imin-3, self%imax+2), j = self%jmin-3, self%jmax+2), k = self%kmin-3, self%kmax+2)
        read(unit_number) ((( self%y_poss_cell_edge       (i,j,k)  , i = self%imin-3, self%imax+2), j = self%jmin-3, self%jmax+2), k = self%kmin-3, self%kmax+2)
        read(unit_number) ((( self%z_poss_cell_edge       (i,j,k)  , i = self%imin-3, self%imax+2), j = self%jmin-3, self%jmax+2), k = self%kmin-3, self%kmax+2)

        close(unit_number)

        ! Read a boundary condition infomation
        call json%initialize()

        call json%load(filename="boundary_condition.json")
        if (json%failed()) then
            call json%check_for_errors(status_ok, error_msg)
            call json%clear_exceptions()
            call json%destroy()
#ifdef _DEBUG
            print *, "DEBUG: JSONFortran error message:"
            print *, error_msg
            print *, "---------------------------------"
#endif
            call call_error("Can not read boundary_condition.json that is writen boundary condition infomation.")
        end if

        call json%get("Xi plus direction", string_value, found)
        if(.not. found) call call_error("Keyword 'Xi plus direction' is not found.")
        self%x_plus_direction = string_to_boundary_type(string_value)
        if(self%x_plus_direction == unknown_boundary_type) call call_error("Unknown boundary type '"//string_value//"' is found.")

        call json%get("Xi minus direction", string_value, found)
        if(.not. found) call call_error("Keyword 'Xi minus direction' is not found.")
        self%x_minus_direction = string_to_boundary_type(string_value)
        if(self%x_minus_direction == unknown_boundary_type) call call_error("Unknown boundary type '"//string_value//"' is found.")

        call json%get("Eta plus direction", string_value, found)
        if(.not. found) call call_error("Keyword 'Eta plus direction' is not found.")
        self%y_plus_direction = string_to_boundary_type(string_value)
        if(self%y_plus_direction == unknown_boundary_type) call call_error("Unknown boundary type '"//string_value//"' is found.")

        call json%get("Eta minus direction", string_value, found)
        if(.not. found) call call_error("Keyword 'Eta minus direction' is not found.")
        self%y_minus_direction = string_to_boundary_type(string_value)
        if(self%y_minus_direction == unknown_boundary_type) call call_error("Unknown boundary type '"//string_value//"' is found.")

        call json%get("Zeta plus direction", string_value, found)
        if(.not. found) call call_error("Keyword 'Zeta plus direction' is not found.")
        self%z_plus_direction = string_to_boundary_type(string_value)
        if(self%z_plus_direction == unknown_boundary_type) call call_error("Unknown boundary type '"//string_value//"' is found.")

        call json%get("Zeta minus direction", string_value, found)
        if(.not. found) call call_error("Keyword 'Zeta minus direction' is not found.")
        self%z_minus_direction = string_to_boundary_type(string_value)
        if(self%z_minus_direction == unknown_boundary_type) call call_error("Unknown boundary type '"//string_value//"' is found.")

        call json%destroy()

        self%parsed = .true.
    end subroutine parse

    subroutine close(self)
        class(nlgrid_parser), intent(inout) :: self

        if(.not. self%parsed)then
            call call_error("'parse' method of nlgrid_parser is not called. But you call 'closs' method.")
        end if

        deallocate(self%x_poss_cell_cent       )
        deallocate(self%y_poss_cell_cent       )
        deallocate(self%z_poss_cell_cent       )
        deallocate(self%cell_vols              )
        deallocate(self%areas_i_drct_face      )
        deallocate(self%areas_j_drct_face      )
        deallocate(self%areas_k_drct_face      )
        deallocate(self%normal_vecs_i_drct_face)
        deallocate(self%normal_vecs_j_drct_face)
        deallocate(self%normal_vecs_k_drct_face)
        deallocate(self%tan1_vecs_i_drct_face  )
        deallocate(self%tan1_vecs_j_drct_face  )
        deallocate(self%tan1_vecs_k_drct_face  )
        deallocate(self%tan2_vecs_i_drct_face  )
        deallocate(self%tan2_vecs_j_drct_face  )
        deallocate(self%tan2_vecs_k_drct_face  )
        deallocate(self%drvts_i_drct           )
        deallocate(self%drvts_j_drct           )
        deallocate(self%drvts_k_drct           )
        deallocate(self%x_poss_cell_edge       )
        deallocate(self%y_poss_cell_edge       )
        deallocate(self%z_poss_cell_edge       )

        self%parsed = .false.
    end subroutine close

    function get_number_of_cells(self) result(n)
        class  (nlgrid_parser), intent(in) :: self
        integer(int_kind     )             :: n

        if(.not. self%parsed)then
            call call_error("'parse' method of nlgrid_parser is not called yet. But you call 'get_number_of_cells' method.")
        end if

        n = (self%imax - self%imin + 1 + 2 * self%num_ghost_cells_) * (self%kmax - self%kmin + 1 + 2 * self%num_ghost_cells_) * (self%jmax - self%jmin + 1 + 2 * self%num_ghost_cells_)
    end function get_number_of_cells

    function get_number_of_faces(self) result(n)
        class  (nlgrid_parser), intent(in) :: self
        integer(int_kind     )             :: n

        if(.not. self%parsed)then
            call call_error("'parse' method of nlgrid_parser is not called yet. But you call 'get_number_of_faces' method.")
        end if

        n = (self%imax - self%imin + 2) * (self%jmax - self%jmin + 1) * (self%kmax - self%kmin + 1) &
          + (self%imax - self%imin + 1) * (self%jmax - self%jmin + 2) * (self%kmax - self%kmin + 1) &
          + (self%imax - self%imin + 1) * (self%jmax - self%jmin + 1) * (self%kmax - self%kmin + 2)
    end function get_number_of_faces

    function get_number_of_ghost_cells(self) result(n)
        class  (nlgrid_parser), intent(in) :: self
        integer(int_kind     )             :: n

        if(.not. self%parsed)then
            call call_error("'parse' method of nlgrid_parser is not called yet. But you call 'get_number_of_ghost_cells' method.")
        end if

        n = self%num_ghost_cells_ ! nlgrid format is followed only 2 ghost cell. But we extend to 3 ghost cell system.
    end function get_number_of_ghost_cells

    function get_number_of_outflow_faces(self) result(n)
        class  (nlgrid_parser), intent(in) :: self
        integer(int_kind     )             :: n

        if(.not. self%parsed)then
            call call_error("'parse' method of nlgrid_parser is not called yet. But you call 'get_number_of_outflow_faces' method.")
        end if

        n = 0

        if(self%x_plus_direction   == outflow_boundary_type) n = n + (self%jmax - self%jmin + 1) * (self%kmax - self%kmin + 1)
        if(self%y_plus_direction   == outflow_boundary_type) n = n + (self%imax - self%imin + 1) * (self%kmax - self%kmin + 1)
        if(self%z_plus_direction   == outflow_boundary_type) n = n + (self%imax - self%imin + 1) * (self%jmax - self%jmin + 1)
        if(self%x_minus_direction == outflow_boundary_type) n = n + (self%jmax - self%jmin + 1) * (self%kmax - self%kmin + 1)
        if(self%y_minus_direction == outflow_boundary_type) n = n + (self%imax - self%imin + 1) * (self%kmax - self%kmin + 1)
        if(self%z_minus_direction == outflow_boundary_type) n = n + (self%imax - self%imin + 1) * (self%jmax - self%jmin + 1)
    end function get_number_of_outflow_faces

    function get_number_of_slipwall_faces(self) result(n)
        class  (nlgrid_parser), intent(in) :: self
        integer(int_kind     )             :: n

        if(.not. self%parsed)then
            call call_error("'parse' method of nlgrid_parser is not called yet. But you call 'get_number_of_slipwall_faces' method.")
        end if

        n = 0

        if(self%x_plus_direction   == slipwall_boundary_type) n = n + (self%jmax - self%jmin + 1) * (self%kmax - self%kmin + 1)
        if(self%y_plus_direction   == slipwall_boundary_type) n = n + (self%imax - self%imin + 1) * (self%kmax - self%kmin + 1)
        if(self%z_plus_direction   == slipwall_boundary_type) n = n + (self%imax - self%imin + 1) * (self%jmax - self%jmin + 1)
        if(self%x_minus_direction == slipwall_boundary_type) n = n + (self%jmax - self%jmin + 1) * (self%kmax - self%kmin + 1)
        if(self%y_minus_direction == slipwall_boundary_type) n = n + (self%imax - self%imin + 1) * (self%kmax - self%kmin + 1)
        if(self%z_minus_direction == slipwall_boundary_type) n = n + (self%imax - self%imin + 1) * (self%jmax - self%jmin + 1)
    end function get_number_of_slipwall_faces

    function get_number_of_symmetric_faces(self) result(n)
        class  (nlgrid_parser), intent(in) :: self
        integer(int_kind     )             :: n

        if(.not. self%parsed)then
            call call_error("'parse' method of nlgrid_parser is not called yet. But you call 'get_number_of_symmetric_faces' method.")
        end if

        n = 0

        if(self%x_plus_direction   == symmetric_boundary_type) n = n + (self%jmax - self%jmin + 1) * (self%kmax - self%kmin + 1)
        if(self%y_plus_direction   == symmetric_boundary_type) n = n + (self%imax - self%imin + 1) * (self%kmax - self%kmin + 1)
        if(self%z_plus_direction   == symmetric_boundary_type) n = n + (self%imax - self%imin + 1) * (self%jmax - self%jmin + 1)
        if(self%x_minus_direction == symmetric_boundary_type) n = n + (self%jmax - self%jmin + 1) * (self%kmax - self%kmin + 1)
        if(self%y_minus_direction == symmetric_boundary_type) n = n + (self%imax - self%imin + 1) * (self%kmax - self%kmin + 1)
        if(self%z_minus_direction == symmetric_boundary_type) n = n + (self%imax - self%imin + 1) * (self%jmax - self%jmin + 1)
    end function get_number_of_symmetric_faces

    function get_number_of_points(self) result(n)
        class  (nlgrid_parser), intent(in) :: self
        integer(int_kind     )             :: n

        if(.not. self%parsed)then
            call call_error("'parse' method of nlgrid_parser is not called yet. But you call 'get_number_of_points' method.")
        end if

        n = (self%imax - self%imin + 2) * (self%jmax - self%jmin + 2) * (self%kmax - self%kmin + 2)
    end function get_number_of_points

    subroutine get_cells(self, centor_positions, volumes, is_real_cell)
        class(nlgrid_parser), intent(in   ) :: self
        real (real_kind    ), intent(inout) :: centor_positions(:,:)
        real (real_kind    ), intent(inout) :: volumes         (:)
        logical             , intent(inout) :: is_real_cell  (:)

        integer(int_kind) :: i, j, k, n
        real (real_kind) :: cell_to_face_vec(3)

        if(.not. self%parsed)then
            call call_error("'parse' method of nlgrid_parser is not called yet. But you call 'get_cells' method.")
        end if

        is_real_cell(:) = .false.

        do k = self%kmin, self%kmax, 1
            do j = self%jmin, self%jmax, 1
                do i = self%imin, self%imax, 1
                    n = convert_structure_index_to_unstructure_index(i, j, k,                            &
                        self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
                        self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
                    centor_positions(n, 1) = self%x_poss_cell_cent(i, j, k)
                    centor_positions(n, 2) = self%y_poss_cell_cent(i, j, k)
                    centor_positions(n, 3) = self%z_poss_cell_cent(i, j, k)
                    volumes         (n)    = self%cell_vols       (i, j, k)
                    is_real_cell    (n)    = .true.
                end do
            end do
        end do

        ! ghosts (x +/- direction)
        do k = self%kmin, self%kmax, 1
            do j = self%jmin, self%jmax, 1
                do i = 1, self%num_ghost_cells_, 1
                    ! x-
                    n = convert_structure_index_to_unstructure_index(self%imin - i, j, k,                &
                        self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
                        self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
                    cell_to_face_vec(1) = 0.5d0 * (self%x_poss_cell_cent(self%imin-1, j  , k  ) &
                                                 + self%x_poss_cell_cent(self%imin-1, j-1, k  ) &
                                                 + self%x_poss_cell_cent(self%imin-1, j  , k-1) &
                                                 + self%x_poss_cell_cent(self%imin-1, j-1, k-1))&
                                        - self%x_poss_cell_cent(self%imin, j, k)
                    cell_to_face_vec(2) = 0.5d0 * (self%y_poss_cell_cent(self%imin-1, j  , k  ) &
                                                 + self%y_poss_cell_cent(self%imin-1, j-1, k  ) &
                                                 + self%y_poss_cell_cent(self%imin-1, j  , k-1) &
                                                 + self%y_poss_cell_cent(self%imin-1, j-1, k-1))&
                                        - self%y_poss_cell_cent(self%imin, j, k)
                    cell_to_face_vec(3) = 0.5d0 * (self%z_poss_cell_cent(self%imin-1, j  , k  ) &
                                                 + self%z_poss_cell_cent(self%imin-1, j-1, k  ) &
                                                 + self%z_poss_cell_cent(self%imin-1, j  , k-1) &
                                                 + self%z_poss_cell_cent(self%imin-1, j-1, k-1))&
                                        - self%z_poss_cell_cent(self%imin, j, k)
                    centor_positions(n, 1) = self%x_poss_cell_cent(self%imin, j, k) + 2.d0 * dble(i) * cell_to_face_vec(1)
                    centor_positions(n, 2) = self%y_poss_cell_cent(self%imin, j, k) + 2.d0 * dble(i) * cell_to_face_vec(2)
                    centor_positions(n, 3) = self%z_poss_cell_cent(self%imin, j, k) + 2.d0 * dble(i) * cell_to_face_vec(3)
                    volumes         (n) = self%cell_vols(self%imin, j, k)
                    ! x+
                    n = convert_structure_index_to_unstructure_index(self%imax + i, j, k,                &
                        self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
                        self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
                    cell_to_face_vec(1) = 0.5d0 * (self%x_poss_cell_cent(self%imax, j  , k  ) &
                                                 + self%x_poss_cell_cent(self%imax, j-1, k  ) &
                                                 + self%x_poss_cell_cent(self%imax, j  , k-1) &
                                                 + self%x_poss_cell_cent(self%imax, j-1, k-1))&
                                        - self%x_poss_cell_cent(self%imax, j, k)
                    cell_to_face_vec(2) = 0.5d0 * (self%y_poss_cell_cent(self%imax, j  , k  ) &
                                                 + self%y_poss_cell_cent(self%imax, j-1, k  ) &
                                                 + self%y_poss_cell_cent(self%imax, j  , k-1) &
                                                 + self%y_poss_cell_cent(self%imax, j-1, k-1))&
                                        - self%y_poss_cell_cent(self%imax, j, k)
                    cell_to_face_vec(3) = 0.5d0 * (self%z_poss_cell_cent(self%imax, j  , k  ) &
                                                 + self%z_poss_cell_cent(self%imax, j-1, k  ) &
                                                 + self%z_poss_cell_cent(self%imax, j  , k-1) &
                                                 + self%z_poss_cell_cent(self%imax, j-1, k-1))&
                                        - self%z_poss_cell_cent(self%imax, j, k)
                    centor_positions(n, 1) = self%x_poss_cell_cent(self%imax, j, k) + 2.d0 * dble(i) * cell_to_face_vec(1)
                    centor_positions(n, 2) = self%y_poss_cell_cent(self%imax, j, k) + 2.d0 * dble(i) * cell_to_face_vec(2)
                    centor_positions(n, 3) = self%z_poss_cell_cent(self%imax, j, k) + 2.d0 * dble(i) * cell_to_face_vec(3)
                    volumes         (n) = self%cell_vols(self%imax, j, k)
                end do
            end do
        end do

        ! ghosts (y +/- direction)
        do k = self%kmin, self%kmax, 1
            do j = 1, self%num_ghost_cells_, 1
                do i = self%imin, self%imax, 1
                    ! y-
                    n = convert_structure_index_to_unstructure_index(i, self%jmin - j, k,                &
                        self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
                        self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
                    cell_to_face_vec(1) = 0.5d0 * (self%x_poss_cell_cent(i  , self%jmin-1, k  ) &
                                                 + self%x_poss_cell_cent(i-1, self%jmin-1, k  ) &
                                                 + self%x_poss_cell_cent(i  , self%jmin-1, k-1) &
                                                 + self%x_poss_cell_cent(i-1, self%jmin-1, k-1))&
                                        - self%x_poss_cell_cent(self%jmin, j, k)
                    cell_to_face_vec(2) = 0.5d0 * (self%y_poss_cell_cent(i  , self%jmin-1, k  ) &
                                                 + self%y_poss_cell_cent(i-1, self%jmin-1, k  ) &
                                                 + self%y_poss_cell_cent(i  , self%jmin-1, k-1) &
                                                 + self%y_poss_cell_cent(i-1, self%jmin-1, k-1))&
                                        - self%y_poss_cell_cent(self%jmin, j, k)
                    cell_to_face_vec(3) = 0.5d0 * (self%z_poss_cell_cent(i  , self%jmin-1, k  ) &
                                                 + self%z_poss_cell_cent(i-1, self%jmin-1, k  ) &
                                                 + self%z_poss_cell_cent(i  , self%jmin-1, k-1) &
                                                 + self%z_poss_cell_cent(i-1, self%jmin-1, k-1))&
                                        - self%z_poss_cell_cent(self%jmin, j, k)
                    centor_positions(n, 1) = self%x_poss_cell_cent(i, self%jmin, k) + 2.d0 * dble(j) * cell_to_face_vec(1)
                    centor_positions(n, 2) = self%y_poss_cell_cent(i, self%jmin, k) + 2.d0 * dble(j) * cell_to_face_vec(2)
                    centor_positions(n, 3) = self%z_poss_cell_cent(i, self%jmin, k) + 2.d0 * dble(j) * cell_to_face_vec(3)
                    volumes         (n)    = self%cell_vols(i, self%jmin, k)
                    ! y+
                    n = convert_structure_index_to_unstructure_index(i, self%jmax + j, k,                &
                        self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
                        self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
                    cell_to_face_vec(1) = 0.5d0 * (self%x_poss_cell_cent(i  , self%jmax, k  ) &
                                                 + self%x_poss_cell_cent(i-1, self%jmax, k  ) &
                                                 + self%x_poss_cell_cent(i  , self%jmax, k-1) &
                                                 + self%x_poss_cell_cent(i-1, self%jmax, k-1))&
                                        - self%x_poss_cell_cent(self%jmax, j, k)
                    cell_to_face_vec(2) = 0.5d0 * (self%y_poss_cell_cent(i  , self%jmax, k  ) &
                                                 + self%y_poss_cell_cent(i-1, self%jmax, k  ) &
                                                 + self%y_poss_cell_cent(i  , self%jmax, k-1) &
                                                 + self%y_poss_cell_cent(i-1, self%jmax, k-1))&
                                        - self%y_poss_cell_cent(self%jmax, j, k)
                    cell_to_face_vec(3) = 0.5d0 * (self%z_poss_cell_cent(i  , self%jmax, k  ) &
                                                 + self%z_poss_cell_cent(i-1, self%jmax, k  ) &
                                                 + self%z_poss_cell_cent(i  , self%jmax, k-1) &
                                                 + self%z_poss_cell_cent(i-1, self%jmax, k-1))&
                                        - self%z_poss_cell_cent(self%jmax, j, k)
                    centor_positions(n, 1) = self%x_poss_cell_cent(i, self%jmax, k) + 2.d0 * dble(j) * cell_to_face_vec(1)
                    centor_positions(n, 2) = self%y_poss_cell_cent(i, self%jmax, k) + 2.d0 * dble(j) * cell_to_face_vec(2)
                    centor_positions(n, 3) = self%z_poss_cell_cent(i, self%jmax, k) + 2.d0 * dble(j) * cell_to_face_vec(3)
                    volumes         (n)    = self%cell_vols(i, self%jmax, k)
                end do
            end do
        end do

        ! ghosts (z +/- direction)
        do k = 1, self%num_ghost_cells_, 1
            do j = self%jmin, self%jmax, 1
                do i = self%imin, self%imax, 1
                    ! z-
                    n = convert_structure_index_to_unstructure_index(i, j, self%kmin - k,                &
                        self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
                        self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
                    cell_to_face_vec(1) = 0.5d0 * (self%x_poss_cell_cent(i  , j  , self%kmin-1) &
                                                 + self%x_poss_cell_cent(i-1, j  , self%kmin-1) &
                                                 + self%x_poss_cell_cent(i  , j-1, self%kmin-1) &
                                                 + self%x_poss_cell_cent(i-1, j-1, self%kmin-1))&
                                        - self%x_poss_cell_cent(self%kmin, j, k)
                    cell_to_face_vec(2) = 0.5d0 * (self%y_poss_cell_cent(i  , j  , self%kmin-1) &
                                                 + self%y_poss_cell_cent(i-1, j  , self%kmin-1) &
                                                 + self%y_poss_cell_cent(i  , j-1, self%kmin-1) &
                                                 + self%y_poss_cell_cent(i-1, j-1, self%kmin-1))&
                                        - self%y_poss_cell_cent(self%kmin, j, k)
                    cell_to_face_vec(3) = 0.5d0 * (self%z_poss_cell_cent(i  , j  , self%kmin-1) &
                                                 + self%z_poss_cell_cent(i-1, j  , self%kmin-1) &
                                                 + self%z_poss_cell_cent(i  , j-1, self%kmin-1) &
                                                 + self%z_poss_cell_cent(i-1, j-1, self%kmin-1))&
                                        - self%z_poss_cell_cent(self%kmin, j, k)
                    centor_positions(n, 1) = self%x_poss_cell_cent(i, j, self%kmin) + 2.d0 * dble(k) * cell_to_face_vec(1)
                    centor_positions(n, 2) = self%y_poss_cell_cent(i, j, self%kmin) + 2.d0 * dble(k) * cell_to_face_vec(2)
                    centor_positions(n, 3) = self%z_poss_cell_cent(i, j, self%kmin) + 2.d0 * dble(k) * cell_to_face_vec(3)
                    volumes         (n)    = self%cell_vols(i, j, self%kmin)
                    ! z+
                    n = convert_structure_index_to_unstructure_index(i, j, self%kmax + k,                &
                        self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
                        self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
                    cell_to_face_vec(1) = 0.5d0 * (self%x_poss_cell_cent(i  , j  , self%kmax) &
                                                 + self%x_poss_cell_cent(i-1, j  , self%kmax) &
                                                 + self%x_poss_cell_cent(i  , j-1, self%kmax) &
                                                 + self%x_poss_cell_cent(i-1, j-1, self%kmax))&
                                        - self%x_poss_cell_cent(self%kmax, j, k)
                    cell_to_face_vec(2) = 0.5d0 * (self%y_poss_cell_cent(i  , j  , self%kmax) &
                                                 + self%y_poss_cell_cent(i-1, j  , self%kmax) &
                                                 + self%y_poss_cell_cent(i  , j-1, self%kmax) &
                                                 + self%y_poss_cell_cent(i-1, j-1, self%kmax))&
                                        - self%y_poss_cell_cent(self%kmax, j, k)
                    cell_to_face_vec(3) = 0.5d0 * (self%z_poss_cell_cent(i  , j  , self%kmax) &
                                                 + self%z_poss_cell_cent(i-1, j  , self%kmax) &
                                                 + self%z_poss_cell_cent(i  , j-1, self%kmax) &
                                                 + self%z_poss_cell_cent(i-1, j-1, self%kmax))&
                                        - self%z_poss_cell_cent(self%kmax, j, k)
                    centor_positions(n, 1) = self%x_poss_cell_cent(i, j, self%kmax) + 2.d0 * dble(j) * cell_to_face_vec(1)
                    centor_positions(n, 2) = self%y_poss_cell_cent(i, j, self%kmax) + 2.d0 * dble(j) * cell_to_face_vec(2)
                    centor_positions(n, 3) = self%z_poss_cell_cent(i, j, self%kmax) + 2.d0 * dble(j) * cell_to_face_vec(3)
                    volumes         (n)    = self%cell_vols(i, j, self%kmax)
                end do
            end do
        end do
    end subroutine get_cells

    subroutine get_faces(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas)
        class(nlgrid_parser), intent(in   ) :: self
        integer(int_kind     ), intent(inout) :: reference_cell_indexs(:,:)
        real   (real_kind    ), intent(inout) :: normal_vectors       (:,:)
        real   (real_kind    ), intent(inout) :: tangential1_vectors  (:,:)
        real   (real_kind    ), intent(inout) :: tangential2_vectors  (:,:)
        real   (real_kind    ), intent(inout) :: positions            (:,:)
        real   (real_kind    ), intent(inout) :: areas                (:)

        integer(int_kind) :: i, j, k
        integer(int_kind) :: face_index

        if(.not. self%parsed)then
            call call_error("'parse' method of nlgrid_parser is not called yet. But you call 'get_faces' method.")
        end if

        face_index = 1
        do k = self%kmin, self%kmax - 1, 1
            do j = self%jmin, self%jmax - 1, 1
                do i = self%imin, self%imax, 1
                    call assign_lower_x_face(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas, face_index, i, j, k)
                    call assign_lower_y_face(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas, face_index, i, j, k)
                    call assign_lower_z_face(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas, face_index, i, j, k)
                end do
                ! # imax boundary cells
                call assign_upper_x_face(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas, face_index, self%imax, j, k)
            end do
            ! # jmax boundary cells
            do i = self%imin, self%imax, 1
                call assign_lower_x_face(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas, face_index, i, self%jmax, k)
                call assign_lower_y_face(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas, face_index, i, self%jmax, k)
                call assign_lower_z_face(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas, face_index, i, self%jmax, k)
                call assign_upper_y_face(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas, face_index, i, self%jmax, k)
            end do
            call assign_upper_x_face(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas, face_index, self%imax, self%jmax, k)
        end do
        ! # kmax boundary cells
        do j = self%jmin, self%jmax - 1, 1
            do i = self%imin, self%imax, 1
                call assign_lower_x_face(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas, face_index, i, j, self%kmax)
                call assign_lower_y_face(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas, face_index, i, j, self%kmax)
                call assign_lower_z_face(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas, face_index, i, j, self%kmax)
                call assign_upper_z_face(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas, face_index, i, j, self%kmax)
            end do
            call assign_upper_x_face(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas, face_index, self%imax, j, self%kmax)
        end do
        do i = self%imin, self%imax, 1
            call assign_lower_x_face(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas, face_index, i, self%jmax, self%kmax)
            call assign_lower_y_face(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas, face_index, i, self%jmax, self%kmax)
            call assign_lower_z_face(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas, face_index, i, self%jmax, self%kmax)
            call assign_upper_y_face(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas, face_index, i, self%jmax, self%kmax)
            call assign_upper_z_face(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas, face_index, i, self%jmax, self%kmax)
        end do
        call assign_upper_x_face(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas, face_index, self%imax, self%jmax, self%kmax)
#ifdef _DEBUG
        print *, "DEBUG: nlgrid parser:"
        print *, "Faces (", face_index-1, "/", self%get_number_of_faces(), ") are assigned."
#endif
    end subroutine get_faces

    subroutine assign_upper_z_face(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas, &
        face_index, i, j, k)
        class  (nlgrid_parser), intent(in   ) :: self
        integer(int_kind     ), intent(inout) :: reference_cell_indexs(:,:)
        real   (real_kind    ), intent(inout) :: normal_vectors       (:,:)
        real   (real_kind    ), intent(inout) :: tangential1_vectors  (:,:)
        real   (real_kind    ), intent(inout) :: tangential2_vectors  (:,:)
        real   (real_kind    ), intent(inout) :: positions            (:,:)
        real   (real_kind    ), intent(inout) :: areas                (:)
        integer(int_kind     ), intent(inout) :: face_index
        integer(int_kind     ), intent(in   ) :: i, j, k

        ! ## cell i,j,k +z face
        normal_vectors     (face_index, 1:3) = self%normal_vecs_k_drct_face(i, j, k, 1:3)
        tangential1_vectors(face_index, 1:3) = self%tan1_vecs_k_drct_face  (i, j, k, 1:3)
        tangential2_vectors(face_index, 1:3) = self%tan2_vecs_k_drct_face  (i, j, k, 1:3)
        areas              (face_index     ) = self%areas_k_drct_face      (i, j, k     )
        positions          (face_index, 1  ) = 0.25d0 *( &
            + self%x_poss_cell_edge(i    , j    , k) &
            + self%x_poss_cell_edge(i - 1, j    , k) &
            + self%x_poss_cell_edge(i    , j - 1, k) &
            + self%x_poss_cell_edge(i - 1, j - 1, k) )
        positions          (face_index, 2  ) = 0.25d0 *( &
            + self%y_poss_cell_edge(i    , j    , k) &
            + self%y_poss_cell_edge(i - 1, j    , k) &
            + self%y_poss_cell_edge(i    , j - 1, k) &
            + self%y_poss_cell_edge(i - 1, j - 1, k) )
        positions          (face_index, 3  ) = 0.25d0 *( &
            + self%z_poss_cell_edge(i    , j    , k) &
            + self%z_poss_cell_edge(i - 1, j    , k) &
            + self%z_poss_cell_edge(i    , j - 1, k) &
            + self%z_poss_cell_edge(i - 1, j - 1, k) )
        reference_cell_indexs(face_index, 1) = convert_structure_index_to_unstructure_index(i, j, k - 2,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 2) = convert_structure_index_to_unstructure_index(i, j, k - 1,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 3) = convert_structure_index_to_unstructure_index(i, j, k,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 4) = convert_structure_index_to_unstructure_index(i, j, k + 1,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 5) = convert_structure_index_to_unstructure_index(i, j, k + 2,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 6) = convert_structure_index_to_unstructure_index(i, j, k + 3,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        face_index = face_index + 1
    end subroutine assign_upper_z_face

    subroutine assign_upper_y_face(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas, &
        face_index, i, j, k)
        class  (nlgrid_parser), intent(in   ) :: self
        integer(int_kind     ), intent(inout) :: reference_cell_indexs(:,:)
        real   (real_kind    ), intent(inout) :: normal_vectors       (:,:)
        real   (real_kind    ), intent(inout) :: tangential1_vectors  (:,:)
        real   (real_kind    ), intent(inout) :: tangential2_vectors  (:,:)
        real   (real_kind    ), intent(inout) :: positions            (:,:)
        real   (real_kind    ), intent(inout) :: areas                (:)
        integer(int_kind     ), intent(inout) :: face_index
        integer(int_kind     ), intent(in   ) :: i, j, k

        ! ## cell i,j,k +y face
        normal_vectors     (face_index, 1:3) = self%normal_vecs_j_drct_face(i, j, k, 1:3)
        tangential1_vectors(face_index, 1:3) = self%tan1_vecs_j_drct_face  (i, j, k, 1:3)
        tangential2_vectors(face_index, 1:3) = self%tan2_vecs_j_drct_face  (i, j, k, 1:3)
        areas              (face_index     ) = self%areas_j_drct_face      (i, j, k     )
        positions          (face_index, 1  ) = 0.25d0 *( &
            + self%x_poss_cell_edge(i    , j, k    ) &
            + self%x_poss_cell_edge(i - 1, j, k    ) &
            + self%x_poss_cell_edge(i    , j, k - 1) &
            + self%x_poss_cell_edge(i - 1, j, k - 1) )
        positions          (face_index, 2  ) = 0.25d0 *( &
            + self%y_poss_cell_edge(i    , j, k    ) &
            + self%y_poss_cell_edge(i - 1, j, k    ) &
            + self%y_poss_cell_edge(i    , j, k - 1) &
            + self%y_poss_cell_edge(i - 1, j, k - 1) )
        positions          (face_index, 3  ) = 0.25d0 *( &
            + self%z_poss_cell_edge(i    , j, k    ) &
            + self%z_poss_cell_edge(i - 1, j, k    ) &
            + self%z_poss_cell_edge(i    , j, k - 1) &
            + self%z_poss_cell_edge(i - 1, j, k - 1) )
        reference_cell_indexs(face_index, 1) = convert_structure_index_to_unstructure_index(i, j - 2, k,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 2) = convert_structure_index_to_unstructure_index(i, j - 1, k,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 3) = convert_structure_index_to_unstructure_index(i, j, k,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 4) = convert_structure_index_to_unstructure_index(i, j + 1, k,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 5) = convert_structure_index_to_unstructure_index(i, j + 2, k,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 6) = convert_structure_index_to_unstructure_index(i, j + 3, k,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        face_index = face_index + 1
    end subroutine assign_upper_y_face

    subroutine assign_upper_x_face(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas, &
        face_index, i, j, k)
        class  (nlgrid_parser), intent(in   ) :: self
        integer(int_kind     ), intent(inout) :: reference_cell_indexs(:,:)
        real   (real_kind    ), intent(inout) :: normal_vectors       (:,:)
        real   (real_kind    ), intent(inout) :: tangential1_vectors  (:,:)
        real   (real_kind    ), intent(inout) :: tangential2_vectors  (:,:)
        real   (real_kind    ), intent(inout) :: positions            (:,:)
        real   (real_kind    ), intent(inout) :: areas                (:)
        integer(int_kind     ), intent(inout) :: face_index
        integer(int_kind     ), intent(in   ) :: i, j, k

        ! ## cell i,j,k +x face
        normal_vectors     (face_index, 1:3) = self%normal_vecs_i_drct_face(i, j, k, 1:3)
        tangential1_vectors(face_index, 1:3) = self%tan1_vecs_i_drct_face  (i, j, k, 1:3)
        tangential2_vectors(face_index, 1:3) = self%tan2_vecs_i_drct_face  (i, j, k, 1:3)
        areas              (     face_index) = self%areas_i_drct_face      (i, j, k     )
        positions          (face_index, 1  ) = 0.25d0 *( &
            + self%x_poss_cell_edge(i, j    , k    ) &
            + self%x_poss_cell_edge(i, j - 1, k    ) &
            + self%x_poss_cell_edge(i, j    , k - 1) &
            + self%x_poss_cell_edge(i, j - 1, k - 1) )
        positions          (face_index, 2  ) = 0.25d0 *( &
            + self%y_poss_cell_edge(i, j    , k    ) &
            + self%y_poss_cell_edge(i, j - 1, k    ) &
            + self%y_poss_cell_edge(i, j    , k - 1) &
            + self%y_poss_cell_edge(i, j - 1, k - 1) )
        positions          (face_index, 3  ) = 0.25d0 *( &
            + self%z_poss_cell_edge(i, j    , k    ) &
            + self%z_poss_cell_edge(i, j - 1, k    ) &
            + self%z_poss_cell_edge(i, j    , k - 1) &
            + self%z_poss_cell_edge(i, j - 1, k - 1) )
        reference_cell_indexs(face_index, 1) = convert_structure_index_to_unstructure_index(i - 2, j, k,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 2) = convert_structure_index_to_unstructure_index(i - 1, j, k,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 3) = convert_structure_index_to_unstructure_index(i, j, k,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 4) = convert_structure_index_to_unstructure_index(i + 1, j, k,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 5) = convert_structure_index_to_unstructure_index(i + 2, j, k,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 6) = convert_structure_index_to_unstructure_index(i + 3, j, k,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        face_index = face_index + 1
    end subroutine assign_upper_x_face

    subroutine assign_lower_x_face(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas, &
        face_index, i, j, k)
        class  (nlgrid_parser), intent(in   ) :: self
        integer(int_kind     ), intent(inout) :: reference_cell_indexs(:,:)
        real   (real_kind    ), intent(inout) :: normal_vectors       (:,:)
        real   (real_kind    ), intent(inout) :: tangential1_vectors  (:,:)
        real   (real_kind    ), intent(inout) :: tangential2_vectors  (:,:)
        real   (real_kind    ), intent(inout) :: positions            (:,:)
        real   (real_kind    ), intent(inout) :: areas                (:)
        integer(int_kind     ), intent(inout) :: face_index
        integer(int_kind     ), intent(in   ) :: i, j, k
        ! # inner cells
        ! ## cell i,j,k -x face
        normal_vectors     (face_index, 1:3) = -1.d0 * self%normal_vecs_i_drct_face(i - 1, j, k, 1:3)
        tangential1_vectors(face_index, 1:3) = self%tan1_vecs_i_drct_face  (i - 1, j, k, 1:3)
        tangential2_vectors(face_index, 1:3) = self%tan2_vecs_i_drct_face  (i - 1, j, k, 1:3)
        areas              (     face_index) = self%areas_i_drct_face      (i - 1, j, k     )
        positions          (face_index, 1  ) = 0.25d0 *( &
            + self%x_poss_cell_edge(i - 1, j    , k    ) &
            + self%x_poss_cell_edge(i - 1, j - 1, k    ) &
            + self%x_poss_cell_edge(i - 1, j    , k - 1) &
            + self%x_poss_cell_edge(i - 1, j - 1, k - 1) )
        positions          (face_index, 2  ) = 0.25d0 *( &
            + self%y_poss_cell_edge(i - 1, j    , k    ) &
            + self%y_poss_cell_edge(i - 1, j - 1, k    ) &
            + self%y_poss_cell_edge(i - 1, j    , k - 1) &
            + self%y_poss_cell_edge(i - 1, j - 1, k - 1) )
        positions          (face_index, 3  ) = 0.25d0 *( &
            + self%z_poss_cell_edge(i - 1, j    , k    ) &
            + self%z_poss_cell_edge(i - 1, j - 1, k    ) &
            + self%z_poss_cell_edge(i - 1, j    , k - 1) &
            + self%z_poss_cell_edge(i - 1, j - 1, k - 1) )
        reference_cell_indexs(face_index, 1) = convert_structure_index_to_unstructure_index(i + 2, j, k,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 2) = convert_structure_index_to_unstructure_index(i + 1, j, k,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 3) = convert_structure_index_to_unstructure_index(i + 0, j, k,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 4) = convert_structure_index_to_unstructure_index(i - 1, j, k,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 5) = convert_structure_index_to_unstructure_index(i - 2, j, k,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 6) = convert_structure_index_to_unstructure_index(i - 3, j, k,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        face_index = face_index + 1
    end subroutine assign_lower_x_face

    subroutine assign_lower_y_face(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas, &
        face_index, i, j, k)
        class  (nlgrid_parser), intent(in   ) :: self
        integer(int_kind     ), intent(inout) :: reference_cell_indexs(:,:)
        real   (real_kind    ), intent(inout) :: normal_vectors       (:,:)
        real   (real_kind    ), intent(inout) :: tangential1_vectors  (:,:)
        real   (real_kind    ), intent(inout) :: tangential2_vectors  (:,:)
        real   (real_kind    ), intent(inout) :: positions            (:,:)
        real   (real_kind    ), intent(inout) :: areas                (:)
        integer(int_kind     ), intent(inout) :: face_index
        integer(int_kind     ), intent(in   ) :: i, j, k

        ! ## cell i,j,k -y face
        normal_vectors     (face_index, 1:3) = -1.d0 * self%normal_vecs_j_drct_face(i, j - 1, k, 1:3)
        tangential1_vectors(face_index, 1:3) = self%tan1_vecs_j_drct_face  (i, j - 1, k, 1:3)
        tangential2_vectors(face_index, 1:3) = self%tan2_vecs_j_drct_face  (i, j - 1, k, 1:3)
        areas              (     face_index) = self%areas_j_drct_face      (i, j - 1, k     )
        positions          (face_index, 1  ) = 0.25d0 *( &
            + self%x_poss_cell_edge(i    , j - 1, k    ) &
            + self%x_poss_cell_edge(i - 1, j - 1, k    ) &
            + self%x_poss_cell_edge(i    , j - 1, k - 1) &
            + self%x_poss_cell_edge(i - 1, j - 1, k - 1) )
        positions          (face_index, 2  ) = 0.25d0 *( &
            + self%y_poss_cell_edge(i    , j - 1, k    ) &
            + self%y_poss_cell_edge(i - 1, j - 1, k    ) &
            + self%y_poss_cell_edge(i    , j - 1, k - 1) &
            + self%y_poss_cell_edge(i - 1, j - 1, k - 1) )
        positions          (face_index, 3  ) = 0.25d0 *( &
            + self%z_poss_cell_edge(i    , j - 1, k    ) &
            + self%z_poss_cell_edge(i - 1, j - 1, k    ) &
            + self%z_poss_cell_edge(i    , j - 1, k - 1) &
            + self%z_poss_cell_edge(i - 1, j - 1, k - 1) )
        reference_cell_indexs(face_index, 1) = convert_structure_index_to_unstructure_index(i, j + 2, k,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 2) = convert_structure_index_to_unstructure_index(i, j + 1, k,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 3) = convert_structure_index_to_unstructure_index(i, j + 0, k,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 4) = convert_structure_index_to_unstructure_index(i, j - 1, k,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 5) = convert_structure_index_to_unstructure_index(i, j - 2, k,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 6) = convert_structure_index_to_unstructure_index(i, j - 3, k,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        face_index = face_index + 1
    end subroutine assign_lower_y_face

    subroutine assign_lower_z_face(self, reference_cell_indexs, normal_vectors, tangential1_vectors, tangential2_vectors, positions, areas, &
        face_index, i, j, k)
        class  (nlgrid_parser), intent(in   ) :: self
        integer(int_kind     ), intent(inout) :: reference_cell_indexs(:,:)
        real   (real_kind    ), intent(inout) :: normal_vectors       (:,:)
        real   (real_kind    ), intent(inout) :: tangential1_vectors  (:,:)
        real   (real_kind    ), intent(inout) :: tangential2_vectors  (:,:)
        real   (real_kind    ), intent(inout) :: positions            (:,:)
        real   (real_kind    ), intent(inout) :: areas                (:)
        integer(int_kind     ), intent(inout) :: face_index
        integer(int_kind     ), intent(in   ) :: i, j, k

        ! ## cell i,j,k -z face
        normal_vectors     (face_index, 1:3) = -1.d0 * self%normal_vecs_k_drct_face(i, j, k - 1, 1:3)
        tangential1_vectors(face_index, 1:3) = self%tan1_vecs_k_drct_face  (i, j, k - 1, 1:3)
        tangential2_vectors(face_index, 1:3) = self%tan2_vecs_k_drct_face  (i, j, k - 1, 1:3)
        areas              (     face_index) = self%areas_k_drct_face      (i, j, k - 1     )
        positions          (face_index, 1  ) = 0.25d0 *( &
            + self%x_poss_cell_edge(i    , j    , k - 1) &
            + self%x_poss_cell_edge(i - 1, j    , k - 1) &
            + self%x_poss_cell_edge(i    , j - 1, k - 1) &
            + self%x_poss_cell_edge(i - 1, j - 1, k - 1) )
        positions          (face_index, 2  ) = 0.25d0 *( &
            + self%y_poss_cell_edge(i    , j    , k - 1) &
            + self%y_poss_cell_edge(i - 1, j    , k - 1) &
            + self%y_poss_cell_edge(i    , j - 1, k - 1) &
            + self%y_poss_cell_edge(i - 1, j - 1, k - 1) )
        positions          (face_index, 3  ) = 0.25d0 *( &
            + self%z_poss_cell_edge(i    , j    , k - 1) &
            + self%z_poss_cell_edge(i - 1, j    , k - 1) &
            + self%z_poss_cell_edge(i    , j - 1, k - 1) &
            + self%z_poss_cell_edge(i - 1, j - 1, k - 1) )
        reference_cell_indexs(face_index, 1) = convert_structure_index_to_unstructure_index(i, j, k + 2,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 2) = convert_structure_index_to_unstructure_index(i, j, k + 1,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 3) = convert_structure_index_to_unstructure_index(i, j, k + 0,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 4) = convert_structure_index_to_unstructure_index(i, j, k - 1,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 5) = convert_structure_index_to_unstructure_index(i, j, k - 2,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        reference_cell_indexs(face_index, 6) = convert_structure_index_to_unstructure_index(i, j, k - 3,&
            self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
            self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
        face_index = face_index + 1
    end subroutine assign_lower_z_face

    subroutine get_boundaries(self, outflow_face_indexs, slipwall_face_indexs, symmetric_face_indexs)
        class  (nlgrid_parser), intent(in   ) :: self
        integer(int_kind     ), intent(inout) :: outflow_face_indexs  (:)
        integer(int_kind     ), intent(inout) :: slipwall_face_indexs (:)
        integer(int_kind     ), intent(inout) :: symmetric_face_indexs(:)

        integer(int_kind) :: i, j, k
        integer(int_kind) :: outflow_index, slipwall_index, symmetric_index, face_index

        if(.not. self%parsed)then
            call call_error("'parse' method of nlgrid_parser is not called yet. But you call 'get_boundaries' method.")
        end if

        outflow_index   = 1
        slipwall_index  = 1
        symmetric_index = 1

        face_index = 1
        do k = self%kmin, self%kmax - 1, 1
            do j = self%jmin, self%jmax - 1, 1
                do i = self%imin, self%imax, 1
                    if(i == self%imin)then
                        if(self%x_minus_direction == outflow_boundary_type  ) call self%assign_boundary(outflow_face_indexs  , outflow_index  , face_index)
                        if(self%x_minus_direction == symmetric_boundary_type) call self%assign_boundary(symmetric_face_indexs, symmetric_index, face_index)
                        if(self%x_minus_direction == slipwall_boundary_type ) call self%assign_boundary(slipwall_face_indexs , slipwall_index , face_index)
                    end if
                    face_index = face_index + 1
                    if(j == self%jmin)then
                        if(self%y_minus_direction == outflow_boundary_type  ) call self%assign_boundary(outflow_face_indexs  , outflow_index  , face_index)
                        if(self%y_minus_direction == symmetric_boundary_type) call self%assign_boundary(symmetric_face_indexs, symmetric_index, face_index)
                        if(self%y_minus_direction == slipwall_boundary_type ) call self%assign_boundary(slipwall_face_indexs , slipwall_index , face_index)
                    end if
                    face_index = face_index + 1
                    if(k == self%kmin)then
                        if(self%z_minus_direction == outflow_boundary_type  ) call self%assign_boundary(outflow_face_indexs  , outflow_index  , face_index)
                        if(self%z_minus_direction == symmetric_boundary_type) call self%assign_boundary(symmetric_face_indexs, symmetric_index, face_index)
                        if(self%z_minus_direction == slipwall_boundary_type ) call self%assign_boundary(slipwall_face_indexs , slipwall_index , face_index)
                    end if
                    face_index = face_index + 1
                end do
                ! # imax boundary cells
                if(self%x_plus_direction == outflow_boundary_type  ) call self%assign_boundary(outflow_face_indexs  , outflow_index  , face_index)
                if(self%x_plus_direction == symmetric_boundary_type) call self%assign_boundary(symmetric_face_indexs, symmetric_index, face_index)
                if(self%x_plus_direction == slipwall_boundary_type ) call self%assign_boundary(slipwall_face_indexs , slipwall_index , face_index)
                face_index = face_index + 1
            end do
            ! # jmax boundary cells
            do i = self%imin, self%imax, 1
                if(i == self%imin)then
                    if(self%x_minus_direction == outflow_boundary_type  ) call self%assign_boundary(outflow_face_indexs  , outflow_index  , face_index)
                    if(self%x_minus_direction == symmetric_boundary_type) call self%assign_boundary(symmetric_face_indexs, symmetric_index, face_index)
                    if(self%x_minus_direction == slipwall_boundary_type ) call self%assign_boundary(slipwall_face_indexs , slipwall_index , face_index)
                end if
                face_index = face_index + 1
                if(self%jmin == self%jmax)then
                    if(self%y_minus_direction == outflow_boundary_type  ) call self%assign_boundary(outflow_face_indexs  , outflow_index  , face_index)
                    if(self%y_minus_direction == symmetric_boundary_type) call self%assign_boundary(symmetric_face_indexs, symmetric_index, face_index)
                    if(self%y_minus_direction == slipwall_boundary_type ) call self%assign_boundary(slipwall_face_indexs , slipwall_index , face_index)
                end if
                face_index = face_index + 1
                if(k == self%kmin)then
                    if(self%z_minus_direction == outflow_boundary_type  ) call self%assign_boundary(outflow_face_indexs  , outflow_index  , face_index)
                    if(self%z_minus_direction == symmetric_boundary_type) call self%assign_boundary(symmetric_face_indexs, symmetric_index, face_index)
                    if(self%z_minus_direction == slipwall_boundary_type ) call self%assign_boundary(slipwall_face_indexs , slipwall_index , face_index)
                end if
                face_index = face_index + 1
                if(self%y_plus_direction == outflow_boundary_type  ) call self%assign_boundary(outflow_face_indexs  , outflow_index  , face_index)
                if(self%y_plus_direction == symmetric_boundary_type) call self%assign_boundary(symmetric_face_indexs, symmetric_index, face_index)
                if(self%y_plus_direction == slipwall_boundary_type ) call self%assign_boundary(slipwall_face_indexs , slipwall_index , face_index)
                face_index = face_index + 1
            end do
            if(self%x_plus_direction == outflow_boundary_type  ) call self%assign_boundary(outflow_face_indexs  , outflow_index, face_index)
            if(self%x_plus_direction == symmetric_boundary_type) call self%assign_boundary(symmetric_face_indexs , symmetric_index, face_index)
            if(self%x_plus_direction == slipwall_boundary_type ) call self%assign_boundary(slipwall_face_indexs, slipwall_index, face_index)
            face_index = face_index + 1
        end do
        ! # kmax boundary cells
        do j = self%jmin, self%jmax - 1, 1
            do i = self%imin, self%imax, 1
                if(i == self%imin)then
                    if(self%x_minus_direction == outflow_boundary_type  ) call self%assign_boundary(outflow_face_indexs  , outflow_index, face_index)
                    if(self%x_minus_direction == symmetric_boundary_type) call self%assign_boundary(symmetric_face_indexs , symmetric_index, face_index)
                    if(self%x_minus_direction == slipwall_boundary_type ) call self%assign_boundary(slipwall_face_indexs, slipwall_index, face_index)
                end if
                face_index = face_index + 1
                if(j == self%jmin)then
                    if(self%y_minus_direction == outflow_boundary_type  ) call self%assign_boundary(outflow_face_indexs  , outflow_index, face_index)
                    if(self%y_minus_direction == symmetric_boundary_type) call self%assign_boundary(symmetric_face_indexs , symmetric_index, face_index)
                    if(self%y_minus_direction == slipwall_boundary_type ) call self%assign_boundary(slipwall_face_indexs, slipwall_index, face_index)
                end if
                face_index = face_index + 1
                if(self%kmin == self%kmax)then
                    if(self%z_minus_direction == outflow_boundary_type  ) call self%assign_boundary(outflow_face_indexs  , outflow_index  , face_index)
                    if(self%z_minus_direction == symmetric_boundary_type) call self%assign_boundary(symmetric_face_indexs, symmetric_index, face_index)
                    if(self%z_minus_direction == slipwall_boundary_type ) call self%assign_boundary(slipwall_face_indexs , slipwall_index , face_index)
                end if
                face_index = face_index + 1
                if(self%z_plus_direction == outflow_boundary_type  ) call self%assign_boundary(outflow_face_indexs  , outflow_index, face_index)
                if(self%z_plus_direction == symmetric_boundary_type) call self%assign_boundary(symmetric_face_indexs , symmetric_index, face_index)
                if(self%z_plus_direction == slipwall_boundary_type ) call self%assign_boundary(slipwall_face_indexs, slipwall_index, face_index)
                face_index = face_index + 1
            end do
            if(self%x_plus_direction == outflow_boundary_type  ) call self%assign_boundary(outflow_face_indexs  , outflow_index, face_index)
            if(self%x_plus_direction == symmetric_boundary_type) call self%assign_boundary(symmetric_face_indexs , symmetric_index, face_index)
            if(self%x_plus_direction == slipwall_boundary_type ) call self%assign_boundary(slipwall_face_indexs, slipwall_index, face_index)
            face_index = face_index + 1
        end do
        do i = self%imin, self%imax, 1
            if(i == self%imin)then
                if(self%x_minus_direction == outflow_boundary_type  ) call self%assign_boundary(outflow_face_indexs  , outflow_index, face_index)
                if(self%x_minus_direction == symmetric_boundary_type) call self%assign_boundary(symmetric_face_indexs , symmetric_index, face_index)
                if(self%x_minus_direction == slipwall_boundary_type ) call self%assign_boundary(slipwall_face_indexs, slipwall_index, face_index)
            end if
            face_index = face_index + 1
            if(self%jmin == self%jmax)then
                if(self%y_minus_direction == outflow_boundary_type  ) call self%assign_boundary(outflow_face_indexs  , outflow_index  , face_index)
                if(self%y_minus_direction == symmetric_boundary_type) call self%assign_boundary(symmetric_face_indexs, symmetric_index, face_index)
                if(self%y_minus_direction == slipwall_boundary_type ) call self%assign_boundary(slipwall_face_indexs , slipwall_index , face_index)
            end if
            face_index = face_index + 1
            if(self%kmin == self%kmax)then
                if(self%z_minus_direction == outflow_boundary_type  ) call self%assign_boundary(outflow_face_indexs  , outflow_index  , face_index)
                if(self%z_minus_direction == symmetric_boundary_type) call self%assign_boundary(symmetric_face_indexs, symmetric_index, face_index)
                if(self%z_minus_direction == slipwall_boundary_type ) call self%assign_boundary(slipwall_face_indexs , slipwall_index , face_index)
            end if
            face_index = face_index + 1
            if(self%y_plus_direction == outflow_boundary_type  ) call self%assign_boundary(outflow_face_indexs  , outflow_index, face_index)
            if(self%y_plus_direction == symmetric_boundary_type) call self%assign_boundary(symmetric_face_indexs , symmetric_index, face_index)
            if(self%y_plus_direction == slipwall_boundary_type ) call self%assign_boundary(slipwall_face_indexs, slipwall_index, face_index)
            face_index = face_index + 1
            if(self%z_plus_direction == outflow_boundary_type  ) call self%assign_boundary(outflow_face_indexs  , outflow_index, face_index)
            if(self%z_plus_direction == symmetric_boundary_type) call self%assign_boundary(symmetric_face_indexs , symmetric_index, face_index)
            if(self%z_plus_direction == slipwall_boundary_type ) call self%assign_boundary(slipwall_face_indexs, slipwall_index, face_index)
            face_index = face_index + 1
        end do
        if(self%x_plus_direction == outflow_boundary_type  ) call self%assign_boundary(outflow_face_indexs  , outflow_index, face_index)
        if(self%x_plus_direction == symmetric_boundary_type) call self%assign_boundary(symmetric_face_indexs , symmetric_index, face_index)
        if(self%x_plus_direction == slipwall_boundary_type ) call self%assign_boundary(slipwall_face_indexs, slipwall_index, face_index)
#ifdef _DEBUG
        print *, "DEBUG: nlgrid parser:"
        print *, "Outflow-faces   (", outflow_index  -1, "/", self%get_number_of_outflow_faces  (), ") are assigned."
        print *, "Symmetric-faces (", symmetric_index-1, "/", self%get_number_of_symmetric_faces(), ") are assigned."
        print *, "Slipwall-faces  (", slipwall_index -1, "/", self%get_number_of_slipwall_faces (), ") are assigned."
#endif
    end subroutine get_boundaries

    subroutine assign_boundary(self, boundary_face_indexs, boundary_index, face_index)
        class  (nlgrid_parser), intent(in   ) :: self
        integer(int_kind     ), intent(inout) :: boundary_face_indexs(:)
        integer(int_kind     ), intent(inout) :: boundary_index
        integer(int_kind     ), intent(in   ) :: face_index

        boundary_face_indexs(boundary_index) = face_index
        boundary_index = boundary_index + 1
    end subroutine assign_boundary

    subroutine get_cell_geometries(self, points, cell_geometries)
        class  (nlgrid_parser), intent(in   ) :: self
        real   (real_kind    ), intent(inout) :: points         (:, :)
        class  (point_id_list), intent(inout) :: cell_geometries(:)

        integer(int_kind) :: i, j, k, n, n_assigned_point
        integer(int_kind) :: p(8)

        if(.not. self%parsed)then
            call call_error("'parse' method of nlgrid_parser is not called yet. But you call 'get_cell_geometries' method.")
        end if

        ! point index loop
        n_assigned_point = 0
        do k = self%kmin - 1, self%kmax, 1
            do j = self%jmin - 1, self%jmax, 1
                do i = self%imin - 1, self%imax, 1
                    n = convert_structure_index_to_unstructure_index(i, j, k, &
                        self%imin - 1, self%jmin - 1, self%kmin - 1,          &
                        self%imax    , self%jmax    , self%kmax                )
                    points(n, 1) = self%x_poss_cell_edge(i, j, k)
                    points(n, 2) = self%y_poss_cell_edge(i, j, k)
                    points(n, 3) = self%z_poss_cell_edge(i, j, k)
                    n_assigned_point = n_assigned_point + 1
                end do
            end do
        end do
#ifdef _DEBUG
        print *, "DEBUG: nlgrid parser:"
        print *, "Points (", n_assigned_point, "/", self%get_number_of_points(), ") are assigned."
#endif

        ! cell index loop
        do k = self%kmin, self%kmax, 1
            do j = self%jmin, self%jmax, 1
                do i = self%imin, self%imax, 1
                    n = convert_structure_index_to_unstructure_index(i, j, k,                                                    &
                        self%imin - self%num_ghost_cells_, self%jmin - self%num_ghost_cells_, self%kmin - self%num_ghost_cells_, &
                        self%imax + self%num_ghost_cells_, self%jmax + self%num_ghost_cells_, self%kmax + self%num_ghost_cells_   )
                    call cell_geometries(n)%initialize(8)

                    ! Numbering is followed VTK format (https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf)
                    p(1) = convert_structure_index_to_unstructure_index(i-1, j-1, k-1, &
                        self%imin - 1, self%jmin - 1, self%kmin - 1, self%imax, self%jmax, self%kmax)
                    p(2) = convert_structure_index_to_unstructure_index(i  , j-1, k-1, &
                        self%imin - 1, self%jmin - 1, self%kmin - 1, self%imax, self%jmax, self%kmax)
                    p(3) = convert_structure_index_to_unstructure_index(i  , j  , k-1, &
                        self%imin - 1, self%jmin - 1, self%kmin - 1, self%imax, self%jmax, self%kmax)
                    p(4) = convert_structure_index_to_unstructure_index(i-1, j  , k-1, &
                        self%imin - 1, self%jmin - 1, self%kmin - 1, self%imax, self%jmax, self%kmax)
                    p(5) = convert_structure_index_to_unstructure_index(i-1, j-1, k  , &
                        self%imin - 1, self%jmin - 1, self%kmin - 1, self%imax, self%jmax, self%kmax)
                    p(6) = convert_structure_index_to_unstructure_index(i  , j-1, k  , &
                        self%imin - 1, self%jmin - 1, self%kmin - 1, self%imax, self%jmax, self%kmax)
                    p(7) = convert_structure_index_to_unstructure_index(i  , j  , k  , &
                        self%imin - 1, self%jmin - 1, self%kmin - 1, self%imax, self%jmax, self%kmax)
                    p(8) = convert_structure_index_to_unstructure_index(i-1, j  , k  , &
                        self%imin - 1, self%jmin - 1, self%kmin - 1, self%imax, self%jmax, self%kmax)

                    call cell_geometries(n)%set_point_id(1, p(1))
                    call cell_geometries(n)%set_point_id(2, p(2))
                    call cell_geometries(n)%set_point_id(3, p(3))
                    call cell_geometries(n)%set_point_id(4, p(4))
                    call cell_geometries(n)%set_point_id(5, p(5))
                    call cell_geometries(n)%set_point_id(6, p(6))
                    call cell_geometries(n)%set_point_id(7, p(7))
                    call cell_geometries(n)%set_point_id(8, p(8))
                end do
            end do
        end do
    end subroutine get_cell_geometries
end module class_nlgrid_parser