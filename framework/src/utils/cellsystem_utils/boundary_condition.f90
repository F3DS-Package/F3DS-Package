module boundary_condition_module
    use typedef_module
    use vector_module
    use matrix_module

    implicit none

    private

    public :: apply_boundary_condition_common_impl
    public :: apply_boundary_condition_empty_impl

    contains

    subroutine apply_boundary_condition_common_impl(    &
            primitive_variables_set                   , &
            face_normal_vector                        , &
            face_tangential1_vector                   , &
            face_tangential2_vector                   , &
            face_to_cell_indexes                      , &
            num_local_cells                           , &
            num_primitive_variables                   , &
            compute_rotate_primitive_variables_function  , &
            compute_unrotate_primitive_variables_function, &
            boundary_condition_function                 &
        )
        real   (real_kind), intent(inout) :: primitive_variables_set    (:,:)
        real   (real_kind), intent(in   ) :: face_normal_vector         (3)
        real   (real_kind), intent(in   ) :: face_tangential1_vector    (3)
        real   (real_kind), intent(in   ) :: face_tangential2_vector    (3)
        integer(int_kind ), intent(in   ) :: face_to_cell_indexes       (:)
        integer(int_kind ), intent(in   ) :: num_local_cells
        integer(int_kind ), intent(in   ) :: num_primitive_variables
        interface
            pure function compute_rotate_primitive_variables_function( &
                primitive_variables                                  , &
                face_normal_vector                                   , &
                face_tangential1_vector                              , &
                face_tangential2_vector                              , &
                num_primitive_variables                                  ) result(rotate_primitive_variables)

                use typedef_module
                real   (real_kind     ), intent(in)  :: primitive_variables     (:)
                real   (real_kind     ), intent(in)  :: face_normal_vector      (3)
                real   (real_kind     ), intent(in)  :: face_tangential1_vector (3)
                real   (real_kind     ), intent(in)  :: face_tangential2_vector (3)
                integer(int_kind      ), intent(in)  :: num_primitive_variables
                real   (real_kind     )              :: rotate_primitive_variables(num_primitive_variables)
            end function compute_rotate_primitive_variables_function

            pure function compute_unrotate_primitive_variables_function( &
                primitive_variables                                    , &
                face_normal_vector                                     , &
                face_tangential1_vector                                , &
                face_tangential2_vector                                , &
                num_primitive_variables                                    ) result(unrotate_primitive_variables)

                use typedef_module
                real   (real_kind     ), intent(in)  :: primitive_variables     (:)
                real   (real_kind     ), intent(in)  :: face_normal_vector      (3)
                real   (real_kind     ), intent(in)  :: face_tangential1_vector (3)
                real   (real_kind     ), intent(in)  :: face_tangential2_vector (3)
                integer(int_kind      ), intent(in)  :: num_primitive_variables
                real   (real_kind     )              :: unrotate_primitive_variables(num_primitive_variables)
            end function compute_unrotate_primitive_variables_function

            pure function boundary_condition_function(inner_primitive_variables, num_primitive_variables) result(ghost_primitive_variables)
                use typedef_module
                real   (real_kind), intent(in) :: inner_primitive_variables(:)
                integer(int_kind ), intent(in) :: num_primitive_variables
                real   (real_kind)             :: ghost_primitive_variables(num_primitive_variables)
            end function boundary_condition_function
        end interface

        integer(int_kind ) :: i


        do i = 1, num_local_cells, 1
            associate(                                                                                      &
                    inner_values => primitive_variables_set(:,face_to_cell_indexes(num_local_cells-(i-1))), &
                    ghost_values => primitive_variables_set(:,face_to_cell_indexes(num_local_cells+(i  )))  &
                )
                ghost_values(:) = compute_unrotate_primitive_variables_function( &
                    boundary_condition_function(                                 &
                        compute_rotate_primitive_variables_function(             &
                            inner_values           ,                             &
                            face_normal_vector     ,                             &
                            face_tangential1_vector,                             &
                            face_tangential2_vector,                             &
                            num_primitive_variables                              &
                        ),                                                       &
                        num_primitive_variables                                  &
                    ),                                                           &
                    face_normal_vector     ,                                     &
                    face_tangential1_vector,                                     &
                    face_tangential2_vector,                                     &
                    num_primitive_variables                                      &
                )
            end associate
        end do
    end subroutine apply_boundary_condition_common_impl

    subroutine apply_boundary_condition_empty_impl(     &
            primitive_variables_set                   , &
            face_normal_vector                        , &
            face_tangential1_vector                   , &
            face_tangential2_vector                   , &
            face_to_cell_indexes                      , &
            num_local_cells                           , &
            num_primitive_variables                   , &
            compute_rotate_primitive_variables_function  , &
            compute_unrotate_primitive_variables_function, &
            boundary_condition_function                 &
        )
        real   (real_kind), intent(inout) :: primitive_variables_set    (:,:)
        real   (real_kind), intent(in   ) :: face_normal_vector         (3)
        real   (real_kind), intent(in   ) :: face_tangential1_vector    (3)
        real   (real_kind), intent(in   ) :: face_tangential2_vector    (3)
        integer(int_kind ), intent(in   ) :: face_to_cell_indexes       (:)
        integer(int_kind ), intent(in   ) :: num_local_cells
        integer(int_kind ), intent(in   ) :: num_primitive_variables
        interface
            pure function compute_rotate_primitive_variables_function( &
                primitive_variables                                  , &
                face_normal_vector                                   , &
                face_tangential1_vector                              , &
                face_tangential2_vector                              , &
                num_primitive_variables                                  ) result(rotate_primitive_variables)

                use typedef_module
                real   (real_kind     ), intent(in)  :: primitive_variables     (:)
                real   (real_kind     ), intent(in)  :: face_normal_vector      (3)
                real   (real_kind     ), intent(in)  :: face_tangential1_vector (3)
                real   (real_kind     ), intent(in)  :: face_tangential2_vector (3)
                integer(int_kind      ), intent(in)  :: num_primitive_variables
                real   (real_kind     )              :: rotate_primitive_variables(num_primitive_variables)
            end function compute_rotate_primitive_variables_function

            pure function compute_unrotate_primitive_variables_function( &
                primitive_variables                                    , &
                face_normal_vector                                     , &
                face_tangential1_vector                                , &
                face_tangential2_vector                                , &
                num_primitive_variables                                    ) result(unrotate_primitive_variables)

                use typedef_module
                real   (real_kind     ), intent(in)  :: primitive_variables     (:)
                real   (real_kind     ), intent(in)  :: face_normal_vector      (3)
                real   (real_kind     ), intent(in)  :: face_tangential1_vector (3)
                real   (real_kind     ), intent(in)  :: face_tangential2_vector (3)
                integer(int_kind      ), intent(in)  :: num_primitive_variables
                real   (real_kind     )              :: unrotate_primitive_variables(num_primitive_variables)
            end function compute_unrotate_primitive_variables_function

            pure function boundary_condition_function(inner_primitive_variables, num_primitive_variables) result(ghost_primitive_variables)
                use typedef_module
                real   (real_kind), intent(in) :: inner_primitive_variables(:)
                integer(int_kind ), intent(in) :: num_primitive_variables
                real   (real_kind)             :: ghost_primitive_variables(num_primitive_variables)
            end function boundary_condition_function
        end interface

        integer(int_kind ) :: i


        do i = 1, num_local_cells, 1
            associate(                                                                                      &
                    inner_values => primitive_variables_set(:,face_to_cell_indexes(num_local_cells-0)), &
                    ghost_values => primitive_variables_set(:,face_to_cell_indexes(num_local_cells+i))  &
                )
                ghost_values(:) = compute_unrotate_primitive_variables_function( &
                    boundary_condition_function(                                 &
                        compute_rotate_primitive_variables_function(             &
                            inner_values           ,                             &
                            face_normal_vector     ,                             &
                            face_tangential1_vector,                             &
                            face_tangential2_vector,                             &
                            num_primitive_variables                              &
                        ),                                                       &
                        num_primitive_variables                                  &
                    ),                                                           &
                    face_normal_vector     ,                                     &
                    face_tangential1_vector,                                     &
                    face_tangential2_vector,                                     &
                    num_primitive_variables                                      &
                )
            end associate
        end do
    end subroutine apply_boundary_condition_empty_impl
end module boundary_condition_module