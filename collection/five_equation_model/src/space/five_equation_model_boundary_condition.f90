module five_equation_model_boundary_condition_module
    use typedef_module
    use vector_module

    implicit none

    private

    public :: five_equation_model_set_boundary_condition

    contains

    function five_equation_model_set_boundary_condition( &
            primitive_variables_set   , &
            face_to_cell_index , &
            face_normal_vectors       , &
            face_tangential1_vectors  , &
            face_tangential2_vectors  , &
            outflow_face_indexs       , &
            slipwall_face_indexs      , &
            symmetric_face_indexs     , &
            num_outflow_faces         , &
            num_slipwall_faces        , &
            num_symmetric_faces       , &
            n_ghost_cells               &
        ) result(error)
        real   (real_kind), intent(inout) :: primitive_variables_set    (:,:)
        integer(int_kind ), intent(in   ) :: face_to_cell_index  (:,:)
        real   (real_kind), intent(in   ) :: face_normal_vectors        (:,:)
        real   (real_kind), intent(in   ) :: face_tangential1_vectors   (:,:)
        real   (real_kind), intent(in   ) :: face_tangential2_vectors   (:,:)
        integer(int_kind ), intent(in   ) :: outflow_face_indexs        (:)
        integer(int_kind ), intent(in   ) :: slipwall_face_indexs       (:)
        integer(int_kind ), intent(in   ) :: symmetric_face_indexs      (:)
        integer(int_kind ), intent(in   ) :: num_outflow_faces
        integer(int_kind ), intent(in   ) :: num_slipwall_faces
        integer(int_kind ), intent(in   ) :: num_symmetric_faces
        integer(int_kind ), intent(in   ) :: n_ghost_cells
        logical                           :: error

        integer(int_kind ) :: bc_face_index

        !$omp parallel do private(bc_face_index)
        do bc_face_index = 1, num_outflow_faces, 1
            associate( &
                face_idx   => outflow_face_indexs(bc_face_index)                                            , &
                ghost_idx1 => face_to_cell_index(outflow_face_indexs(bc_face_index), n_ghost_cells+1), &
                ghost_idx2 => face_to_cell_index(outflow_face_indexs(bc_face_index), n_ghost_cells+2), &
                ghost_idx3 => face_to_cell_index(outflow_face_indexs(bc_face_index), n_ghost_cells+3), &
                inner_idx1 => face_to_cell_index(outflow_face_indexs(bc_face_index), n_ghost_cells+0), &
                inner_idx2 => face_to_cell_index(outflow_face_indexs(bc_face_index), n_ghost_cells-1), &
                inner_idx3 => face_to_cell_index(outflow_face_indexs(bc_face_index), n_ghost_cells-2)  &
            )
                primitive_variables_set(ghost_idx1, :) = make_ghost_primitive_variables_outflow( &
                    primitive_variables_set , &
                    face_normal_vectors     , &
                    face_tangential1_vectors, &
                    face_tangential2_vectors, &
                    face_idx                , &
                    inner_idx1                &
                )
                primitive_variables_set(ghost_idx2, :) = make_ghost_primitive_variables_outflow( &
                    primitive_variables_set , &
                    face_normal_vectors     , &
                    face_tangential1_vectors, &
                    face_tangential2_vectors, &
                    face_idx                , &
                    inner_idx1                &
                )
                primitive_variables_set(ghost_idx3, :) = make_ghost_primitive_variables_outflow( &
                    primitive_variables_set , &
                    face_normal_vectors     , &
                    face_tangential1_vectors, &
                    face_tangential2_vectors, &
                    face_idx                , &
                    inner_idx1                &
                )
            end associate
        end do

        !$omp parallel do private(bc_face_index)
        do bc_face_index = 1, num_slipwall_faces, 1
            associate( &
                face_idx   => slipwall_face_indexs(bc_face_index)                                            , &
                ghost_idx1 => face_to_cell_index(slipwall_face_indexs(bc_face_index), n_ghost_cells+1), &
                ghost_idx2 => face_to_cell_index(slipwall_face_indexs(bc_face_index), n_ghost_cells+2), &
                ghost_idx3 => face_to_cell_index(slipwall_face_indexs(bc_face_index), n_ghost_cells+3), &
                inner_idx1 => face_to_cell_index(slipwall_face_indexs(bc_face_index), n_ghost_cells+0), &
                inner_idx2 => face_to_cell_index(slipwall_face_indexs(bc_face_index), n_ghost_cells-1), &
                inner_idx3 => face_to_cell_index(slipwall_face_indexs(bc_face_index), n_ghost_cells-2)  &
            )
                primitive_variables_set(ghost_idx1, :) = make_ghost_primitive_variables_slipwall( &
                    primitive_variables_set , &
                    face_normal_vectors     , &
                    face_tangential1_vectors, &
                    face_tangential2_vectors, &
                    face_idx                , &
                    inner_idx1                &
                )
                primitive_variables_set(ghost_idx2, :) = make_ghost_primitive_variables_slipwall( &
                    primitive_variables_set , &
                    face_normal_vectors     , &
                    face_tangential1_vectors, &
                    face_tangential2_vectors, &
                    face_idx                , &
                    inner_idx2                &
                )
                primitive_variables_set(ghost_idx3, :) = make_ghost_primitive_variables_slipwall( &
                    primitive_variables_set , &
                    face_normal_vectors     , &
                    face_tangential1_vectors, &
                    face_tangential2_vectors, &
                    face_idx                , &
                    inner_idx3                &
                )
            end associate
        end do

        !$omp parallel do private(bc_face_index)
        do bc_face_index = 1, num_symmetric_faces, 1
            associate( &
                face_idx   => symmetric_face_indexs(bc_face_index)                                            , &
                ghost_idx1 => face_to_cell_index(symmetric_face_indexs(bc_face_index), n_ghost_cells+1), &
                ghost_idx2 => face_to_cell_index(symmetric_face_indexs(bc_face_index), n_ghost_cells+2), &
                ghost_idx3 => face_to_cell_index(symmetric_face_indexs(bc_face_index), n_ghost_cells+3), &
                inner_idx1 => face_to_cell_index(symmetric_face_indexs(bc_face_index), n_ghost_cells+0), &
                inner_idx2 => face_to_cell_index(symmetric_face_indexs(bc_face_index), n_ghost_cells-1), &
                inner_idx3 => face_to_cell_index(symmetric_face_indexs(bc_face_index), n_ghost_cells-2)  &
            )
                primitive_variables_set(ghost_idx1, :) = make_ghost_primitive_variables_symmetric( &
                    primitive_variables_set , &
                    face_normal_vectors     , &
                    face_tangential1_vectors, &
                    face_tangential2_vectors, &
                    face_idx                , &
                    inner_idx1                &
                )
                primitive_variables_set(ghost_idx2, :) = make_ghost_primitive_variables_symmetric( &
                    primitive_variables_set , &
                    face_normal_vectors     , &
                    face_tangential1_vectors, &
                    face_tangential2_vectors, &
                    face_idx                , &
                    inner_idx1                &
                )
                primitive_variables_set(ghost_idx3, :) = make_ghost_primitive_variables_symmetric( &
                    primitive_variables_set , &
                    face_normal_vectors     , &
                    face_tangential1_vectors, &
                    face_tangential2_vectors, &
                    face_idx                , &
                    inner_idx1                &
                )
            end associate
        end do

        error = .false.
    end function five_equation_model_set_boundary_condition

    pure function make_ghost_primitive_variables_outflow( &
            primitive_variables_set   , &
            face_normal_vectors       , &
            face_tangential1_vectors  , &
            face_tangential2_vectors  , &
            face_index                , &
            inner_cell_index              ) result(ghost_primitive_variables)
        real   (real_kind), intent(in   ) :: primitive_variables_set    (:,:)
        real   (real_kind), intent(in   ) :: face_normal_vectors        (:,:)
        real   (real_kind), intent(in   ) :: face_tangential1_vectors   (:,:)
        real   (real_kind), intent(in   ) :: face_tangential2_vectors   (:,:)
        integer(int_kind ), intent(in   ) :: face_index
        integer(int_kind ), intent(in   ) :: inner_cell_index
        real   (real_kind) :: ghost_primitive_variables(7)
        real   (real_kind) :: local_inner_vector(3), local_ghost_vector(3)
        ghost_primitive_variables(1:2) = primitive_variables_set(inner_cell_index, 1:2)
        local_inner_vector(:) = vector_rotate(           &
            primitive_variables_set   (inner_cell_index, 3:5), &
            face_normal_vectors       (face_index      ,  : ), &
            face_tangential1_vectors  (face_index      ,  : ), &
            face_tangential2_vectors  (face_index      ,  : )  &
        )
        local_ghost_vector(1) = local_inner_vector(1)
        local_ghost_vector(2) = local_inner_vector(2)
        local_ghost_vector(3) = local_inner_vector(3)
        ghost_primitive_variables(3:5) = vector_unrotate( &
            local_ghost_vector        (:)             , &
            face_normal_vectors       (face_index, : ), &
            face_tangential1_vectors  (face_index, : ), &
            face_tangential2_vectors  (face_index, : )  &
        )
        ghost_primitive_variables(6:7) = primitive_variables_set(inner_cell_index, 6:7)
    end function make_ghost_primitive_variables_outflow

    pure function make_ghost_primitive_variables_slipwall( &
            primitive_variables_set   , &
            face_normal_vectors       , &
            face_tangential1_vectors  , &
            face_tangential2_vectors  , &
            face_index                , &
            inner_cell_index              ) result(ghost_primitive_variables)
        real   (real_kind), intent(in   ) :: primitive_variables_set    (:,:)
        real   (real_kind), intent(in   ) :: face_normal_vectors        (:,:)
        real   (real_kind), intent(in   ) :: face_tangential1_vectors   (:,:)
        real   (real_kind), intent(in   ) :: face_tangential2_vectors   (:,:)
        integer(int_kind ), intent(in   ) :: face_index
        integer(int_kind ), intent(in   ) :: inner_cell_index
        real   (real_kind) :: ghost_primitive_variables(7)
        real   (real_kind) :: local_inner_vector(3), local_ghost_vector(3)
        ghost_primitive_variables(1:2) = primitive_variables_set(inner_cell_index, 1:2)
        local_inner_vector(:) = vector_rotate(           &
            primitive_variables_set   (inner_cell_index, 3:5), &
            face_normal_vectors       (face_index      ,  : ), &
            face_tangential1_vectors  (face_index      ,  : ), &
            face_tangential2_vectors  (face_index      ,  : )  &
        )
        local_ghost_vector(1) = -1.d0 * local_inner_vector(1)
        local_ghost_vector(2) = local_inner_vector(2)
        local_ghost_vector(3) = local_inner_vector(3)
        ghost_primitive_variables(3:5) = vector_unrotate( &
            local_ghost_vector        (:)             , &
            face_normal_vectors       (face_index, : ), &
            face_tangential1_vectors  (face_index, : ), &
            face_tangential2_vectors  (face_index, : )  &
        )
        ghost_primitive_variables(6:7) = primitive_variables_set(inner_cell_index, 6:7)
    end function make_ghost_primitive_variables_slipwall

    pure function make_ghost_primitive_variables_symmetric( &
            primitive_variables_set   , &
            face_normal_vectors       , &
            face_tangential1_vectors  , &
            face_tangential2_vectors  , &
            face_index                , &
            inner_cell_index              ) result(ghost_primitive_variables)
        real   (real_kind), intent(in   ) :: primitive_variables_set    (:,:)
        real   (real_kind), intent(in   ) :: face_normal_vectors        (:,:)
        real   (real_kind), intent(in   ) :: face_tangential1_vectors   (:,:)
        real   (real_kind), intent(in   ) :: face_tangential2_vectors   (:,:)
        integer(int_kind ), intent(in   ) :: face_index
        integer(int_kind ), intent(in   ) :: inner_cell_index
        real   (real_kind) :: ghost_primitive_variables(7)
        real   (real_kind) :: local_inner_vector(3), local_ghost_vector(3)
        ghost_primitive_variables(1:2) = primitive_variables_set(inner_cell_index, 1:2)
        local_inner_vector(:) = vector_rotate(           &
            primitive_variables_set   (inner_cell_index, 3:5), &
            face_normal_vectors       (face_index      ,  : ), &
            face_tangential1_vectors  (face_index      ,  : ), &
            face_tangential2_vectors  (face_index      ,  : )  &
        )
        local_ghost_vector(1) = -1.d0 * local_inner_vector(1)
        local_ghost_vector(2) = local_inner_vector(2)
        local_ghost_vector(3) = local_inner_vector(3)
        ghost_primitive_variables(3:5) = vector_unrotate( &
            local_ghost_vector        (:)             , &
            face_normal_vectors       (face_index, : ), &
            face_tangential1_vectors  (face_index, : ), &
            face_tangential2_vectors  (face_index, : )  &
        )
        ghost_primitive_variables(6:7) = primitive_variables_set(inner_cell_index, 6:7)
    end function make_ghost_primitive_variables_symmetric

end module five_equation_model_boundary_condition_module