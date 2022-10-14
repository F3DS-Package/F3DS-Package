module matrix_module
    use typedef_module

    implicit none

    private

    public :: matrix_multiply

    contains

    pure function matrix_multiply(matrix, vector) result(v)
        real(real_kind), intent(in) :: matrix(:,:) ! 1:n_row, 1:n_col
        real(real_kind), intent(in) :: vector(:)
        real(real_kind)             :: v     (size(vector))

        v(1) = matrix(1,1) * vector(1) + matrix(2,1) * vector(2) + matrix(3,1) * vector(3)
        v(2) = matrix(1,2) * vector(1) + matrix(2,2) * vector(2) + matrix(3,2) * vector(3)
        v(3) = matrix(1,3) * vector(1) + matrix(2,3) * vector(2) + matrix(3,3) * vector(3)
    end function matrix_multiply
end module matrix_module