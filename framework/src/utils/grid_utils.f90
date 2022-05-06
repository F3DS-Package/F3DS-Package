module grid_utils_module
    use typedef_module

    implicit none

    private

    public :: convert_structure_index_to_unstructure_index

    contains

    pure function convert_structure_index_to_unstructure_index(i, j, k, imin, jmin, kmin, imax, jmax, kmax) result(n)
        integer(int_kind), intent(in) :: i, j, k
        integer(int_kind), intent(in) :: imin, jmin, kmin
        integer(int_kind), intent(in) :: imax, jmax, kmax
        integer(int_kind)             :: n

        n = (i - imin) + (j - jmin) * (imax - imin) + (k - kmin) * (imax - imin) * (jmax - jmin) + 1
    end function convert_structure_index_to_unstructure_index
end module grid_utils_module