module face_type_module
    use typedef_module

    implicit none

    enum, bind(c)
        enumerator :: boundary_face_type = 0
        enumerator :: outflow_face_type
        enumerator :: slip_and_symmetric_face_type
        enumerator :: nonslip_wall_face_type
        enumerator :: empty_face_type
        enumerator :: unknown_face_type
    end enum

    integer, public, parameter :: boundary_face_type_kind = kind(boundary_face_type)

    integer(int_kind), public :: number_of_boundary_face_type = 5

    public :: string_to_boundary_face_type

contains

    pure function string_to_boundary_face_type(string) result(type)
        character(len=*), intent(in) :: string
        integer(boundary_face_type_kind) :: type
        if (string == "outflow") then
            type = outflow_face_type
        else if (string == "symmetric") then
            type = slip_and_symmetric_face_type
        else if (string == "slip wall") then
            type = slip_and_symmetric_face_type
        else if (string == "non-slip wall") then
            type = nonslip_wall_face_type
        else if (string == "empty") then
            type = empty_face_type
        else
            type = unknown_face_type
        end if
    end function string_to_boundary_face_type
end module face_type_module