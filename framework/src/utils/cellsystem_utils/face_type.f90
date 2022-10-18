module face_type_module
    enum, bind(c)
        enumerator :: face_type = 0
        enumerator :: outflow_face_type
        enumerator :: symmetric_face_type
        enumerator :: wall_face_type
        enumerator :: empty_face_type
        enumerator :: unknown_face_type
    end enum

    public :: string_to_face_type

contains

    pure function string_to_face_type(string) result(type)
        character(len=*), intent(in) :: string
        integer(kind(face_type)) :: type
        if (string == "outflow") then
            type = outflow_face_type
        else if (string == "symmetric") then
            type = symmetric_face_type
        else if (string == "wall") then
            type = wall_face_type
        else if (string == "empty") then
            type = empty_face_type
        else
            type = unknown_face_type
        end if
    end function string_to_face_type
end module face_type_module