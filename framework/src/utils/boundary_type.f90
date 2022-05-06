module boundary_type_module
    enum, bind(c)
        enumerator :: boundary_type = 0
        enumerator :: outflow_boundary_type
        enumerator :: symmetric_boundary_type
        enumerator :: slipwall_boundary_type
        enumerator :: unknown_boundary_type
    end enum

    public :: string_to_boundary_type

contains

    pure function string_to_boundary_type(string) result(type)
        character(len=*), intent(in) :: string
        integer(kind(boundary_type)) :: type
        if (string == "outflow") then
            type = outflow_boundary_type
        else if (string == "symmetric") then
            type = symmetric_boundary_type
        else if (string == "slipwall") then
            type = slipwall_boundary_type
        else
            type = unknown_boundary_type
        end if
    end function string_to_boundary_type
end module boundary_type_module