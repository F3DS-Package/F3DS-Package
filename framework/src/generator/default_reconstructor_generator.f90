module default_reconstructor_generator_module
    use abstract_configuration
    use stdio_module
    use abstract_reconstructor
    use class_mp
    use class_minmod_muscl3
    use class_weno5_js
    use class_weno5_z
    implicit none

    private

    type, public, extends(reconstructor_generator) :: default_reconstructor_generator
        contains
        procedure, pass(self) :: generate_from_configuration
        procedure, pass(self) :: generate_from_name
    end type default_reconstructor_generator

    contains

    subroutine generate_from_configuration(self, a_reconstructor, a_config)
        class(default_reconstructor_generator),          intent(inout) :: self
        class(reconstructor                  ), pointer, intent(inout) :: a_reconstructor
        class(configuration                  ),          intent(inout) :: a_config

        logical          :: found
        character(len=:), allocatable :: name

        call a_config%get_char("Reconstructor.Name", name, found)
        if(.not. found) call call_error("'Reconstructor.Name' is not found in configuration file you set.")

        call self%generate_from_name(a_reconstructor, name)
    end subroutine generate_from_configuration

    subroutine generate_from_name(self, a_reconstructor, name)
        class    (default_reconstructor_generator),              intent(inout) :: self
        class    (reconstructor                  ), pointer    , intent(inout) :: a_reconstructor
        character(len=:                          ), allocatable, intent(in   ) :: name

        if(name == "Minmod MUSCL3") then
            allocate(minmod_muscl3 :: a_reconstructor)
        else if(name == "WENO5-JS") then
            allocate(weno5_js :: a_reconstructor)
        else if(name == "WENO5-Z") then
            allocate(weno5_z :: a_reconstructor)
        else if(name == "MP") then
            allocate(mp :: a_reconstructor)
        else
            call call_error("Unknow reconstructor name '"//name//"' is found.")
        end if
    end subroutine generate_from_name
end module default_reconstructor_generator_module