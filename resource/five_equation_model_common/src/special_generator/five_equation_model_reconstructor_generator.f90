module five_equation_model_reconstructor_generator_module
    use abstract_configuration
    use stdio_module
    use abstract_reconstructor
    use default_reconstructor_generator_module
    use class_rho_thinc
    implicit none

    private

    type, public, extends(reconstructor_generator) :: five_equation_model_reconstructor_generator
        contains
        procedure, pass(self) :: generate_from_configuration
        procedure, pass(self) :: generate_from_name
    end type five_equation_model_reconstructor_generator

    contains

    subroutine generate_from_configuration(self, a_reconstructor, a_config)
        class(five_equation_model_reconstructor_generator),          intent(inout) :: self
        class(reconstructor                              ), pointer, intent(inout) :: a_reconstructor
        class(configuration                              ),          intent(inout) :: a_config

        logical          :: found
        character(len=:), allocatable :: name

        call a_config%get_char("Reconstructor.Name", name, found)
        if(.not. found) call call_error("'Reconstructor.Name' is not found in configuration file you set.")

        call self%generate_from_name(a_reconstructor, name)
    end subroutine generate_from_configuration

    subroutine generate_from_name(self, a_reconstructor, name)
        class    (five_equation_model_reconstructor_generator),              intent(inout) :: self
        class    (reconstructor                              ), pointer    , intent(inout) :: a_reconstructor
        character(len=:                                      ), allocatable, intent(in   ) :: name

        type(default_reconstructor_generator) :: default_generator

        if(name == "rho-THINC") then
            allocate(rho_thinc :: a_reconstructor)
        else
            call default_generator%generate(a_reconstructor, name)
        end if
    end subroutine generate_from_name
end module five_equation_model_reconstructor_generator_module