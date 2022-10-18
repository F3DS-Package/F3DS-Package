module viscous_five_equation_model_reconstructor_generator_module
    use abstract_configuration
    use stdio_module
    use abstract_reconstructor
    use class_minmod_muscl3
    use class_weno5_js
    use class_mp_weno5_js
    use class_rho_thinc
    implicit none

    private

    public :: viscous_five_equation_model_reconstructor_generator

    interface viscous_five_equation_model_reconstructor_generator
        module procedure with_config, with_name
    end interface viscous_five_equation_model_reconstructor_generator

    contains

    subroutine with_config(config, method)
        class(configuration), intent(inout) :: config
        class(reconstructor), pointer, intent(inout) :: method

        logical          :: found
        character(len=:), allocatable :: name

        call config%get_char("Reconstructor.Name", name, found)
        if(.not. found) call call_error("'Reconstructor.Name' is not found in configuration file you set.")

        call with_name(name, method)
    end subroutine with_config

    subroutine with_name(name, method)
        character(len=:), allocatable, intent(in) :: name
        class(reconstructor), pointer, intent(inout) :: method

        if(name == "Minmod MUSCL3") then
            allocate(minmod_muscl3 :: method)
        else if(name == "WENO5-JS") then
            allocate(weno5_js :: method)
        else if(name == "MP-WENO5-JS") then
            allocate(mp_weno5_js :: method)
        else if(name == "rho-THINC") then
            allocate(rho_thinc :: method)
        else
            call call_error("Unknow reconstructor name '"//name//"' is found.")
        end if
    end subroutine with_name
end module viscous_five_equation_model_reconstructor_generator_module