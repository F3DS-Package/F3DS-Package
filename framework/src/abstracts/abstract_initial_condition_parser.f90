module abstract_initial_condition_parser
    implicit none

    private

    type, public, abstract :: initial_condition_parser
        contains
        procedure(parse_interface), pass(self), deferred :: parse
        procedure(close_interface), pass(self), deferred :: close

        procedure(get_conservative_variables_set_interface), pass(self), deferred :: get_conservative_variables_set
    end type initial_condition_parser

    abstract interface
        subroutine parse_interface(self, filepath)
            import initial_condition_parser
            class    (initial_condition_parser), intent(inout) :: self
            character(len=*)                   , intent(in   ) :: filepath
        end subroutine parse_interface

        subroutine close_interface(self)
            import initial_condition_parser
            class(initial_condition_parser), intent(inout) :: self
        end subroutine close_interface

        subroutine get_conservative_variables_set_interface(self, conservative_variables_set)
            use typedef_module
            import initial_condition_parser
            class(initial_condition_parser), intent(inout) :: self
            real (real_kind               ), intent(inout) :: conservative_variables_set(:,:)
        end subroutine get_conservative_variables_set_interface
    end interface
end module abstract_initial_condition_parser