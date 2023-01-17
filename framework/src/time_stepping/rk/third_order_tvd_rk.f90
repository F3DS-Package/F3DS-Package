module class_third_order_tvd_rk
    use typedef_module
    use abstract_time_stepping
    use abstract_configuration
    use abstract_eos

    implicit none

    private

    integer(int_kind ), parameter :: nmu_stage_ = 3
    real   (real_kind), parameter :: alpha1_(1:3) = [real(real_kind) :: 1.d0, 3.d0/4.d0, 1.d0/3.d0]
    real   (real_kind), parameter :: alpha2_(1:3) = [real(real_kind) :: 0.d0, 1.d0/4.d0, 2.d0/3.d0]
    real   (real_kind), parameter :: beta_  (1:3) = [real(real_kind) :: 1.d0, 1.d0/4.d0, 2.d0/3.d0]

    type, public, extends(time_stepping) :: third_order_tvd_rk
        private

        ! Stored init conservative-variables of each in cells
        ! Elm. 1) 1:number of conservative variables
        ! Elm. 2) 1:number of cells
        real(real_kind), allocatable :: init_conservative_variables_set(:,:)

        contains

        procedure, public, pass(self) :: initialize
        procedure, public, pass(self) :: compute_next_stage
        procedure, public, pass(self) :: prepare_time_stepping
        procedure, public, pass(self) :: get_number_of_stages
    end type third_order_tvd_rk

    contains

    subroutine initialize(self, config, num_cells, num_conservative_variables)
        class  (third_order_tvd_rk), intent(inout) :: self
        class  (configuration     ), intent(in   ) :: config
        integer(int_kind          ), intent(in   ) :: num_cells
        integer(int_kind          ), intent(in   ) :: num_conservative_variables

        allocate(self%init_conservative_variables_set(1:num_conservative_variables, 1:num_cells))
    end subroutine initialize

    subroutine compute_next_stage(   &
        self                       , &
        cell_index                 , &
        stage_num                  , &
        time_increment             , &
        conservative_variables     , &
        residuals                      )

        class  (third_order_tvd_rk), intent(inout) :: self
        integer(int_kind          ), intent(in   ) :: cell_index
        integer(int_kind          ), intent(in   ) :: stage_num
        real   (real_kind         ), intent(in   ) :: time_increment
        real   (real_kind         ), intent(inout) :: conservative_variables(:)
        real   (real_kind         ), intent(inout) :: residuals             (:)

        conservative_variables(:) = alpha1_(stage_num) * self%init_conservative_variables_set(:,cell_index) &
                                  + alpha2_(stage_num) * conservative_variables    (:)                      &
                                  + beta_  (stage_num) * time_increment * residuals(:)
        residuals             (:) = 0.d0
    end subroutine compute_next_stage

    subroutine prepare_time_stepping(   &
        self                     , &
        cell_index               , &
        conservative_variables   , &
        residuals                    )

        class   (third_order_tvd_rk), intent(inout) :: self
        integer(int_kind           ), intent(in   ) :: cell_index
        real   (real_kind          ), intent(inout) :: conservative_variables(:)
        real   (real_kind          ), intent(inout) :: residuals             (:)

        self%init_conservative_variables_set(:,cell_index) = conservative_variables(:)
    end subroutine prepare_time_stepping

    pure function get_number_of_stages(self) result(n)
        class  (third_order_tvd_rk), intent(in) :: self
        integer(int_kind          )             :: n
        n = nmu_stage_
    end function get_number_of_stages
end module class_third_order_tvd_rk