module five_equation_model_variables_module
    use typedef_module
    use stdio_module

    implicit none

    private

    ! Elm. 1) following variables are saved
    ! conservative_variables_set(:, 1  )   = Z1*rho1   : Z1*density of fluid1
    ! conservative_variables_set(:, 2  )   = Z2*rho2   : Z2*density of fluid2
    ! conservative_variables_set(:, 3:5) = rho*v     : momentum vector
    ! conservative_variables_set(:, 6  )   = e         : energy density
    ! conservative_variables_set(:, 7  )   = Z1        : volume fraction of fluid1
    ! Elm. 2) 1 : {@code num_cells}, cell index
    real(real_kind), public, allocatable :: conservative_variables_set(:,:)

    ! Elm. 1) following variables are saved
    ! primitive_variables_set(:, 1  )   : volume fraction * density of fluid1
    ! primitive_variables_set(:, 2  )   : volume fraction * density of fluid2
    ! primitive_variables_set(:, 3:5)   : velocity vector (u,v,w)
    ! primitive_variables_set(:, 6  )   : specific internal energy (not energy density)
    ! primitive_variables_set(:, 7  )   : volume fraction of fluid1
    ! Elm. 2) 1 : {@code num_cells}, cell index
    real(real_kind), public, allocatable :: primitive_variables_set(:,:)

    public :: conservative_to_primitive
    public :: primitive_to_conservative
    public :: initialise_variables
    public :: finalize_variables

    contains

    pure function primitive_to_conservative(primitive_variables_set) result(conservative_variables_set)
        real(real_kind), intent(in)  :: primitive_variables_set   (:)
        real(real_kind), allocatable :: conservative_variables_set(:)

        real(8) :: rho

        allocate(conservative_variables_set(7))

        associate(                              &
                rho1_z1 => primitive_variables_set(1), &
                rho2_z2 => primitive_variables_set(2), &
                u       => primitive_variables_set(3), &
                v       => primitive_variables_set(4), &
                w       => primitive_variables_set(5), &
                ie      => primitive_variables_set(6), &
                z1      => primitive_variables_set(7)  &
            )
            conservative_variables_set(1) = rho1_z1
            conservative_variables_set(2) = rho2_z2
            rho = conservative_variables_set(1) + conservative_variables_set(2)
            conservative_variables_set(3) = u  * rho
            conservative_variables_set(4) = v  * rho
            conservative_variables_set(5) = w  * rho
            conservative_variables_set(6) = ie * rho + 0.5d0 * (u**2.d0 + v**2.d0 + w**2.d0) / rho
            conservative_variables_set(7) = z1
        end associate
    end function primitive_to_conservative

    pure function conservative_to_primitive(conservative_variables_set) result(primitive_variables_set)
        real(real_kind), intent(in)  :: conservative_variables_set(:)
        real(real_kind), allocatable :: primitive_variables_set   (:)

        real(real_kind) :: rho, u, v, w

        allocate(primitive_variables_set(7))

        associate(                                    &
                rho1_z1 => conservative_variables_set(1), &
                rho2_z2 => conservative_variables_set(2), &
                rho_u   => conservative_variables_set(3), &
                rho_v   => conservative_variables_set(4), &
                rho_w   => conservative_variables_set(5), &
                e       => conservative_variables_set(6), &
                z1      => conservative_variables_set(7)  &
            )
            rho = rho1_z1 + rho2_z2
            primitive_variables_set(1) = rho1_z1
            primitive_variables_set(2) = rho2_z2
            if (rho == 0.d0)then
                primitive_variables_set(3:6) = 0.d0
            else
                primitive_variables_set(3) = rho_u / rho
                primitive_variables_set(4) = rho_v / rho
                primitive_variables_set(5) = rho_w / rho
                primitive_variables_set(6) = (e - 0.5d0 * (rho_u**2.d0 + rho_v**2.d0 + rho_w**2.d0) / rho) / rho
            end if
            primitive_variables_set(7) = z1
        end associate
    end function conservative_to_primitive

    subroutine initialise_variables(num_cells)
        integer(int_kind), intent(in) :: num_cells

        if(allocated(conservative_variables_set))then
            call call_error("Array conservative_variables_set is already allocated. But you call the initialiser for variables module.")
        end if
        allocate(conservative_variables_set(num_cells, 7))

        if(allocated(primitive_variables_set))then
            call call_error("Array primitive_variables_set is already allocated. But you call the initialiser for variables module.")
        end if
        allocate(primitive_variables_set(num_cells, 7))
    end subroutine initialise_variables

    subroutine finalize_variables()
        if(.not. allocated(conservative_variables_set))then
            call call_error("Array conservative_variables_set is not allocated. But you call the finalizer for variables module.")
        end if
        deallocate(conservative_variables_set)

        if(.not. allocated(primitive_variables_set))then
            call call_error("Array primitive_variables_set is not allocated. But you call the finalizer for variables module.")
        end if
        deallocate(primitive_variables_set)
    end subroutine finalize_variables
end module five_equation_model_variables_module