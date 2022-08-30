module five_equation_model_variables_module
    use typedef_module
    use stdio_module
    use abstract_eos
    use math_constant_module
    use vector_module

    implicit none

    private

    ! surface tension coef. = 72e-3
    real(real_kind), public, parameter :: weber_number = 1.446d0

    integer(int_kind), public, parameter :: num_conservative_variables = 7
    integer(int_kind), public, parameter :: num_primitive_variables    = 8
    integer(int_kind), public, parameter :: num_surface_tension_variables = 4

    ! Elm. 1) following variables are saved
    ! conservative_variables_set(1  , :)   = Z1*rho1   : Z1*density of fluid1
    ! conservative_variables_set(2  , :)   = Z2*rho2   : Z2*density of fluid2
    ! conservative_variables_set(3:5, :)   = rho*v     : momentum vector
    ! conservative_variables_set(6  , :)   = e         : energy density
    ! conservative_variables_set(7  , :)   = Z1        : volume fraction of fluid1
    ! Elm. 2) 1 : {@code num_cells}, cell index
    real(real_kind), public, allocatable :: conservative_variables_set(:,:)

    ! Elm. 1) following variables are saved
    ! primitive_variables_set(1  , :)   : density of fluid1
    ! primitive_variables_set(2  , :)   : density of fluid2
    ! primitive_variables_set(3:5, :)   : velocity vector (u,v,w)
    ! primitive_variables_set(6  , :)   : pressre
    ! primitive_variables_set(7  , :)   : volume fraction of fluid1
    ! primitive_variables_set(8  , :)   : pressure jump σκ (based on the surface tension coeficient σ and interface curvature κ)
    ! Elm. 2) 1 : {@code num_cells}, cell index
    real(real_kind), public, allocatable :: primitive_variables_set(:,:)

    ! Elm. 1) following variables are saved
    ! conservative_variables_set(1  , :)   = Z1*rho1   : Z1*density of fluid1
    ! conservative_variables_set(2  , :)   = Z2*rho2   : Z2*density of fluid2
    ! conservative_variables_set(3:5, :)   = rho*v     : momentum vector
    ! conservative_variables_set(6  , :)   = e         : energy density
    ! conservative_variables_set(7  , :)   = Z1        : volume fraction of fluid1
    ! Elm. 2) 1 : {@code num_cells}, cell index
    real(real_kind), public, allocatable :: residual_set(:,:)

    ! Elm. 1) following variables are saved
    ! surface_tension_variables_set(1:3, :) = gradient volume fraction (to be normalized)
    ! surface_tension_variables_set(4  , :) = curvature
    ! Elm. 2) 1 : {@code num_cells}, cell index
    real(real_kind), public, allocatable :: surface_tension_variables_set(:,:)

    public :: conservative_to_primitive
    public :: primitive_to_conservative
    public :: spectral_radius
    public :: rotate_primitive
    public :: unrotate_primitive
    public :: normarize_gradient_volume_fraction
    public :: rotate_gradient_value
    public :: unrotate_gradient_value
    public :: compute_pressure_jump

    real(real_kind), parameter :: interface_threshold = 1e-2

    contains

    pure function primitive_to_conservative(an_eos, primitive_variables, num_conservatives) result(conservative_variables)
        class  (eos         ), intent(in)  :: an_eos
        real   (real_kind   ), intent(in)  :: primitive_variables   (:)
        integer(int_kind    ), intent(in)  :: num_conservatives
        real   (real_kind   )              :: conservative_variables(num_conservatives)

        real(8) :: rho

        associate(                                 &
                rho1    => primitive_variables(1), &
                rho2    => primitive_variables(2), &
                u       => primitive_variables(3), &
                v       => primitive_variables(4), &
                w       => primitive_variables(5), &
                p       => primitive_variables(6), &
                z1      => primitive_variables(7)  &
            )
            rho = rho1 * z1 + rho2 * (1.d0 - z1)
            conservative_variables(1) = rho1 * z1
            conservative_variables(2) = rho2 * (1.d0 - z1)
            conservative_variables(3) = u  * rho
            conservative_variables(4) = v  * rho
            conservative_variables(5) = w  * rho
            conservative_variables(6) = an_eos%compute_internal_energy_density(p, rho, z1) + 0.5d0 * (u**2.d0 + v**2.d0 + w**2.d0) * rho
            conservative_variables(7) = z1
        end associate
    end function primitive_to_conservative

    pure function conservative_to_primitive(an_eos, conservative_variables, num_primitives) result(primitive_variables)
        class  (eos      ), intent(in)  :: an_eos
        real   (real_kind), intent(in)  :: conservative_variables(:)
        integer(int_kind ), intent(in)  :: num_primitives
        real   (real_kind)              :: primitive_variables   (num_primitives)

        real(real_kind) :: rho, ie

        associate(                                    &
                rho1_z1 => conservative_variables(1), &
                rho2_z2 => conservative_variables(2), &
                rho_u   => conservative_variables(3), &
                rho_v   => conservative_variables(4), &
                rho_w   => conservative_variables(5), &
                e       => conservative_variables(6), &
                z1      => conservative_variables(7)  &
            )
            rho = rho1_z1 + rho2_z2
            if(z1 < machine_epsilon)then
                primitive_variables(1) = 0.d0
                primitive_variables(2) = rho2_z2
                primitive_variables(7) = 0.d0
            else if(z1 > 1.d0 - machine_epsilon)then
                primitive_variables(1) = rho1_z1
                primitive_variables(2) = 0.d0
                primitive_variables(7) = 1.d0
            else
                primitive_variables(1) = rho1_z1 / z1
                primitive_variables(2) = rho2_z2 / (1.d0 - z1)
                primitive_variables(7) = z1
            end if

            if(rho < machine_epsilon)then
                primitive_variables(3:6) = 0.d0
            else
                ie  = e / rho - 0.5d0 * (rho_u**2.d0 + rho_v**2.d0 + rho_w**2.d0) / rho**2.d0
                primitive_variables(3) = rho_u / rho
                primitive_variables(4) = rho_v / rho
                primitive_variables(5) = rho_w / rho
                primitive_variables(6) = an_eos%compute_pressure(ie, rho, z1)
            endif
        end associate
    end function conservative_to_primitive

    pure function spectral_radius(an_eos, primitive_variables) result(r)
        class(eos      ), intent(in) :: an_eos
        real (real_kind), intent(in) :: primitive_variables(:)
        real (real_kind) :: r

        real (real_kind) :: c, rho

        associate(                            &
            rho1 => primitive_variables(1)  , &
            rho2 => primitive_variables(2)  , &
            vel  => primitive_variables(3:5), &
            p    => primitive_variables(6)  , &
            z1   => primitive_variables(7)    &
        )
            rho = z1 * rho1 + (1.d0 - z1) * rho2
            c   = an_eos%compute_soundspeed(p, rho, z1)
            r = c + vector_magnitude(vel)
        end associate
    end function spectral_radius

    pure function rotate_primitive( &
        primitives                , &
        face_normal_vector        , &
        face_tangential1_vector   , &
        face_tangential2_vector   , &
        num_primitives                ) result(p)

        real   (real_kind     ), intent(in)  :: primitives              (:)
        real   (real_kind     ), intent(in)  :: face_normal_vector      (3)
        real   (real_kind     ), intent(in)  :: face_tangential1_vector (3)
        real   (real_kind     ), intent(in)  :: face_tangential2_vector (3)
        integer(int_kind      ), intent(in)  :: num_primitives
        real   (real_kind     )              :: p(num_primitives)

        p(1  ) = primitives(1)
        p(2  ) = primitives(2)
        p(3:5) = vector_rotate(primitives(3:5), face_normal_vector, face_tangential1_vector, face_tangential2_vector)
        p(6  ) = primitives(6)
        p(7  ) = primitives(7)
        p(8  ) = primitives(8)
    end function rotate_primitive

    pure function unrotate_primitive( &
        primitives                  , &
        face_normal_vector          , &
        face_tangential1_vector     , &
        face_tangential2_vector     , &
        num_primitives                  ) result(p)

        real   (real_kind     ), intent(in)  :: primitives              (:)
        real   (real_kind     ), intent(in)  :: face_normal_vector      (3)
        real   (real_kind     ), intent(in)  :: face_tangential1_vector (3)
        real   (real_kind     ), intent(in)  :: face_tangential2_vector (3)
        integer(int_kind      ), intent(in)  :: num_primitives
        real   (real_kind     )              :: p(num_primitives)

        p(1  ) = primitives(1)
        p(2  ) = primitives(2)
        p(3:5) = vector_unrotate(primitives(3:5), face_normal_vector, face_tangential1_vector, face_tangential2_vector)
        p(6  ) = primitives(6)
        p(7  ) = primitives(7)
        p(8  ) = primitives(8)
    end function unrotate_primitive

    pure function normarize_gradient_volume_fraction(surface_tension_variables, primitives, num_variables) result(dst_surface_tension_variables)
        real   (real_kind ), intent(in) :: surface_tension_variables(:), primitives(:)
        integer(int_kind  ), intent(in) :: num_variables
        real   (real_kind )             :: dst_surface_tension_variables(num_variables)

        associate(                                                   &
            mag => vector_magnitude(surface_tension_variables(1:3)), &
            z   => primitives(7)                                     &
        )
            if((machine_epsilon < mag))then
                ! {@code normarize_gradient_volume_fraction} must always be oriented toward the heaviest fluid.
                dst_surface_tension_variables(1:3) = surface_tension_variables(1:3) / mag
            else
                dst_surface_tension_variables(1:3) = 0.d0
            endif
            dst_surface_tension_variables(4) = surface_tension_variables(4)
        end associate
    end function normarize_gradient_volume_fraction

    pure function rotate_gradient_value( &
        gradient_value                 , &
        face_normal_vector             , &
        face_tangential1_vector        , &
        face_tangential2_vector        , &
        num_values                         ) result(p)

        real   (real_kind     ), intent(in)  :: gradient_value          (:)
        real   (real_kind     ), intent(in)  :: face_normal_vector      (3)
        real   (real_kind     ), intent(in)  :: face_tangential1_vector (3)
        real   (real_kind     ), intent(in)  :: face_tangential2_vector (3)
        integer(int_kind      ), intent(in)  :: num_values
        real   (real_kind     )              :: p(num_values)

        p(1:3) = vector_rotate(gradient_value(1:3), face_normal_vector, face_tangential1_vector, face_tangential2_vector)
    end function rotate_gradient_value

    pure function unrotate_gradient_value( &
        gradient_value              , &
        face_normal_vector          , &
        face_tangential1_vector     , &
        face_tangential2_vector     , &
        num_values                  ) result(p)

        real   (real_kind     ), intent(in)  :: gradient_value          (:)
        real   (real_kind     ), intent(in)  :: face_normal_vector      (3)
        real   (real_kind     ), intent(in)  :: face_tangential1_vector (3)
        real   (real_kind     ), intent(in)  :: face_tangential2_vector (3)
        integer(int_kind      ), intent(in)  :: num_values
        real   (real_kind     )              :: p(num_values)

        p(1:3) = vector_unrotate(gradient_value(1:3), face_normal_vector, face_tangential1_vector, face_tangential2_vector)
    end function unrotate_gradient_value

    pure function compute_pressure_jump(primitives, surface_tension_variables, num_variables) result(dst_primitives)
        real   (real_kind ), intent(in) :: surface_tension_variables(:), primitives(:)
        integer(int_kind  ), intent(in) :: num_variables
        real   (real_kind )             :: dst_primitives(num_variables)

        dst_primitives(1:7) = primitives(1:7)
        associate(                                         &
            z     => primitives(7)                       , &
            kappa => -1.d0 * surface_tension_variables(4)  &
        )
            if(z < interface_threshold)then
                dst_primitives(8) = 0.d0
            else if(z > 1.d0 - interface_threshold)then
                dst_primitives(8) = 0.d0
            else
                dst_primitives(8) = (1.d0 / weber_number) * kappa
            endif
        end associate
    end function compute_pressure_jump
end module five_equation_model_variables_module