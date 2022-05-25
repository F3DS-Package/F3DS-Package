module vector_module
    use typedef_module

    implicit none

    private

    public :: vector_angle
    public :: vector_magnitude
    public :: vector_distance
    public :: multiply_vector
    public :: rotate_vector
    public :: reverse_vector

    contains

    pure function vector_angle(vector1, vector2) result(angle)
        real(real_kind), intent(in) :: vector1 (3)
        real(real_kind), intent(in) :: vector2 (3)
        real(real_kind)             :: angle
        angle = acos( &
            (vector1(1) * vector2(1) + vector1(2) * vector2(2) + vector1(3) * vector2(3)) &
            / (sqrt(vector1(1)**2.d0 + vector1(2)**2.d0 + vector1(3)**2.d0) + sqrt(vector2(1)**2.d0 + vector2(2)**2.d0 + vector2(3)**2.d0)) &
        )
    end function vector_angle

    pure function vector_magnitude(vector) result(magnitude)
        real(real_kind), intent(in) :: vector(3)
        real(real_kind)             :: magnitude
        magnitude = sqrt(vector(1)**2.d0 + vector(2)**2.d0 + vector(3)**2.d0)
    end function vector_magnitude

    pure function vector_distance(vector1, vector2) result(distance)
        real(real_kind), intent(in) :: vector1 (3)
        real(real_kind), intent(in) :: vector2 (3)
        real(real_kind)             :: v(3)
        real(real_kind)             :: distance
        v(1) = vector1(1) - vector2(1)
        v(2) = vector1(2) - vector2(2)
        v(3) = vector1(3) - vector2(3)
        distance = vector_magnitude(v)
    end function vector_distance

    pure function multiply_vector(vector1, vector2) result(scalar)
        real(real_kind), intent(in) :: vector1 (3)
        real(real_kind), intent(in) :: vector2 (3)
        real(real_kind)             :: scalar
        scalar = vector1(1) * vector2(1) + vector1(2) * vector2(2) + vector1(3) * vector2(3)
    end function

    pure function rotate_vector(vector, normal_vector, tangential1_vector, tangential2_vector) result(rotated_vector)
        real(real_kind), intent(in) :: vector             (3)
        real(real_kind), intent(in) :: normal_vector      (3)
        real(real_kind), intent(in) :: tangential1_vector (3)
        real(real_kind), intent(in) :: tangential2_vector (3)

        real(real_kind) :: rotated_vector (3)

        rotated_vector(1) = multiply_vector(vector, normal_vector     )
        rotated_vector(2) = multiply_vector(vector, tangential1_vector)
        rotated_vector(3) = multiply_vector(vector, tangential2_vector)
    end function

    pure function reverse_vector(vector, normal_vector, tangential1_vector, tangential2_vector) result(reversed_vector)
        real(real_kind), intent(in) :: vector             (3)
        real(real_kind), intent(in) :: normal_vector      (3)
        real(real_kind), intent(in) :: tangential1_vector (3)
        real(real_kind), intent(in) :: tangential2_vector (3)

        real(real_kind) :: reversed_vector (3)

        real(real_kind) :: xs(3), ys(3), zs(3)

        xs(1) = normal_vector     (1)
        xs(2) = tangential1_vector(1)
        xs(3) = tangential2_vector(1)

        ys(1) = normal_vector     (2)
        ys(2) = tangential1_vector(2)
        ys(3) = tangential2_vector(2)

        zs(1) = normal_vector     (3)
        zs(2) = tangential1_vector(3)
        zs(3) = tangential2_vector(3)

        reversed_vector(1) = multiply_vector(vector, xs)
        reversed_vector(2) = multiply_vector(vector, ys)
        reversed_vector(3) = multiply_vector(vector, zs)
    end function
end module