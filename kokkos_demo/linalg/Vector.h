#pragma once

#include <cmath>
#include <ostream>

using namespace std;

template<typename S>
using enable_if_scalar = enable_if_t<is_arithmetic<remove_reference_t<S>>::value>;

template<typename T,
    size_t N>
class Vector
{
    T elements[N];

public:
    // constructors
    KOKKOS_INLINE_FUNCTION
    Vector():
        elements{}
    {}

    template<typename ...U,
        typename = enable_if_t<sizeof...(U) == N>>
    KOKKOS_INLINE_FUNCTION
    Vector(U ...elements):
        elements{ static_cast<T>(elements)... }
    {}

    template<typename U,
      typename = enable_if_scalar<U>>
    KOKKOS_INLINE_FUNCTION
    Vector(U value)
    {
        for (int i = 0; i < N; i++)
        {
            (*this)[i] = static_cast<T>(value);
        }
    }

    template<typename U>
    KOKKOS_INLINE_FUNCTION
    explicit Vector(const Vector<U, N> &vector)
    {
        for (int i = 0; i < N; i++) (*this)[i] = static_cast<T>(vector[i]);
    }

    Vector(const Vector<T, N> &vector) = default;

    Vector(Vector<T, N> &&vector) = default;

    Vector<T, N> &
    operator=(const Vector<T, N> &vector) = default;

    Vector<T, N> &
    operator=(Vector<T, N> &&vector) = default;

    // accessors
    KOKKOS_INLINE_FUNCTION T
    operator[](size_t i) const
    {
        return elements[i];
    }

    KOKKOS_INLINE_FUNCTION T &
    operator[](size_t i)
    {
        return elements[i];
    }

    template<size_t M = N,
        typename = enable_if_t<(M > 0)>>
    KOKKOS_INLINE_FUNCTION T
    x() const
    { 
        return (*this)[0]; 
    }

    template<size_t M = N,
        typename = enable_if_t<(M > 0)>>
    KOKKOS_INLINE_FUNCTION T &
    x()
    { 
        return (*this)[0]; 
    }

    template<size_t M = N,
        typename = enable_if_t<(M > 1)>>
    KOKKOS_INLINE_FUNCTION T
    y() const
    { 
        return (*this)[1];
    }

    template<size_t M = N,
        typename = enable_if_t<(M > 1)>>
    KOKKOS_INLINE_FUNCTION T &
    y()
    { 
        return (*this)[1];
    }

    template<size_t M = N, 
        typename = enable_if_t<(M > 2)>>
    KOKKOS_INLINE_FUNCTION T 
    z() const
    {
        return (*this)[2];
    }

    template<size_t M = N, 
        typename = enable_if_t<(M > 2)>>
    KOKKOS_INLINE_FUNCTION T &
    z()
    {
        return (*this)[2];
    }

    template<size_t M = N,
        typename = enable_if_t<(M > 3)>>
    KOKKOS_INLINE_FUNCTION T
    w() const
    {
        return (*this)[3];
    }

    template<size_t M = N,
        typename = enable_if_t<(M > 3)>>
    KOKKOS_INLINE_FUNCTION T &
    w()
    {
        return (*this)[3];
    }

    // addition/subtraction
    template<typename U>
    KOKKOS_INLINE_FUNCTION Vector<T, N> &
    operator+=(const Vector<U, N> &vector)
    {
        for (int i = 0; i < N; i++)
            (*this)[i] += vector[i];
        return *this;
    }

    KOKKOS_INLINE_FUNCTION Vector<T, N> &
    operator-=(const Vector<T, N> &vector)
    {
        for (int i = 0; i < N; i++)
            (*this)[i] -= vector[i];
        return *this;
    }

    KOKKOS_INLINE_FUNCTION Vector<T, N>
    operator+(const Vector<T, N> &vector) const
    {
        return Vector<T, N>(*this) += vector;
    }

    KOKKOS_INLINE_FUNCTION Vector<T, N>
    operator-(const Vector<T, N> &vector) const
    {
        return Vector<T, N>(*this) -= vector;
    }

    // scalar multiplication/division
    template<typename S,
        typename = enable_if_scalar<S>>
    KOKKOS_INLINE_FUNCTION Vector<T, N> &
    operator*=(S &&scalar)
    {
        for (int i = 0; i < N; i++)
            (*this)[i] *= scalar;
        return *this;
    }

    template<typename S,
        typename = enable_if_scalar<S>>
    KOKKOS_INLINE_FUNCTION Vector<T, N> &
    operator/=(S &&scalar)
    {
        for (int i = 0; i < N; i++)
            (*this)[i] /= scalar;
        return *this;
    }

    template<typename S,
        typename = enable_if_scalar<S>>
    KOKKOS_INLINE_FUNCTION Vector<T, N>
    operator*(S &&scalar) const
    {
        return Vector<T, N>(*this) *= scalar;
    }

    template<typename S,
        typename = enable_if_scalar<S>>
    KOKKOS_INLINE_FUNCTION Vector<T, N>
    operator/(S &&scalar) const
    {
        return Vector<T, N>(*this) /= scalar;
    }

    KOKKOS_INLINE_FUNCTION Vector<T, N>
    operator-() const
    {
        return -1*Vector<T, N>(*this);
    }

    // dot product
    template<typename S = T,
        typename = enable_if_scalar<S>>
    KOKKOS_INLINE_FUNCTION T 
    operator*(const Vector<T, N> &vector) const
    {
        T result{};
        for (int i = 0; i < N; i++)
            result += (*this)[i]*vector[i];
        return result;
    }

    // norms and units
    template<typename S = T,
        typename = enable_if_scalar<S>>
    KOKKOS_INLINE_FUNCTION T
    norm2() const
    {
        return (*this)**this;
    }

    template<typename S = T,
        typename = enable_if_scalar<S>>
    KOKKOS_INLINE_FUNCTION T
    norm() const
    { 
        return sqrt(norm2());
    }

    template<typename S = T,
        typename = enable_if_scalar<S>>
    KOKKOS_INLINE_FUNCTION Vector<T, N>
    unit() const
    { 
        return (*this)/norm();
    }

    template<typename S = T,
        typename = enable_if_scalar<S>>
    KOKKOS_INLINE_FUNCTION Vector<T, N> &
    normalize()
    {
        return (*this) /= norm();
    }

    KOKKOS_INLINE_FUNCTION Vector<T, N>
    square()
    {
        Vector<T, N> vector(*this);

        for (int i = 0; i < N; i++)
        {
            vector[i] *= vector[i];
        }

        return vector;
    }
};

// scale
template<typename S,
    typename T,
    size_t N,
    typename = enable_if_scalar<S>>
KOKKOS_INLINE_FUNCTION Vector<T, N>
operator*(S &&scalar,
    const Vector<T, N> &vector)
{
    return vector*scalar;
}

// printing
template<typename T,
    size_t N>
KOKKOS_INLINE_FUNCTION ostream &
operator<<(ostream &os,
    const Vector<T, N> &vector)
{
    os << "(" << vector[0];

    for (int i = 1; i < N; i++)
    {
        os << ", " << vector[i];
    }

    return os << ")";
}

// typedefs
typedef Vector<double, 2> Vector2d;
typedef Vector<double, 3> Vector3d;
typedef Vector<double, 4> Vector4d;
typedef Vector<float,  2> Vector2f;
typedef Vector<float,  3> Vector3f;
typedef Vector<float,  4> Vector4f;
typedef Vector<int,    2> Vector2i;
typedef Vector<int,    3> Vector3i;
typedef Vector<int,    4> Vector4i;