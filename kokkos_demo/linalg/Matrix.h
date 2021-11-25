#pragma once

#include <Vector.h>

#include <ostream>

using namespace std;

template<typename T,
    size_t M, 
    size_t N>
class Matrix:
    public Vector<Vector<T, N>, M>
{
public:
    // constructors
    Matrix() = default;

    template<typename ...U,
        typename = enable_if_t<sizeof...(U) == M>>
    KOKKOS_INLINE_FUNCTION
    Matrix(U ...rows):
        Vector<Vector<T, N>, M>{ (rows)... }
    {}

    // implement properly later
    template<typename U00,
        typename U01,
        typename U10,
        typename U11,
        size_t O = M,
        typename = enable_if_t<O == 2 && O == N>>
    KOKKOS_INLINE_FUNCTION
    Matrix(U00 m00, U01 m01, U10 m10, U11 m11):
        Vector<Vector<T, N>, M>{ Vector<T, N>{ m00, m01 }, Vector<T, M>{ m00, m01 } }
    {}

    Matrix(const Matrix<T, M, N> &matrix) = default;

    Matrix(Matrix<T, M, N> &&matrix) = default;

    Matrix<T, M, N> &
    operator=(const Matrix<T, M, N> &matrix) = default;

    Matrix<T, M, N> &
    operator=(Matrix<T, M, N> &&matrix) = default;

    // addition/subtraction REMOVE THIS
    KOKKOS_INLINE_FUNCTION Matrix<T, M, N> &
    operator+=(const Matrix<T, M, N> &vector)
    {
        for (int i = 0; i < N; i++)
            (*this)[i] += vector[i];
        return *this;
    }

    KOKKOS_INLINE_FUNCTION Matrix<T, M, N> &
    operator-=(const Matrix<T, M, N> &vector)
    {
        for (int i = 0; i < N; i++)
            (*this)[i] -= vector[i];
        return *this;
    }

    KOKKOS_INLINE_FUNCTION Matrix<T, M, N>
    operator+(const Matrix<T, M, N> &vector) const
    {
        return Matrix<T, M, N>(*this) += vector;
    }

    KOKKOS_INLINE_FUNCTION Matrix<T, M, N>
    operator-(const Matrix<T, M, N> &vector) const
    {
        return Matrix<T, M, N>(*this) -= vector;
    }

    // scalar multiplication/division REMOVE THIS
    template<typename S,
        typename = enable_if_scalar<S>>
    KOKKOS_INLINE_FUNCTION Matrix<T, M, N> &
    operator*=(S &&scalar)
    {
        for (int i = 0; i < N; i++)
            (*this)[i] *= scalar;
        return *this;
    }

    template<typename S,
        typename = enable_if_scalar<S>>
    KOKKOS_INLINE_FUNCTION Matrix<T, M, N> &
    operator/=(S &&scalar)
    {
        for (int i = 0; i < N; i++)
            (*this)[i] /= scalar;
        return *this;
    }

    template<typename S,
        typename = enable_if_scalar<S>>
    KOKKOS_INLINE_FUNCTION Matrix<T, M, N>
    operator*(S &&scalar) const
    {
        return Matrix<T, M, N>(*this) *= scalar;
    }

    template<typename S,
        typename = enable_if_scalar<S>>
    KOKKOS_INLINE_FUNCTION Matrix<T, M, N>
    operator/(S &&scalar) const
    {
        return Matrix<T, M, N>(*this) /= scalar;
    }

    // vector pre-multiplication
    KOKKOS_INLINE_FUNCTION Vector<T, M>
    operator*(const Vector<T, N> &vector) const
    {
        Vector<T, M> result;

        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                result[i] += (*this)[i][j]*vector[j];

        return result;
    }

    // matrix product
    template<size_t P>
    KOKKOS_INLINE_FUNCTION Matrix<T, M, P>
    operator*(const Matrix<T, N, P> &matrix) const
    {
        Matrix<T, M, P> result;

        for(int i = 0; i < M; i++)
            for(int j = 0; j < N; j++) 
                for(int k = 0; k < P; k++) 
                    result[i][k] += (*this)[i][j]*matrix[j][k];
                
        return result;
    }

    // transpose
    KOKKOS_INLINE_FUNCTION Matrix<T, N, M>
    transpose() const
    {
        Matrix<T, N, M> result;

        for(int i = 0; i < M; i++)
            for(int j = 0; j < N; j++) 
                result[j][i] = (*this)[i][j];

        return result;
    }

};

// vector post-multiplication
template<typename T,
    size_t M,
    size_t N>
KOKKOS_INLINE_FUNCTION Vector<T, M>
operator*(const Vector<T, N> &vector,
    const Matrix<T, M, N> &matrix)
{
    Vector<T, M> result;

    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
            result[i] += matrix[i][j]*vector[j];

    return result;
}

// scale REMOVE THIS ALSO
template<typename S,
    typename T,
    size_t M,
    size_t N,
    typename = enable_if_scalar<S>>
KOKKOS_INLINE_FUNCTION Matrix<T, M, N>
operator*(S &&scalar,
    const Matrix<T, M, N> &matrix)
{
    return matrix*scalar;
}

// typedefs
typedef Matrix<double, 2, 2> Matrix2d;
typedef Matrix<double, 3, 3> Matrix3d;
typedef Matrix<double, 4, 4> Matrix4d;
typedef Matrix<float,  2, 2> Matrix2f;
typedef Matrix<float,  3, 3> Matrix3f;
typedef Matrix<float,  4, 4> Matrix4f;