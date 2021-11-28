#pragma once

#include <Kokkos_Core.hpp>

template<typename T>
class LazyReference {
    T *source;
    T local;

public:
    LazyReference() = default;

    KOKKOS_INLINE_FUNCTION
    LazyReference(T *source):
        source(source),
        local(*source)
    {}
    
    KOKKOS_INLINE_FUNCTION
    ~LazyReference() {
        if (*source != local) {
            *source = local;
        }
    }

    KOKKOS_INLINE_FUNCTION
    operator const T &() const {
      return local;
    }

    KOKKOS_INLINE_FUNCTION
    operator T &() {
        return local;
    }

    KOKKOS_INLINE_FUNCTION T &
    operator=(const T &value)
    {
      return local = value;
    }

    KOKKOS_INLINE_FUNCTION T &
    operator+=(const T &value)
    {
      return local += value;
    }

    KOKKOS_INLINE_FUNCTION T &
    operator-=(const T &value)
    {
      return local -= value;
    }

    KOKKOS_INLINE_FUNCTION T &
    operator*=(const T &value)
    {
      return local *= value;
    }

    KOKKOS_INLINE_FUNCTION T &
    operator/=(const T &value)
    {
      return local /= value;
    }
};