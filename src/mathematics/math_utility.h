#pragma once

template<typename T>
KOKKOS_INLINE_FUNCTION T
sign(T value)
{
  return value < 0? -1: 1;
}