#pragma once

template<typename T>
T
constexpr sign(T value)
{
  return value < 0? -1: 1;
}