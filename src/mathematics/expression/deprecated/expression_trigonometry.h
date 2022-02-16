#pragma once

#include <expression_arithmetic.h>

#include <cmath>

template<typename T>
class ExpressionSin: public Expression<T>
{
  Expression<T> angle;

  KOKKOS_INLINE_FUNCTION void
  updateValueCache(const Solid *solid, int ip, const Grid *grid, int in) const override
  {
    return std::sin(angle.getValue(solid, ip, grid, in));
  }

  void
  print(std::ostream &stream) const override
  {
    stream << "sin(" << angle << ")";
  }

  constexpr int
  precedence() const override
  {
    return 2;
  }

public:
  ExpressionSum(const Expression<T> &angle):
    Expression(angle),
    angle(angle)
  {}
  
  KOKKOS_INLINE_FUNCTION Expression<T>
  timeDerivative() const override;
};

template<typename T>
class ExpressionCos: public Expression<T>
{
  Expression<T> angle;

  KOKKOS_INLINE_FUNCTION void
  updateValueCache(const Solid *solid, int ip, const Grid *grid, int in) const override
  {
    return std::cos(angle.getValue(solid, ip, grid, in));
  }

  void
  print(std::ostream &stream) const override
  {
    stream << "cos(" << angle << ")";
  }

  constexpr int
  precedence() const override
  {
    return 2;
  }

public:
  ExpressionCos(const Expression<T> &angle):
    Expression(angle),
    angle(angle)
  {}
  
  KOKKOS_INLINE_FUNCTION Expression<T>
  timeDerivative() const override;
};

template<typename T>
Expression<T>
ExpressionSin<T>::timeDerivative() const
{
  return ExpressionCos<T>(angle)*angle.timeDerivative();
}

template<typename T>
Expression<T>
ExpressionCos<T>::timeDerivative() const
{
  return -ExpressionSin<T>(angle)*angle.timeDerivative();
}
