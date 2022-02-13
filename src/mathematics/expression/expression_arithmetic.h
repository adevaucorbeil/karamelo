#pragma once

#include <expression.h>

template<typename T>
class ExpressionSum: public Expression<T>
{
  Expression<T> first;
  Expression<T> second;

  KOKKOS_INLINE_FUNCTION void
  updateValueCache(const Solid *solid, int ip, const Grid *grid, int in) override
  {
    return first.getValue(solid, ip, grid, in) + second.getValue(solid, ip, grid, in);
  }

  void
  print(std::ostream &stream) override
  {
    stream << first << " + " << second;
  }

  constexpr int
  precedence() override
  {
    return 6;
  }

public:
  ExpressionSum(const Expression<T> &first, const Expression<T> &second):
    Expression(first, second),
    first(first),
    second(second)
  {}
  
  KOKKOS_INLINE_FUNCTION Expression<T>
  timeDerivative() override
  {
    return first.timeDerivative() + second.timeDerivative();
  }
};

template<typename T>
ExpressionSum<T>
operator+(const Expression<T> &first, const Expression<T> &second)
{
  return ExpressionSum<T>(first, second);
}

template<typename T>
class ExpressionDifference: public Expression<T>
{
  Expression<T> first;
  Expression<T> second;

  KOKKOS_INLINE_FUNCTION void
  updateValueCache(const Solid *solid, int ip, const Grid *grid, int in) override
  {
    return first.getValue(solid, ip, grid, in) - second.getValue(solid, ip, grid, in);
  }

  void
  print(std::ostream &stream) override
  {
    stream << first << " - " << second;
  }

  constexpr int
  precedence() override
  {
    return 6;
  }

public:
  ExpressionDifference(const Expression<T> &first, const Expression<T> &second):
    Expression(first, second),
    first(first),
    second(second)
  {}
  
  KOKKOS_INLINE_FUNCTION Expression<T>
  timeDerivative() override
  {
    return first.timeDerivative() - second.timeDerivative();
  }
};

template<typename T>
ExpressionDifference<T>
operator-(const Expression<T> &first, const Expression<T> &second)
{
  return ExpressionDifference<T>(first, second);
}

template<typename T>
class ExpressionProduct: public Expression<T>
{
  Expression<T> first;
  Expression<T> second;

  KOKKOS_INLINE_FUNCTION void
  updateValueCache(const Solid *solid, int ip, const Grid *grid, int in) override
  {
    return first.getValue(solid, ip, grid, in)*second.getValue(solid, ip, grid, in);
  }

  void
  print(std::ostream &stream) override
  {
    stream << first << "/" << second;
  }

  constexpr int
  precedence() override
  {
    return 5;
  }

public:
  ExpressionProduct(const Expression<T> &first, const Expression<T> &second):
    Expression(first, second),
    first(first),
    second(second)
  {}
  
  KOKKOS_INLINE_FUNCTION Expression<T>
  timeDerivative() override
  {
    return first*second.timeDerivative() + first.timeDerivative()*second;
  }
};

template<typename T>
ExpressionProduct<T>
operator*(const Expression<T> &first, const Expression<T> &second)
{
  return ExpressionProduct<T>(first, second);
}

template<typename T>
class ExpressionQuotient: public Expression<T>
{
  Expression<T> numerator;
  Expression<T> denominator;

  KOKKOS_INLINE_FUNCTION void
  updateValueCache(const Solid *solid, int ip, const Grid *grid, int in) override
  {
    return numerator.getValue(solid, ip, grid, in)/denominator.getValue(solid, ip, grid, in);
  }

  void
  print(std::ostream &stream) override
  {
    stream << first << "/" << second;
  }

  constexpr int
  precedence() override
  {
    return 5;
  }

public:
  ExpressionQuotient(const Expression<T> &numerator, const Expression<T> &denominator):
    Expression(numerator, denominator),
    numerator(numerator),
    denominator(denominator)
  {}
  
  KOKKOS_INLINE_FUNCTION Expression<T>
  timeDerivative() override
  {
    return (first.timeDerivative()*second - first*second.timeDerivative())/second/second;
  }
};

template<typename T>
ExpressionQuotient<T>
operator/(const Expression<T> &numerator, const Expression<T> &denominator)
{
  return ExpressionQuotient(numerator, denominator);
}
