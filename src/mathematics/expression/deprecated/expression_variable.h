#pragma once

#include <expression.h>

#include <string>

template<typename T>
class ExpressionConstant: public Expression<T>
{
  std::string name;

  KOKKOS_INLINE_FUNCTION void
  updateValueCache(const Solid *solid, int ip, const Grid *grid, int in) override
  {}

  void
  print(std::ostream &stream) override
  {
    stream << name;
  }

  constexpr int
  precedence() override
  {
    return 0;
  }

public:
  ExpressionConstant(const T &constant, const std::string &name):
    name(name)
  {
    value_cache = constant;
  }
  
  KOKKOS_INLINE_FUNCTION Expression<T>
  timeDerivative() override
  {
    return ExpressionConstant(T{});
  }
};

class ExpressionTime: public Expression<double>
{
  KOKKOS_INLINE_FUNCTION void
  updateValueCache(const Solid *solid, int ip, const Grid *grid, int in) override
  {
    value_cache = solid->update->atime;
  }

  void
  print(std::ostream &stream) override
  {
    stream << "t";
  }

  constexpr int
  precedence() override
  {
    return 0;
  }

public:
  ExpressionTime()
  {}
  
  KOKKOS_INLINE_FUNCTION Expression<double>
  timeDerivative() override
  {
    return ExpressionConstant(1);
  }
};

template<typename T, vector<T> Solid::*VARIABLE>
class ExpressionParticle: public Expression<T>
{
  KOKKOS_INLINE_FUNCTION void
  updateValueCache(const Solid *solid, int ip, const Grid *grid, int in) override
  {
    value_cache = solid->*VARIABLE.at(ip);
  }

  void
  print(std::ostream &stream) override
  {
    // EB: revisit
    stream << "var";
  }

  constexpr int
  precedence() override
  {
    return 0;
  }

public:
  ExpressionParticle()
  {}
  
  KOKKOS_INLINE_FUNCTION Expression<T>
  timeDerivative() override
  {
    return ExpressionConstant(T{});
  }
};
