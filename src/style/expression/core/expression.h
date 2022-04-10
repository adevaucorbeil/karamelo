#pragma once

#include <style_factory.h>

#include <deque>
#include <memory>
#include <unordered_set>

#include <Kokkos_Core.hpp>

class Solid;
class Grid;

class Expression
{
  Kokkos::View<double**> registers;
  int index = 0;

  class Operation;
  std::deque<std::unique_ptr<Operation>> operations;
  std::unordered_set<Expression *> expression_dependencies;

  template<typename, bool, bool, int, int>
  friend class ExpressionOperation;
  friend class Input;

public:
  ~Expression();

  void evaluate();
  void evaluate(Solid &solid);
  void evaluate(Grid &grid);

  KOKKOS_INLINE_FUNCTION double
  operator[](int i) const
  {
    return registers(0, registers.extent(1) > 1? i: 0);
  }
};
