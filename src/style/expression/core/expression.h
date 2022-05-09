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
public:
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

  double
  getConstant()
  {
    evaluate();
    auto element_device = Kokkos::subview(registers, 0, 0);
    double element_host;
    Kokkos::deep_copy(element_host, element_device);
    return element_host;
  }
};
