#pragma once

#include <style_factory.h>

#include <deque>
#include <memory>

#include <Kokkos_Core.hpp>

class Expression
{
  Kokkos::View<double**> registers;
  int index = 0;

  class Operation;
  std::deque<std::unique_ptr<Operation>> operations;
  static StyleFactory<Operation> operation_factory;

  template<typename, bool, bool, int, int>
  friend class ExpressionOperation;

public:
  static void initialize();

  Expression(const std::string &expression);
  ~Expression();

  void evaluate() const;
  void print() const;
};
