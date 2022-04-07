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

  static std::map<std::string, std::unique_ptr<Expression>> named_expressions;

  template<typename, bool, bool, int, int>
  friend class ExpressionOperation;

public:
  static void initialize();
  static void finalize();

  static Expression &make_named_expression(const std::string &name,
                                           const std::string &expression);

  Expression(const std::string &expression);
  ~Expression();

  void evaluate() const;
  void print() const;

  KOKKOS_INLINE_FUNCTION double
  operator[](int i)
  {
    return registers(0, i);
  }
};
