#pragma once

#include <expression.h>

class Expression::Operation
{
  Expression *expression;
  Kokkos::View<double**> registers;

  void setExpression(Expression *expression)
  {
    this->expression = expression;
    this->registers = expression->registers;
  }

  virtual bool isOperand() const = 0;
  virtual bool isFunction() const = 0;
  virtual bool isParenthesis() const = 0;
  virtual int precedence() const = 0;
  virtual void apply() = 0;

  friend class Expression;
  template<typename, bool, bool, int, int>
  friend class ExpressionOperation;
};

template<typename DERIVED, bool IS_OPERAND, bool IS_FUNCTION, int PRECEDENCE, int ARITY>
class ExpressionOperation:
  public Expression::Operation
{
  int index;

  bool
  isOperand() const override
  {
    return IS_OPERAND;
  }

  bool
  isFunction() const override
  {
    return IS_FUNCTION;
  }

  bool
  isParenthesis() const override
  {
    return false;
  }

  int
  precedence() const override
  {
    return PRECEDENCE;
  }

public:
  void
  apply() override
  {
    index = expression->index;

    DERIVED derived(*static_cast<const DERIVED *>(this));
    int index = this->index;
    Kokkos::View<double**> registers = this->registers;

    Kokkos::parallel_for("", expression->registers.extent(1),
    KOKKOS_LAMBDA (int i)
    {
      registers(index - ARITY, i) = derived.evaluate(i);
    });

    expression->index -= ARITY - 1;
  }

protected:
  KOKKOS_INLINE_FUNCTION double
  get_value(int iregister, int i) const
  {
    return registers(index - ARITY + iregister, i);
  }
};
