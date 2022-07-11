#pragma once

#include <expression.h>

class Expression::Operation
{
  virtual bool isOperand() const = 0;
  virtual bool isFunction() const = 0;
  virtual int precedence() const = 0;
  virtual int arity() const = 0;
  virtual void set_solid(const Solid &solid) {};
  virtual void set_grid(const Grid &grid) {};
  virtual void apply(Expression *expression) = 0;

  friend class Expression;
  template<typename, bool, bool, int, int>
  friend class ExpressionOperation;
  friend class Input;
};

template<typename DERIVED, bool IS_OPERAND, bool IS_FUNCTION, int PRECEDENCE, int ARITY>
class ExpressionOperation:
  public Expression::Operation
{
  int index;
  Kokkos::View<float **> registers;

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

  int
  precedence() const override
  {
    return PRECEDENCE;
  }

  int
  arity() const override
  {
    return ARITY;
  }

public:
  void
  apply(Expression *expression) override
  {
    int index = this->index = expression->index;
    Kokkos::View<float **> registers = this->registers = expression->registers;
    DERIVED derived(*static_cast<const DERIVED *>(this));

    Kokkos::parallel_for(__PRETTY_FUNCTION__, expression->registers.extent(1),
    KOKKOS_LAMBDA (int i)
    {
      registers(index - ARITY, i) = derived.evaluate(i);
    });

    expression->index -= ARITY - 1;
  }

protected:
  KOKKOS_INLINE_FUNCTION float
  get_value(int iregister, int i) const
  {
    return registers(index - ARITY + iregister, i);
  }
};
