#pragma once

#include <expression_operand.h>

#include <solid.h>

class ExpressionOperandTime:
  public ExpressionOperand<ExpressionOperandTime>
{
  double t;

public:
  void
  set_solid(const Solid &solid) override
  {
    t = solid.update->atime;
  }

  void
  set_grid(const Grid &grid) override
  {
    t = grid.update->atime;
  }

  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    return t;
  }
};
