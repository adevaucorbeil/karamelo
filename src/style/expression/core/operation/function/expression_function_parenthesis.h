#pragma once

#include <expression_function.h>

class ExpressionFunctionParenthesis:
  public ExpressionFunction<ExpressionFunctionParenthesis, 1>
{
public:
  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    return 0;
  }

  bool
  isParenthesis() const override
  {
    return true;
  }
};
