#pragma once

#include <expression_function.h>

class ExpressionFunctionEvaluate:
  public ExpressionFunction<ExpressionFunctionEvaluate, 0>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return 0;
  }
};
