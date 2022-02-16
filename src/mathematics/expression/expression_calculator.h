#pragma once

#include <Kokkos_Macros.hpp>

class ExpressionCalculator
{
  enum class CalculatorType: char
  {
    PRODUCT,
    QUOTIENT,
    SUM,
    DIFFERENCE
  };

  static KOKKOS_INLINE_FUNCTION double
  calculate(CalculatorType calculator_type)
  {
    switch (calculator_type)
    {
    case CalculatorType::PRODUCT:

    case CalculatorType::QUOTIENT:

    case CalculatorType::SUM:

    case CalculatorType::DIFFERENCE
 
    }
  }
};
