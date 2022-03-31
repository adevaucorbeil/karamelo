#ifdef EXPRESSION_CLASS

ExpressionStyle(inv_norm, ExpressionFunctionInverseNormal)

#else

#ifndef EXPRESSION_FUNCTION_INVERSE_NORMAL_H
#define EXPRESSION_FUNCTION_INVERSE_NORMAL_H

class ExpressionFunctionInverseNormal:
  public ExpressionFunction<ExpressionFunctionInverseNormal, 1>
{
  static KOKKOS_INLINE_FUNCTION double
  normalCDF(double x)
  {
    return Kokkos::Experimental::erfc(-x/Kokkos::Experimental::sqrt(2))/2;
  }

public:
  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    constexpr double eps = 1e-5;

    double value = get_value(0, i);

    if (value < eps || value > 1 - eps)
      return Kokkos::Experimental::nan("");

    double x_min = -1, x_max = 1;

    while (normalCDF(x_min) > value)
    {
      x_max = x_min;
      x_min *= 2;
    }

    while (normalCDF(x_max) < value)
    {
      x_min = x_max;
      x_max *= 2;
    }

    while (x_max - x_min > eps)
    {
      double x = (x_min + x_max)/2;

      (normalCDF(x) < value? x_min: x_max) = x;
    }

    return x_min;
  }
};

#endif
#endif