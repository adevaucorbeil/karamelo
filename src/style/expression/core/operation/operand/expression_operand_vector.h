#pragma once

#include <expression_operand.h>

#include <solid.h>

template <Kokkos::View<Vector3d *> Solid::*SOLID_TO_VECTOR,
  Kokkos::View<Vector3d *> Grid::*GRID_TO_VECTOR,
  int N>
class ExpressionOperandVector:
  public ExpressionOperand<ExpressionOperandVector<SOLID_TO_VECTOR, GRID_TO_VECTOR, N>>
{
  Kokkos::View<Vector3d*> vector;

public:
  void
  set_solid(const Solid &solid) override
  {
    vector = solid.*SOLID_TO_VECTOR;
  }

  void
  set_grid(const Grid &grid) override
  {
    vector = grid.*GRID_TO_VECTOR;
  }

  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    return vector[i][N];
  }
};
