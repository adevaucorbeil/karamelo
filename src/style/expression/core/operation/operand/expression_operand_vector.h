#pragma once

#include <expression_operand.h>

#include <solid.h>

template <Kokkos::View<Vector3d *, MemorySpace> Solid::*SOLID_TO_VECTOR,
  Kokkos::View<Vector3d *, MemorySpace> Grid::*GRID_TO_VECTOR,
  int N>
class ExpressionOperandVector:
  public ExpressionOperand<ExpressionOperandVector<SOLID_TO_VECTOR, GRID_TO_VECTOR, N>>
{
  Kokkos::View<Vector3d*, MemorySpace> vector;

public:
  void
  set_solid(Solid &solid) override
  {
    vector = solid.*SOLID_TO_VECTOR;
  }

  void
  set_grid(Grid &grid) override
  {
    vector = grid.*GRID_TO_VECTOR;
  }

  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    return vector[i][N];
  }
};
