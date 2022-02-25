#include <grid.h>
#include <update.h>

KOKKOS_INLINE_FUNCTION void
Method::reset_mass_nodes(Grid &grid, int in) const
{
  if (!is_TL || !update->atimestep)
    grid.mass[in] = 0;
}

KOKKOS_INLINE_FUNCTION void
Method::reset_velocity_nodes(Grid &grid, int in) const
{
  grid.v[in] = Vector3d();
  grid.mb[in] = Vector3d();
  if (temp)
    grid.T[in] = 0;
}

KOKKOS_INLINE_FUNCTION void
Method::reset_force_nodes(Grid &grid, int in) const
{
  grid.f[in] = Vector3d();
  grid.mb[in] = Vector3d();
  if (temp)
  {
    grid.Qint[in] = 0;
    grid.Qext[in] = 0;
  }
}