#include <expression_base.h>

#include <solid.h>
#include <update.h>

bool
ExpressionBase::cacheNeedsUpdating(const Solid *solid, int ip, const Grid *grid, int in)
{
  if (constant)
    return false;

  bigint t = (solid? solid->update: grid->update)->ntimestep;

  bool cache_needs_updating = false;

  if (cache_needs_updating |= solid_cache != solid)
    solid_cache = solid;
  if (cache_needs_updating |= grid_cache != grid)
    grid_cache = grid;
  if (cache_needs_updating |= t_cache != t)
    t_cache = t;

  if (constant_if_constant_time && !cache_needs_updating)
    return false;

  if (cache_needs_updating |= ip_cache != ip)
    ip_cache = ip;
  if (cache_needs_updating |= in_cache != in)
    in_cache = in;

  return cache_needs_updating;
}
