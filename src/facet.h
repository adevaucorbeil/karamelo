#pragma once

#include "constants.h"
#include "bounding_box.h"

#include <array>

using namespace std;

class Facet : public array<Vector3d, 3> {
 public:
  Vector3d normal;

  bool intersects(const Vector3d &origin, const Vector3d &direction) const {
    Vector3d edge0 = at(1) - at(0);
    Vector3d edge1 = at(2) - at(0);
    Vector3d h = direction.cross(edge1);
    double a = edge0.dot(h);
    if (abs(a) < eps)
      return false;    // This ray is parallel to this triangle.
    double f = 1/a;
    Vector3d s = origin - at(0);
    double u = f*s.dot(h);
    if (u < 0.0 || u > 1.0)
      return false;
    Vector3d q = s.cross(edge0);
    double v = f*direction.dot(q);
    if (v < 0.0 || u + v > 1.0)
      return false;
    // At this stage we can compute t to find out where the intersection point is on the line.
    double t = f*edge1.dot(q);
    if (t > eps) // ray intersection
    {
      return true;
    }
    // This means that there is a line intersection but not a ray intersection.
    return false;
  }

  bool intersects(const BoundingBox &bounding_box) const {
    // REVISIT: use AABB vs. triangle intersection test for better performance
    BoundingBox this_bounding_box;

    for (const Vector3d &vertex: *this) {
      this_bounding_box.add(vertex);
    }

    return this_bounding_box.intercepts(bounding_box);
  }
};