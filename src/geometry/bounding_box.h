#pragma once

#include <interval.h>

#include <Eigen/Dense>

using namespace Eigen;

class BoundingBox {
public:
  Interval interval_x;
  Interval interval_y;
  Interval interval_z;

  BoundingBox() = default;

  BoundingBox(double x0, double x1,
              double y0, double y1,
              double z0, double z1) :
    interval_x(x0, x1),
    interval_y(y0, y1),
    interval_z(z0, z1)
  {}

  BoundingBox(const BoundingBox &bounding_box) = default;

  void add(const Vector3d &p) {
    interval_x.add(p.x());
    interval_y.add(p.y());
    interval_z.add(p.z());
  }

  bool contains(const Vector3d &p) const {
    return interval_x.contains(p.x()) &&
           interval_y.contains(p.y()) &&
           interval_z.contains(p.z());
  }

  bool intercepts(const BoundingBox &bounding_box) const {
    return interval_x.intercepts(bounding_box.interval_x) &&
           interval_y.intercepts(bounding_box.interval_y) &&
           interval_z.intercepts(bounding_box.interval_z);  
  }

  bool intercepts(const Vector3d &origin, const Vector3d &direction) const {
     double t[9];
     t[1] = (interval_x.x0 - origin.x())/direction.x();
     t[2] = (interval_x.x1 - origin.x())/direction.x();
     t[3] = (interval_y.x0 - origin.y())/direction.y();
     t[4] = (interval_y.x1 - origin.y())/direction.y();
     t[5] = (interval_z.x0 - origin.z())/direction.z();
     t[6] = (interval_z.x1 - origin.z())/direction.z();
     t[7] = max(max(min(t[1], t[2]), min(t[3], t[4])), min(t[5], t[6]));
     t[8] = min(min(max(t[1], t[2]), max(t[3], t[4])), max(t[5], t[6]));
     return t[8] > 0 && t[7] < t[8];
  }
};