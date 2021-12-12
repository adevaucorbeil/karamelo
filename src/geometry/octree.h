#pragma once

#include <bounding_box.h>
#include <facet.h>

#include <deque>
#include <unordered_set>

class OctreeNode : public BoundingBox {
public:
  bool children_are_leaves = true;
  union {
    OctreeNode *nodes[8];
    deque<const Facet *> leaves;
  };

  OctreeNode(BoundingBox bounding_box) :
    BoundingBox(bounding_box) {
    new(&leaves) deque<const Facet*>();
  }

  OctreeNode(const OctreeNode &node):
    BoundingBox(node) {
    if (children_are_leaves = node.children_are_leaves) {
      new(&leaves) deque<const Facet*>(node.leaves);
    }
    else {
      for (int i = 0; i < 8; i++) {
        nodes[i] = node.nodes[i];
      }
    }
  }

  OctreeNode(OctreeNode &&node):
    BoundingBox(node) {
    if (children_are_leaves = node.children_are_leaves) {
      new(&leaves) deque<const Facet*>(move(node.leaves));
    }
    else {
      for (int i = 0; i < 8; i++) {
        nodes[i] = node.nodes[i];
      }
    }
  }

  ~OctreeNode()
  {
    if (children_are_leaves)
    {
      leaves.~deque();
    }
  }

  void add(const Facet *facet, deque<OctreeNode> &nodes_deque,
           int max_leaves_per_node, int max_depth, int depth) {
    if (children_are_leaves) {
      leaves.push_back(facet);

      if (depth < max_depth && leaves.size() > max_leaves_per_node) {
        // turn into intermediate node
        children_are_leaves = false;

        deque<const Facet *> leaves_copy(move(leaves));

        leaves.~deque();

        // create 8 new children nodes
        for (int x = 0; x < 2; x++)
          for (int y = 0; y < 2; y++)
            for (int z = 0; z < 2; z++) {
              nodes_deque.emplace_back(BoundingBox(((2 - x)*interval_x.x0 +      x *interval_x.x1)/2,
                                                   ((1 - x)*interval_x.x0 + (1 + x)*interval_x.x1)/2,
                                                   ((2 - y)*interval_y.x0 +      y *interval_y.x1)/2,
                                                   ((1 - y)*interval_y.x0 + (1 + y)*interval_y.x1)/2,
                                                   ((2 - z)*interval_z.x0 +      z *interval_z.x1)/2,
                                                   ((1 - z)*interval_z.x0 + (1 + z)*interval_z.x1)/2));

              nodes[2*(2*x + y) + z] = &nodes_deque.back();
            }

        // add back leaves
        for (const Facet *leaf: leaves_copy) {
          add(leaf, nodes_deque, max_leaves_per_node, max_depth, depth);
        }
      }
    }
    else {
      for (OctreeNode *node: nodes) {
        if (facet->intersects(*node)) {
          node->add(facet, nodes_deque, max_leaves_per_node, max_depth, depth + 1);
        }
      }
    }
  }

  void intersections(const Vector3d &origin, const Vector3d &direction, unordered_set<const Facet *> &facets) const {
    int intersections = 0;

    if (children_are_leaves) {
      for (const Facet *facet: leaves)
        if (facet->intersects(origin, direction))
          facets.insert(facet);
    }
    else {
      for (const OctreeNode *node: nodes)
        if (node->intercepts(origin, direction))
          node->intersections(origin, direction, facets);
    }
  }
};

class Octree {
public:
  deque<OctreeNode> nodes;
  int max_leaves_per_node;
  int max_depth;
  
  Octree() = default;

  Octree(BoundingBox bounding_box,
         int max_leaves_per_node,
         int max_depth) :
    nodes{ bounding_box },
    max_leaves_per_node(max_leaves_per_node),
    max_depth(max_depth)
  {}

  Octree(const Octree &) = default;
  Octree(Octree &&) = default;
  Octree &operator=(const Octree &) = default;
  Octree &operator=(Octree &&) = default;

  void add(const Facet &facet) {
    nodes.front().add(&facet, nodes, max_leaves_per_node, max_depth, 0);
  }

  unordered_set<const Facet *> intersections(const Vector3d &origin, const Vector3d &direction) const {
    unordered_set<const Facet *> facets;

    nodes.front().intersections(origin, direction, facets);

    return facets;
  }
};