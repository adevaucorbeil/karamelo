#include "mpm.h"
#include "grid.h"
#include "material.h"
#include "input.h"
#include <vector>
#include "memory.h"

using namespace std;


Grid::Grid(MPM *mpm) :
  Pointers(mpm)
{
  cout << "Creating new grid" << endl;

  x= NULL;
  v = v_update = NULL;
  b = f = NULL;
}

Grid::~Grid()
{
  memory->destroy(x);
  memory->destroy(v);
  memory->destroy(v_update);
  memory->destroy(b);
  memory->destroy(f);
}

void Grid::init(string cs){
  cellsize = input->parse(cs);
  cout << "Set grid cellsize to " << cellsize << endl;
}

void Grid::grow(int nn){
  nnodes = nn;

  if (x == NULL) x = new Eigen::Vector3d[nn];
  else {
    cout << "Error: x already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (v == NULL) v = new Eigen::Vector3d[nn];
  else {
    cout << "Error: v already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (v_update == NULL) v_update = new Eigen::Vector3d[nn];
  else {
    cout << "Error: v_update already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (b == NULL) b = new Eigen::Vector3d[nn];
  else {
    cout << "Error: b already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (f == NULL) f = new Eigen::Vector3d[nn];
  else {
    cout << "Error: f already exists, I don't know how to grow it!\n";
    exit(1);
  }
}
