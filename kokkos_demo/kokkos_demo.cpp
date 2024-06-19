#include <Kokkos_Core.hpp>

#include <mls_mpm.h>

#include <cstdio>

using namespace std;
using namespace std::chrono;

class garbage {
  int an_int = 1;
  int padding[8];
  public:
  KOKKOS_INLINE_FUNCTION garbage operator+(const garbage other) const {
    garbage g;
    g.an_int = an_int + other.an_int; 
    g.padding[2] = 5;
    return g;
  }
};

int main(int argc,
         char *argv[]) {
  Kokkos::initialize(argc, argv);

  for (int i = 0; i < 55; i++)
    mls_mpm(i + 1);

  Kokkos::finalize();
}
