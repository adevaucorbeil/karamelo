#include <Kokkos_Core.hpp>

#include <kokkos_test.h>
#include <mls_mpm.h>

int main(int argc,
         char *argv[])
{
  Kokkos::initialize(argc, argv);

  kokkos_test();
  mls_mpm();

  Kokkos::finalize();
}
