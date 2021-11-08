#include <Kokkos_Core.hpp>

#include <mls_mpm.h>

#include <cstdio>

using namespace std;
using namespace Kokkos;
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
  initialize(argc, argv);

  /*FILE *out = fopen("gputest.csv", "w");
  for (int i = 0; i < 300; i++) {
    int n = 50*(i + 1);

    // View<int*> results("results", n);

    // View<int[1]> dummy("dummy");
    // parallel_for(1, KOKKOS_LAMBDA(int j) { dummy(0) = 0; });

    View<garbage[1]> dummy("dummy");

    time_point<steady_clock> start_time = steady_clock::now();
    parallel_for(n, KOKKOS_LAMBDA(int j) {
      for (int k = 0;; k++)
      {
        for (int l = 0; l < k; l++) {
          if (l*k == 5000000) {
            // results(j) = k;

            // if (dummy(0) == 6)
            //   dummy(0) = 123;

            atomic_add(&dummy(0), garbage());

            return;
          }
        }
      }
    });
    fence();
    time_point<steady_clock> end_time = steady_clock::now();

    fprintf(out, "%d,%d\n", n, (int)duration_cast<milliseconds>(end_time - start_time).count());
    cout << i << endl;
  }
  fclose(out);

  for (int i = 0; i < 40; i++)
    mls_mpm(i + 1);*/
  
  mls_mpm(20);

  finalize();
}
