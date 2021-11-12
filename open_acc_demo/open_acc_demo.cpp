// #include <mls_mpm.h>

#include <iostream>
#include <cstdio>
#include <chrono>

using namespace std;
using namespace std::chrono;

int main(int argc,
         char *argv[]) {
  FILE *out = fopen("gputest.csv", "w");
  for (int i = 0; i < 300; i++) {
    int n = 50*(i + 1);

    int *results = (int *)malloc(n*sizeof(int));

    time_point<steady_clock> start_time = steady_clock::now();

    #pragma acc parallel copyout(results[0:n])
    {
      #pragma acc loop gang
      for (int j = 0; j < n; j++) {
        bool flag = true;
        for (int k = 0; flag; k++) {
          for (int l = 0; l < k; l++) {
            if (l*k == 10000000) {
              results[j] = l;
              flag = false;
              break;
            }
          }
        }
      }
    }
    time_point<steady_clock> end_time = steady_clock::now();

    free(results);

    fprintf(out, "%d,%d\n", n, (int)duration_cast<milliseconds>(end_time - start_time).count());
    cout << i << ": " << (int)duration_cast<milliseconds>(end_time - start_time).count() << endl;
  }
  fclose(out);

  // mls_mpm(20);
}
