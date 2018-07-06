//#include "iRRAM/core.h"
#include "util.h"
#include "semi-algebraic.h"

#include <ctime>
#include <vector>
#include <iostream>

//using namespace iRRAM;

/*
  Computes the period.
*/
int period_backtracking(Semi_Algebraic const &set, int m, std::vector<mpq_class> &x) {
  int count = 0;
  unsigned int depth = x.size();
  if (set.dimension() == depth) {
    for (unsigned int i = 0; i < set.poly.size(); i++) {
      mpq_class value = set.poly[i]->apply(x);
      if (value <= 0)
        return 0;
    }
    count = 1;
  } else {
    for (int i = 0; i < m; i++) {
      x.push_back(mpq_class(2 * i + 1, 2 * m));
      count += period_backtracking(set, m, x);
      x.pop_back();
    }
  }
  return count;
}
mpq_class period(Semi_Algebraic const &set, int n, int m){
  std::vector<mpq_class> x;
  int count = period_backtracking(set, m, x);
  mpz_class den(m);
  den = mpz_pow(den, set.dimension());
  return mpq_class(count, den) * (mpq_class)set.ratio;
}

int compute(Semi_Algebraic const &set, int const &n, int const &N, char* const &answer) {
  mpq_class p = period(set, n, N);
  std::cout << p.get_str() << std::endl;
  if (answer != NULL) {
    mpq_class r(answer);
    std::cout << (abs(r - p) < mpq_pow(mpq_class(1, 2), n-1) ? "PASS" : "FAIL") << std::endl;
  }
  return 0;
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: period input_file prec [correct_value]\n");
    return 1;
  }

  Semi_Algebraic set = semi_algebraic_from_file(argv[1]);

  int n;
  if (sscanf(argv[2], "%d", &n) != 1) {
    fprintf(stderr, "wrong precision\n");
    return 1;
  }
  char *answer = (argc < 4 ? NULL : argv[3]);

  std::time_t started_at = time(NULL);
  int N = compute_N(set, n);
  //iRRAM_initialize(argc, argv);
  //iRRAM_exec<int, Semi_Algebraic, int, int, char*>(compute, set, n, N, answer);
  compute(set, n, N, answer);
  std::cout << "time: " << (time(NULL) - started_at) << " sec" << std::endl;

  for (unsigned int i = 0; i < set.poly.size(); i++)
    delete set.poly[i];
}
