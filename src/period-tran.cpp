#include "iRRAM/core.h"
#include "util.h"
#include "semi-algebraic.h"

#include <ctime>
#include <vector>

using namespace iRRAM;

/*
  Computes the period.
*/
int period_backtracking(Semi_Algebraic const &set, int m, std::vector<REAL> &x, std::vector<REAL> &tr) {
  int count = 0;
  unsigned int depth = x.size();
  if (set.dimension() == depth) {
    std::vector<REAL> xp;
    for (unsigned int i = 0; i < x.size(); i++)
      xp.push_back(x[i] + tr[i]);
    LAZY_BOOLEAN lzbl = true;
    for (unsigned int i = 0; i < set.poly.size(); i++)
      lzbl = lzbl && (set.poly[i]->apply(xp) > 0);
		if (lzbl)
      count++;
  } else {
    for (int i = 0; i < m; i++) {
      x.push_back((REAL) i / m);
      count += period_backtracking(set, m, x, tr);
      x.pop_back();
    }
  }
  return count;
}
REAL period(Semi_Algebraic const &set, int n, int m){
  std::vector<int> p = {2, 3, 5, 7};
  std::vector<REAL> tr;
  REAL r = 2;
  for (unsigned int i = 0; i < set.dimension(); i++) {
    tr.push_back(exp(power(p[i], REAL(1) / 2)));
    while (tr[i] > r)
      r *= 2;
  }
  for (unsigned int i = 0; i < set.dimension(); i++)
    tr[i] /= r * m;

  std::vector<REAL> x;
  int count = period_backtracking(set, m, x, tr);
  return count / power(m, set.dimension()) * set.ratio;
}

int compute(Semi_Algebraic const &set, int const &n, int const &N, char* const &answer) {
  REAL p = period(set, n, N);
  cout << p << std::endl;
  if (answer != NULL) {
    REAL r = answer;
    cout << (abs(r - p) < power(2, -n+1) ? "PASS" : "FAIL") << std::endl;
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
  iRRAM_initialize(argc, argv);
  iRRAM_exec<int, Semi_Algebraic, int, int, char*>(compute, set, n, N, answer);
  std::cout << "time: " << (time(NULL) - started_at) << " sec" << std::endl;

  for (unsigned int i = 0; i < set.poly.size(); i++)
    delete set.poly[i];
}
