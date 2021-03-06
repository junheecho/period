#include "iRRAM/core.h"
#include "util.h"
#include "semi-algebraic.h"

#include <ctime>
#include <vector>

using namespace iRRAM;

/*
  Computes the period.
*/
void period_more_backtracking(Semi_Algebraic const &set, int m, std::vector<REAL> &x, unsigned int depth, std::vector<LAZY_BOOLEAN> &choice) {
	if (set.dimension() == depth) {
    LAZY_BOOLEAN lzbl = true;
    for (unsigned int i = 0; i < set.poly.size(); i++)
      lzbl = lzbl && (set.poly[i]->apply(x) > 0);
    choice.push_back(lzbl);
    choice.push_back(!lzbl);
	} else {
		REAL x0 = x[depth];
    int k = set.max_degree();
    REAL inc = (REAL) 1 / m / (k + 1);
		for (int i = 0; i <= k; i++) {
			period_more_backtracking(set, m, x, depth + 1, choice);
			x[depth] += inc;
		}
		x[depth] = x0;
	}
}
int period_backtracking(Semi_Algebraic const &set, int m, std::vector<REAL> &x) {
  int count = 0;
  unsigned int depth = x.size();
  if (set.dimension() == depth) {
		std::vector<LAZY_BOOLEAN> choice;
		period_more_backtracking(set, m, x, 0, choice);
		if (choose(choice) & 1)
      count++;
  } else {
    for (int i = 0; i < m; i++) {
      x.push_back((REAL) i / m);
      count += period_backtracking(set, m, x);
      x.pop_back();
    }
  }
  return count;
}
REAL period(Semi_Algebraic const &set, int n, int m) {
  std::vector<REAL> x;
  int count = period_backtracking(set, m, x);
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
