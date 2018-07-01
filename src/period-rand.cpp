#include "iRRAM/core.h"
#include "util.h"
#include "random.h"
#include "semi-algebraic.h"

#include <ctime>
#include <vector>

using namespace iRRAM;

int RANDOM_POINTS = 1;

/*
  Computes the period.
*/
int period_backtracking(Semi_Algebraic const &set, int m, std::vector<std::vector<REAL> > &x) {
  int count = 0;
  unsigned int depth = x[0].size();
  if (set.dimension() == depth) {
    std::vector<LAZY_BOOLEAN> choice;
    for (int j = 0; j < RANDOM_POINTS; j++) {
      LAZY_BOOLEAN lzbl = true;
      for (unsigned int i = 0; i < set.poly.size(); i++)
        lzbl = lzbl && (set.poly[i]->apply(x[j]) > 0);
      choice.push_back(lzbl);
      choice.push_back(!lzbl);
    }
    if (choose(choice) & 1)
      count++;
  } else {
    std::vector<REAL> random;
    for (int j = 0; j < RANDOM_POINTS; j++)
      random.push_back(uniform_real());
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < RANDOM_POINTS; j++)
        x[j].push_back((i + random[j]) / m);
      count += period_backtracking(set, m, x);
      for (int j = 0; j < RANDOM_POINTS; j++)
        x[j].pop_back();
    }
  }
  return count;
}
REAL period(Semi_Algebraic const &set, int n, int m) {
  std::vector<std::vector<REAL> > x(RANDOM_POINTS);
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

/*
 * the entry point
 */
int main(int argc, char **argv) {
  if (argc < 4) {
    fprintf(stderr, "usage: period input_file n_random_points prec [correct_value]\n");
    return 1;
  }

  Semi_Algebraic set = semi_algebraic_from_file(argv[1]);

  int n;
  sscanf(argv[2], "%d", &RANDOM_POINTS);
  if (sscanf(argv[3], "%d", &n) != 1) {
    fprintf(stderr, "wrong precision\n");
    return 1;
  }
  char *answer = (argc < 5 ? NULL : argv[4]);

  std::time_t started_at = time(NULL);
  int N = compute_N(set, n);
  iRRAM_initialize(argc, argv);
  iRRAM_exec<int, Semi_Algebraic, int, int, char*>(compute, set, n, N, answer);
  std::cout << "time: " << (time(NULL) - started_at) << " sec" << std::endl;

  for (unsigned int i = 0; i < set.poly.size(); i++)
    delete set.poly[i];
}
