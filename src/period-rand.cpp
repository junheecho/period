#include "iRRAM/core.h"
#include "random.h"

#include <ctime>
#include <vector>

#define max(x, y) ((x) > (y) ? (x) : (y))

using namespace iRRAM;

struct  Monomial {
  std::vector<int> exp;
  int p, q;
};

struct Polynomial {
  int nm;
  std::vector<Monomial> term;
};

struct Semi_Algebraic {
  int d, p, q, np;
  std::vector<Polynomial> poly;
};

/*
  Computes the period.

  d: the dimension
  p: the polynomial
  n: the precision 2^-n
*/
int period_backtracking(Semi_Algebraic const &set, int m, std::vector<REAL> &x) {
  int count = 0;
  int depth = x.size();
  if (set.d == depth) {
    int check = 1;
    for (int k = 0; k < set.np; k++) {
      REAL v = 0;
      for (int i = 0; i < set.poly[k].nm; i++) {
        REAL w = (REAL) set.poly[k].term[i].p / set.poly[k].term[i].q;
        for (int j = 0; j < set.d; j++)
          w *= power(x[j], set.poly[k].term[i].exp[j]);
        v += w;
      }
      if (!(v > 0)) {
        check = 0;
        break;
      }
    }
    count += check;
  } else {
    REAL random = uniform_real();
    REAL size = (REAL) 1 / m;
    for (int i = 0; i < m; i++) {
      x.push_back((i + random) / m);
      count += period_backtracking(set, m, x);
      x.pop_back();
    }
  }
  return count;
}
REAL period(Semi_Algebraic const &set, int n){
  int k = 0;
  for (int l = 0; l < set.np; l++) {
    for (int i = 0; i < set.poly[l].nm; i++) {
      int deg = 0;
      for (int j = 0; j < set.d; j++)
        deg += set.poly[l].term[i].exp[j];
      k = max(k, deg);
    }
  }

  int m = round(k * (REAL) power(2, set.d + n)) + 1;

  std::vector<REAL> x;
  int count = period_backtracking(set, m, x);
  return count / power(m, set.d) * set.p / set.q;
}

int compute(Semi_Algebraic const &set, int const &n, char* const &answer) {
  REAL p = period(set, n);
  cout << p << std::endl;
  if (answer != NULL) {
    REAL r = answer;
    cout << (abs(r - p) < power(2, -n+1) ? "PASS" : "FAIL") << std::endl;
  }
  return 0;
}

int error(FILE *f) {
  fclose(f);
  fprintf(stderr, "wrong format\n");
  return 1;
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: period input_file prec [correct_value]\n");
    return 1;
  }

  FILE *f = fopen(argv[1], "r");
  if (f == NULL) {
    fprintf(stderr, "cannot open the file\n");
    return 1;
  }

  Semi_Algebraic set;
  if (fscanf(f, "%d%d%d%d", &set.d, &set.p, &set.q, &set.np) != 4) return error(f);
  for (int i = 0; i < set.np; i++) {
    Polynomial p;
    if (fscanf(f, "%d", &p.nm) != 1) return error(f);
    for (int j = 0; j < p.nm; j++) {
      Monomial m;
      if (fscanf(f, "%d%d", &m.p, &m.q) != 2) return error(f);
      m.exp.resize(set.d);
      for (int k = 0; k < set.d; k++)
        if (fscanf(f, "%d", &m.exp[k]) != 1) return error(f);
      p.term.push_back(m);
    }
    set.poly.push_back(p);
  }

  fclose(f);

  int n;
  if (sscanf(argv[2], "%d", &n) != 1) {
    fprintf(stderr, "wrong precision\n");
    return 1;
  }
  char *answer = (argc < 4 ? NULL : argv[3]);

  std::time_t started_at = time(NULL);
  iRRAM_initialize(argc, argv);
  iRRAM_exec<int, Semi_Algebraic, int, char*>(compute, set, n, answer);
  std::cout << "time: " << (time(NULL) - started_at) << " sec" << std::endl;
}
