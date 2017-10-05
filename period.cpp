#include "iRRAM.h"

#include <ctime>
#include <vector>

#define max(x, y) ((x) > (y) ? (x) : (y))

using namespace iRRAM;

/*
  Computes the period.

  d: the dimension
  p: the polynomial
  box_min, box_max: the bounding box
  n: the precision 2^-n
*/
int period_backtracking(int d, std::vector<std::vector<int> > &p, std::vector<REAL> &box_min, std::vector<REAL> &box_max, int m, std::vector<REAL> &x) {
  int count = 0;
  int depth = x.size();
  if (d == depth) {
    REAL v = 0;
    for (int i = 0; (unsigned) i < p.size(); i++) {
      REAL w = p[i][d];
      for (int j = 0; j < d; j++)
        w *= power(x[j], p[i][j]);
      v += w;
    }
    count = (v > 0 ? 1 : 0);
  } else {
    REAL random = (REAL) 1 / (REAL) 2;
    REAL size = (box_max[depth] - box_min[depth]) / m;
    for (int i = 0; i < m; i++) {
      x.push_back(box_min[depth] + size * (i + random));
      count += period_backtracking(d, p, box_min, box_max, m, x);
      x.pop_back();
    }
  }
  return count;
}
REAL period(int d, std::vector<std::vector<int> > p, std::vector<REAL> box_min, std::vector<REAL> box_max, int n){
  int k = 0;
  for (int i = 0; (unsigned) i < p.size(); i++) {
    int deg = 0;
    for (int j = 0; j < d; j++)
      deg += p[i][j];
    k = max(k, deg);
  }

  REAL S = 1;
  for (int i = 0; i < d; i++)
    S *= box_max[i] - box_min[i];

  int m = round(S * k * (REAL) power(2, d + n)) + 1;

  std::vector<REAL> x;
  int count = period_backtracking(d, p, box_min, box_max, m, x);
  return count * S / power(m, d);
}

// pi : 3.141592...
void run_pi(int n) {
  static std::time_t started_at = time(NULL);
  int d = 2;
  std::vector<std::vector<int> > p = {
    { 2, 0, -1 },
    { 0, 2, -1 },
    { 0, 0,  1 }
  };
  std::vector<REAL> box_min = { -1, -1 };
  std::vector<REAL> box_max = {  1,  1 };
  REAL pr = period(d, p, box_min, box_max, n);
  REAL error = pi() - pr;
  bool pass = abs(error) < power(2, -n);

  cout << (pass ? "PASS" : "FAIL");
  cout << " pi = " << pr << " " << error;
  cout << " [" << time(NULL) - started_at << " sec]\n";
}

// log 2 : 0.693147...
void run_log2(int n) {
  static std::time_t started_at = time(NULL);
  int d = 2;
  std::vector<std::vector<int> > p = {
    { 1, 1, -1 },
    { 0, 0,  1 }
  };
  std::vector<REAL> box_min = { 1, 0 };
  std::vector<REAL> box_max = { 2, 1 };
  REAL pr = period(d, p, box_min, box_max, n);
  REAL error = (REAL) log(2) - pr;
  bool pass = abs(error) < power(2, -n);

  cout << (pass ? "PASS" : "FAIL");
  cout << " log 2 = " << pr << " " << error;
  cout << " [" << time(NULL) - started_at << " sec]\n";
}

void compute() {
  int n;

  cout << "Input n for precision 2^-n: ";
  cin >> n;

  cout << "\n";

  //run_pi(n);
  run_log2(n);
}
