#ifndef UTIL_H
#define UTIL_H

#include "semi-algebraic.h"

/*
 * compute N
 */
inline bool cmp_2(long long N, int n, int k) {
  return N * N >= pow(2LL, n) * (1 + (k - 1) * (k - 2) / 2 + 2 * (N + 1) * k);
}

int compute_N(Semi_Algebraic set, int n) {
  unsigned int d = set.dimension();
  unsigned int k = set.max_degree();
  unsigned int K = set.degree();

  if (d == 1) {
    return k;
  } else if (d == 2) {
    int l, r, m;
    for (r = 2; !cmp_2(r, n, k); r <<= 1);
    l = r >> 1;
    while (l <= r) {
      m = (l + r) >> 1;
      if (cmp_2(m, n, k))
        r = m - 1;
      else
        l = m + 1;
    }
    return l;
  } else {
    return K * pow(2LL, d + n) + 1;
  }
}
#endif
