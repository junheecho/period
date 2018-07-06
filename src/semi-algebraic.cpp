#include "semi-algebraic.h"

#include <cstring>
#include <stdexcept>
#include <vector>

#define MAX_BUF 20

template<typename T>
T mp_pow(T q, unsigned int exp) {
  T p = 1;
  for (unsigned int e = 1; e <= exp; e <<= 1) {
    if (e & exp)
      p *= q;
    q *= q;
  }
  return p;
}

mpz_class mpz_pow(mpz_class base, unsigned int exp) {
  return mp_pow<mpz_class>(base, exp);
}
mpq_class mpq_pow(mpq_class base, unsigned int exp) {
  return mp_pow<mpq_class>(base, exp);
}

Polynomial::~Polynomial() {}

std::invalid_argument exc("invalid input");

std::invalid_argument clean(FILE *f, std::vector<Polynomial*> p) {
  fclose(f);
  for (unsigned int i = 0; i < p.size(); i++)
    delete p[i];
  return exc;
}

Semi_Algebraic semi_algebraic_from_file(char *fn) {
  FILE *f = fopen(fn, "r");
  if (f == NULL)
    throw std::invalid_argument("cannot open the file");
  unsigned d, np;
  int p, q;
  std::vector<Polynomial*> poly;
  if (fscanf(f, "%u%d%d%u", &d, &p, &q, &np) != 4)
    throw clean(f, poly);
  for (unsigned int i = 0; i < np; i++) {
    char type[MAX_BUF];
    if (fscanf(f, "%s", type) != 1)
      throw clean(f, poly);
    if (!strcmp("canonical", type)) {
      unsigned int nm;
      if (fscanf(f, "%u", &nm) != 1)
        throw clean(f, poly);
      std::vector<Monomial<Rational> > terms;
      for (unsigned int j = 0; j < nm; j++) {
        int p, q;
        if (fscanf(f, "%d%d", &p, &q) != 2)
          throw clean(f, poly);
        std::vector<unsigned int> exp(d);
        for (unsigned int k = 0; k < d; k++) {
          if (fscanf(f, "%u", &exp[k]) != 1)
            throw clean(f, poly);
        }
        terms.push_back(Monomial<Rational>(Rational(p, q), exp));
      }
      poly.push_back(new CanonicalPolynomial<Rational>(terms));
    } else if (!strcmp("rpn", type)) {
      std::vector<RPNLiteral> rpn;
      char buf[MAX_BUF];
      while (1) {
        if (fscanf(f, "%s", buf) != 1)
          throw clean(f, poly);
        RPNLiteral lit;
        if (buf[0] == '.') {
          break;
        } else if (buf[0] == '+') {
          lit.op = PLUS;
        } else if (buf[0] == '-') {
          lit.op = MINUS;
        } else if (buf[0] == '*') {
          lit.op = TIMES;
        } else if (buf[0] == '/') {
          lit.op = DIV;
        } else if (buf[0] == '^') {
          lit.op = POWER;
        } else if (buf[0] == 'x') {
          lit.var = new RPNVariable();
          sscanf(buf + 1, "%u", &lit.var->n);
        } else {
          lit.con = new RPNConstant();
          sscanf(buf, "%u", &lit.con->n);
        }
        rpn.push_back(lit);
      }
      poly.push_back(new RPNPolynomial(rpn));
    }
  }
  fclose(f);
  return Semi_Algebraic(poly, Rational(p, q));
}
