#ifndef SEMI_ALGEBRAIC_H
#define SEMI_ALGEBRAIC_H

#include "iRRAM/core.h"

#include <vector>
#include <stack>
#include <algorithm>

#include <gmpxx.h>

using namespace iRRAM;

mpz_class mpz_pow(mpz_class base, unsigned int exp);
mpq_class mpq_pow(mpq_class base, unsigned int exp);

class Rational {
private:
  int p, q;

public:
  Rational() : p(0), q(1) {}
  Rational(int p, int q) : p(p), q(q) {}

  operator REAL() const {
    return (REAL) p / q;
  }
  operator mpq_class() const {
    return mpq_class(p, q);
  }
};

class Polynomial {
public:
  virtual ~Polynomial() = 0;
  virtual unsigned int dimension() = 0;
  virtual unsigned int degree() = 0;
  virtual unsigned int max_degree() = 0;
  virtual REAL apply(std::vector<REAL> x) = 0;
  virtual mpq_class apply(std::vector<mpq_class> x) = 0;
};

template<typename T>
class Monomial {
private:
  T coeff;
  std::vector<unsigned int> exp;
  unsigned int _degree;
  unsigned int _max_degree;

public:
  Monomial<T>(T coeff, std::vector<unsigned int> exp) : coeff(coeff), exp(exp) {
    _degree = _max_degree = 0;
    for (unsigned int i = 0; i < exp.size(); i++) {
      _degree += exp[i];
      _max_degree = std::max(_max_degree, exp[i]);
    }
  }

  unsigned int dimension() { return exp.size(); }
  unsigned int degree() { return _degree; }
  unsigned int max_degree() { return _max_degree; }

  REAL apply(std::vector<REAL> x) {
    REAL v = (REAL) coeff;
    for (unsigned int i = 0; i < exp.size(); i++)
      v *= power(x[i], exp[i]);
    return v;
  }

  mpq_class apply(std::vector<mpq_class> x) {
    mpq_class v = (mpq_class) coeff;
    for (unsigned int i = 0; i < exp.size(); i++)
      v *= mpq_pow(x[i], exp[i]);
    return v;
  }
};

template<typename T>
class CanonicalPolynomial : public Polynomial {
private:
  std::vector<Monomial<T> > terms;
  unsigned int _dimension;
  unsigned int _degree;
  unsigned int _max_degree;

public:
  CanonicalPolynomial(std::vector<Monomial<T> > terms) : terms(terms) {
    _dimension = _degree = _max_degree = 0;
    for (unsigned int i = 0; i < terms.size(); i++) {
      _dimension = std::max(_dimension, terms[i].dimension());
      _degree = std::max(_degree, terms[i].degree());
      _max_degree = std::max(_max_degree, terms[i].max_degree());
    }
  }
  ~CanonicalPolynomial() {}

  unsigned int dimension() { return _dimension; }
  unsigned int degree() { return _degree; }
  unsigned int max_degree() { return _max_degree; }

  REAL apply(std::vector<REAL> x) {
    REAL v = 0;
    for (unsigned int i = 0; i < terms.size(); i++)
      v += terms[i].apply(x);
    return v;
  }

  mpq_class apply(std::vector<mpq_class> x) {
    mpq_class v = 0;
    for (unsigned int i = 0; i < terms.size(); i++)
      v += terms[i].apply(x);
    return v;
  }
};

class RPNVariable {
public:
  unsigned int n;
  RPNVariable() : n(0) {}
  RPNVariable(unsigned int n) : n(n) {}
};
class RPNConstant {
public:
  int n;
  RPNConstant() : n(0) {}
  RPNConstant(int n) : n(n) {}
};
enum RPNOperator { NA, PLUS, MINUS, TIMES, DIV, POWER };
class RPNLiteral {
public:
  RPNVariable *var = 0;
  RPNConstant *con = 0;
  RPNOperator op = NA;
};

class RPNPolynomial : public Polynomial {
private:
  std::vector<RPNLiteral> rpn;
  unsigned int _dimension;
  unsigned int _degree;
  unsigned int _max_degree;

  std::vector<int> vector_init() {
    std::vector<int> z(_dimension);
    return z;
  }
  std::vector<int> vector_init(int k) {
    std::vector<int> z(_dimension);
    z[k] = 1;
    return z;
  }
  std::vector<int> vector_add(std::vector<int> &x, std::vector<int> &y) {
    std::vector<int> z(_dimension);
    for (unsigned int i = 0; i < _dimension; i++)
      z[i] = x[i] + y[i];
    return z;
  }
  std::vector<int> vector_mult(std::vector<int> &x, int k) {
    std::vector<int> z(_dimension);
    for (unsigned int i = 0; i < _dimension; i++)
      z[i] = x[i] * k;
    return z;
  }
  std::vector<int> vector_max(std::vector<int> &x, std::vector<int> &y) {
    std::vector<int> z(_dimension);
    for (unsigned int i = 0; i < _dimension; i++)
      z[i] = std::max(x[i], y[i]);
    return z;
  }

public:
  RPNPolynomial(std::vector<RPNLiteral> rpn) : rpn(rpn) {
    _dimension = _degree = _max_degree = 0;
    for (unsigned int i = 0; i < rpn.size(); i++)
      if (rpn[i].var)
        _dimension = std::max(_dimension, rpn[i].var->n + 1);
    std::stack<int> deg, val;
    std::stack<std::vector<int> > max_deg;
    for (unsigned int i = 0; i < rpn.size(); i++) {
      if (rpn[i].var) {
        deg.push(1);
        max_deg.push(vector_init(rpn[i].var->n));
        val.push(0);
      } else if (rpn[i].con) {
        deg.push(0);
        max_deg.push(vector_init());
        val.push(rpn[i].con->n);
      } else {
        int d2 = deg.top(); deg.pop();
        int d1 = deg.top(); deg.pop();
        std::vector<int> md2 = max_deg.top(); max_deg.pop();
        std::vector<int> md1 = max_deg.top(); max_deg.pop();
        int v2 = val.top(); val.pop();
        int v1 = val.top(); val.pop();
        if (rpn[i].op == PLUS) {
          deg.push(std::max(d1, d2));
          max_deg.push(vector_max(md1, md2));
          val.push(v1 + v2);
        } else if (rpn[i].op == MINUS) {
          deg.push(std::max(d1, d2));
          max_deg.push(vector_max(md1, md2));
          val.push(v1 - v2);
        } else if (rpn[i].op == TIMES) {
          deg.push(d1 + d2);
          max_deg.push(vector_add(md1, md2));
          val.push(v1 * v2);
        } else if (rpn[i].op == DIV) {
          if (d2) throw std::invalid_argument("not a polynomial");
          deg.push(d1);
          max_deg.push(md1);
          val.push(v1 / v2);
        } else if (rpn[i].op == POWER) {
          if (d2) throw std::invalid_argument("not a polynomial");
          deg.push(d1 * v2);
          max_deg.push(vector_mult(md1, v2));
          val.push(pow(v1, v2));
        } else throw std::invalid_argument("not a polynomial");
      }
    }
    _degree = deg.top(); deg.pop();
    if (!deg.empty()) throw std::invalid_argument("not a polynomial");
    _max_degree = 0;
    for (unsigned int i = 0; i < _dimension; i++)
      _max_degree = std::max(_max_degree, (unsigned int)max_deg.top()[i]);
  }
  ~RPNPolynomial() {
    for (unsigned int i = 0; i < rpn.size(); i++) {
      delete rpn[i].var;
      delete rpn[i].con;
    }
  }

  unsigned int dimension() { return _dimension; }
  unsigned int degree() { return _degree; }
  unsigned int max_degree() { return _max_degree; }

  REAL apply(std::vector<REAL> x) {
    std::stack<REAL> val;
    for (unsigned int i = 0; i < rpn.size(); i++) {
      if (rpn[i].var) {
        val.push(x[rpn[i].var->n]);
      } else if (rpn[i].con) {
        val.push(rpn[i].con->n);
      } else {
        REAL v2 = val.top(); val.pop();
        REAL v1 = val.top(); val.pop();
        if (rpn[i].op == PLUS) {
          val.push(v1 + v2);
        } else if (rpn[i].op == MINUS) {
          val.push(v1 - v2);
        } else if (rpn[i].op == TIMES) {
          val.push(v1 * v2);
        } else if (rpn[i].op == DIV) {
          val.push(v1 / v2);
        } else if (rpn[i].op == POWER) {
          val.push(power(v1, v2));
        }
      }
    }
    return val.top();
  }

  mpq_class apply(std::vector<mpq_class> x) {
    std::stack<mpq_class> val;
    for (unsigned int i = 0; i < rpn.size(); i++) {
      if (rpn[i].var) {
        val.push(x[rpn[i].var->n]);
      } else if (rpn[i].con) {
        val.push(rpn[i].con->n);
      } else {
        mpq_class v2 = val.top(); val.pop();
        mpq_class v1 = val.top(); val.pop();
        if (rpn[i].op == PLUS) {
          val.push(v1 + v2);
        } else if (rpn[i].op == MINUS) {
          val.push(v1 - v2);
        } else if (rpn[i].op == TIMES) {
          val.push(v1 * v2);
        } else if (rpn[i].op == DIV) {
          val.push(v1 / v2);
        } else if (rpn[i].op == POWER) {
          val.push(mpq_pow(v1, v2.get_num().get_ui()));
        }
      }
    }
    return val.top();
  }
};

#include <stdio.h>

class Semi_Algebraic {
private:
  unsigned int _dimension;
  unsigned int _degree;
  unsigned int _max_degree;

public:
  std::vector<Polynomial*> poly;
  Rational ratio;

  Semi_Algebraic(std::vector<Polynomial*> poly, Rational ratio) : poly(poly), ratio(ratio) {
    _dimension = _degree = _max_degree = 0;
    for (unsigned int i = 0; i < poly.size(); i++) {
      _dimension = std::max(_dimension, poly[i]->dimension());
      _degree = std::max(_degree, poly[i]->degree());
      _max_degree = std::max(_max_degree, poly[i]->max_degree());
    }
  }

  unsigned int dimension() const { return _dimension; }
  unsigned int degree() const { return _degree; }
  unsigned int max_degree() const { return _max_degree; }
};

Semi_Algebraic semi_algebraic_from_file(char *fn);

#endif
