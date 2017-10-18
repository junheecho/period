#ifndef SEMI_ALGEBRAIC_H
#define SEMI_ALGEBRAIC_H

#include "iRRAM/core.h"

#include <vector>
#include <stack>

using namespace iRRAM;

class Rational {
private:
  int p, q;

public:
  Rational() : p(0), q(1) {}
  Rational(int p, int q) : p(p), q(q) {}

  operator REAL() const {
    return (REAL) p / q;
  }
};

class Polynomial {
public:
  virtual ~Polynomial() = 0;
  virtual unsigned int dimension() = 0;
  virtual unsigned int degree() = 0;
  virtual REAL apply(std::vector<REAL> x) = 0;
};

template<typename T>
class Monomial {
private:
  T coeff;
  std::vector<unsigned int> exp;
  unsigned int _degree;

public:
  Monomial<T>(T coeff, std::vector<unsigned int> exp) : coeff(coeff), exp(exp) {
    _degree = 0;
    for (unsigned int i = 0; i < exp.size(); i++)
      _degree += exp[i];
  }

  unsigned int dimension() { return exp.size(); }
  unsigned int degree() { return _degree; }

  REAL apply(std::vector<REAL> x) {
    REAL v = (REAL) coeff;
    for (unsigned int i = 0; i < exp.size(); i++)
      v *= power(x[i], exp[i]);
    return v;
  }
};

template<typename T>
class CanonicalPolynomial : public Polynomial {
private:
  std::vector<Monomial<T> > terms;
  unsigned int _dimension;
  unsigned int _degree;

public:
  CanonicalPolynomial(std::vector<Monomial<T> > terms) : terms(terms) {
    for (unsigned int i = 0; i < terms.size(); i++) {
      _dimension = std::max(_dimension, terms[i].dimension());
      _degree = std::max(_degree, terms[i].degree());
    }
  }
  ~CanonicalPolynomial() {}

  unsigned int dimension() { return _dimension; }
  unsigned int degree() { return _degree; }

  REAL apply(std::vector<REAL> x) {
    REAL v = 0;
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

public:
  RPNPolynomial(std::vector<RPNLiteral> rpn) : rpn(rpn) {
    _dimension = _degree = 0;
    for (unsigned int i = 0; i < rpn.size(); i++)
      if (rpn[i].var)
        _dimension = std::max(_dimension, rpn[i].var->n + 1);
    std::stack<int> deg, val;
    for (unsigned int i = 0; i < rpn.size(); i++) {
      if (rpn[i].var) {
        deg.push(1);
        val.push(0);
      } else if (rpn[i].con) {
        deg.push(0);
        val.push(rpn[i].con->n);
      } else {
        int d2 = deg.top(); deg.pop();
        int d1 = deg.top(); deg.pop();
        int v2 = val.top(); val.pop();
        int v1 = val.top(); val.pop();
        if (rpn[i].op == PLUS) {
          deg.push(std::max(d1, d2));
          val.push(v1 + v2);
        } else if (rpn[i].op == MINUS) {
          deg.push(std::max(d1, d2));
          val.push(v1 - v2);
        } else if (rpn[i].op == TIMES) {
          deg.push(d1 + d2);
          val.push(v1 * v2);
        } else if (rpn[i].op == DIV) {
          if (d2) throw std::invalid_argument("not a polynomial");
          deg.push(d1);
          val.push(v1 / v2);
        } else if (rpn[i].op == POWER) {
          if (d2) throw std::invalid_argument("not a polynomial");
          deg.push(d1 * v2);
          val.push(pow(v1, v2));
        } else throw std::invalid_argument("not a polynomial");
      }
    }
    _degree = deg.top(); deg.pop();
    if (!deg.empty()) throw std::invalid_argument("not a polynomial");
  }
  ~RPNPolynomial() {
    for (unsigned int i = 0; i < rpn.size(); i++) {
      delete rpn[i].var;
      delete rpn[i].con;
    }
  }

  unsigned int dimension() { return _dimension; }
  unsigned int degree() { return _degree; }

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
};

class Semi_Algebraic {
private:
  unsigned int _dimension;
  unsigned int _degree;

public:
  std::vector<Polynomial*> poly;
  Rational ratio;

  Semi_Algebraic(std::vector<Polynomial*> poly, Rational ratio) : poly(poly), ratio(ratio) {
    _dimension = _degree = 0;
    for (unsigned int i = 0; i < poly.size(); i++) {
      _dimension = std::max(_dimension, poly[i]->dimension());
      _degree = std::max(_degree, poly[i]->degree());
    }
  }

  unsigned int dimension() const { return _dimension; }
  unsigned int degree() const { return _degree; }
};

Semi_Algebraic semi_algebraic_from_file(char *fn);

#endif
