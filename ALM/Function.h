// Copyright https://github.com/PatWie/CppNumericalSolvers, MIT license
#ifndef INCLUDE_CPPOPTLIB_FUNCTION_H_
#define INCLUDE_CPPOPTLIB_FUNCTION_H_
#define EIGEN_USE_MKL_ALL
#include "Eigen/Core"
#include "Eigen/SparseCore"
#include "derivatives.h"

namespace cppoptlib::function {

// Specifies a current function state.
template <class scalar_t, class vector_t, class matrix_t>
struct State {
  int dim;
  int order;

  scalar_t value = 0;  // The objective value.
  vector_t x;          // The current input value in x.
  vector_t gradient;   // The gradient in x.
  matrix_t hessian;    // The Hessian in x;

  matrix_t H;
  vector_t f;
  vector_t lb;         // lower bound constraint
  vector_t ub;         // upper bound
  matrix_t Aeq;        // equality constraint Aeq * x = beq
  vector_t beq;

  vector_t lambda;
  scalar_t c;

  // TODO(patwie): There is probably a better way.
  State() : dim(-1), order(-1) {}

  State(const int dim, const int order)
      : dim(dim),
        order(order),
        x(dim),
        gradient(dim),
        hessian(dim, dim) {}

  // (x,H, f, Aeq, beq, lambda,c);
  State operator=(const State<scalar_t, vector_t, matrix_t> &rhs) {
    assert(rhs.order > -1);
    dim = rhs.dim;
    order = rhs.order;
    value = rhs.value;
    x = rhs.x.eval();
    if (order >= 1) {
      gradient = rhs.gradient.eval();
    }
    if (order >= 2) {
      hessian = rhs.hessian.eval();
    }
    return *this;
  }
};

template <class TScalar, int TDim = Eigen::Dynamic>
class Function {

public:
  using scalar_t  = TScalar;
  using vector_t  = Eigen::SparseVector<TScalar, Eigen::ColMajor, int>;
  using hessian_t = Eigen::SparseMatrix<TScalar, Eigen::ColMajor, int>;
  using matrix_t  = Eigen::SparseMatrix<TScalar, Eigen::ColMajor, int>;

    static const int Dim = TDim;
public:
  Function() = default;
  virtual ~Function() = default;

  // Computes the value of a function.
  virtual scalar_t operator()(const vector_t &x,  const matrix_t &H, const vector_t &f, const matrix_t &Aeq, const vector_t &beq, const vector_t &lambda, const scalar_t c) const = 0;

  // Computes the gradient of a function.
  virtual void Gradient(const vector_t &x,  const matrix_t &H, const vector_t &f, const matrix_t &Aeq, const vector_t &beq, const vector_t &lambda, const scalar_t c, vector_t *grad) const {};

  // Computes the Hessian of a function.
  virtual void Hessian(const vector_t &x,  const matrix_t &H, const vector_t &f, const matrix_t &Aeq, const vector_t &beq, const vector_t &lambda, const scalar_t c, hessian_t *hessian) const {}

  virtual int Order() const { return 1; }

    virtual State <scalar_t, vector_t, hessian_t>
    Eval(const vector_t &x, const matrix_t &H, const vector_t &f, const matrix_t &Aeq, const vector_t &beq, const vector_t &lb, const vector_t &ub, const vector_t &lambda, const scalar_t c ,const int order = 2) const {
        State<scalar_t, vector_t, hessian_t> state(x.rows(), order);
        state.value = this->operator()(x,H, f, Aeq, beq, lambda,c);
        state.x      = x;
        state.lambda = lambda;
        state.H      = H;
        state.f      = f;
        state.Aeq    = Aeq;
        state.beq    = beq;
        state.c      = c;
        state.lb     = lb;
        state.ub     = ub;
        if (order >= 1) {
            this->Gradient(x, H, f, Aeq, beq, lambda, c, &state.gradient);
        }
        if (order >= 2) {
            this->Hessian(x, H, f, Aeq, beq, lambda, c, &state.hessian);
        }
        return state;
    }
};

}  // namespace cppoptlib::function

#endif  // INCLUDE_CPPOPTLIB_FUNCTION_H_
