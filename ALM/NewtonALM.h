// Copyright 2020, https://github.com/PatWie/CppNumericalSolvers
#ifndef INCLUDE_CPPOPTLIB_SOLVER_ALM_BOUND_H_
#define INCLUDE_CPPOPTLIB_SOLVER_ALM_BOUND_H_

#define EIGEN_USE_MKL_ALL
#include "Armijo.h"
#include "Solver.h"  // NOLINT
#include "speye.h"
#include "MUMPSSupport"
#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>

namespace cppoptlib::solver {

template <typename function_t>
class NewtonBound : public Solver<function_t> {
 private:
  using Superclass  = Solver<function_t>;
  using state_t     = typename Superclass::state_t;

  using scalar_t   = typename function_t::scalar_t;
  using hessian_t  = typename function_t::hessian_t;
  using matrix_t   = typename function_t::matrix_t;
  using vector_t   = typename function_t::vector_t;
  using function_state_t = typename function_t::state_t;

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  int Order() const override { return 2; }

  void InitializeSolver(const function_state_t &initial_state) override {
    dim_ = initial_state.x.rows();
  }

  function_state_t OptimizationStep(const function_t &function,
                                    const function_state_t &current,
                                    const state_t & /*state*/) override {

    // create k = 0 for initial and next for the rest
    // check constraint violation
    // update Lagrange multiplier and constant
    // update penalty and constraint violation

    function_state_t next = current;

    constexpr scalar_t safe_guard = 1e-5;
    hessian_t identity;
    speye(dim_, identity);
    const hessian_t hessian = next.hessian + safe_guard * identity;

    //Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    //solver.compute(hessian);

    Eigen::MUMPSLDLT<Eigen::SparseMatrix<double>, Eigen::Upper|Eigen::Lower> solver;
    //Eigen::MUMPSLDLT<Eigen::SparseMatrix<double>, Eigen::Upper> solver;

    //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(hessian);

    const vector_t delta_x = solver.solve(std::move(-next.gradient));
    //const scalar_t rate = linesearch::Armijo<function_t, 2>::Search(next.x, next.H, next.f, next.Aeq, next.beq, next.lambda, next.c, delta_x, function);
    const scalar_t rate = 1.0;

    return function.Eval(boundProj(next.x + rate * delta_x, next.ub, next.lb), next.H, next.f, next.Aeq, next.beq, next.lb, next.ub, next.lambda, next.c, 2);
  }

 private:
  int dim_;

  vector_t boundProj(const vector_t &x, const vector_t &upper, const vector_t &lower){
      Eigen::SparseVector<double> px(x.size());
      px.setZero();
      px = x.cwiseMin(upper);
      px = px.cwiseMax(lower);
      return px;
  }

};

}  // namespace cppoptlib::solver

#endif  // INCLUDE_CPPOPTLIB_SOLVER_ALM_BOUND_H_
