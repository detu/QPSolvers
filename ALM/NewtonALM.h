// Copyright 2020, https://github.com/PatWie/CppNumericalSolvers
#ifndef INCLUDE_CPPOPTLIB_SOLVER_ALM_BOUND_H_
#define INCLUDE_CPPOPTLIB_SOLVER_ALM_BOUND_H_

#include "Armijo.h"
#include "Eigen/Dense"
#include "Solver.h"  // NOLINT

// Hardcore constraint Problem 19.3 of book Michael Bierlaier

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
    const hessian_t hessian =
        next.hessian + safe_guard * hessian_t::Identity(dim_, dim_);

    const vector_t delta_x = hessian.lu().solve(-next.gradient);
    //const scalar_t rate = linesearch::Armijo<function_t, 2>::Search(next.x, next.lambda, next.c, delta_x, function);
      const scalar_t rate = linesearch::Armijo<function_t, 2>::Search(next.x, next.H, next.f, next.Aeq, next.beq, next.lambda, next.c, delta_x, function);

    // project (next.x + rate * delta_x to bound constraint using KKT_boundProjection method
    //return function.Eval(next.x + rate * delta_x, next.lambda, next.c, next.lb, next.ub, 2);
      return function.Eval(next.x + rate * delta_x, next.H, next.f, next.Aeq, next.beq, next.lb, next.ub, next.lambda, next.c, 2);
  }

  // Add a method for gradient projection for bound constraint
  vector_t KKT_boundProjection(const function_state_t &current){
      // see gradproj.m or GradProj.h
  }

 private:
  int dim_;

};

}  // namespace cppoptlib::solver

#endif  // INCLUDE_CPPOPTLIB_SOLVER_ALM_BOUND_H_
