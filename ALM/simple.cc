// Copyright 2020, https://github.com/PatWie/CppNumericalSolvers

#include <iostream>

#include "function.h"
#include "bfgs.h"
#include "conjugated_gradient_descent.h"
#include "gradient_descent.h"
#include "lbfgs.h"
#include "lbfgsb.h"
#include "newton_descent.h"
#include "alm_bound.h"

using FunctionXd = cppoptlib::function::Function<double>;

//class Function : public FunctionXd {
// public:
//  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
//
//  using FunctionXd::hessian_t;
//  using FunctionXd::vector_t;
//
//  scalar_t operator()(const vector_t &x) const override {
//    return 5 * x[0] * x[0] + 100 * x[1] * x[1] + 5;
//  }
//
//  void Gradient(const vector_t &x, vector_t *grad) const override {
//    (*grad)[0] = 2 * 5 * x[0];
//    (*grad)[1] = 2 * 100 * x[1];
//  }
//
//  void Hessian(const vector_t &x, hessian_t *hessian) const override {
//    (*hessian)(0, 0) = 1200 * x[0] * x[0] - 400 * x[1] + 2;
//    (*hessian)(0, 1) = -400 * x[0];
//    (*hessian)(1, 0) = -400 * x[0];
//    (*hessian)(1, 1) = 200;
//  }
//};

class Function : public FunctionXd {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using FunctionXd::hessian_t;
    using FunctionXd::vector_t;

    scalar_t operator()(const vector_t &x, const scalar_t lambda, const scalar_t p) const override {
        return 2 * ( x[0]*x[0] + x[1]*x[1] - 1 ) - x[0] + lambda * (x[0]*x[0] + x[1]*x[1] - 1) + p/2 * ((x[0]*x[0] + x[1]*x[1] - 1)*(x[0]*x[0] + x[1]*x[1] - 1));
    }

    void Gradient(const vector_t &x, const scalar_t lambda, const scalar_t c, vector_t *grad) const override {
        (*grad)[0] = 2 * x[0] - 1 + 2*lambda*x[0] + c/2 * (2*(x[0]*x[0] + x[1]*x[1] - 1)*2*x[0]);
        (*grad)[1] = 2 * x[1] + 2*lambda*x[1] + c/2* (2*(x[0]*x[0] + x[1]*x[1] - 1)*2*x[1]);
    }

    void Hessian(const vector_t &x, const scalar_t lambda, const scalar_t c, hessian_t *hessian) const override {
        (*hessian)(0, 0) = 2 + 2*lambda + 2*c*(3*x[0]*x[0] + x[1]*x[1] - 1);
        (*hessian)(0, 1) = x[1]*x[1];
        (*hessian)(1, 0) = x[0]*x[0];
        (*hessian)(1, 1) = 2 + 2*lambda + 2*c*(x[0]*x[0] + 3*x[1]*x[1] - 1);
    }
};

int main(int argc, char const *argv[]) {
  // using Solver = cppoptlib::solver::NewtonDescent<Function>;
  // using Solver = cppoptlib::solver::GradientDescent<Function>;
  // using Solver = cppoptlib::solver::ConjugatedGradientDescent<Function>;
  // using Solver = cppoptlib::solver::Bfgs<Function>;
  // using Solver = cppoptlib::solver::Lbfgs<Function>;
  // using Solver = cppoptlib::solver::Lbfgsb<Function>;
    using Solver = cppoptlib::solver::ALMBound<Function>;

  Function f;
  Function::vector_t x(2);
  x << -1.0, 0.1;
  Function::scalar_t lambda(0);
  Function::scalar_t c(1);

  Function::vector_t lb(2);
  Function::vector_t ub(2);
  Function::matrix_t Aeq(2,2);
  Function::vector_t beq(2);

  auto state = f.Eval(x,lambda, c);
  std::cout << "this" << std::endl;

  std::cout << f(x,lambda,c) << std::endl;
  std::cout << state.gradient << std::endl;
  std::cout << state.hessian << std::endl;

  // std::cout << cppoptlib::utils::IsGradientCorrect(f, x) << std::endl;
  // std::cout << cppoptlib::utils::IsHessianCorrect(f, x) << std::endl;

  Solver solver;

  //auto [solution, solver_state] = solver.Minimize(f, x);

  // bikin loop disini buat outer loop dari ALM!
  // dari function class bikin method buat:
  // 1. compute objective function value dengan Augmented Lagrangian terms
  // 2. modified gradient and hessian values.

  //auto [solution, solver_state] = solver.Minimize(f, Aeq,beq,ub,lb,x);
    auto [solution, solver_state] = solver.Minimize(f,x,lambda,c);

  std::cout << "argmin " << solution.x.transpose() << std::endl;
  std::cout << "f in argmin " << solution.value << std::endl;
  std::cout << "iterations " << solver_state.num_iterations << std::endl;
  std::cout << "status " << solver_state.status << std::endl;

  return 0;
}
