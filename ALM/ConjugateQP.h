#include <Eigen/Dense>

struct ConjugateQP{

  ConjugateQP(Eigen::MatrixXf Q, Eigen::VectorXf b, Eigen::VectorXf x0): Q_(Q), b_(b), x0_(x0){
    // code Conjugate Gradient algorithm here!
  }

private:
  Eigen::MatrixXf Q_;
  Eigen::VectorXf b_;
  Eigen::VectorXf x0_;
  Eigen::VectorXf x_;
};
