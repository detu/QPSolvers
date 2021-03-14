#include <Eigen/Dense>

struct DirectQP{

  DirectQP(Eigen::MatrixXf Q, Eigen::VectorXf b): Q_(Q), b_(b){
    x_ = Q_.colPivHouseholderQr().solve(b_);
  }

  Eigen::VectorXf solve() const noexcept{
    return x_;
  }

private:
  Eigen::MatrixXf Q_;
  Eigen::VectorXf b_;
  Eigen::VectorXf x_;
};

