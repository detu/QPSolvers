// Implementation of Newton Method with Line Search 
#include <Eigen/Dense>

struct NewtonLS{

  NewtonLS(Eigen::MatrixXd Q, Eigen::VectorXd b, int numIter=100): Q_(Q), b_(b), numIter_(numIter){
      double alpha0 = 1.0;
      double beta1  = 1.0e-4;
      double beta2  = 0.99;
      double lambda = 2;
      while ( norm(g) >= eps || k <= numIter_){
          //[L,tau] = modifiedCholesky(Q_);
          //z = L \ g;
          //d = - L' \ z;
          //alpha = lineSearch();
          //xk = xk + alpha * d;
          // k = k + 1;
      }
  }

  Eigen::VectorXd solve() const noexcept{
    return x_;
  }
  
private:
  Eigen::MatrixXd Q_;
  Eigen::VectorXd b_;
  Eigen::VectorXd x_;
  int numIter_;

  // modifiedCholesky()

  // lineSearch()
};
