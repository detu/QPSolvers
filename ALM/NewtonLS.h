// Implementation of Newton Method with Line Search 
#include <Eigen/Dense>

struct NewtonLS{

  NewtonLS(Eigen::MatrixXd Q, Eigen::VectorXd b, Eigen::VectorXd x0, int maxIter=100): Q_(Q), b_(b), x0_(x0), maxIter_(maxIter){
      double alpha0 = 1.0;
      double beta1  = 1.0e-4;
      double beta2  = 0.99;
      double lambda = 2;
      double eps    = 1.0e-8;
      int    k      = 0;
      Eigen::VectorXd xk = x0_;
      while ( b.norm() >= eps || k <= maxIter_){
        Eigen::VectorXd d = Q_.llt().solve(b_);
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
  Eigen::VectorXd x0_;
  int maxIter_;

  // modifiedCholesky()

  double lineSearch(Eigen::MatrixXd &Q, Eigen::VectorXd &b, double alpha0, double beta1, double beta2, double lambda){

  }
};
