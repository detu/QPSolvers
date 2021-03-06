// Implementation of Newton Method with Line Search 
#include <Eigen/Dense>

struct NewtonLS{

  NewtonLS(Eigen::MatrixXd Q, Eigen::VectorXd b, Eigen::VectorXd x0, int maxIter=100): Q_(Q), b_(b), x0_(x0), maxIter_(maxIter){
      double alpha0 = 1.0;
      double beta1  = 1.0e-4;
      double beta2  = 0.99;
      double lambda = 2;
      double eps    = 1e-8;
      int    k      = 0;
      x_ = x0_;
      Eigen::VectorXd g  = Q_*x_ + b_;
      double gNorm = g.norm();
      while ( gNorm > eps && k <= maxIter_){
        Eigen::VectorXd d = Q_.llt().solve(b_);
//        double alpha = lineSearch(Q_,b_,g,x_,d,alpha0,beta1,beta2,lambda);
//        x_ -= alpha * d;
        x_ -= d;
        g  = Q_*x_ + b_;
        k++;
        gNorm = g.norm();
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

  double lineSearch(Eigen::MatrixXd &Q, Eigen::VectorXd &b,Eigen::VectorXd &g, Eigen::VectorXd &x, Eigen::VectorXd &d, double alpha0, double beta1, double beta2, double lambda){
      double alpha  = alpha0;
      double alphal = 0;
      double alphar = 999999;
      double deriv  = g.transpose()*d;
      double f        = (x.transpose()*Q*x + x.transpose()*b).value();
      int ok = 0;
      while (ok == 0){
          Eigen::VectorXd xnew = x + alpha*d;
          double fnew          = (xnew.transpose()*Q*xnew + xnew.transpose()*b).value();
          Eigen::VectorXd gnew = Q*xnew + b;

          if ( fnew > (f + alpha*beta1*deriv)){
              alphar = alpha;
              alpha  = (alphal + alphar) / 2.0;
          }
          else if (gnew.transpose()*d < beta2*deriv){
              alphal = alpha;
              if(alphar == 999999)
                  alpha = lambda*alpha;
              else
                  alpha = (alphal + alphar) / 2.0;
          } else{
              ok = 1;
          }
      }

      return alpha;
  }
};
