#include <Eigen/Sparse>

template <typename T>
void speye(const int n, const int m, Eigen::SparseMatrix<T> &I){
  
  int d = (m < n ? m : n );
  I = Eigen::SparseMatrix<T>(m,n);
  I.reserve(d);
  for(int i=0; i < d; i++){
    I.insert(i,i) = 1.0;
  }
  I.finalize();
}

template <typename T>
void speye(const int n, Eigen::SparseMatrix<T> &I){
  return speye(n,n,I);
}
