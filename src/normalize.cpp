#include <RcppEigen.h>
#include <progress.hpp>
#include <cmath>

using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]

// [[Rcpp::export]]
Eigen::SparseMatrix<double> CosineNormSparse(Eigen::SparseMatrix<double> data, bool display_progress = true){
  Progress p(data.outerSize(), display_progress);
  Eigen::SparseMatrix<double> data2 = data.cwiseProduct(data);
  Eigen::ArrayXd col2Sums = data2.transpose() * Eigen::VectorXd::Ones(data2.cols());
  col2Sums = col2Sums.sqrt();
  for (int k=0; k < data.outerSize(); ++k){
    p.increment();
    for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
      data.coeffRef(it.row(), it.col()) = double(it.value()) / col2Sums[k];
    }
  }
  return data;
}

// [[Rcpp::export]]
Eigen::MatrixXd CosineNorm(Eigen::MatrixXd data, bool display_progress = true){
  Progress p(data.cols(), display_progress);
  Eigen::MatrixXd data2 = data.cwiseProduct(data);
  Eigen::ArrayXd col2Sums = data2.colwise().sum();
  col2Sums = col2Sums.sqrt();
  for (int k=0; k < data.cols(); ++k){
    p.increment();
    Eigen::ArrayXd c = data.col(k).array();
    data.col(k) = c / col2Sums[k];
  }
  return data;
}
