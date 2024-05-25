//df_P11.cpp
#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
 arma::mat df_P11(const arma::mat& df_P11_term, const arma::vec& P11){
   
	arma::mat out = df_P11_term.each_row() % P11.t();
  
	return out;
}
