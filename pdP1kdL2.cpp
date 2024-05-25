//pdP1kdL2.cpp
#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
 arma::mat pdP1kdL2(const arma::mat& ind_mat, const arma::vec& p1k){
   
  int d = p1k.n_elem;
  
  arma::vec j(d, arma::fill::ones);
   
	arma::mat pdterm2 = cumsum(j * p1k.t() % ind_mat, 1);
  
	return pdterm2;
}
