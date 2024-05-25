//rowcumsum.cpp
#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
 arma::mat rowcumsum(const arma::mat& inmat){
   
	arma::mat out = cumsum(inmat, 1);
  
	return out;
}
