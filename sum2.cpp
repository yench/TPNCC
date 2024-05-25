//sum2.cpp
#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
arma::mat sum2 (const arma::mat& df){
  int ncols = df.n_cols;
  int nrows = df.n_rows/2;
  arma::mat out(nrows, ncols);
  arma::uvec id(nrows, arma::fill::ones);
  arma::uvec odd_id = 2*arma::cumsum(id)-1;
  arma::uvec eve_id = odd_id-1;
	out = df.rows(odd_id) + df.rows(eve_id);
	return out;
}