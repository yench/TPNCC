//Rho.cpp
#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
 arma::mat Rho(const arma::vec& H, int control_num, const arma::vec& control_times, const arma::vec& case_times){
   int nct = case_times.size();
	 arma::mat Rho(control_num, control_num);

	 for(int i=0; i<(control_num-1); i++){
		  for(int j=i+1; j<(control_num); j++){
		    arma::uvec ind(nct);
        ind = arma::find((case_times<control_times(i)) && (case_times<control_times(j)));
  			Rho(i,j) = arma::prod(H(ind))-1;
		}
	}
	return Rho;
}