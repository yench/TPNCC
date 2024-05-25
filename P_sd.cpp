//P_sd.cpp
#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
List P_sd (const arma::mat& df_beta12, const arma::mat& df_bint12, const arma::mat& df0_beta12, const arma::mat& df0_bint12,
           const arma::mat& df_beta13, const arma::mat& df_bint13, const arma::mat& df0_beta13, const arma::mat& df0_bint13,
           const arma::mat& df_beta23, const arma::mat& df_bint23, const arma::mat& df0_beta23, const arma::mat& df0_bint23,
           const arma::mat& tsubtable, const arma::uvec& subset2, const arma::uvec& subset02, const arma::vec& zstar, 
           const arma::vec& beta12, const arma::vec& dLambda120, 
           const arma::vec& beta13, const arma::vec& dLambda130,
           const arma::vec& beta23, const arma::vec& dLambda230, int& d)
{
  arma::cube termary(3,3,d); 
  arma::uword n_df = df_bint12.n_rows, n_df0 = df0_bint12.n_rows;
  arma::mat diag3, term(3,3), pdbeta12sum(3,3), pdbeta13sum(3,3), pdbeta23sum(3,3), 
            sum_lambda(n_df,4), sum0_lambda(n_df0,4);
  diag3.eye(3,3);
  arma::mat P = diag3, pdbeta12, pdbeta13, pdbeta23, pdlambda;
  
  termary.each_slice() = diag3;
  arma::vec df_binthk(n_df), df0_binthk(n_df0);
  arma::uvec idxP = {0,3,6,7};
  double dLhk = 0;
  int hk, idxj;
  bool omg = (n_df0>2);

  for(int j=0; j<d; j++){
    hk = tsubtable(j,0);
    idxj = tsubtable(j,2)-1;
    if(hk==12) {
      dLhk = exp(dot(zstar, beta12))*dLambda120(idxj);
      termary(0,0,j) = 1-dLhk;
      termary(0,1,j) = dLhk;
    }
    if(hk==13) {
      dLhk = exp(dot(zstar, beta13))*dLambda130(idxj);
      termary(0,0,j) = 1-dLhk;
      termary(0,2,j) = dLhk;
    }
    if(hk==23) {
      dLhk = exp(dot(zstar, beta23))*dLambda230(idxj);
      termary(1,1,j) = 1-dLhk;
      termary(1,2,j) = dLhk;
    }
    P = P*termary.slice(j);
  }

  //compute sum of pd(beta) without zstar and sum of pd(lambda) * IF(lambda)
  for(int j=0; j<d; j++){
    hk = tsubtable(j,0);
    idxj = tsubtable(j, 2)-1;
    pdbeta12.eye(3,3);
    pdbeta13.eye(3,3);
    pdbeta23.eye(3,3);
    pdlambda.eye(3,3);

    for(int l=0; l<d; l++){
      if(l==j){
        term = termary.slice(j)-diag3;
        if(hk==12) {
          pdbeta12 = pdbeta12 * term;
          pdlambda = pdlambda * (term/dLambda120(idxj));
        } else if(hk==13){
          pdbeta13 = pdbeta13 * term;
          pdlambda = pdlambda * (term/dLambda130(idxj));
        } else if(hk==23){
          pdbeta23 = pdbeta23 * term;
          pdlambda = pdlambda * (term/dLambda230(idxj));
        }
      } else {
        if(hk==12) pdbeta12 = pdbeta12 * termary.slice(l);
        if(hk==13) pdbeta13 = pdbeta13 * termary.slice(l);
        if(hk==23) pdbeta23 = pdbeta23 * termary.slice(l);
        pdlambda = pdlambda * termary.slice(l);
      }
    }
    if(hk==12) {
      pdbeta12sum += pdbeta12;
      df_binthk = df_bint12.col(idxj);
      df0_binthk = df0_bint12.col(std::min<int>(df0_bint12.n_cols-1, idxj));
    } else if(hk==13) {
      pdbeta13sum += pdbeta13;
      df_binthk = df_bint13.col(idxj);
      df0_binthk = df0_bint13.col(std::min<int>(df0_bint13.n_cols-1, idxj));
    } else if(hk==23) {
      pdbeta23sum += pdbeta23;
      df_binthk.zeros();
      df0_binthk.zeros();
      df_binthk(subset2) = df_bint23.col(idxj);
      df0_binthk(subset02) = df0_bint23.col(std::min<int>(df0_bint23.n_cols-1, idxj));
    }
    sum_lambda += df_binthk * (pdlambda(idxP).t());
    sum0_lambda += df0_binthk * (pdlambda(idxP).t());
  }
  termary.clear();
  arma::vec df_beta23z(n_df);
  df_beta23z(subset2) = df_beta23 * zstar;
  arma::mat df_P = (df_beta12 * zstar) * (pdbeta12sum(idxP).t()) +
                   (df_beta13 * zstar) * (pdbeta13sum(idxP).t()) +
                   df_beta23z * (pdbeta23sum(idxP).t()) + sum_lambda;
  if(omg){
    df_beta23z.zeros();
    df_beta23z(subset2) = df0_beta23 * zstar;
    arma::mat df_P0 = (df0_beta12 * zstar) * (pdbeta12sum(idxP).t()) +
                      (df0_beta13 * zstar) * (pdbeta13sum(idxP).t()) +
                      df_beta23z * (pdbeta23sum(idxP).t()) + sum0_lambda;
    return List::create(Named("P") = P, Named("df_P") = df_P, Named("df_P0") = df_P0);
  }
  return List::create(Named("P") = P, Named("df_P") = df_P);
}