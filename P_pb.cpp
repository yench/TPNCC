//P_pb.cpp
#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
List P_pb (const arma::mat& df_beta1, const arma::mat& df_bint1, const arma::mat& df0_beta1, const arma::mat& df0_bint1,
		   const arma::mat& df_beta2, const arma::mat& df_bint2, const arma::mat& df0_beta2, const arma::mat& df0_bint2,
		   const arma::mat& tsubtable, const arma::uvec& subset2, const arma::uvec& subset02, const arma::vec& zstar, 
           const arma::vec& beta1, const arma::vec& dLambda10, 
           const arma::vec& beta2, const arma::vec& dLambda20, int& d, int& d1)
{
  arma::cube termary(3,3,d); 
  arma::uword n_df = df_bint1.n_rows, n_df0 = df0_bint1.n_rows;
  arma::mat diag3, term12(3,3), term13(3,3), pdbeta12sum(3,3), pdbeta13sum(3,3), pdbeta23sum(3,3), 
            sum_lambda(n_df,4), sum0_lambda(n_df0,4);
  diag3.eye(3,3);
  arma::mat P = diag3, pdbeta12, pdbeta13, pdbeta23, pdlambda;

  termary.each_slice() = diag3;
  arma::vec dL12v(d1), dL13v(d1), df_binth(n_df), df0_binth(n_df0);
  arma::uvec idxP = {0,3,6,7};
  double dL12 = 0, dL13 = 0, dL23 = 0;
  int j1 = 0, hk, h, idxj;
  bool omg = (n_df0>2);

  for(int j=0; j<d; j++){
    hk = tsubtable(j,0);
    h = hk/10;
    idxj = tsubtable(j,3)-1;
    
    if(h==1){
      dL12 = exp(dot(zstar, beta1(arma::span(0,1))))*dLambda10(idxj);
      dL13 = exp(beta1(2) + dot(zstar, beta1(arma::span(3,4))))*dLambda10(idxj);
      termary(0,0,j) = 1-dL12-dL13;
      termary(0,1,j) = dL12;
      termary(0,2,j) = dL13;
      dL12v(j1) = dL12;
      dL13v(j1) = dL13;
      j1++;
      
    } else {
      dL23 = exp(dot(zstar, beta2))*dLambda20(idxj);
      termary(1,1,j) = 1-dL23;
      termary(1,2,j) = dL23;
    }
    P = P*termary.slice(j);
  }
  //compute sum of pd(beta) without zstar and sum of pd(lambda) * IF(lambda)
  j1 = 0;
  for(int j=0; j<d; j++){
    h = tsubtable(j,0)/10;
    idxj = tsubtable(j, 3)-1;
    pdbeta12.eye(3,3);
    pdbeta13.eye(3,3);
    pdbeta23.eye(3,3);
    pdlambda.eye(3,3);
    for(int l=0; l<d; l++){
      if(h==1){
        if(l==j){
          term12(0,0)=-dL12v(j1); 
          term12(0,1)= dL12v(j1);
          term13(0,0)=-dL13v(j1);
          term13(0,2)= dL13v(j1);
          pdbeta12 = pdbeta12 * term12;
          pdbeta13 = pdbeta13 * term13;
          pdlambda = pdlambda * ((termary.slice(j)-diag3)/dLambda10(idxj));
        } else {
          pdbeta12 = pdbeta12 * termary.slice(l);
          pdbeta13 = pdbeta13 * termary.slice(l);
          pdlambda = pdlambda * termary.slice(l);
        }
      } else if(h==2){
        if(l==j){
          pdbeta23 = pdbeta23 * (termary.slice(j)-diag3);
          pdlambda = pdlambda * ((termary.slice(j)-diag3)/dLambda20(idxj));
        } else {
          pdbeta23 = pdbeta23 * termary.slice(l);
          pdlambda = pdlambda * termary.slice(l);
        }
      }
    }

    if(h==1) {
      pdbeta12sum += pdbeta12;
      pdbeta13sum += pdbeta13;
      df_binth = df_bint1.col(idxj);
      df0_binth = df0_bint1.col(std::min<int>(df0_bint1.n_cols-1, idxj));
      j1++;
    } else if(h==2) {
      pdbeta23sum += pdbeta23;
      df_binth.zeros();
      df0_binth.zeros();
      df_binth(subset2) = df_bint2.col(idxj);
      df0_binth(subset02) = df0_bint2.col(std::min<int>(df0_bint2.n_cols-1, idxj));
    }

    sum_lambda += df_binth * (pdlambda(idxP).t());
    sum0_lambda += df0_binth * (pdlambda(idxP).t());
  }
  termary.clear();
  arma::vec df_beta2z(n_df);
  df_beta2z(subset2) = df_beta2 * zstar;
  arma::mat df_P = (df_beta1.cols(0,1) * zstar) * (pdbeta12sum(idxP).t()) +
                   (df_beta1.col(2) + df_beta1.cols(3,4) * zstar) * (pdbeta13sum(idxP).t()) +
                   df_beta2z * (pdbeta23sum(idxP).t()) + sum_lambda;
  if(omg){
    df_beta2z.zeros();
    df_beta2z(subset2) = df0_beta2 * zstar;
    arma::mat df_P0 = (df0_beta1.cols(0,1) * zstar) * (pdbeta12sum(idxP).t()) +
                      (df0_beta1.col(2) + df0_beta1.cols(3,4) * zstar) * (pdbeta13sum(idxP).t()) +
                       df_beta2z * (pdbeta23sum(idxP).t()) + sum0_lambda;
    return List::create(Named("P") = P, Named("df_P") = df_P, Named("df_P0") = df_P0);
  }
  return List::create(Named("P") = P, Named("df_P") = df_P);
}