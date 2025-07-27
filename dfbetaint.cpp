//dfbetaint.cpp
#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
List dfbetaint (const arma::mat& Z, const arma::mat& fullA, const arma::mat& naive_var,
                const arma::mat& eqbeta, const arma::mat& df_beta_no_eta,
                const arma::vec& wexpbZ, const arma::vec& trantime,  
                const arma::vec& lefttime, const arma::vec& eventime, 
                const arma::vec& wgt, const arma::vec& nevent_wt,
                const arma::uvec& subset, const arma::uvec& trantime_order,
                int& nsub, int& ntime, int& nbeta,
                bool& eta, bool& omg){
  
  arma::mat S1(ntime, nbeta);
  arma::vec S0(ntime);
  arma::uvec idxR;
  
  for(int j=0; j<ntime; j++){
    idxR = arma::find(lefttime<trantime(j) && eventime>=trantime(j));
    S0(j) = sum(wexpbZ(idxR));
    S1.row(j) = wexpbZ(idxR).t() * Z.rows(idxR);
  }
  
  arma::vec dL0t = nevent_wt/S0;
  arma::mat eqlambda(nsub, ntime), df_lambda, df_lambda_temp;
  arma::uvec id(nsub, arma::fill::ones);
  id = arma::cumsum(id);
  for(int j=0; j<ntime; j++){
    eqlambda.col(j) = wgt % (id==trantime_order(j)) - dL0t(j) * (wexpbZ % (lefttime<trantime(j) && eventime>=trantime(j)));
  }
  if(eta){
    int naux = fullA.n_cols;       // number of auxiliary variables in A
    arma::cube AZ(naux, nbeta, nsub);
    arma::mat subA, eqeta, I11(naux, naux), df_eta, T1, Psi2(nsub, nbeta), Ascore_sum(naux, nbeta), wAdNt_sum, df_beta;
    
    // eta: estimating equation and influence function ####
    subA = fullA.rows(subset);
    eqeta = -fullA * (!omg);
    eqeta.rows(subset) =  subA.each_col() % (wgt-(!omg)*arma::ones(nsub)); // each term in the vector is -IF_i(Psi_1) in (14)

    for(int i=0; i<nsub; i++){I11 = I11 + wgt(i)*subA.row(i).t()*subA.row(i);} // negative Psi_11 in appendix
    df_eta = eqeta * (-arma::inv(I11));  // transpose of IF_i(eta\hat) in (15)
    
    // beta: estimating equation and influence function ####
    T1.zeros(ntime, naux);
    for(int j=0; j<ntime; j++){
      idxR = arma::find(lefttime<trantime(j) && eventime>=trantime(j));
      T1.row(j) = wexpbZ(idxR).t() * subA.rows(idxR);
    }
    Psi2.rows(trantime_order-1) = (Z.rows(trantime_order-1) - S1.each_col()/S0);
    Psi2 = Psi2.each_col() % wgt;
    wAdNt_sum.zeros(naux, ntime);
    for(int i=0; i<nsub; i++){
      AZ.slice(i) = subA.row(i).t() * Z.row(i);
      Ascore_sum = Ascore_sum + subA.row(i).t() * Psi2.row(i);
      wAdNt_sum = wAdNt_sum + wgt(i) * subA.row(i).t() * ((i+1)==trantime_order).t();
    }
    Psi2.clear();

    // I12 is the transpose of Psi_21   
    arma::mat S1T1(naux, nbeta), ST(naux, nbeta), I12(naux, nbeta);
    for(int j=0; j<ntime; j++){
      S1T1 = T1.row(j).t() * S1.row(j); 
      ST.zeros();
      for(int i=0; i<nsub; i++){
        ST = ST + wexpbZ(i) * (lefttime(i)<trantime(j)) * (trantime(j)<=eventime(i)) * AZ.slice(i);
      }
      I12 = I12 - (nevent_wt(j)/S0(j)) * ST + (nevent_wt(j)/pow(S0(j), 2)) * S1T1;
    }
    I12 = I12 + Ascore_sum;
    
    arma::mat eqbeta_full(df_eta.n_rows, nbeta);
    eqbeta_full.rows(subset) = eqbeta;
    df_beta = (eqbeta_full + df_eta * I12) * naive_var;
    
    df_lambda = -(df_beta * (S1.t()));
    df_lambda = df_lambda.each_row() % ((dL0t/S0).t());
    df_lambda.rows(subset) += eqlambda.each_row()/(S0.t());
    df_lambda_temp = df_eta * (wAdNt_sum - ((T1.each_col() % dL0t).t()));
    df_lambda = df_lambda + df_lambda_temp.each_row()/(S0.t());
    return List::create(Named("bint") = dL0t, Named("df_beta") = df_beta, Named("df_lambda") = df_lambda);
    
  } else {
    df_lambda = -(df_beta_no_eta * (S1.t()));
    df_lambda = df_lambda.each_row() % ((dL0t/S0).t());
    df_lambda += eqlambda.each_row() / (S0.t());
    return List::create(Named("bint") = dL0t, Named("df_beta") = df_beta_no_eta, Named("df_lambda") = df_lambda);
  }
}