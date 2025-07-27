# Weibull illness-death model data
simul.idm = function(n,               # cohort size
                     mctrl=1,         # number of controls per case
                     max.diag.t=6,    # administrative censoring time of the NCC study
                     max.death.t=10,  # administrative censoring time of cohort follow-up
                     max.enter=2,     # maximum late entry time into the cohort
                     Smax.censor=0.9, # baseline survival fraction from random censoring
                     Smax12=0.98, beta121=0.2, beta122=0.5, alpha12=2, # parameters for 1->2 transition   
                     Smax13=0.90, beta131=0.1, beta132=0.4, alpha13=2, # parameters for 1->3 transition
                     Smax23=0.70, beta231=0.3, beta232=0.7, alpha23=2, # parameters for 2->3 transition 
                     mu=c(z1=0, z2=0, u=0),                 # mean of Z1, Z2, and U
                     sig=c(sdz1=1, sdz2=1, sdu=1),          # marginal variance of Z1, Z2, and U 
                     cov12=0, cov1u=0, cov2u=0.75){         # cov(Z1, Z2), cov(Z1, U), cov(Z2, U)
  
  Sig = matrix(0, nrow=3, ncol=3); diag(Sig) = sig^2
  Sig[upper.tri(Sig)] = c(cov12, cov1u, cov2u)
  Sig[lower.tri(Sig)] = c(cov12, cov1u, cov2u)
  Covmat = data.table(MASS::mvrnorm(n=n, mu=mu, Sigma=Sig))
  
  # coefficients
  gamma12 = log(-log(Smax12)/(max.death.t^alpha12))
  gamma13 = log(-log(Smax13)/(max.death.t^alpha13))
  gamma23 = log(-log(Smax23)/(max.death.t^alpha23))
  
  eta12 = with(Covmat, exp((gamma12 + beta121*z1 + beta122*z2)))
  eta13 = with(Covmat, exp((gamma13 + beta131*z1 + beta132*z2)))
  eta23 = with(Covmat, exp((gamma23 + beta231*z1 + beta232*z2)))
  
  # generate transition times
  entrtime = runif(n, 0, max.enter)
  droptime = rexp(n, rate=-log(Smax.censor)/max.death.t) # drop out for any reason
  censtime = pmin(max.death.t-entrtime, droptime) # end of study
  
  evntime12 = rweibull(n, shape=alpha12, scale=eta12^(-1/alpha12))
  evntime13 = rweibull(n, shape=alpha13, scale=eta13^(-1/alpha13))
  evntime1 = pmin(evntime12, evntime13)
  evntime23 = (evntime12^alpha23-log(1-runif(n))/eta23)^(1/alpha23)
  
  nonttime = pmin(evntime1, censtime)
  ind.jmp1 = as.integer(evntime1<censtime)
  ind.nont = ifelse(ind.jmp1==1, as.integer(evntime12<evntime13), 0)
  termtime = ifelse(ind.jmp1==1 & ind.nont==0, evntime1, pmin(evntime23, censtime))
  ind.term = as.integer(termtime<censtime)
  
  # generate diagnosis time subject to the administrative censoring time of the incidence study
  diagtime = pmin(max.diag.t-entrtime, droptime, evntime1)
  ind.diag = as.integer(ind.nont==1 & evntime1<pmin(droptime, max.diag.t-entrtime))
  
  DATA = data.table(ID=1:n, Covmat, diagtime=diagtime, ind.diag=ind.diag, 
                    nonttime=nonttime, ind.12=ind.nont, termtime=termtime, 
                    ind.13=as.integer(ind.nont==0 & ind.term==1), ind.23=ind.nont*ind.term)
  
  # NCC sampling
  DATA[ , sample := .(ifelse(ind.diag==1, 
                             lapply(ID, 
                                    function(i){
                                      riskset = which(diagtime>=diagtime[i] & ID!=ID[i])
                                      out = NA
                                      if (length(riskset)>mctrl) out = sample(riskset, mctrl)
                                      else if (length(riskset)==mctrl) out = riskset
                                      return(out)}), NA))]
  DATA[ , nrisk := rank(-diagtime)] 
  
  # ID and indicator for inclusion in the NCC study
  ph2id = DATA[ind.diag==1, .(ID, sample)]
  
  for(h in 1:mctrl){
    ph2id[ , paste0("ctrl", h) := as.integer(lapply(.I, function(x){sample[[x]][h]}))]
  }
  setnames(ph2id, "ID", "case")
  ph2id[ , matched.id := paste0(1:.N)]
  ph2id[ , sample := NULL]
  
  # compute Design weights
  ph2dat = melt(ph2id, id.vars="matched.id",
                measure.vars=c("case", paste0("ctrl", 1:mctrl)))
  ph2dat[ , `:=`(case = ifelse(variable=="case", 1, 0), ind.ncc = TRUE)]
  setnames(ph2dat, "value", "ID")
  nccwt = DATA[, .(ID, diagtime, ind.diag, nrisk, ind.ncc = FALSE)]
  nccwt[ph2dat, on=.(ID), ind.ncc := i.ind.ncc]
  nccwt[ , term := ifelse(ind.diag==1, 1- mctrl/(nrisk-1), 1)]
  setorder(nccwt, diagtime)
  nccwt[ , incl.prob := ifelse(ind.diag==1, 1, 1 - cumprod(term))]
  nccwt[ , `:=`(wgt = 1/incl.prob, diagtime = NULL, ind.diag = NULL, term = NULL)]
  
  DATA[ , sample := NULL]
  DATA[nccwt, `:=` (incl.prob = incl.prob, wgt = wgt, ind.ncc = ind.ncc), on=.(ID)]
  return(DATA)
}


# Breslow-type calibration
calib_bre = function (datph1,     # full cohort data (replicated or not) 
                      trantime,   # name of the transition time
                      status.var, # name of the transition status
                      covars){    # names of covariates used in step 2 of weight calibration
                                  # include fully observed covariates and predicted values of partially observed covariates
  subset = datph1[ , ind.ncc]
  wgts = datph1[ , wgt]
  
  form = as.formula(paste0("Surv(", trantime, ",", status.var, ") ~ ", paste(covars, collapse="+")))
  fit = coxph(form, data=datph1)
  DFBETA = residuals(fit, type="dfbeta")
  A = cbind(1, DFBETA)
  g = rakecal(A[subset, ], wgts[subset], colSums(A))
  if(is.null(g)) {print(paste("Brewgt", substr(status.var, 5, 6), "did not converge")); return(NULL)}
  
  return(list(wgt = g, A = A))    # wgt = Breslow weight, A = matrix of auxiliary variables for the Breslow weight calibration procedure
}


# raking
rakecal = function(A,            # matrix of auxiliary variables among those in the NCC sample
                   w0,           # design weights among those in the NCC sample
                   total,        # sum of auxiliary variables in the full cohort
                   max_iter=500, # maximum number of iterations
                   EPS=1e-11){   # threshold for convergence{
  
  lambda = as.matrix(rep(0, ncol(A)))
  w1 = as.vector(w0 * exp(A %*% lambda))
  conv = 0
  
  for (l in 1:max_iter) {
    phi = crossprod(A, w1) - total
    phiprim = crossprod(A * w1, A)
    lambda1 = lambda - MASS::ginv(phiprim, tol=.Machine$double.eps) %*% phi
    w1 = as.vector(w0 * exp(A %*% lambda1))    
    if (any(is.na(w1)) | any(is.infinite(w1))) {
      warning("No convergence")
      g = NULL
      break
    }
    tr = crossprod(A, w1)
    
    if (max(abs(tr - total)) < EPS) {
      conv = 1
      break
    } else lambda = lambda1
  }
  
  if (!conv){
    warning(paste0("Calibration does not converge after ", max_iter, "iterations."))
    g = NULL
  } else g = w1
  
  attributes(g) <- list(lambda = lambda1)
  
  return(g) # calibrated weights
} 

## fits the standard and proportional baselines model, estimate baseline intensities and predict transition probabilities
coxModel = function (datph1,        # full cohort data
                     datph1.repli,  # replicated full cohort data for the proportional baselines model
                     tau,           # tau = (tau0, tau1], the interval to obtain prediction for
                     method,        # "ful" = full cohort,    "des" = design weight
                                    # "bre" = Breslow weight, "shi" = Shin weight
                     jcov.byp=NULL, # matrix for joint inclusion probabilities (required for all weighted analyses)
                     cal=NULL){     # list of weights and matrices of auxiliary variables (required for calibrated weights)
                                    # list(wgt12=, A12=, wgt13=, A13=, wgt23=, A23=, wgt1=, A1=)
                                    #      standard model                            proportional baselines model
  nfull = nrow(datph1)          # cohort size
  subset2 = (datph1$ind.12==1)  # indicator of experiencing a 1->2 transition
  n2 = sum(subset2)             # number of 1->2 transitions in the cohort
  
  #### beta: relative hazards ####
  if(method=="ful") { # full cohort analysis
    subset1 = rep(TRUE, nfull)  # include all subjects in the analysis
    wgt12 = wgt13 = rep(1, nfull)
    wgt23 = rep(1, n2)
    wgt1 = rep(1, 2*nfull)
    dfeta = FALSE

  } else { # analysis with NCC data
    subset1 = datph1$ind.ncc   # include thoe in the NCC sample in the analysis
    
    if(method=="des"){ # Design weights
      wgt12 = wgt13 = datph1[subset1, wgt]
      wgt23 = datph1[subset1 & subset2, wgt]
      wgt1 = rep(wgt12, each=2)
      dfeta = FALSE

    } else if(method %in% c("bre", "shi")){  # calibrated weights
      wgt12 = cal$wgt12
      wgt13 = cal$wgt13 
      wgt23 = cal$wgt23
      wgt1  = cal$wgt1
      dfeta = TRUE} 
  }
  #### fitting standard models ####
  fit12 = coxph(Surv(nonttime, ind.12) ~ z1 + z2, data=datph1[subset1], weights=wgt12, timefix=FALSE) 
  fit13 = coxph(Surv(nonttime, ind.13) ~ z1 + z2, data=datph1[subset1], weights=wgt13, timefix=FALSE)
  fit23 = coxph(Surv(nonttime, termtime, ind.23) ~ z1 + z2, data=datph1[subset1 & subset2, ], weights=wgt23, timefix=FALSE)
  #### fitting proportional baselines model to transitions out of the healthy state ####
  fit1 =  coxph(Surv(nonttime, ind.1) ~ z1 + z2 + z3 + z4 + z5, data=datph1.repli[rep(subset1, each=2)], weights=wgt1, timefix=FALSE)
  
  #### dLambdas: baseline intensities ####
  #### compute IF_1i of eta, beta, and dLambda ####
  df2.12 = df.coxph(fit12, dfeta=dfeta, fullA=cal$A12, subset=subset1)
  df2.13 = df.coxph(fit13, dfeta=dfeta, fullA=cal$A13, subset=subset1)
  df2.23 = df.coxph(fit23, dfeta=dfeta, fullA=cal$A23, subset=datph1[subset2, ind.ncc])
  df2.1  = df.coxph(fit1,  dfeta=dfeta, fullA=cal$A1,  subset=rep(subset1, each=2), pb=TRUE)
  df20.12 = df20.13 = df20.23 = df20.1 = NULL
  
  #### compute IF_2i of eta, beta, and dLambda ####
  if (method %in% c("bre", "shi")){ 
    df20.12 = df.coxph(fit12, dfeta=dfeta, fullA=cal$A12, subset=subset1, omg=TRUE)
    df20.13 = df.coxph(fit13, dfeta=dfeta, fullA=cal$A13, subset=subset1, omg=TRUE)
    df20.23 = df.coxph(fit23, dfeta=dfeta, fullA=cal$A23, subset=datph1[subset2, ind.ncc], omg=TRUE)
    df20.1 =  df.coxph(fit1,  dfeta=dfeta, fullA=cal$A1,  subset=rep(subset1, each=2), omg=TRUE, pb=TRUE)
  }
 
  #### objects for computation of variance ####
  wgts = datph1[ , wgt]
  wgt0 = rep(1, nfull) 
  wgt0[subset1] = wgts[subset1]
  ind.undg = (datph1$ind.diag==0)
  se.comp.ls = list(wgts = wgts, wgt0=wgt0, ind.undg=ind.undg, jcov.byp=jcov.byp, subset1=subset1, method=method, nh=nfull, n2=n2)
  se.comp.ls2 = list(wgts = wgts[subset2], wgt0=wgt0[subset2], ind.undg=ind.undg[subset2], 
                     jcov.byp=jcov.byp[datph1[subset1 & ind.diag==0, ind.12==1], datph1[subset1 & ind.diag==0, ind.12==1]], 
                     subset1=subset1[subset2], method=method, nh=n2)
  
  #### SE of beta ####
  beta12.se = se.comp(df2.12$df_beta, df20.12$df_beta, se.comp.ls)
  beta13.se = se.comp(df2.13$df_beta, df20.13$df_beta, se.comp.ls)
  beta23.se = se.comp(df2.23$df_beta, df20.23$df_beta, se.comp.ls2)
  beta1.se  = se.comp(df2.1$df_beta,  df20.1$df_beta,  se.comp.ls)
 
  #### Lambda: cumulative baseline intensities ####
  t.table = data.table(hk = rep(c(12, 13, 23), c(length(df2.12$bint), length(df2.13$bint), length(df2.23$bint))),
                       tj = c(attr(df2.12$bint, "time"), attr(df2.13$bint, "time"), attr(df2.23$bint, "time")),
                       tj.sd.idx = c(1:length(df2.12$bint), 1:length(df2.13$bint), 1:length(df2.23$bint)), key = "tj")
  t.table[ , tj.pb.idx := rank(tj), by=.(hk %/% 10)]
  bint12.idx = t.table[hk==12, tau[1]<tj & tj<=tau[2]]
  bint13.idx = t.table[hk==13, tau[1]<tj & tj<=tau[2]]
  bint23.idx = t.table[hk==23, tau[1]<tj & tj<=tau[2]]
  bint1.idx =  t.table[hk %in% c(12, 13), tau[1]<tj & tj<=tau[2]]
  
  bInt12 = sum(df2.12$bint[bint12.idx])
  bInt13 = sum(df2.13$bint[bint13.idx])
  bInt23 = sum(df2.23$bint[bint23.idx])
  bInt1  = sum(df2.1$bint[bint1.idx])
  bInt13.pb = bInt1*exp(fit1$coefficients["z3"])
  
  #### SE of Lambda ####
  se.bInt12 = se.comp(as.matrix(rowSums(df2.12$df_bint[ , bint12.idx])), as.matrix(rowSums(df20.12$df_bint[ , bint12.idx])), se.comp.ls)
  se.bInt13 = se.comp(as.matrix(rowSums(df2.13$df_bint[ , bint13.idx])), as.matrix(rowSums(df20.13$df_bint[ , bint13.idx])), se.comp.ls)
  se.bInt23 = se.comp(as.matrix(rowSums(df2.23$df_bint[ , bint23.idx])), as.matrix(rowSums(df20.23$df_bint[ , bint23.idx])), se.comp.ls2)
  se.bInt1  = se.comp(as.matrix(rowSums(df2.1$df_bint[, bint1.idx])), as.matrix(rowSums(df20.1$df_bint[, bint1.idx])), se.comp.ls)
  se.bInt13.pb = se.comp(as.matrix(rowSums(df2.1$df_bint[ , bint1.idx]) + df2.1$df_beta[ , "z3"]*bInt1)*exp(fit1$coefficients["z3"]), 
                         as.matrix(rowSums(df20.1$df_bint[, bint1.idx]) + df20.1$df_beta[, "z3"]*bInt1)*exp(fit1$coefficients["z3"]), se.comp.ls)

  #### transition probabilities ####
  tp.sd = tp.sd(t.table, zstar, tau, se.comp.ls, subset2,
                fit12$coefficients, df2.12, df20.12, 
                fit13$coefficients, df2.13, df20.13,
                fit23$coefficients, df2.23, df20.23)
  tp.pb = tp.pb(t.table, zstar, tau, se.comp.ls, subset2,
                fit1$coefficients,  df2.1,  df20.1, 
                fit23$coefficients, df2.23, df20.23)
 
  #### summarize result ####
  est.sd = c(fit12$coefficients, fit13$coefficients, fit23$coefficients, bInt12, bInt13, bInt23, 
             tp.sd[[1]]$P, tp.sd[[2]]$P, tp.sd[[3]]$P, tp.sd[[4]]$P, use.names=TRUE)
  est.pb = c(fit1$coefficients, fit23$coefficients, bInt1, bInt13.pb, bInt23, tp.pb[[1]]$P, tp.pb[[2]]$P, tp.pb[[3]]$P, tp.pb[[4]]$P, use.names=TRUE)
  
  se.sd = c(beta12.se, beta13.se, beta23.se, se.bInt12, se.bInt13, se.bInt23, tp.sd[[1]]$se, tp.sd[[2]]$se, tp.sd[[3]]$se, tp.sd[[4]]$se)
  se.pb = c(beta1.se, beta23.se, se.bInt1, se.bInt13.pb, se.bInt23, tp.pb[[1]]$se, tp.pb[[2]]$se, tp.pb[[3]]$se, tp.pb[[4]]$se)
  
  out.sd = rbind(est.sd, se.sd)
  out.pb = rbind(est.pb, se.pb)
  
  colnames(out.sd) = c("beta121", "beta122", "beta131", "beta132", "beta231", "beta232", "bInt12", "bInt13", "bInt23",
                       paste0(c("P11", "P12", "P13", "P23"), ".", rep(1:4, each=4)))
  colnames(out.pb) = c("beta121", "beta122", "gamma", "beta131", "beta132", "beta231", "beta232", "bInt1", "bInt13", "bInt2",
                       paste0(c("P11", "P12", "P13", "P23"), "." ,rep(1:4, each=4)))
  
  out.sd = data.table(out.sd, method = method, type = c("est", "se"), model = "sd")
  out.pb = data.table(out.pb, method = method, type = c("est", "se"), model = "pb")
  out = rbind(out.pb, out.sd, fill=TRUE) # estimates/predictions and SE 
  return(list(out = out, bint12 = df2.12$bint, bint13 = df2.13$bint))
}


# predict transition probabilities using the AJ estimator 
tp.sd = function(t.table,    # lookup table for indexing transition times within tau
                 zstar, tau, # covariates and time intervals to obtain predictions for    
                 se.comp.ls, # list of objects needed for computing SE
                 subset2,    # TRUE/FALSE vector of indicators for subjects ever in state 2
                 beta12, df12, df012,  # betas and baseline intensity estimates, and their influence functions for 1->2 transitions
                 beta13, df13, df013,  # betas and baseline intensity estimates, and their influence functions for 1->3 transitions
                 beta23, df23, df023){ # betas and baseline intensity estimates, and their influence functions for 2->3 transitions
  d = t.table[tau[1]<tj & tj<=tau[2], .N]
  nzstar = nrow(zstar)
  
  if(is.null(df012)){ # to avoid P_sd issueing errors for full-cohort and design weight analyses
    dummyM = matrix(0, 2, 2, dimnames = list(NULL, paste0("z", 1:2)));
    df012 = df013 = df023 = list(df_beta = dummyM, df_bint = dummyM)
    subset02 = rep(TRUE, 2)
  } else {subset02 = subset2}
  
  # design weights
  if(nrow(df12$df_beta)<se.comp.ls$nh) {subset2 = subset2[se.comp.ls$subset1]}
  
  # Aalen-Johansen estimator
  P.ret = lapply(1:nzstar, function(x){
    P_sd(df12$df_beta, df12$df_bint, df012$df_beta, df012$df_bint,
         df13$df_beta, df13$df_bint, df013$df_beta, df013$df_bint,
         df23$df_beta, df23$df_bint, df023$df_beta, df023$df_bint,
         as.matrix(t.table[tau[1]<tj & tj<=tau[2], ]), which(subset2)-1, which(subset02)-1, zstar[x, ], 
         beta12, df12$bint, beta13, df13$bint, beta23, df23$bint, d)})
  P = lapply(1:nzstar, function(x){c(P11 = P.ret[[x]]$P[1, 1], P12 = P.ret[[x]]$P[1, 2], P13 = P.ret[[x]]$P[1, 3], P23 = P.ret[[x]]$P[2, 3])})
  se = lapply(1:nzstar, function(x){se.comp(P.ret[[x]]$df_P, P.ret[[x]]$df_P0, se.comp.ls, c(rep(se.comp.ls$nh, 3), se.comp.ls$n2))})
  out = lapply(1:nzstar, function(x){list(P = P[[x]], se = se[[x]])})
  return(out) # list of predicted transition probabilities and their SE
}


tp.pb = function(t.table, zstar, tau, se.comp.ls, subset2, # see comments for tp.sd()
                 beta1, df1, df01,  # betas and baseline intensity estimates and their influence functions for transition out of state 1
                 beta2, df2, df02){ # betas and baseline intensity estimates and their influence functions for transition out of state 2
  d = t.table[tau[1]<tj & tj<=tau[2], .N]
  d1 = t.table[tau[1]<tj & tj<=tau[2] & hk %in% c(12,13), .N]
  nzstar = nrow(zstar)
  
  if(is.null(df01)){
    dummyM = matrix(0, 2, 5, dimnames = list(NULL, paste0("z", 1:5)));
    df01 = df02 = list(df_beta = dummyM, df_bint = dummyM)
    subset02 = rep(TRUE, 2)
  } else {subset02 = subset2}
  
  # design weights
  if(nrow(df1$df_beta)<se.comp.ls$nh) {subset2 = subset2[se.comp.ls$subset1]}
  
  P.ret = lapply(1:nzstar, function(x){
    P_pb(df1$df_beta, df1$df_bint, df01$df_beta, df01$df_bint,
         df2$df_beta, df2$df_bint, df02$df_beta, df02$df_bint,
         as.matrix(t.table[tau[1]<tj & tj<=tau[2], ]), which(subset2)-1, which(subset02)-1, zstar[x, ], 
         beta1, df1$bint, beta2, df2$bint, d, d1)})
  P = lapply(1:nzstar, function(x){c(P11 = P.ret[[x]]$P[1, 1], P12 = P.ret[[x]]$P[1, 2], P13 = P.ret[[x]]$P[1, 3], P23 = P.ret[[x]]$P[2, 3])})
  se = lapply(1:nzstar, function(x){se.comp(P.ret[[x]]$df_P, P.ret[[x]]$df_P0, se.comp.ls, c(rep(se.comp.ls$nh, 3), se.comp.ls$n2))})
  out = lapply(1:nzstar, function(x){list(P = P[[x]], se = se[[x]])})
  return(out) # list of predicted transition probabilities and their SE
}


# influence functions of log-relative hazards and baseline intensities
df.coxph = function(coxph_object, # a coxph object from fitting the standard or the proportional baselines model
                    dfeta=TRUE,   # compute IF(eta)? TRUE for calibrated weights
                    fullA=NULL,   # NULL for full cohort and design weights, the matrix of auxiliary variables for calibrated weights 
                    subset=NULL,  # a vector of TRUE/FALSE indicators for inclusion in the analysis
                    omg=FALSE,    # FALSE to compute IF_i, TRUE to compute IF_2i
                    pb=FALSE){    # TRUE if the coxph object comes from fitting a proportional baselines model
  
  coxph_detail = coxph.detail(coxph_object)
  if(is.null(coxph_detail$nevent.wt)) {
    nevent.wt = coxph_detail$nevent            # dN_i(t)
  } else {nevent.wt = coxph_detail$nevent.wt}  # w_i*dN_i(t)

  Z = model.matrix(coxph_object)
  bcoef = coxph_object$coefficients
  if(is.null(coxph_object$weights)) wgt = rep(1, nrow(Z)) else wgt = coxph_object$weights
  wexpbZ = wgt*exp(Z %*% bcoef)
  
  ncoly = ncol(coxph_object$y)
  eventime = coxph_object$y[ , ncoly-1]
  trantime = sort(eventime[coxph_object$y[ , ncoly]==1])
  trantime.order = as.integer(names(trantime))
  names(trantime) = NULL
  if(ncoly==2){lefttime = 0} else {lefttime = coxph_object$y[ , 1]}
  
  ntime = length(trantime) # number of unique transition times
  nfull = length(subset)   # number of rows in the full cohort data
  nsub = coxph_object$n    # number of subjects in the phase 2 sample
  nbeta = length(bcoef)
  
  eqbeta  = residuals(coxph_object, type="score",  weighted=TRUE)
  df_beta = residuals(coxph_object, type="dfbeta", weighted=TRUE)
  naive.var = coxph_object$naive.var

  if(is.null(fullA)) fullA = matrix(0, 2, 2)
  if(is.null(naive.var)) naive.var = matrix(0, 2, 2)
  
  ret = dfbetaint(Z, fullA, naive.var, eqbeta, df_beta, wexpbZ, trantime, rep(lefttime, ifelse(length(lefttime)==1, nsub, 1)), eventime,
                  wgt, nevent.wt, which(subset)-1, trantime.order, nsub, ntime, nbeta, dfeta, omg)
  
  if(pb==TRUE){
    ret$df_beta = sum2(ret$df_beta)
    ret$df_lambda = sum2(ret$df_lambda)
  }
  
  colnames(ret$df_beta) = names(bcoef)
  bint = ret$bint
  attr(bint, "time") = trantime
  
  return(list(bint = bint,              # estimated increments in the baseline intensity
              df_beta = ret$df_beta,    # influence functions of betas
              df_bint = ret$df_lambda)) # influence functions of increments in the baseline intensity
   
}


# joint inclusion probabilities
jointVP = function(datph2,   # the NCC data 
                   m){       # number of controls per case
  control.num = sum(datph2$ind.diag==0)
  case.times = datph2[ind.diag==1, diagtime, drop=TRUE]
  control.times = datph2[ind.diag==0, diagtime, drop=TRUE]
  H = datph2[ind.diag==1, (1-2*..m/(nrisk-1)+..m*(..m-1)/((nrisk-1)*(nrisk-2)))/(1-..m/(nrisk-1))^2]
  
  Rho = Rho(H, control.num, control.times, case.times)
  Rhony = Rho + t(Rho)
  p0 = datph2[ind.diag==0, incl.prob]
  
  Vij = Rhony * tcrossprod(1-p0) + diag((1-p0)*(p0)) 
  Pij = Rhony * tcrossprod(1-p0) + tcrossprod(p0) 
  diag(Pij) <- p0 
  
  return(Vij/Pij)
}


# standard error
se.comp = function(df, df0, se.comp.ls, nh=NULL){
  if(is.null(nh)) nh = se.comp.ls$nh
  method = se.comp.ls$method
  
  if(method=="ful") {
    se = sqdg(crossprod(df)*nh/(nh-1))
  } else if(method=="des") {
    omega = t(df[with(se.comp.ls, ind.undg[subset1]), ]) %*% se.comp.ls$jcov.byp %*% 
      df[with(se.comp.ls, ind.undg[subset1]), ]
    se = sqdg(crossprod(df/sqrt(se.comp.ls$wgts[se.comp.ls$subset1]))*nh/(nh-1) + omega)
  } else if(method %in% c("bre", "shi")) {
    omega = t(df0[se.comp.ls$subset1 & se.comp.ls$ind.undg,]) %*% se.comp.ls$jcov.byp %*% 
      df0[se.comp.ls$subset1 & se.comp.ls$ind.undg,]
    se = sqdg((crossprod(df0/sqrt(se.comp.ls$wgt0)) + 2*crossprod(df0, df - df0) + crossprod(df - df0))*nh/(nh-1) + omega)
  } 
  return(se)
}


# general functions
sqdg = function (x) {sqrt(diag(x))}
btn = function(t, intl, intr){intl<t & t<=intr}
thi = function(int1l, int1r, int2l, int2r){(!(int1r<int2l|int1l>int2r))*(pmin(int1r, int2r) - pmax(int1l, int2l))}


# CPP functions
library(Rcpp)
library(RcppArmadillo)
sourceCpp(paste0(workdir, 'dfbetaint.cpp'))
sourceCpp(paste0(workdir, 'sum2.cpp'))
sourceCpp(paste0(workdir, 'P_sd.cpp'))
sourceCpp(paste0(workdir, 'P_pb.cpp'))
sourceCpp(paste0(workdir, 'Rho.cpp'))