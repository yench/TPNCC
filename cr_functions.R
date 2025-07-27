# exponential competing risks data
simul.cr = function(n,               # cohort size 
                    mctrl=1,         # number of controls per case
                    max.diag.t=6,    # administrative censoring time of the NCC study
                    max.cohort.t=10, # administrative censoring time of cohort follow-up
                    max.enter=2,     # maximum late entry time into the cohort
                    Smax.censor=0.9, # baseline survival fraction from random censoring
                    Smax12=0.98, beta121=1.2, beta122=2,   # parameters for 1->2 transition   
                    Smax13=0.90, beta131=1.5, beta132=2.5, # parameters for 1->3 transition   
                    mu=c(z1=0, z2=0, u=0),                 # mean of Z1, Z2, and U
                    sig=c(sdz1=1, sdz2=1, sdu=1),          # marginal variance of Z1, Z2, and U
                    cov12=0, cov1u=0, cov2u=0.75){         # cov(Z1, Z2), cov(Z1, U), cov(Z2, U)
  
  Sig = matrix(0, nrow = 3, ncol = 3); diag(Sig) = sig^2
  Sig[upper.tri(Sig)] = c(cov12, cov1u, cov2u)
  Sig[lower.tri(Sig)] = c(cov12, cov1u, cov2u)
  Covmat = data.frame(MASS::mvrnorm(n=n, mu=mu, Sigma=Sig))
  
  # coefficients
  gamma12 = log(-log(Smax12)/max.cohort.t)
  gamma13 = log(-log(Smax13)/max.cohort.t)
  
  eta12 = with(Covmat, exp((gamma12 + beta121*z1 + beta122*z2)))
  eta13 = with(Covmat, exp((gamma13 + beta131*z1 + beta132*z2)))
  
  # survival data
  entrtime = runif(n, 0, max.enter)
  droptime = rexp(n, rate=-log(Smax.censor)/max.cohort.t) # drop out for any reason
  censtime = pmin(max.cohort.t - entrtime, droptime)      # end of study
  
  survtime = rexp(n, rate=eta12+eta13)
  causetype = 2 + rbinom(n, 1, prob=eta13/(eta12+eta13))
  
  diagtime = pmin(survtime, droptime, max.diag.t - entrtime)
  ind.diag = as.integer(causetype==2 & (survtime<pmin(droptime, max.diag.t - entrtime)))
  
  eventime = pmin(survtime, censtime)
  ind.fail = ifelse(survtime<=censtime, causetype, 0)

  DATA = data.table(cbind(ID = 1:n, Covmat, diagtime, ind.diag, eventime, ind.fail))
  
  # NCC sampling
  DATA[ , sample := .(ifelse(ind.diag==1, 
                             lapply(ID, 
                                    function(i){
                                      riskset = which(diagtime>=diagtime[i] & ID!=ID[i])
                                      out = NA
                                      if (length(riskset)>mctrl) out = sample(riskset, mctrl)
                                      else if (length(riskset)==mctrl) out = riskset
                                      return(out)
                                    }), NA))]
  DATA[ , nrisk := rank(-diagtime)] 
  
  # ID and indicator for inclusion in NCC
  ph2id = DATA[ind.diag==1, .(ID, sample)]
  
  for(h in 1:mctrl){
    ph2id[ , paste0("ctrl", h) := as.integer(lapply(.I, function(x){sample[[x]][h]}))]
  }
  setnames(ph2id, "ID", "case")
  ph2id[ , matched.id := paste0(1:.N)]
  ph2id[ , sample := NULL]
  
  # compute Design weights
  ph2dat = melt(ph2id, id.vars="matched.id", 
                measure.vars = c("case", paste0("ctrl", 1:mctrl)))
  ph2dat[ , `:=`(case = ifelse(variable=="case", 1, 0), ind.ncc = TRUE)]
  setnames(ph2dat, "value", "ID")
  nccwt = DATA[, .(ID, diagtime, ind.diag, nrisk, ind.ncc = FALSE)]
  nccwt[ph2dat, on=.(ID), ind.ncc := i.ind.ncc]
  nccwt[ , term := ifelse(ind.diag==1, 1- mctrl/(nrisk-1), 1)]
  setorder(nccwt, diagtime)
  nccwt[ , incl.prob := ifelse(ind.diag==1, 1, 1 - cumprod(term))]
  nccwt[ , `:=`(wgt = 1/incl.prob, diagtime = NULL, ind.diag = NULL, term = NULL)]
  
  DATA[ , sample := NULL]
  DATA[nccwt, `:=` (incl.prob = incl.prob, wgt = wgt, ind.ncc = ind.ncc), on = .(ID)]
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
  
  return(list(wgt = g, A = A)) # wgt = Breslow weight, A = matrix of auxiliary variables for the Breslow weight calibration procedure
}


# raking
rakecal = function(A,            # matrix of auxiliary variables among those in the NCC sample
                   w0,           # design weights among those in the NCC sample
                   total,        # sum of auxiliary variables in the full cohort
                   max_iter=500, # maximum number of iterations
                   EPS=1e-11){   # threshold for convergence
  
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


# fit standard models and predict transition probabilities using plug-in estimators
coxModel.sd = function (method,       # "ful" = full cohort,    "des" = design weight
                                      # "bre" = Breslow weight, "shi" = Shin weight
                        datph1,       # full cohort data
                        aux12,         # matrix of auxiliary variables for the 1->2 transition
                        aux13,         # matrix of auxiliary variables for the 1->3 transition
                        covars.final, # name of covariates included in the standard model
                        tau1, zstar, 
                        jcov.byp=NULL){ # matrix of joint inclusion probabilities
  subset = datph1[ , ind.ncc]
  wgts = datph1[ , wgt]
  ntau1 = length(tau1)
  
  form12 = as.formula(paste0("Surv(eventime, ind.12) ~ ", paste(covars.final, collapse="+")))
  form13 = as.formula(paste0("Surv(eventime, ind.13) ~ ", paste(covars.final, collapse="+")))
  
  if (method=="ful") { # full cohort
    fit12 = coxph(form12, data=datph1, timefix=FALSE)
    fit13 = coxph(form13, data=datph1, timefix=FALSE)
    df12 = df.coxph(fit12, dfeta=FALSE, subset=rep(TRUE, nrow(datph1)))
    df13 = df.coxph(fit13, dfeta=FALSE, subset=rep(TRUE, nrow(datph1)))
    tp = tp.sd.method(fit12$coefficients, df12, fit13$coefficients, df13, method, zstar, tau1)
    out = sd.out.method(fit12$coefficients, df12, NULL, fit13$coefficients, df13, NULL, tp, datph1, jcov.byp, method)
    
  } else if (method=="des") { # design weights
    fit12 = coxph(form12, data=datph1[subset, , drop=FALSE], weights=wgts[subset], timefix=FALSE)
    fit13 = coxph(form13, data=datph1[subset, , drop=FALSE], weights=wgts[subset], timefix=FALSE)
    df12 = df.coxph(fit12, dfeta=FALSE, subset=subset)
    df13 = df.coxph(fit13, dfeta=FALSE, subset=subset)
    tp = tp.sd.method(fit12$coefficients, df12, fit13$coefficients, df13, method, zstar, tau1)
    out = sd.out.method(fit12$coefficients, df12, NULL, fit13$coefficients, df13, NULL, tp, datph1, jcov.byp, method)
    
  } else if (method %in% c("bre", "shi")) { # calibrated weights
    fit12 = coxph(form12, data=datph1[subset, , drop=FALSE], weights=aux12$wgt, timefix=FALSE)
    fit13 = coxph(form13, data=datph1[subset, , drop=FALSE], weights=aux13$wgt, timefix=FALSE)
    df12  = df.coxph(fit12, dfeta=TRUE, fullA=aux12$A, subset=subset)
    df13  = df.coxph(fit13, dfeta=TRUE, fullA=aux13$A, subset=subset)
    df120 = df.coxph(fit12, dfeta=TRUE, fullA=aux12$A, subset=subset, omg=TRUE)
    df130 = df.coxph(fit13, dfeta=TRUE, fullA=aux13$A, subset=subset, omg=TRUE)
    tp  = tp.sd.method(fit12$coefficients, df12, fit13$coefficients, df13, method, zstar, tau1)
    tp0 = tp.sd.method(fit12$coefficients, df12, fit13$coefficients, df13, method, zstar, tau1)
    out = sd.out.method(fit12$coefficients, df12, df120, fit13$coefficients, df13, df130, tp, datph1, jcov.byp, method, tp0)
  } 
  
  if(method %in% c("ful", "bre")) {return(list(out = out, bint12 = df12$bint, bint13 = df13$bint, tp.lookup=tp$tp.lookup))
  } else {return(list(out = out, bint12 = df12$bint, bint13 = df13$bint))}
}


# fit proportional baselines models
coxModel.pb = function (method,         # "ful" = full cohort,    "des" = design weight
                                        # "bre" = Breslow weight, "shi" = Shin weight
                        datph1.repli,   # replicated full cohort data
                        aux1,            # matrix of auxiliary variables based on the proportional baselines model
                        covars.final.repli,  # name of covariates included in the proportional baselines model
                        tau1, zstar, 
                        tp.lookup,      # lookup table for indexing transition times within (0, tau1]
                        jcov.byp=NULL){ # matrix of joint inclusion probabilities
  subset = datph1.repli[ , ind.ncc]
  wgts = datph1.repli[ , wgt]
  ntau1 = length(tau1)
  
  form1 = as.formula(paste0("Surv(eventime, ind.1) ~ ", paste(covars.final.repli, collapse="+")))
  
  if (method=="ful") { # full cohort
    fit1 = coxph(form1, data=datph1.repli, timefix=FALSE)
    df1 = df.coxph(fit1, dfeta=FALSE, subset=rep(TRUE, nrow(datph1.repli)), pb=TRUE)
    tp = tp.pb.method(fit1$coefficients, df1, method, zstar, tau1, tp.lookup)
    out = pb.out.method(fit1$coefficients, df1, NULL, tp, datph1.repli, jcov.byp, method)
    
  } else if (method=="des") { # design weights
    fit1 = coxph(form1, data=datph1.repli[subset, , drop=FALSE], weights=wgts[subset], timefix=FALSE)
    df1 = df.coxph(fit1, dfeta=FALSE, subset=subset, pb=TRUE)
    tp = tp.pb.method(fit1$coefficients, df1, method, zstar, tau1, tp.lookup)
    out = pb.out.method(fit1$coefficients, df1, NULL, tp, datph1.repli, jcov.byp, method)
    
  } else if (method %in% c("bre", "shi")) { # calibrated weights
    fit1 = coxph(form1, data=datph1.repli[subset, , drop=FALSE], weights=aux1$wgt, timefix=FALSE)
    df1  = df.coxph(fit1, dfeta=TRUE, fullA=aux1$A, subset=subset, pb=TRUE)
    df10 = df.coxph(fit1, dfeta=TRUE, fullA=aux1$A, subset=subset, pb=TRUE, omg=TRUE)
    tp  = tp.pb.method(fit1$coefficients, df1, method, zstar, tau1, tp.lookup)
    tp0 = tp.pb.method(fit1$coefficients, df1, method, zstar, tau1, tp.lookup)
    out = pb.out.method(fit1$coefficients, df1, df10, tp, datph1.repli, jcov.byp, method, tp0)
  } 

  return(list(out = out, bint = df1$bint)) 
}


# predicted transition probabilities based on plug-in estimators
tp.sd.method = function(beta12, df12,  # beta and baseline intensity estimates and their influence functions for 1->2 transitions
                        beta13, df13,  # beta and baseline intensity estimates and their influence functions for 1->3 transitions
                        method, zstar, tau1, tp.lookup=NULL){ 

  bint12 = df12$bint
  bint13 = df13$bint
  
  S12 = S_1k(beta12, cumsum(bint12), zstar)
  S13 = S_1k(beta13, cumsum(bint13), zstar)
  
  haz.zst2 = haz.zst(beta12, bint12, zstar)
  haz.zst3 = haz.zst(beta13, bint13, zstar)
  
  df_P11.term.2 = df.P11.term.k(beta12, df12, zstar)
  df_P11.term.3 = df.P11.term.k(beta13, df13, zstar)
  
  d12 = nrow(S12); d13 = nrow(S13); d1 = d12 + d13
  nDF = nrow(df12$df_beta); nzstar = nrow(zstar)
  
  if(is.null(tp.lookup)){
    tp.lookup = data.table(cause = rep(2:3, c(d12, d13)),
                           tj = c(attr(bint12, "time"), attr(bint13, "time")),
                           tj.ind = c(1:d12, 1:d13), key="tj")
    tp.lookup[ , `:=` (tj.ind12 = cummax(ifelse(cause==2, tj.ind, 0)),  # carry forward the last index
                       tj.ind13 = cummax(ifelse(cause==3, tj.ind, 0)))]
  }
  
  cause2 = tp.lookup[ , cause==2]
  tj.ind12 = tp.lookup[ , tj.ind12];   tj.ind12.start = min(which(tj.ind12==1))
  tj.ind13 = tp.lookup[ , tj.ind13];   tj.ind13.start = min(which(tj.ind13==1))
  
  bint12.tj = bint13.tj = rep(0, d1)
  bint12.tj[ cause2] = bint12 
  bint13.tj[!cause2] = bint13
  bInt12.tj = cumsum(bint12.tj)
  bInt13.tj = cumsum(bint13.tj)
  bInt12.tj.lag = c(0, bInt12.tj[-d1])
  bInt13.tj.lag = c(0, bInt13.tj[-d1])
  
  tau1.loc  = sapply(tau1, function(x){tp.lookup[tj<=x, .N]})
  tau1.loc2 = sapply(tau1, function(x){tp.lookup[cause==2, ][tj<=x, .N]})
  tau1.loc3 = sapply(tau1, function(x){tp.lookup[cause==3, ][tj<=x, .N]})
  ntau1 = length(tau1)

  bInt1.out = cbind(bInt12.tj[tau1.loc], bInt13.tj[tau1.loc])
  dimnames(bInt1.out) = list(tau1, c("bInt12", "bInt13"))
  df_bInt1.out =  cbind(df12$df_bInt[ , tau1.loc2], df13$df_bInt[ , tau1.loc3])
  dimnames(df_bInt1.out) = list(NULL, paste0("bInt1", 2:3))
  tp.out.cols = paste0(rep(c("P11.zst", "P12.zst", "P13.zst"), nzstar), rep(1:nzstar, each=3))
  df.out.cols = paste0("zst", rep(1:nzstar))
  tp.out = matrix(NA, 1, 3*nzstar, dimnames=list(tau1, tp.out.cols))
  df_P11.out = df_P12.out = df_P13.out = matrix(NA, nDF, nzstar, dimnames=list(NULL, df.out.cols))
  
  for (i in 1:nzstar){
    S12.tj = S13.tj = rep(1, d1)
    S12.tj[tj.ind12.start:d1] = S12[tj.ind12, i]
    S13.tj[tj.ind13.start:d1] = S13[tj.ind13, i]
    
    df_P11.term.2.tj = df_P11.term.3.tj = matrix(0, nDF, d1)
    df_P11.term.2.tj[ , tj.ind12.start:d1] = df_P11.term.2[ , tj.ind12, i]
    df_P11.term.3.tj[ , tj.ind13.start:d1] = df_P11.term.3[ , tj.ind13, i]
    
    tp.i = matrix(NA, d1, 6, 
                  dimnames = list(rownames(tp.lookup), 
                                  c("P11.tj", "P11.tj.lag", "p12.tj", "p13.tj", "P12.tj", "P13.tj")))
    tp.i[ , "P11.tj"] = S12.tj*S13.tj
    tp.i[ , "P11.tj.lag"] = c(1, tp.i[-d1, "P11.tj"])
    tp.i[ , c("p12.tj", "p13.tj")] = 0
    tp.i[ cause2, "p12.tj"] = tp.i[ cause2, "P11.tj.lag"]*haz.zst2[ , i]
    tp.i[!cause2, "p13.tj"] = tp.i[!cause2, "P11.tj.lag"]*haz.zst3[ , i]
    tp.i[ , "P12.tj"] = cumsum(tp.i[ , "p12.tj"])
    tp.i[ , "P13.tj"] = cumsum(tp.i[ , "p13.tj"])
    
    df_P11.tj.i = df_P11(df_P11.term.2.tj + df_P11.term.3.tj, tp.i[ , "P11.tj"])
    df_P12.tj.i = df.P1k.tj.l(2, beta12, 2, beta12, df12, tp.i, zstar[i, ], bInt12.tj.lag, bint12.tj,  cause2) +
      df.P1k.tj.l(2, beta12, 3, beta13, df13, tp.i, zstar[i, ], bInt13.tj.lag, bint13.tj, !cause2) 
    df_P13.tj.i = df.P1k.tj.l(3, beta13, 2, beta12, df12, tp.i, zstar[i, ], bInt12.tj.lag, bint12.tj,  cause2) +
      df.P1k.tj.l(3, beta13, 3, beta13, df13, tp.i, zstar[i, ], bInt13.tj.lag, bint13.tj, !cause2)
    
    tp.col.loc = ((i-1)*3+1):(i*3)
    df.col.loc = i
    tp.out[ , tp.col.loc] = tp.i[tau1.loc, c("P11.tj", "P12.tj", "P13.tj")]
    df_P11.out[ , df.col.loc] = df_P11.tj.i[ , tau1.loc]
    df_P12.out[ , df.col.loc] = df_P12.tj.i[ , tau1.loc]
    df_P13.out[ , df.col.loc] = df_P13.tj.i[ , tau1.loc]
  }
  DF = cbind(df_bInt1.out, df_P11.out, df_P12.out, df_P13.out)
  colnames(DF)[-(1:ncol(df_bInt1.out))] = paste0(rep(c("P11", "P12", "P13"), each=ncol(df_P11.out)), 
                                                ".", colnames(DF)[-(1:ncol(df_bInt1.out))])
  if(method %in% c("ful", "bre")){return(list(est = cbind(bInt1.out, tp.out), DF = DF, tp.lookup = tp.lookup))
  } else {return(list(est = cbind(bInt1.out, tp.out), DF = DF))}
}


tp.pb.method = function(beta1, df1, # betas and baseline intensity estimates and their influence functions for transition out of state 1
                        method, zstar, tau1, tp.lookup){

  bint1 = df1$bint
  nDF = nrow(df1$df_beta); nzstar = nrow(zstar); d1 = length(bint1)
  
  cause2 = tp.lookup[ , cause==2]
  zstar2 = cbind(zstar, matrix(0, nzstar, ncol(zstar)+1))
  zstar3 = cbind(matrix(0, nzstar, ncol(zstar)), rep(1, nzstar), zstar)
  colnames(zstar2) = names(beta1)
  colnames(zstar3) = names(beta1)

  S12 = S_1k(beta1, cumsum(bint1), zstar2)
  S13 = S_1k(beta1, cumsum(bint1), zstar3)
  
  haz.zst2 = haz.zst(beta1, bint1, zstar2)  
  haz.zst3 = haz.zst(beta1, bint1, zstar3) 
  
  df_P11.term.2 = df.P11.term.k(beta1, df1, zstar2) 
  df_P11.term.3 = df.P11.term.k(beta1, df1, zstar3)
 
  bint12.tj = bint1    
  bint13.tj = bint1 * exp(beta1["z3"]) 
  bInt12.tj = cumsum(bint12.tj)
  bInt13.tj = cumsum(bint13.tj)
  bInt12.tj.lag = c(0, bInt12.tj[-d1])
  bInt13.tj.lag = c(0, bInt13.tj[-d1])
  
  tau1.loc  = sapply(tau1, function(x){tp.lookup[tj<=x, .N]})
  bInt1.out = cbind(bInt12.tj[tau1.loc], bInt13.tj[tau1.loc])
  dimnames(bInt1.out) = list(tau1, c("bInt12", "bInt13"))
  df_bInt1.out =  cbind(df1$df_bInt[ , tau1.loc], 
                       (df1$df_bInt[ , tau1.loc] + df1$df_beta[ , "z3"] %o% cumsum(bint1)[tau1.loc])*exp(beta1["z3"]))
  dimnames(df_bInt1.out) = list(NULL, paste0("bInt1", 2:3))
  tp.out.cols = paste0(rep(c("P11.zst", "P12.zst", "P13.zst"), nzstar), rep(1:nzstar, each=3))
  df.out.cols = paste0("zst", 1:nzstar)
  tp.out = matrix(NA, 1, 3*nzstar, dimnames=list(tau1, tp.out.cols))
  df_P11.out = df_P12.out = df_P13.out = matrix(NA, nDF, nzstar, dimnames=list(NULL, df.out.cols))
  
  for (i in 1:nzstar){
    S12.tj = S12[ , i]
    S13.tj = S13[ , i]
    df_P11.term.2.tj = df_P11.term.2[ , , i]
    df_P11.term.3.tj = df_P11.term.3[ , , i]
    tp.i = matrix(NA, d1, 6, dimnames=list(rownames(tp.lookup), 
                                          c("P11.tj", "P11.tj.lag", "p12.tj", "p13.tj", "P12.tj", "P13.tj")))
    tp.i[ , "P11.tj"] = S12.tj*S13.tj
    tp.i[ , "P11.tj.lag"] = c(1, tp.i[-d1, "P11.tj"])
    tp.i[ , c("p12.tj", "p13.tj")] = 0
    tp.i[ , "p12.tj"] = tp.i[ , "P11.tj.lag"]*haz.zst2[ , i]
    tp.i[ , "p13.tj"] = tp.i[ , "P11.tj.lag"]*haz.zst3[ , i]
    tp.i[ , "P12.tj"] = cumsum(tp.i[ , "p12.tj"])
    tp.i[ , "P13.tj"] = cumsum(tp.i[ , "p13.tj"])
    df_P11.tj.i = df_P11(df_P11.term.2.tj + df_P11.term.3.tj, tp.i[ , "P11.tj"])
    df_P12.tj.i = df.P1k.tj.l(2, beta1, 2, beta1, df1, tp.i, zstar2[i, ], bInt12.tj.lag, bint12.tj, rep(TRUE, d1)) +
                  df.P1k.tj.l(2, beta1, 3, beta1, df1, tp.i, zstar3[i, ], rep(0, d1),    bint13.tj, rep(TRUE, d1)) 
    df_P13.tj.i = df.P1k.tj.l(3, beta1, 2, beta1, df1, tp.i, zstar2[i, ], bInt12.tj.lag, bint12.tj, rep(TRUE, d1)) +
                  df.P1k.tj.l(3, beta1, 3, beta1, df1, tp.i, zstar3[i, ], rep(0, d1),    bint13.tj, rep(TRUE, d1))
    
    tp.col.loc = ((i-1)*3+1):(i*3)
    df.col.loc = i
    tp.out[ , tp.col.loc] = tp.i[tau1.loc, c("P11.tj", "P12.tj", "P13.tj")]
    df_P11.out[ , df.col.loc] = df_P11.tj.i[ , tau1.loc]
    df_P12.out[ , df.col.loc] = df_P12.tj.i[ , tau1.loc]
    df_P13.out[ , df.col.loc] = df_P13.tj.i[ , tau1.loc]
  }
  DF = cbind(df_bInt1.out, df_P11.out, df_P12.out, df_P13.out)
  colnames(DF)[-(1:ncol(df_bInt1.out))] = paste0(rep(c("P11", "P12", "P13"), each=ncol(df_P11.out)), 
                                                ".", colnames(DF)[-(1:ncol(df_bInt1.out))])
  return(list(est = cbind(bInt1.out, tp.out), DF = DF))
}


# cause-specific survival at zstar given estimates of betas and the cumulative baseline intensity
S_1k = function(beta, bInt, zstar){
  apply(zstar, 1, function(x){exp(-bInt*as.vector(exp(beta %*% x)))})
}


# hazard given zstar and estimates of betas and the baseline intensity
haz.zst = function(beta, bint, zstar){
  apply(zstar, 1, function(x){bint*as.vector(exp(beta %*% x))})
}


# compute increments in the baseline intensity and its influence functions and influence functions of betas
df.coxph = function(coxph_object,  # a coxph object from fitting the standard or the proportional baselines model
                    dfeta=TRUE,    # compute IF(eta)? TRUE for calibrated weights
                    fullA=NULL,    # NULL for full cohort and design weights, the matrix of auxiliary variables for calibrated weights
                    subset=NULL,   # a vector of TRUE/FALSE indicators for inclusion in the analysis
                    omg=FALSE,     # FALSE to compute IF_i, TRUE to compute IF_2i
                    pb=FALSE){     # TRUE if the coxph object comes from fitting a proportional baselines model
  
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
  bint = as.vector(ret$bint)
  attr(bint, "time") = trantime
  
  return(list(bint = bint,             # estimates of increments in cumulative baseline intensity
              df_beta = ret$df_beta,   # influence functions of beta estimates
              df_bint = ret$df_lambda, # influence functions of estimates of increments in cumulative baseline intensity
              df_bInt = rowcumsum(ret$df_lambda))) # influence functions of cumulative baseline intensity estimates
}


# functions for computing parts of influence functions of predicted transition probabilities
df.P11.term.k = function(beta, df, zstar){
  term1 = aperm((df$df_beta %*% t(zstar * as.vector(-exp(zstar %*% beta)))) %o% cumsum(df$bint), c(1, 3, 2))
  term2 = sapply(-exp(zstar %*% beta), FUN=function(x){x*df$df_bInt}, simplify="array")
  return(term1 + term2)
}


df.P1k.tj.l = function(k, betak, l, betal, dfl, tp.i, zstar.i, bInt1l.tj.lag, bint1l.tj, causel){
  P1kobj = paste0('P1', k, '.tj')
  p1kobj = paste0('p1', k, '.tj')

  df_bintl = dfl$df_bint
  pd.P1k.tj.betal = tp.i[ , P1kobj] * as.integer(k==l) - as.vector(exp(betal %*% zstar.i)) * cumsum(tp.i[ , p1kobj]*bInt1l.tj.lag)
  
  d = nrow(tp.i)
  ind.mat = matrix(0, d, d)
  ind.mat[upper.tri(ind.mat, diag=TRUE)] = 1
  if(k!=l){pd.P1k.tj.dLaml0.term1 = matrix(0, sum(causel), d)
  } else {pd.P1k.tj.dLaml0.term1 = (matrix(tp.i[ , "P11.tj.lag"], d, d) * ind.mat)[causel, ] * as.vector(exp(betak %*% zstar.i))}
  diag(ind.mat) = 0
  pd.P1k.tj.dLaml0.term2 = pdP1kdL2(ind.mat, tp.i[ , p1kobj])[causel, ] * as.vector(exp(betal %*% zstar.i))
  pd.P1k.tj.dLaml0 = pd.P1k.tj.dLaml0.term1 - pd.P1k.tj.dLaml0.term2
  return(dfl$df_beta %*% (zstar.i %o% pd.P1k.tj.betal) + df_bintl %*% pd.P1k.tj.dLaml0)
}


# compute joint inclusion probabilities
jointVP = function(datph2, m){
  
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


# organize output and compute standard error
sd.out.method = function(beta12, df12, df120=NULL, 
                         beta13, df13, df130=NULL, 
                         tp_obj,        # output of tp.sd.method() with omg=FALSE
                         datph1, jcov.byp, method, 
                         tp0_obj=NULL){ # output of tp.sd.method() with omg=TRUE
  tau1 = as.numeric(rownames(tp_obj$est))

  subset = datph1[ , ind.ncc]
  wgts = datph1[ , wgt]
  nfull = nrow(datph1)
  wgt0 = rep(1, nfull) 
  wgt0[subset] = wgts[subset]
  
  if(method %in% c("bre", "shi")){
    df7  = cbind(df12$df_beta,  df13$df_beta,  tp_obj$DF)
    df70 = cbind(df120$df_beta, df130$df_beta, tp0_obj$DF)
    ind.diag = datph1[ , ind.diag]
    omega = t(df70[subset & ind.diag==0,]) %*% jcov.byp %*% df70[subset & ind.diag==0,]
    se = sqdg((crossprod(df70/sqrt(wgt0)) + 2*crossprod(df70, df7 - df70) + crossprod(df7 - df70))*nfull/(nfull-1) + omega)
    
  } else {
    df7 = cbind(df12$df_beta, df13$df_beta, tp_obj$DF)
    df70 = NULL
    if(method=="ful") {se = sqdg(crossprod(df7)*nfull/(nfull-1))}
    if(method=="des") {
      ind.diag = datph1[subset, ind.diag]
      omega = t(df7[ind.diag==0,]) %*% jcov.byp %*% df7[ind.diag==0,]
      se = sqdg(crossprod(df7/sqrt(wgts[subset]))*nfull/(nfull-1) + omega)}
  }
  se.cols = unique(gsub("[.]\\d+", "", colnames(df7)[-(1:(length(beta12) + length(beta13)))]))
  est = c(beta12, beta13, c(tp_obj$est[ , se.cols])) 
  out = rbind(est, se)
  colnames(out) = names(se)
  colnames(out)[1:(length(beta12) + length(beta13))] = c(paste0("beta12", 1:2), paste0("beta13", 1:2))
  out = data.table(out, method = method, type = c("est", "se"), model = "sd")
  return(out)
}


pb.out.method = function(beta1, df, df0=NULL, tp_obj, datph1.repli, jcov.byp, method, tp0_obj=NULL){
  tau1 = as.numeric(rownames(tp_obj$est))

  subset = datph1.repli[cause==2, ind.ncc]
  wgts = datph1.repli[cause==2 , wgt]
  nfull = length(subset)
  wgt0 = rep(1, nfull) 
  wgt0[subset] = wgts[subset]
  
  if(method %in% c("bre", "shi")){
    df7  = cbind(df$df_beta,  tp_obj$DF)
    df70 = cbind(df0$df_beta, tp0_obj$DF)
    ind.diag = datph1.repli[cause==2, ind.diag]
    omega = t(df70[subset & ind.diag==0,]) %*% jcov.byp %*% df70[subset & ind.diag==0,]
    se = sqdg((crossprod(df70/sqrt(wgt0)) + 2*crossprod(df70, df7 - df70) + crossprod(df7 - df70))*nfull/(nfull-1) + omega)

  } else {
    df7 = cbind(df$df_beta, tp_obj$DF)
    df70 = NULL
    if(method=="ful") {se = sqdg(crossprod(df7)*nfull/(nfull-1))}
    if(method=="des") {
      ind.diag = datph1.repli[cause==2 & ind.ncc==TRUE, ind.diag]
      omega = t(df7[ind.diag==0,]) %*% jcov.byp %*% df7[ind.diag==0,]
      se = sqdg(crossprod(df7/sqrt(wgts[subset]))*nfull/(nfull-1) + omega)}
  }
  se.cols = unique(gsub("[.]\\d+", "", colnames(df7)[-(1:length(beta1))]))
  beta_len = length(beta1)
  beta_col = c(paste0("beta12", 1:2), "gamma", paste0("beta13", 1:2))
  est = c(beta1, c(tp_obj$est[ , se.cols])) 
  out = rbind(est, se)
  colnames(out) = names(se)
  colnames(out)[1:beta_len] = beta_col
  out = data.table(out, method = method, type = c("est", "se"), model = "pb")
  return(out)
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
sourceCpp(paste0(workdir, 'rowcumsum.cpp'))
sourceCpp(paste0(workdir, 'df_P11.cpp'))
sourceCpp(paste0(workdir, 'pdP1kdL2.cpp'))
sourceCpp(paste0(workdir, 'Rho.cpp'))