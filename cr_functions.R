# exponential competing risks data
simul.cr = function(n, mctrl=1, 
                    max.diag.t=6, max.cohort.t=10,         # administrative censoring times
                    max.enter=2, Smax.censor=0.9, 
                    Smax12=0.98, beta121=1.2, beta122=2,   # hazard 0->2   
                    Smax13=0.90, beta131=1.5, beta132=2.5, # hazard 0->3   
                    mu=c(z1=0, z2=0, u=0),
                    sig=c(sdz1=1, sdz2=1, sdu=1),
                    cov12=0, cov1u=0, cov2u=0.75){      
  
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
calib_bre = function (datph1, trantime, status.var, covars){
  subset = datph1[ , ind.ncc]
  wgts = datph1[ , wgt]
  
  form = as.formula(paste0("Surv(", trantime, ",", status.var, ") ~ ", paste(covars, collapse="+")))
  fit = coxph(form, data=datph1)
  DFBETA = residuals(fit, type="dfbeta")
  A = cbind(1, DFBETA)
  g = rakecal(A[subset, ], wgts[subset], colSums(A))
  if(is.null(g)) {print(paste("Brewgt", substr(status.var, 5, 6), "did not converge")); return(NULL)}
  
  return(list(wgt = g, A = A))
}


# raking
rakecal = function(A, w0, total, max_iter=500, EPS=1e-11, eta=TRUE){
  
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
  
  if(eta== T) attributes(g) <- list(lambda = lambda1)
  
  return(g)
} 


# fit standard models
coxModel.sd = function (method, datph1, aux2, aux3, covars.final, tau1, zstar, jcov.byp=NULL){
  subset = datph1[ , ind.ncc]
  wgts = datph1[ , wgt]
  ntau1 = length(tau1)
  
  form2 = as.formula(paste0("Surv(eventime, ind.12) ~ ", paste(covars.final, collapse="+")))
  form3 = as.formula(paste0("Surv(eventime, ind.13) ~ ", paste(covars.final, collapse="+")))
  
  if (method=="ful") { # full cohort
    fit2 = coxph(form2, data=datph1, timefix=FALSE)
    fit3 = coxph(form3, data=datph1, timefix=FALSE)
    df2 = df.coxph(fit2, type="bint", subset=rep(TRUE, nrow(datph1)))
    df3 = df.coxph(fit3, type="bint", subset=rep(TRUE, nrow(datph1)))
    tp = tp.sd.method(fit2$coefficients, df2, fit3$coefficients, df3, method, zstar, tau1)
    out = sd.out.method(fit2$coefficients, df2, NULL, fit3$coefficients, df3, NULL, tp, datph1, jcov.byp, method)
    
  } else if (method=="des") { # Design
    fit2 = coxph(form2, data=datph1[subset, , drop=FALSE], weights=wgts[subset], timefix=FALSE)
    fit3 = coxph(form3, data=datph1[subset, , drop=FALSE], weights=wgts[subset], timefix=FALSE)
    df2 = df.coxph(fit2, type="bint", subset=subset)
    df3 = df.coxph(fit3, type="bint", subset=subset)
    tp = tp.sd.method(fit2$coefficients, df2, fit3$coefficients, df3, method, zstar, tau1)
    out = sd.out.method(fit2$coefficients, df2, NULL, fit3$coefficients, df3, NULL, tp, datph1, jcov.byp, method)
    
  } else if (method %in% c("bre", "shi")) { # Calibrated
    fit2 = coxph(form2, data=datph1[subset, , drop=FALSE], weights=aux2$wgt, timefix=FALSE)
    fit3 = coxph(form3, data=datph1[subset, , drop=FALSE], weights=aux3$wgt, timefix=FALSE)
    df2  = df.coxph(fit2, type=c("eta", "bint"), fullA=aux2$A, subset=subset)
    df3  = df.coxph(fit3, type=c("eta", "bint"), fullA=aux3$A, subset=subset)
    df20 = df.coxph(fit2, type=c("eta", "bint"), fullA=aux2$A, subset=subset, omg=TRUE)
    df30 = df.coxph(fit3, type=c("eta", "bint"), fullA=aux3$A, subset=subset, omg=TRUE)
    tp  = tp.sd.method(fit2$coefficients, df2, fit3$coefficients, df3, method, zstar, tau1)
    tp0 = tp.sd.method(fit2$coefficients, df2, fit3$coefficients, df3, method, zstar, tau1)
    out = sd.out.method(fit2$coefficients, df2, df20, fit3$coefficients, df3, df30, tp, datph1, jcov.byp, method, tp0)
  } 

  zph2 = cox.zph(fit2, terms=FALSE)
  zph3 = cox.zph(fit3, terms=FALSE)
  
  if(method %in% c("ful", "bre")) {return(list(out = out, bint2 = df2$bint, bint3 = df3$bint, zph2 = zph2, zph3 = zph3, tp.lookup=tp$tp.lookup))
  } else {return(list(out = out, bint2 = df2$bint, bint3 = df3$bint, zph2 = zph2, zph3 = zph3))}
}


# fit proportional baselines models
coxModel.pb = function (method, datph1.repli, aux, covars.final.repli, tau1, zstar, tp.lookup, jcov.byp=NULL){
  subset = datph1.repli[ , ind.ncc]
  wgts = datph1.repli[ , wgt]
  ntau1 = length(tau1)
  
  form = as.formula(paste0("Surv(eventime, ind.1) ~ ", paste(covars.final.repli, collapse="+")))
  
  if (method=="ful") { # full cohort
    fit = coxph(form, data=datph1.repli, timefix=FALSE)
    df = df.coxph(fit, type="bint", subset=rep(TRUE, nrow(datph1.repli)), pb=TRUE)
    tp = tp.pb.method(fit$coefficients, df, method, zstar, tau1, tp.lookup)
    out = pb.out.method(fit$coefficients, df, NULL, tp, datph1.repli, jcov.byp, method)
    
  } else if (method=="des") { # Design
    fit = coxph(form, data=datph1.repli[subset, , drop=FALSE], weights=wgts[subset], timefix=FALSE)
    df = df.coxph(fit, type="bint", subset=subset, pb=TRUE)
    tp = tp.pb.method(fit$coefficients, df, method, zstar, tau1, tp.lookup)
    out = pb.out.method(fit$coefficients, df, NULL, tp, datph1.repli, jcov.byp, method)
    
  } else if (method %in% c("bre", "shi")) { # Calibrated
    fit = coxph(form, data=datph1.repli[subset, , drop=FALSE], weights=aux$wgt, timefix=FALSE)
    df  = df.coxph(fit, type=c("eta", "bint"), fullA=aux$A, subset=subset, pb=TRUE)
    df0 = df.coxph(fit, type=c("eta", "bint"), fullA=aux$A, subset=subset, pb=TRUE, omg=TRUE)
    tp  = tp.pb.method(fit$coefficients, df, method, zstar, tau1, tp.lookup)
    tp0 = tp.pb.method(fit$coefficients, df, method, zstar, tau1, tp.lookup)
    out = pb.out.method(fit$coefficients, df, df0, tp, datph1.repli, jcov.byp, method, tp0)
  } 
  zph = cox.zph(fit, terms=FALSE)

  return(list(out = out, bint = df$bint, zph = zph))
}


# plug-in predicted transition probabilities
tp.sd.method = function(est2, df2, est3, df3, method, zstar, tau1, tp.lookup=NULL){

  bint2 = df2$bint
  bint3 = df3$bint
  
  S2 = S_k(est2, cumsum(bint2), zstar)
  S3 = S_k(est3, cumsum(bint3), zstar)
  
  haz.zst2 = haz.zst(est2, bint2, zstar)
  haz.zst3 = haz.zst(est3, bint3, zstar)
  
  df_P11.term.2 = df.P11.term.k(est2, df2, zstar)
  df_P11.term.3 = df.P11.term.k(est3, df3, zstar)
  
  d2 = nrow(S2); d3 = nrow(S3); d = d2 + d3
  nDF = nrow(df2$df_beta); nzstar = nrow(zstar)
  
  if(is.null(tp.lookup)){
    tp.lookup = data.table(cause = rep(2:3, c(d2, d3)),
                           tj = c(attr(bint2, "time"), attr(bint3, "time")),
                           tj.ind = c(1:d2, 1:d3), key="tj")
    tp.lookup[ , `:=` (tj.ind2 = cummax(ifelse(cause==2, tj.ind, 0)),  # carry forward the last index
                       tj.ind3 = cummax(ifelse(cause==3, tj.ind, 0)))]
  }
  
  cause2 = tp.lookup[ , cause==2]
  tj.ind2 = tp.lookup[ , tj.ind2];   tj.ind2.start = min(which(tj.ind2==1))
  tj.ind3 = tp.lookup[ , tj.ind3];   tj.ind3.start = min(which(tj.ind3==1))
  
  bint2.tj = bint3.tj = rep(0, d)
  bint2.tj[ cause2] = bint2 
  bint3.tj[!cause2] = bint3
  bInt2.tj = cumsum(bint2.tj)
  bInt3.tj = cumsum(bint3.tj)
  bInt2.tj.lag = c(0, bInt2.tj[-d])
  bInt3.tj.lag = c(0, bInt3.tj[-d])
  
  tau1.loc  = sapply(tau1, function(x){tp.lookup[tj<=x, .N]})
  tau1.loc2 = sapply(tau1, function(x){tp.lookup[cause==2, ][tj<=x, .N]})
  tau1.loc3 = sapply(tau1, function(x){tp.lookup[cause==3, ][tj<=x, .N]})
  ntau1 = length(tau1)

  bInt.out = cbind(bInt2.tj[tau1.loc], bInt3.tj[tau1.loc])
  dimnames(bInt.out) = list(tau1, c("bInt2", "bInt3"))
  df_bInt.out =  cbind(df2$df_bInt[ , tau1.loc2], df3$df_bInt[ , tau1.loc3])
  dimnames(df_bInt.out) = list(NULL, paste0("bInt", rep(2:3, each=ntau1),".", rep(tau1, 2)))
  tp.out.cols = paste0(rep(c("P11.zst", "P12.zst", "P13.zst"), nzstar), rep(1:nzstar, each=3))
  df.out.cols = paste0("zst", rep(1:nzstar, each=ntau1), ".", rep(tau1, nzstar))
  tp.out = matrix(NA, ntau1, 3*nzstar, dimnames=list(tau1, tp.out.cols))
  df_P11.out = df_P12.out = df_P13.out = matrix(NA, nDF, ntau1*nzstar, dimnames=list(NULL, df.out.cols))
  
  for (i in 1:nzstar){
    S2.tj = S3.tj = rep(1, d)
    S2.tj[tj.ind2.start:d] = S2[tj.ind2, i]
    S3.tj[tj.ind3.start:d] = S3[tj.ind3, i]
    
    df_P11.term.2.tj = df_P11.term.3.tj = matrix(0, nDF, d)
    df_P11.term.2.tj[ , tj.ind2.start:d] = df_P11.term.2[ , tj.ind2, i]
    df_P11.term.3.tj[ , tj.ind3.start:d] = df_P11.term.3[ , tj.ind3, i]
    
    tp.i = matrix(NA, d, 6, 
                  dimnames = list(rownames(tp.lookup), 
                                  c("P11.tj", "P11.tj.lag", "p12.tj", "p13.tj", "P12.tj", "P13.tj")))
    tp.i[ , "P11.tj"] = S2.tj*S3.tj
    tp.i[ , "P11.tj.lag"] = c(1, tp.i[-d, "P11.tj"])
    tp.i[ , c("p12.tj", "p13.tj")] = 0
    tp.i[ cause2, "p12.tj"] = tp.i[ cause2, "P11.tj.lag"]*haz.zst2[ , i]
    tp.i[!cause2, "p13.tj"] = tp.i[!cause2, "P11.tj.lag"]*haz.zst3[ , i]
    tp.i[ , "P12.tj"] = cumsum(tp.i[ , "p12.tj"])
    tp.i[ , "P13.tj"] = cumsum(tp.i[ , "p13.tj"])
    
    df_P11.tj.i = df_P11(df_P11.term.2.tj + df_P11.term.3.tj, tp.i[ , "P11.tj"])
    df_P12.tj.i = df.P1k.tj.l(2, est2, 2, est2, df2, tp.i, zstar[i, ], bInt2.tj.lag, bint2.tj,  cause2) +
      df.P1k.tj.l(2, est2, 3, est3, df3, tp.i, zstar[i, ], bInt3.tj.lag, bint3.tj, !cause2) 
    df_P13.tj.i = df.P1k.tj.l(3, est3, 2, est2, df2, tp.i, zstar[i, ], bInt2.tj.lag, bint2.tj,  cause2) +
      df.P1k.tj.l(3, est3, 3, est3, df3, tp.i, zstar[i, ], bInt3.tj.lag, bint3.tj, !cause2)
    
    tp.col.loc = ((i-1)*3+1):(i*3)
    df.col.loc = ((i-1)*ntau1+1):(i*ntau1)
    tp.out[ , tp.col.loc] = tp.i[tau1.loc, c("P11.tj", "P12.tj", "P13.tj")]
    df_P11.out[ , df.col.loc] = df_P11.tj.i[ , tau1.loc]
    df_P12.out[ , df.col.loc] = df_P12.tj.i[ , tau1.loc]
    df_P13.out[ , df.col.loc] = df_P13.tj.i[ , tau1.loc]
  }
  DF = cbind(df_bInt.out, df_P11.out, df_P12.out, df_P13.out)
  colnames(DF)[-(1:ncol(df_bInt.out))] = paste0(rep(c("P11", "P12", "P13"), each=ncol(df_P11.out)), 
                                                ".", colnames(DF)[-(1:ncol(df_bInt.out))])
  if(method %in% c("ful", "bre")){return(list(est = cbind(bInt.out, tp.out), DF = DF, tp.lookup = tp.lookup))
  } else {return(list(est = cbind(bInt.out, tp.out), DF = DF))}
}


tp.pb.method = function(est, df, method, zstar, tau1, tp.lookup){

  bint = df$bint
  nDF = nrow(df$df_beta); nzstar = nrow(zstar); d = length(bint)
  
  cause2 = tp.lookup[ , cause==2]
  zstar2 = cbind(zstar, matrix(0, nzstar, ncol(zstar)+1))
  zstar3 = cbind(matrix(0, nzstar, ncol(zstar)), rep(1, nzstar), zstar)
  colnames(zstar2) = names(est)
  colnames(zstar3) = names(est)

  S2 = S_k(est, cumsum(bint), zstar2)
  S3 = S_k(est, cumsum(bint), zstar3)
  
  haz.zst2 = haz.zst(est, bint, zstar2)  
  haz.zst3 = haz.zst(est, bint, zstar3) 
  
  df_P11.term.2 = df.P11.term.k(est, df, zstar2) 
  df_P11.term.3 = df.P11.term.k(est, df, zstar3)
 
  bint2.tj = bint    
  bint3.tj = bint * exp(est["z3"]) 
  bInt2.tj = cumsum(bint2.tj)
  bInt3.tj = cumsum(bint3.tj)
  bInt2.tj.lag = c(0, bInt2.tj[-d])
  bInt3.tj.lag = c(0, bInt3.tj[-d])
  
  tau1.loc  = sapply(tau1, function(x){tp.lookup[tj<=x, .N]})
  ntau1 = length(tau1)
  bInt.out = cbind(bInt2.tj[tau1.loc], bInt3.tj[tau1.loc])
  dimnames(bInt.out) = list(tau1, c("bInt2", "bInt3"))
  df_bInt.out =  cbind(df$df_bInt[ , tau1.loc], 
                       (df$df_bInt[ , tau1.loc] + df$df_beta[ , "z3"] %o% cumsum(bint)[tau1.loc])*exp(est["z3"]))
  dimnames(df_bInt.out) = list(NULL, paste0("bInt", rep(2:3, each=ntau1),".", rep(tau1, 2)))
  tp.out.cols = paste0(rep(c("P11.zst", "P12.zst", "P13.zst"), nzstar), rep(1:nzstar, each=3))
  df.out.cols = paste0("zst", rep(1:nzstar, each=ntau1), ".", rep(tau1, nzstar))
  tp.out = matrix(NA, ntau1, 3*nzstar, dimnames=list(tau1, tp.out.cols))
  df_P11.out = df_P12.out = df_P13.out = matrix(NA, nDF, ntau1*nzstar, dimnames=list(NULL, df.out.cols))
  
  for (i in 1:nzstar){
    S2.tj = S2[ , i]
    S3.tj = S3[ , i]
    df_P11.term.2.tj = df_P11.term.2[ , , i]
    df_P11.term.3.tj = df_P11.term.3[ , , i]
    tp.i = matrix(NA, d, 6, dimnames=list(rownames(tp.lookup), 
                                          c("P11.tj", "P11.tj.lag", "p12.tj", "p13.tj", "P12.tj", "P13.tj")))
    tp.i[ , "P11.tj"] = S2.tj*S3.tj
    tp.i[ , "P11.tj.lag"] = c(1, tp.i[-d, "P11.tj"])
    tp.i[ , c("p12.tj", "p13.tj")] = 0
    tp.i[ , "p12.tj"] = tp.i[ , "P11.tj.lag"]*haz.zst2[ , i]
    tp.i[ , "p13.tj"] = tp.i[ , "P11.tj.lag"]*haz.zst3[ , i]
    tp.i[ , "P12.tj"] = cumsum(tp.i[ , "p12.tj"])
    tp.i[ , "P13.tj"] = cumsum(tp.i[ , "p13.tj"])
    # df_P11.tj.i = sweep(df_P11.term.2.tj + df_P11.term.3.tj, MARGIN=2, tp.i[ , "P11.tj"], '*')
    df_P11.tj.i = df_P11(df_P11.term.2.tj + df_P11.term.3.tj, tp.i[ , "P11.tj"])
    df_P12.tj.i = df.P1k.tj.l(2, est, 2, est, df, tp.i, zstar2[i, ], bInt2.tj.lag, bint2.tj, rep(TRUE, d)) +
                  df.P1k.tj.l(2, est, 3, est, df, tp.i, zstar3[i, ], rep(0, d),    bint3.tj, rep(TRUE, d)) 
    df_P13.tj.i = df.P1k.tj.l(3, est, 2, est, df, tp.i, zstar2[i, ], bInt2.tj.lag, bint2.tj, rep(TRUE, d)) +
                  df.P1k.tj.l(3, est, 3, est, df, tp.i, zstar3[i, ], rep(0, d),    bint3.tj, rep(TRUE, d))
    
    tp.col.loc = ((i-1)*3+1):(i*3)
    df.col.loc = ((i-1)*ntau1+1):(i*ntau1)
    tp.out[ , tp.col.loc] = tp.i[tau1.loc, c("P11.tj", "P12.tj", "P13.tj")]
    df_P11.out[ , df.col.loc] = df_P11.tj.i[ , tau1.loc]
    df_P12.out[ , df.col.loc] = df_P12.tj.i[ , tau1.loc]
    df_P13.out[ , df.col.loc] = df_P13.tj.i[ , tau1.loc]
  }
  DF = cbind(df_bInt.out, df_P11.out, df_P12.out, df_P13.out)
  colnames(DF)[-(1:ncol(df_bInt.out))] = paste0(rep(c("P11", "P12", "P13"), each=ncol(df_P11.out)), 
                                                ".", colnames(DF)[-(1:ncol(df_bInt.out))])
  return(list(est = cbind(bInt.out, tp.out), DF = DF))
}


# cause-specific survival (at zstar)
S_k = function(bcoef, bInt, zstar){
  apply(zstar, 1, function(x){exp(-bInt*as.vector(exp(bcoef %*% x)))})
}


# hazard given beta, baseline estimates, and zstar
haz.zst = function(bcoef, bint, zstar){
  apply(zstar, 1, function(x){bint*as.vector(exp(bcoef %*% x))})
}


# influence functions of log-relative hazards and cumulative baseline hazards
df.coxph = function(coxph_object, type=c("eta", "bint"), fullA=NULL, subset=NULL, omg=FALSE, pb=FALSE){
  
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
                  wgt, nevent.wt, which(subset)-1, trantime.order, nsub, ntime, nbeta, "eta" %in% type, "bint" %in% type, omg)
  
  if(pb==TRUE){
    ret$df_beta = sum2(ret$df_beta)
    ret$df_lambda = sum2(ret$df_lambda)
  }
  
  colnames(ret$df_beta) = names(bcoef)
  bint = as.vector(ret$bint)
  attr(bint, "time") = trantime
  
  return(list(bint = bint, df_beta = ret$df_beta, df_bint = ret$df_lambda, df_bInt = rowcumsum(ret$df_lambda)))
}


# influence functions of transition probabilities
df.P11.term.k = function(est, df, zstar){
  term1 = aperm((df$df_beta %*% t(zstar * as.vector(-exp(zstar %*% est)))) %o% cumsum(df$bint), c(1, 3, 2))
  term2 = sapply(-exp(zstar %*% est), FUN=function(x){x*df$df_bInt}, simplify="array")
  return(term1 + term2)
}


df.P1k.tj.l = function(k, estk, l, estl, dfl, tp.i, zstar.i, bIntl.tj.lag, bintl.tj, causel){
  P1kobj = paste0('P1', k, '.tj')
  p1kobj = paste0('p1', k, '.tj')

  df_bintl = dfl$df_bint
  pd.P1k.tj.betal = tp.i[ , P1kobj] * as.integer(k==l) - as.vector(exp(estl %*% zstar.i)) * cumsum(tp.i[ , p1kobj]*bIntl.tj.lag)
  
  d = nrow(tp.i)
  ind.mat = matrix(0, d, d)
  ind.mat[upper.tri(ind.mat, diag=TRUE)] = 1
  if(k!=l){pd.P1k.tj.dLaml0.term1 = matrix(0, sum(causel), d)
  } else {pd.P1k.tj.dLaml0.term1 = (matrix(tp.i[ , "P11.tj.lag"], d, d) * ind.mat)[causel, ] * as.vector(exp(estk %*% zstar.i))}
  diag(ind.mat) = 0
  pd.P1k.tj.dLaml0.term2 = pdP1kdL2(ind.mat, tp.i[ , p1kobj])[causel, ] * as.vector(exp(estl %*% zstar.i))
  pd.P1k.tj.dLaml0 = pd.P1k.tj.dLaml0.term1 - pd.P1k.tj.dLaml0.term2
  return(dfl$df_beta %*% (zstar.i %o% pd.P1k.tj.betal) + df_bintl %*% pd.P1k.tj.dLaml0)
}


# joint inclusion probabilities
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


# standard error
sd.out.method = function(est2, df2, df20=NULL, est3, df3, df30=NULL, tp_obj, datph1, jcov.byp, method, tp0_obj = NULL){
  tau1 = as.numeric(rownames(tp_obj$est))
  ntau1 = length(tau1)

  subset = datph1[ , ind.ncc]
  wgts = datph1[ , wgt]
  nfull = nrow(datph1)
  wgt0 = rep(1, nfull) 
  wgt0[subset] = wgts[subset]
  
  if(method %in% c("bre", "shi")){
    df7  = cbind(df2$df_beta,  df3$df_beta,  tp_obj$DF)
    df70 = cbind(df20$df_beta, df30$df_beta, tp0_obj$DF)
    ind.diag = datph1[ , ind.diag]
    omega = t(df70[subset & ind.diag==0,]) %*% jcov.byp %*% df70[subset & ind.diag==0,]
    se = sqdg((crossprod(df70/sqrt(wgt0)) + 2*crossprod(df70, df7 - df70) + crossprod(df7 - df70))*nfull/(nfull-1) + omega)
    
  } else {
    df7 = cbind(df2$df_beta, df3$df_beta, tp_obj$DF)
    df70 = NULL
    if(method=="ful") {se = sqdg(crossprod(df7)*nfull/(nfull-1))}
    if(method=="des") {
      ind.diag = datph1[subset, ind.diag]
      omega = t(df7[ind.diag==0,]) %*% jcov.byp %*% df7[ind.diag==0,]
      se = sqdg(crossprod(df7/sqrt(wgts[subset]))*nfull/(nfull-1) + omega)}
  }
  se.cols = unique(gsub("[.]\\d+", "", colnames(df7)[-(1:(length(est2) + length(est3)))]))
  est = c(est2, est3, c(tp_obj$est[ , se.cols])) 
  out = rbind(est, se)
  colnames(out) = names(se)
  colnames(out)[1:(length(est2) + length(est3))] = c(paste0("beta12", 1:2), paste0("beta13", 1:2))
  out = data.table(out, method = method, type = c("est", "se"), model = "sd")
  return(out)
}


pb.out.method = function(est, df, df0=NULL, tp_obj, datph1.repli, jcov.byp, method, tp0_obj=NULL){
  tau1 = as.numeric(rownames(tp_obj$est))
  ntau1 = length(tau1)
  
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
  se.cols = unique(gsub("[.]\\d+", "", colnames(df7)[-(1:length(est))]))
  beta_len = length(est)
  beta_col = c(paste0("beta12", 1:2), "gamma", paste0("beta13", 1:2))
  est = c(est, c(tp_obj$est[ , se.cols])) 
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