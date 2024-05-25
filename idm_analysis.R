library(data.table)
library(survival)
workdir = 'C:/Users/Yen/Documents/GitHub/TPNCC/'
source(paste0(workdir, 'TPNCC_idm_functions.R'))

n = 10000
mctrl = 5
max.diag.t = 6
max.death.t = 10
cov2u = 0.75
Smax12 = 0.90
Smax13 = 0.85
Smax23 = 0.6
beta121 = 0.2; beta122 = 0.5
beta131 = 0.1; beta132 = 0.4
beta231 = 0.3; beta232 = 0.7
zstar = rbind(c(0, 0), c(1, 0), c(0, 1), c(1, 1))
tau = c(2, 5) # (tau0, tau1]

# simulate illness-death model data
datph1 = setDT(simul.idm(n, mctrl=mctrl, max.diag.t=max.diag.t, max.death.t=max.death.t, Smax12=Smax12, Smax13=Smax13, Smax23=Smax23,
                         beta121=beta121, beta122=beta122, beta131=beta131, beta132=beta132, beta231=beta231, beta232=beta232, cov2u=cov2u))

# exclude subjects with 0 prob to be included in the NCC sample 
# (subjects with smaller event time than the smallest 1-2 transition time)
datph1 = datph1[incl.prob>0, ]


##################################################################
#                      Weight calibration                        #
##################################################################
# step 1: predict Z2
imp.model = glm(z2 ~ z1 + u, weights=wgt, data=datph1[ind.ncc==TRUE, ])
datph1[ , z2.pred := predict(imp.model, datph1)]

# step 2: Breslow-type calibration
## step 2.1
covars_calib = c("z1", "z2.pred")
aux12 = calib_bre(datph1, "nonttime", "ind.12", covars_calib)  
aux13 = calib_bre(datph1, "nonttime", "ind.13", covars_calib)
aux23 = calib_bre(datph1[ind.12==1], "nonttime, termtime", "ind.23", covars_calib)

## step 2.2
### replicate data for the proportional baselines model
datph1.repli = datph1[rep(1:.N, each=2), ]
datph1.repli[ , end := rep(2:3, .N/2)]
datph1.repli[ , `:=`(ind.1 = as.integer((end==2 & ind.12==1)|(end==3 & ind.13==1)),
                     z3=as.integer(end==3), z4=z1, z5=z2, z5.pred=z2.pred)] # cell-mean coding
datph1.repli[end==2, `:=` (z4=0, z5=0, z5.pred=0)]
datph1.repli[end==3, `:=` (z1=0, z2=0, z2.pred=0)]
covars_calib.repli = c("z1", "z2.pred", "z3", "z4", "z5.pred")
aux1  = calib_bre(datph1.repli, "nonttime", "ind.1", covars_calib.repli)

# step 3: estimation with Breslow weights
jcov.byp = jointVP(datph1[ind.ncc==TRUE, ], m=mctrl)
outbre = coxModel(datph1, datph1.repli, tau, "bre", jcov.byp, 
                  list(wgt12=aux12$wgt, A12=aux12$A, wgt13=aux13$wgt, A13=aux13$A,
                       wgt23=aux23$wgt, A23=aux23$A, wgt1 =aux1$wgt,  A1 =aux1$A))

# step 4: Shin-type calibration
## step 4.1
shiA12 = cbind(aux12$A, scale(btn(datph1$nonttime, tau[1], tau[2])*datph1$ind.12),
               scale(thi(0, datph1$nonttime, tau[1], tau[2]) *
                     exp(as.matrix(datph1[ , ..covars_calib]) %*% t(outbre[type=="est" & model=="sd", paste0("beta12", 1:2)]))))

shiA13 = cbind(aux13$A, scale(btn(datph1$nonttime, tau[1], tau[2])*datph1$ind.13),
               scale(thi(0, datph1$nonttime, tau[1], tau[2]) * 
                     exp(as.matrix(datph1[ , ..covars_calib]) %*% t(outbre[type=="est" & model=="sd", paste0("beta13", 1:2)]))))

shiA23 = cbind(aux23$A, scale(btn(datph1[ind.12==1, termtime], tau[1], tau[2])*datph1[ind.12==1, ind.23]),
               scale(thi(datph1[ind.12==1, nonttime], datph1[ind.12==1, termtime], tau[1], tau[2]) *
                     exp(as.matrix(datph1[ind.12==1 , ..covars_calib]) %*% t(outbre[type=="est" & model=="sd", paste0("beta23", 1:2)]))))

subset = datph1[ , ind.ncc]
shiwgt12 = rakecal(shiA12[subset, ], datph1[subset, wgt], colSums(shiA12), eta=TRUE) 
shiwgt13 = rakecal(shiA13[subset, ], datph1[subset, wgt], colSums(shiA13), eta=TRUE) 
shiwgt23 = rakecal(shiA23[datph1[ind.12==1, ind.ncc], ], datph1[ind.12==1 & subset, wgt], colSums(shiA23), eta=TRUE) 

## step 4.2
shiA1 = cbind(aux1$breA, scale(btn(datph1.repli$nonttime, tau[1], tau[2])*datph1.repli$ind.1),
              scale(thi(0, datph1.repli$nonttime, tau[1], tau[2]) * 
                    exp(as.matrix(datph1.repli[ , ..covars_calib.repli]) %*% 
                        t(outbre[type=="est" & model=="pb", c(paste0("beta12", 1:2), "gamma", paste0("beta13", 1:2))]))))

shiwgt1 = rakecal(shiA1[rep(subset, each=2), ], datph1.repli[rep(subset, each=2), wgt], colSums(shiA1),  eta=TRUE)

# step 5: estimation with Shin weights
outshi = coxModel(datph1, datph1.repli, tau, "shi", jcov.byp,
                  list(wgt12=shiwgt12, A12=shiA12, wgt13=shiwgt13, A13=shiA13,
                       wgt23=shiwgt23, A23=shiA23, wgt1 =shiwgt1,  A1 =shiA1))


##################################################################
#                           Benchmark                            #
##################################################################
# full cohort analysis
outful = coxModel(datph1, datph1.repli, tau, "ful")

# design weights
outdes = coxModel(datph1, datph1.repli, tau, "des", jcov.byp)

out = rbind(outful, outdes, outbre, outshi)

# calculate true values
beta12 = c(beta121, beta122); gamma12 = log(-log(Smax12)/max.death.t)
beta13 = c(beta131, beta132); gamma13 = log(-log(Smax13)/max.death.t) 
beta23 = c(beta231, beta232); gamma23 = log(-log(Smax23)/max.death.t) 
gamma = gamma13 - gamma12

bInt12 = bInt1 = exp(gamma12)*(tau[2]-tau[1])
bInt13 = exp(gamma13)*(tau[2]-tau[1])
bInt23 = bInt2 = exp(gamma23)*(tau[2]-tau[1])

P11.zst.t = function(z, tau) {exp(-c(exp(gamma12 + z %*% beta12) + exp(gamma13 + z %*% beta13))*(tau[2]-tau[1]))}
P12.zst.t = function(z, tau){
  rr12 = exp(gamma12 + z %*% beta12)
  rr23 = exp(gamma23 + z %*% beta23)
  rr1 = rr12 + exp(gamma13 + z %*% beta13)
  out = rr12*exp(rr1*tau[1]-rr23*tau[2])*(exp(-tau[1]*(rr1-rr23))-exp(-tau[2]*(rr1-rr23)))/(rr1-rr23)
  return(out)}
P22.zst.t = function(z, tau) {exp(-c(exp(gamma23 + z %*% beta23))*(tau[2]-tau[1]))}

P11 = apply(zstar, 1, function(z){P11.zst.t(z, tau)})
P12 = apply(zstar, 1, function(z){P12.zst.t(z, tau)})
P13 = 1-P11-P12
P23 = apply(zstar, 1, function(z){1-P22.zst.t(z, tau)})

beta.param = c("beta121", "beta122", "gamma", "beta131", "beta132", "beta231", "beta232")
beta.true = unlist(mget(beta.param))
bInt.param = c("bInt12", "bInt13", "bInt23", "bInt1", "bInt2")
P.param = paste0("P", rep(c(11, 12, 13, 23), each=nrow(zstar)), ".", rep(1:nrow(zstar)))
param.lookup = data.table(param = as.factor(c(beta.param, bInt.param, P.param)), 
                          true.value = c(beta.true, bInt12, bInt13, bInt23, bInt1, bInt2, P11, P12, P13, P23))

out2 = melt(out, id.vars = c("model", "method", "type"),
            measure.vars = c(beta.param, bInt.param, paste0(rep(c("P11.", "P12.", "P13.", "P23."), each=4), 1:4)), variable.name = "param")
out2 = dcast(out2, model + method + param ~ type)
out2[param.lookup, true.value := true.value, on = .(param)]
