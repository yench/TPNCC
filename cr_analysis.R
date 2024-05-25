library(data.table)
library(survival)
workdir = 'C:/Users/Yen/Documents/GitHub/TPNCC/'
source(paste0(workdir, 'TPNCC_cr_functions.R'))

n = 10000
mctrl = 5
max.diag.t = 6    # the incidence NCC study has shorter follow-up (6) than the cohort follow-up (10)
max.cohort.t = 10
cov2u = 0.75
Smax12 = 0.90
Smax13 = 0.85
beta121 = 0.2; beta122 = 0.5
beta131 = 0.1; beta132 = 0.4
zstar = rbind(c(0, 0), c(1, 0), c(0, 1), c(1, 1))
tau1 = 5 # (tau0, tau1] = (0, 5]

# simulate competing risks data
datph1 = setDT(simul.cr(n, mctrl=mctrl, max.diag.t=max.diag.t, max.cohort.t=max.cohort.t, Smax12=Smax12, Smax13=Smax13,
                        beta121=beta121, beta122=beta122, beta131=beta131, beta132=beta132, cov2u=cov2u))

# exclude subjects with 0 prob to be included in the NCC sample 
# (subjects with smaller event time than the smallest 1-2 transition time)
datph1 = datph1[incl.prob>0, ]

datph1[ , `:=` (ind.12 = as.integer(ind.fail==2), ind.13 = as.integer(ind.fail==3))]


##################################################################
#                      Weight calibration                        #
##################################################################
# step 1: predict Z2
imp.model = glm(z2 ~ z1 + u, weights=wgt, data=datph1[ind.ncc==TRUE, ])
datph1[ , z2.pred := predict(imp.model, datph1)]

# step 2: Breslow-type calibration
## step 2.1
covars_calib = c("z1", "z2.pred")
aux2 = calib_bre(datph1, "eventime", "ind.12", covars_calib)
aux3 = calib_bre(datph1, "eventime", "ind.13", covars_calib)

## step 2.2
### replicate data for the proportional baselines model
datph1.repli = datph1[rep(1:.N, each=2), ]
datph1.repli[ , cause := rep(2:3, .N/2)]
datph1.repli[ , `:=`(ind.1 = as.integer((cause==2 & ind.12==1)|(cause==3 & ind.13==1)),
                     z3=as.integer(cause==3), z4=z1, z5=z2, z5.pred=z2.pred)] # cell-mean coding
datph1.repli[cause==2, `:=` (z4=0, z5=0, z5.pred=0)]
datph1.repli[cause==3, `:=` (z1=0, z2=0, z2.pred=0)]
covars_calib.repli = c("z1", "z2.pred", "z3", "z4", "z5.pred")
aux  = calib_bre(datph1.repli, "eventime", "ind.1", covars_calib.repli)

# step 3: estimation with Breslow weights
jcov.byp = jointVP(datph1[ind.ncc==TRUE, ], m=mctrl)
covars_final = c("z1", "z2")
covars_final.repli = c("z1", "z2", "z3", "z4", "z5")
sdbre = coxModel.sd("bre", datph1, aux2, aux3, covars_final, tau1, zstar, jcov.byp)
pbbre = coxModel.pb("bre", datph1.repli, aux, covars_final.repli, tau1, zstar, sdbre$tp.lookup, jcov.byp)

# step 4: Shin-type calibration
## step 4.1
shiA12 = cbind(aux2$A, scale(btn(datph1$eventime, 0, tau1)*datph1$ind.12),
               scale(thi(0, datph1$eventime, 0, tau1) *
               exp(as.matrix(datph1[ , ..covars_calib]) %*% t(sdbre$out[type=="est", paste0("beta12", 1:2)]))))

shiA13 = cbind(aux3$A, scale(btn(datph1$eventime, 0, tau1)*datph1$ind.13),
               scale(thi(0, datph1$eventime, 0, tau1) * 
               exp(as.matrix(datph1[ , ..covars_calib]) %*% t(sdbre$out[type=="est", paste0("beta13", 1:2)]))))

subset = datph1[ , ind.ncc]
shiwgt12 = rakecal(shiA12[subset, ], datph1[subset, wgt], colSums(shiA12), eta=TRUE) 
shiwgt13 = rakecal(shiA13[subset, ], datph1[subset, wgt], colSums(shiA13), eta=TRUE) 

## step 4.2
shiA1 = cbind(aux$A, scale(btn(datph1.repli$eventime, 0, tau1)*datph1.repli$ind.1),
              scale(thi(0, datph1.repli$eventime, 0, tau1) * 
                    exp(as.matrix(datph1.repli[ , ..covars_calib.repli]) %*% 
                          t(pbbre$out[type=="est", c(paste0("beta12", 1:2), "gamma", paste0("beta13", 1:2))]))))
shiwgt1 = rakecal(shiA1[rep(subset, each=2), ], datph1.repli[rep(subset, each=2), wgt], colSums(shiA1),  eta=TRUE)

# step 5: estimation with Shin weights
sdshi = coxModel.sd("shi", datph1, list(A=shiA12, wgt=shiwgt12), list(A=shiA13, wgt=shiwgt13), covars_final, tau1, zstar, jcov.byp)
pbshi = coxModel.pb("shi", datph1.repli, list(A=shiA1, wgt=shiwgt1), covars_final.repli, tau1, zstar, sdbre$tp.lookup, jcov.byp)
  

##################################################################
#                           Benchmark                            #
##################################################################
# full cohort analysis
sdful = coxModel.sd("ful", datph1, NULL, NULL, covars_final, tau1, zstar)
pbful = coxModel.pb("ful", datph1.repli, NULL, covars_final.repli, tau1, zstar, sdful$tp.lookup)

# design weights
sddes = coxModel.sd("des", datph1, NULL, NULL, covars_final, tau1, zstar, jcov.byp)
pbdes = coxModel.pb("des", datph1.repli, NULL, covars_final.repli, tau1, zstar, sdbre$tp.lookup, jcov.byp)

out = rbind(sdful$out, sddes$out, sdbre$out, sdshi$out, pbful$out, pbdes$out, pbbre$out, pbshi$out, fill = TRUE)

# calculate true values
beta12 = c(beta121, beta122); gamma12 = log(-log(Smax12)/max.cohort.t)
beta13 = c(beta131, beta132); gamma13 = log(-log(Smax13)/max.cohort.t) 
gamma = gamma13 - gamma12

bInt2 = exp(gamma12)*tau1
bInt3 = exp(gamma13)*tau1

P11.zst.t = function(z, tau1) {exp(-c(exp(gamma12 + z %*% beta12) + exp(gamma13 + z %*% beta13))*tau1)}
P12.zst.t = function(z, tau1){
  rr12 = exp(gamma12 + z %*% beta12)
  rr1 = rr12 + exp(gamma13 + z %*% beta13)
  out = rr12*(1 - exp(-rr1*tau1))/rr1
  return(out)
}

P11 = apply(zstar, 1, function(z){P11.zst.t(z, tau1)})
P12 = apply(zstar, 1, function(z){P12.zst.t(z, tau1)})
P13 = 1-P11-P12

beta.param = c("beta121", "beta122", "gamma", "beta131", "beta132")
beta.true = unlist(mget(beta.param))
bInt.param = paste0(c("bInt2", "bInt3"), '.', tau1)
P.param = paste0("P", rep(c(11, 12, 13), each=nrow(zstar)), ".zst", rep(1:nrow(zstar)), '.', tau1)
param.lookup = data.table(param = as.factor(c(beta.param, bInt.param, P.param)), 
                          true.value = c(beta.true, bInt2, bInt3, P11, P12, P13))

out2 = melt(out, id.vars = c("model", "method", "type"),
            measure.vars = c(beta.param, bInt.param, P.param), variable.name = "param")
out2 = dcast(out2, model + method + param ~ type)
out2[param.lookup, true.value := true.value, on = .(param)]
