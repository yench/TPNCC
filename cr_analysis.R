library(data.table)
library(survival)
library(kableExtra)
library(ggplot2)
workdir = 'C:/Users/Yen/Documents/GitHub/TPNCC/'
source(paste0(workdir, 'cr_functions.R'))
set.seed(442)

n = 10000         # cohort size
mctrl = 5         # number of controls per case
max.diag.t = 6    # administrative censoring time of the original NCC study
max.cohort.t = 10 # administrative censoring time of the cohort study
cov2u = 0.75      # cov(Z2, U): covariance between the partially observed Z2 and the ancillary variable U
Smax12 = 0.90     # marginal baseline survival fraction of the 1->2 transition at max.death.t
Smax13 = 0.85     # marginal baseline survival fraction of the 1->3 transition at max.death.t
beta121 = 0.2; beta122 = 0.5
beta131 = 0.1; beta132 = 0.4
zstar = rbind(c(0, 0), c(1, 0), c(0, 1), c(1, 1)) # covariates to obtain predictions for
tau1 = 5          # (tau0, tau1] = (0, 5], where tau0 is assumed to be 0

# simulate competing risks data in Figure 1a
datph1 = setDT(simul.cr(n, mctrl=mctrl, max.diag.t=max.diag.t, max.cohort.t=max.cohort.t, Smax12=Smax12, Smax13=Smax13,
                        beta121=beta121, beta122=beta122, beta131=beta131, beta132=beta132, cov2u=cov2u))
## the simulated data have the following columns
### ID
### z1: the covariate observed in all cohort members
### z2: the covariate only collected among members in the NCC study
###     here all members have non-missing z2 in order to perform a full cohort analysis
### u: the ancillary variable correlated with z2
### diagtime: time to disease diagnosis (incidence) in the original NCC study
### ind.diag: indicator of diagnosis (incidence) in the original NCC study
### eventime: time to death, subject to administrative and random censoring
### ind.fail: failure type indicator
###           (0=censoring, 2=death with the disease, 3=death without the disease)
### nrisk: number of subjects at risk in the cohort at diagtime
### incl.prob: probability to be included in the NCC study
### wgt: design weight = 1/incl.prob
### ind.ncc: TRUE if included in the NCC sample


# exclude subjects with 0 prob to be included in the NCC sample 
# (non-cases (ind.diag=0) with diagtime smaller than the smallest case (ind.diag=1) diagtime)
datph1 = datph1[incl.prob>0, ]

datph1[ , `:=` (ind.12 = as.integer(ind.fail==2), ind.13 = as.integer(ind.fail==3))]


##################################################################
#                      Weight calibration                        #
##################################################################
# step 1: predict Z2 using fully-observed variables
imp.model = glm(z2 ~ z1 + u, weights=wgt, data=datph1[ind.ncc==TRUE, ])
datph1[ , z2.pred := predict(imp.model, datph1)]

# step 2: Breslow-type calibration
## step 2.1
covars_calib = c("z1", "z2.pred")
aux12 = calib_bre(datph1, "eventime", "ind.12", covars_calib)
aux13 = calib_bre(datph1, "eventime", "ind.13", covars_calib)

## step 2.2
### replicate data for the proportional baselines model
datph1.repli = datph1[rep(1:.N, each=2), ]
datph1.repli[ , cause := rep(2:3, .N/2)]
datph1.repli[ , `:=`(ind.1 = as.integer((cause==2 & ind.12==1)|(cause==3 & ind.13==1)),
                     z3=as.integer(cause==3), z4=z1, z5=z2, z5.pred=z2.pred)] # cell-mean coding
datph1.repli[cause==2, `:=` (z4=0, z5=0, z5.pred=0)]
datph1.repli[cause==3, `:=` (z1=0, z2=0, z2.pred=0)]
covars_calib.repli = c("z1", "z2.pred", "z3", "z4", "z5.pred")
aux1  = calib_bre(datph1.repli, "eventime", "ind.1", covars_calib.repli)

# step 3: estimation with Breslow weights
jcov.byp = jointVP(datph1[ind.ncc==TRUE, ], m=mctrl)
covars_final = c("z1", "z2")
covars_final.repli = c("z1", "z2", "z3", "z4", "z5")
sdbre = coxModel.sd("bre", datph1, aux12, aux13, covars_final, tau1, zstar, jcov.byp)
pbbre = coxModel.pb("bre", datph1.repli, aux1, covars_final.repli, tau1, zstar, sdbre$tp.lookup, jcov.byp)

# step 4: Shin-type calibration
## step 4.1
shiA12 = cbind(aux12$A, scale(btn(datph1$eventime, 0, tau1)*datph1$ind.12),
               scale(thi(0, datph1$eventime, 0, tau1) *
               exp(as.matrix(datph1[ , ..covars_calib]) %*% t(sdbre$out[type=="est", paste0("beta12", 1:2)]))))

shiA13 = cbind(aux13$A, scale(btn(datph1$eventime, 0, tau1)*datph1$ind.13),
               scale(thi(0, datph1$eventime, 0, tau1) * 
               exp(as.matrix(datph1[ , ..covars_calib]) %*% t(sdbre$out[type=="est", paste0("beta13", 1:2)]))))

subset = datph1[ , ind.ncc]
shiwgt12 = rakecal(shiA12[subset, ], datph1[subset, wgt], colSums(shiA12)) 
shiwgt13 = rakecal(shiA13[subset, ], datph1[subset, wgt], colSums(shiA13)) 

## step 4.2
shiA1 = cbind(aux1$A, scale(btn(datph1.repli$eventime, 0, tau1)*datph1.repli$ind.1),
              scale(thi(0, datph1.repli$eventime, 0, tau1) * 
                    exp(as.matrix(datph1.repli[ , ..covars_calib.repli]) %*% 
                          t(pbbre$out[type=="est", c(paste0("beta12", 1:2), "gamma", paste0("beta13", 1:2))]))))
shiwgt1 = rakecal(shiA1[rep(subset, each=2), ], datph1.repli[rep(subset, each=2), wgt], colSums(shiA1))

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

out = rbind(sdful$out, sddes$out, sdbre$out, sdshi$out, pbful$out, pbdes$out, pbbre$out, pbshi$out, fill=TRUE)

# calculate true values
beta12 = c(beta121, beta122); gamma12 = log(-log(Smax12)/max.cohort.t)
beta13 = c(beta131, beta132); gamma13 = log(-log(Smax13)/max.cohort.t) 
gamma = gamma13 - gamma12

bInt12 = exp(gamma12)*tau1
bInt13 = exp(gamma13)*tau1

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
bInt1.param = c("bInt12", "bInt13")
P.param = paste0("P", rep(c(11, 12, 13), each=nrow(zstar)), ".zst", rep(1:nrow(zstar)))
param.lookup = data.table(param = as.factor(c(beta.param, bInt1.param, P.param)), 
                          true.value = c(beta.true, bInt12, bInt13, P11, P12, P13))

out2 = melt(out, id.vars = c("model", "method", "type"),
            measure.vars = c(beta.param, bInt1.param, P.param), variable.name="param")
out2 = dcast(out2, model + method + param ~ type)
out2[param.lookup, true.value := true.value, on=.(param)] 
out2[ , order := 1:.N, by=.(model, method)]
out.wide = dcast(out2, model + param + order + true.value ~ method, value.var=c("est", "se"))
setcolorder(out.wide, c("model", "param", "true.value", "est_ful", "se_ful", "est_des", "se_des", "est_bre", "se_bre", "est_shi", "se_shi"))

# output results
kable(out.wide[ , -"order"], col.names=c("Model", "Parameter", "True value", rep(c("Estimate", "SE"), 4))) %>% 
  add_header_above(c(" "=3, "Full"=2, "Design"=2, "Breslow"=2, "Shin"=2)) %>% 
  kable_classic()

# graphically evaluate the proportional baselines assumption
baseline.wide = data.table(hk = c(rep(12, length(sddes$bint12)), rep(13, length(sddes$bint13))),
                           tj = c(attr(sddes$bint12, "time"), attr(sddes$bint13, "time")), 
                           Design  = c(sddes$bint12, sddes$bint13), 
                           Breslow = c(sdbre$bint12, sdbre$bint13),
                           Shin    = c(sdshi$bint12, sdshi$bint13))
baseline = melt(baseline.wide, id.vars=c("hk", "tj"), measure.vars=c("Design", "Breslow", "Shin"), 
                value.name="dLambda", variable.name="method")
baseline[ , Lambda := cumsum(dLambda), by=.(hk, method)]

# the proportional baselines model appears valid
ggplot(data=baseline) +
  facet_wrap(vars(method)) +
  geom_line(aes(x=tj, y=Lambda, group=hk, color=factor(hk)))
