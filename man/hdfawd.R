library(foreign)
library(rdrobust)
library(rdd)
library(rddtools)
library(RItools)
library(rddensity)
hdfaw = readRDS("~/Downloads/GD_award.rds")
#load("/Users/hayyen/Downloads/GD_award.rds")
head(hdfaw)
hdfawdat<-hdfaw[,c(1,3,5,11:20,34)]
head(hdfawdat)
hdfawdata<-hdfawdat[,-10]
summary(hdfawdata)
write.csv(hdfawdata, "~/Downloads/GD_award.cvs")
#####
library(rdd)
plot(hdfawdata$yviews[hdfawdata$ypatient >= 1 & hdfawdata$ypatient <= 30], 
     hdfawdata$ypatient[hdfawdata$ypatient >= 1 & hdfawdata$ypatient <= 30])

hdfawdata_rd = RDestimate(yviews~ypatient, data=hdfawdata[hdfawdata$yviews >= 100 
                                                          & hdfawdata$yviews <= 30000,])
hhs_rd
IKbandwidth(X=hhs$rv[hhs$rv >= -10 & hhs$rv <= 10], Y=hhs$dv[hhs$rv >= -10 & hhs$rv <= 10], cutpoint=0)
summary(hhs_rd)
plot(hhs_rd)



##### 
#Different packages with more built in options
library(rddtools)
library(RItools)
hhs_rd2 = rdd_data(x=rv, y=dv, covar= c(rv*treat, treat), cutpoint=0,data=hhs[hhs$rv >= -10 & hhs$rv <= 10,])
head(hhs_rd2)
summary(hhs_rd2)
bw_hhs <- rdd_bw_ik(hhs_rd2) #Imbens-Kalyanaraman Optimal Bandwidth Calculation
plot(hhs_rd2, nbins=100)
#plot(hhs_rd2, nbins=10)
hhs_rd2_reg = rdd_reg_lm(hhs_rd2) #linear RDD
hhs_rd2_np = rdd_reg_np(hhs_rd2) #non-parametric RDD

print(hhs_rd2_reg)
summary(hhs_rd2_reg)
print(hhs_rd2_np)

plot(hhs_rd2_reg, nbins=100)
plot(hhs_rd2_np, nbins=100)



###########
hdfawcsv<-hdfawdata[,c(-1,-2)]
hdfawcsv<-as.numeric(hdfawcsv)
hdfawcsv[2,1]
library(sqldf)
hdfaw<-sqldf("select * from hdfawcsv where award like '%16%'")
hdfaw2<-sqldf("select * from hdfaw where award like '%;%'")
hdfaw1<-sqldf("select * from hdfaw where award not like '%;%'")
hdfawn<-sqldf("select * from hdfawcsv where award not like '%16%'")
#hdfawcsv$treat<- ifelse(!is.na(hdfawcsv$award), 0, 1)
hdfaw$treat<- 1
hdfawn$treat<- 0
hdfaw1$once<-1
hdfaw2$once<-0
hdfaw1num<-hdfaw1[,-1]
hdfaw2num<-hdfaw2[,-1]
write.csv(hdfawcsv, "~/Downloads/hdfawcsv.cvs")

hdfawdat1<-sqldf("select * from hdfawdata where award not like '%;%'")
hdfawdat2<-sqldf("select * from hdfawdata where award  like '%;%'")
 hdfawdat1$once<-1
 hdfawdat2$once<-0
 hdfonce<-rbind(hdfawdat1,hdfawdat2)
 hdfoncenum<-hdfonce[,c(-2,-3)]
 write.csv(hdfoncenum, "~/Downloads/hdfoncenum.cvs")
 
 
 
 ###################################
 ## Load packages and data
 
 library(ggplot2)
 library(granova)
 library(granovaGG)
 library(Matching)
 library(MatchIt)
 library(party)
 library(PSAgraphics)
 library(rbounds)
 library(rpart)
 library(tree)
 library(rgenoud)
 library(reshape)
 
# data(lalonde, package='Matching')
 
 haodf_psm<-read.csv("file:///D:/icthealth/web/award1807/haodf14w_t.csv",header = TRUE)
 
 ################################################################################
 ## Phase I
 #rename(haodf_psm, c("hot" = "DRR"))
 ## Using logistic regression for estimating propensity scores 
 hdf.formu <- treat ~ hot+ thank+gift+contr+article
 hdf.glm <- glm(hdf.formu, family=binomial, data=haodf_psm)
 summary(hdf.glm)
 
 hdf.formu2 <- treat ~  hot + thank + gift
 hdf.glm2 <- glm(hdf.formu2, data=haodf_psm, family=binomial)
 summary(hdf.glm2)
 
 df <- data.frame(ps1=fitted(hdf.glm), ps2=fitted(hdf.glm2), Tr=haodf_psm$treat)
 
 ggplot(df, aes(x=ps1)) + geom_histogram() + facet_wrap(~ Tr, ncol=1)
 
 ggplot(df, aes(x=ps2)) + geom_histogram() + facet_wrap(~ Tr, ncol=1)
 
 summary(hdf.glm)
 
 
 hdf_ps <- fitted(hdf.glm)  # Propensity scores
 hdf_Y  <- haodf_psm$patient1  # Dependent variable, patient1
 hdf_Tr <- haodf_psm$treat # Treatment indicator
 
 hdf_rr  <- Match(Y=hdf_Y,Tr=hdf_Tr,X=hdf_ps,M=1)         #M=1, 1:1
 ?Match
 
 summary(hdf_rr)
 #hdf_1tc<-hdf_rr$mdata
 hdf_1tc<-rbind(haodf_psm[hdf_rr$index.treated,],haodf_psm[hdf_rr$index.control,])
write.csv(hdf_1tc,"file:///C:/Users/Administrator/Desktop/hdf_1tc.csv")  
#hdf_1tc<-hdf_rr$index.treated

#########
# Rosenbaum sensitivity Test
haodf_psm_m <- Match(Y=hdf_Y,Tr=hdf_Tr,X = hdf_ps, replace = FALSE)
?Match
# for sensitivity to p-value
psens(haodf_psm_m, Gamma = 2, GammaInc = 0.1) 
hdf_p_psen<-psens(haodf_psm_m, Gamma = 2, GammaInc = 0.25) 
hdf_p_psen
# for sensitivity to ATT
#hdf_att_psen<-hlsens(haodf_psm_m, Gamma = 2, GammaInc = 0.1, 0.1) 
hdf_att_psen<-hlsens(haodf_psm_m, pr=0.1, Gamma = 2, GammaInc = 0.25) 
hdf_att_psen

# Rosenbaum sensitivity Test
hdf_ps <- fitted(hdf.glm)  # Propensity scores
hdf_Yv  <- haodf_psm$views  # Dependent variable, patient1
hdf_Tr <- haodf_psm$treat # Treatment indicator
haodf_psm_mv <- Match(Y=hdf_Yv,Tr=hdf_Tr,X = hdf_ps, replace = FALSE)
?Match
# for sensitivity to p-value
#hdfv_p_sen<-psens(haodf_psm_mv, Gamma = 2, GammaInc = 0.1) 
hdfv_p_sen<-psens(haodf_psm_mv, Gamma = 2, GammaInc = 0.25) 
hdfv_p_sen
# for sensitivity to ATT
#hlsens(haodf_psm_mv, Gamma = 2, GammaInc = 0.1, 0.1) 
?hlsens
#hlsens(haodf_psm_mv, pr=0.1, Gamma = 2, GammaInc=0.25) 
hdfv_att_sen<-hlsens(haodf_psm_mv, pr=0.2, Gamma = 2, GammaInc=0.25)
#PSM Matching                          #
############################################################

############################################################
#               Analyses before Matching                   #
############################################################
#hdf_awnoaw1prof<-hdf_awnoaw[hdf_awnoaw$period==1,]
# T-test on the original data
t.test(views ~ treat, data =haodf_psm)
t.test(patient1 ~ treat, data =haodf_psm)

# Regression on the original data
haodf_psm.formu <- treat ~ hot+ thank+gift+contr+article
haodf_psm.glm <- glm(haodf_psm.formu, family=binomial(), data=haodf_psm)
summary(haodf_psm.glm)


############################################################
#                    PSM Matching                          #
############################################################

# Run optimal matching
haodf_psm.out <- matchit(haodf_psm.formu, 
                          data = haodf_psm, method = "optimal", distance = "logit", ratio = 1)

# Checking Balance 
summary(haodf_psm.out)###############difference of balance
plot(haodf_psm.out, type="jitter")
plot(haodf_psm.out, type="hist",ylab="density(%)")
plot(haodf_psm.out, type="QQ")

plot(summary(haodf_psm.out, standardize = T), interactive = T,col=c(1:5))#
#legend(1, 95, legend=c("DRR","thank","gift","contr","article"),  col=c(1:5))

################

#Optimal_Match
#match_on
library("optmatch", lib.loc="C:/Program Files/R/R-3.5.0/library")
match_on.examples <- list()
haodf_psm<-read.csv("file:///D:/icthealth/web/award1807/haodf2w_t.csv",header = TRUE)

################################################################################
## Phase I

## Using logistic regression for estimating propensity scores 
(hdf.formu <- treat ~ hot+ thank+gift+contr+article)
hdf.glm <- glm(hdf.formu, family=binomial(), data=haodf_psm)
summary(hdf.glm)
############match_on, pS (t,c), dist(t,c)
match_on.examples$ps1 <- match_on(hdf.glm)
write.csv(match_on.examples$ps1,"file:///C:/Users/Administrator/Desktop/hdf_ps.csv")  

#hdfPS <- predict(hdf.glm)
hdf_fm1 <- fullmatch(match_on.examples$ps1, data =  haodf_psm)
write.csv(hdf_fm1,"file:///C:/Users/Administrator/Desktop/hdf_fm.csv")  
hdf_fm1_ma<-as.matrix(hdf_fm1)


### Pair matching on the variable `t1`:  2-->20

hdf_mhd <- match_on(hdf.formu, data = haodf_psm) + caliper(match_on(hdf.glm), 20)
#write.csv(hdf_mhd,"file:///C:/Users/Administrator/Desktop/hdf_mhd.csv")  

( hdf_pm2 <- pairmatch(hdf_mhd, data = haodf_psm) )
summary(hdf_pm2)

### Propensity balance assessment. Requires RItools package.
if(require(RItools)) summary(hdf_pm2, hdf.glm)

### to make sure observations are in the proper order:
all.equal(names(hdf_pm2), row.names(haodf_psm))
### So our data frame including the matched sets is just
hdf_pm<-cbind(haodf_psm, matches=hdf_pm2)
write.csv(hdf_pm,"file:///C:/Users/Administrator/Desktop/hdf_pm.csv")  


#(absdist <- match_on(hdf.glm))
#hdf_pm<-pairmatch(absdist, data = haodf_psm)

 ## Matching
 # one-to-one matching with replacement (the "M=1" option).
 # Estimating the treatment effect on the treated (default is ATT).
 hdf_rr.att <- Match(Y=hdf_Y, Tr=hdf_Tr, X=hdf_ps, M=1, estimand='ATT')
 summary(hdf_rr.att) # The default estimate is ATT here
 ls(hdf_rr.att)
 hdf_rr.att$mdata$Y
 
 
 hdf_rr.ate <- Match(Y=hdf_Y, Tr=hdf_Tr, X=hdf_ps, M=1, estimand='ATE')#Average treatment effect
 summary(hdf_rr.ate)
 
 hdf_rr.atc <- Match(Y=hdf_Y, Tr=hdf_Tr, X=hdf_ps, M=1, estimand='ATC')#Average treatment effect of the controls
 summary(hdf_rr.atc)
 
 
 hdf_rr.att$est
 hdf_rr.ate$est
 hdf_rr.att$estimand
 hdf_rr.ate$estimand
 
 
 
 #Solutions to fundamental Matching problem: do we achieve balanced sets between T and C among the observed covariates?
 
 ## Genetic Matching
 hdf_rr.gen <- GenMatch(Tr=hdf_Tr, X=hdf_ps, 
                    BalanceMatrix=haodf_psm[,all.vars(hdf.formu)[-1]],
                    estimand='ATE', M=1, pop.size=16)
 hdf_rr.gen$matches
 hdf_rr.gen.mout <- Match(Y=hdf_Y, Tr=hdf_Tr, X=hdf_ps, estimand='ATE', Weight.matrix=hdf_rr.gen)
 summary(hdf_rr.gen.mout)
 
 
 
 
 
 library(CBPS)
 
 ##Load the LaLonde data
 data(LaLonde, package="CBPS")
 ## Estimate CBPS via logistic regression
 fit <- CBPS(formula=lalonde.formu, data = lalonde, ATT = TRUE, twostep = FALSE, standardize = TRUE)
 
 summary(fit)
 
 plot(fit) #expand window before plotting
 balance(fit)
 
 m.out <- Match(treat ~ fitted(fit), method = "nearest", data = lalonde, 
                replace = TRUE)
 rr.att.CBPS <- Match(Y=Y, Tr=Tr, X=fitted(fit), M=1, ties=FALSE, replace=FALSE, estimand='ATT')
 
 summary(rr.att.CBPS)
 
 
 # Can we use Machine Learning to estimate Better Propensity Scores?
 
 #RANDOM FOREST
 lalonde.formu
 rf_model = randomForest(formula = as.factor(treat) ~ age + educ + black + hisp + married + nodegr + re74 + 
                           re75, data=lalonde, type="classification")
 
 varImpPlot(rf_model)
 
 rf_ps =   as.numeric(predict(rf_model, type="prob")[,2])
 
 rr.att.RF <- Match(Y=Y, Tr=Tr, X=rf_ps, M=1, ties=FALSE, replace=FALSE, estimand='ATT')
 summary(rr.att.RF)
 
 
hdf_c355=readRDS("D:/icthealth/web/award1807/match/GD_noaward.rds")
 #hdf_noaward<-hdf_c355[,c(hot, thank, gift, contr, views, article, patient1)]
# which( colnames(hdf_c355)=="article" )
hdf_noaward<-hdf_c355[,c(1,11, 12, 13, 14, 15, 18, 19,34)]
hdf_noaward$treat<-0
write.csv(hdf_noaward,"file:///C:/Users/Administrator/Desktop/hdf_noaward.csv")  
summary(hdf_noaward[hdf_noaward$period==1,c(2:9)])


hdf_t355=read.table("clipboard",header = TRUE)
#hdf_noaward<-hdf_c355[,c(hot, thank, gift, contr, views, article, patient1)]
# which( colnames(hdf_c355)=="article" )
summary(hdf_t355[hdf_t355$period==1,c(2:8)])

hdf_awnoaw<-rbind(hdf_t355,hdf_noaward)



## Using logistic regression for estimating propensity scores 
hdf_awnoaw.formu <- treat ~ hot+ thank+gift+contr+article
hdf_awnoaw.glm <- glm(hdf_awnoaw.formu, family=binomial(), data=hdf_awnoaw)
summary(hdf_awnoaw.glm)

hdf_awnoaw.formu2 <- treat ~  hot + thank + gift
hdf_awnoaw.glm2 <- glm(hdf_awnoaw.formu2, data=hdf_awnoaw, family=binomial)
summary(hdf_awnoaw.glm2)

hdf_awnoaw_df <- data.frame(hdf_awnoaw_ps1=fitted(hdf_awnoaw.glm), hdf_awnoaw_ps2=fitted(hdf_awnoaw.glm2), hdf_awnoaw_Tr=hdf_awnoaw$treat)


plot_ps1<-ggplot(hdf_awnoaw_df, aes(x=hdf_awnoaw_ps1)) + geom_histogram() + facet_wrap(~ hdf_awnoaw_Tr, ncol=1)
plot_ps1+ labs(colour = "Cylinders")+ labs(x = "Propensity Score")+ xlim(0,1)


plot_ps2<-ggplot(hdf_awnoaw_df, aes(x=hdf_awnoaw_ps2)) + geom_histogram() + facet_wrap(~ hdf_awnoaw_Tr, ncol=1)
plot_ps2+ labs(colour = "Cylinders")+ labs(x = "Propensity Score")+ xlim(0,1)

summary(hdf_awnoaw.glm)
#summary(hdf_awnoaw_ps1)


#########################Profiles

#hdf_awnoaw_prof<-read.table("clipboard",header = TRUE)
## Using logistic regression for estimating propensity scores 
hdf_awnoaw_prof.formu <- treat ~ hot+ thank+gift+contr+article
hdf_awnoaw_prof.glm <- glm(hdf_awnoaw_prof.formu, family=binomial(), data=hdf_awnoaw_prof)
summary(hdf_awnoaw_prof.glm)

hdf_awnoaw_prof.formu2 <- treat ~  hot + thank + gift
hdf_awnoaw_prof.glm2 <- glm(hdf_awnoaw_prof.formu2, data=hdf_awnoaw_prof, family=binomial)
summary(hdf_awnoaw_prof.glm2)

hdf_awnoaw_prof_df <- data.frame(hdf_awnoaw_prof_ps1=fitted(hdf_awnoaw_prof.glm), hdf_awnoaw_prof_ps2=fitted(hdf_awnoaw_prof.glm2), hdf_awnoaw_prof_Tr=hdf_awnoaw_prof$treat)

plot_prof_ps1<-ggplot(hdf_awnoaw_prof_df, aes(x=hdf_awnoaw_prof_ps1)) + geom_histogram() + facet_wrap(~ hdf_awnoaw_prof_Tr, ncol=1)
plot_prof_ps1+ labs(colour = "Cylinders")+ labs(x = "Propensity Score")+ xlim(0,1)

plot_prof_ps2<-ggplot(hdf_awnoaw_prof_df, aes(x=hdf_awnoaw_prof_ps2)) + geom_histogram() + facet_wrap(~ hdf_awnoaw_prof_Tr, ncol=1)
plot_prof_ps2+ labs(colour = "Cylinders")+ labs(x = "Propensity Score")+ xlim(0,1)

summary(hdf_awnoaw_prof.glm)


par(mfrow = c(1, 2))
plot(plot_ps1+ labs(x = "Propensity Score")+ xlim(0,1))
plot(plot_prof_ps1+  labs(x = "Propensity Score")+ xlim(0,1))



# T-test on the original data
t.test(views ~ treat, data =hdf_awnoaw_prof )
t.test(patient1 ~ treat, data =hdf_awnoaw_prof)

# Regression on the original data
hdf_awnoaw_prof.formu <- treat ~ hot+ thank+gift+contr+article
hdf_awnoaw_prof.glm <- glm(hdf_awnoaw.formu, family=binomial(), data=hdf_awnoaw_prof)
summary(hdf_awnoaw_prof.glm)


# Run optimal matching
hdf_awnoaw_prof.out <- matchit(hdf_awnoaw_prof.formu, 
                          data = hdf_awnoaw_prof, method = "optimal", distance = "logit", ratio = 1)

# Checking Balance 
summary(hdf_awnoaw_prof.out)
plot(hdf_awnoaw_prof.out, type="jitter")
plot(hdf_awnoaw_prof.out, type="hist")
plot(hdf_awnoaw_prof.out, type="QQ")
plot(summary(hdf_awnoaw_prof.out, standardize = T), interactive = T)


############
hdf_tcpanel<-read.table("clipboard",header = TRUE)
head(hdf_tcpanel)
#hdf_tcpanel_lmv<-lm(Views~DRR+Thank+Gift+Contri+Article+D, data=hdf_tcpanel)
hdf_tcpanel_lmv<-lm(Views~Thank+Gift+Contri+D, data=hdf_tcpanel)
summary(hdf_tcpanel_lmv)

#hdf_tcpanel<-data.frame(hdf_tcpanel)
hdf_tcpanel_lmp<-lm(Patient ~DRR+Thank+Gift+Contri+Article+D, data=hdf_tcpanel)
#hdf_tcpanel_lmp<-lm(Patient ~Thank+Gift+Contri+D, data=hdf_tcpanel)
summary(hdf_tcpanel_lmp)
