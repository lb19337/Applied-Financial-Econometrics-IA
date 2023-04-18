#set working directory
setwd("")
library("xts")
library("zoo")
library("urca")
library("ggplot2")
library("vars")
library("tsDyn")
library("coefplot")
library("forecast")
library("rugarch")
library("rmgarch")
library("ARDL")

EP <- readRDS("Coursework.rds")

EPpre <- EP[0:670,]
EPpost <- EP[671:1254,]

dtrdOMXS30pre <- as.xts(residuals(lm(EPpre$Close ~ c(1:length(EPpre$Close)))))
names(dtrdOMXS30pre)[1]<-"dtrdClose"
EPpre<-cbind(EPpre,dtrdOMXS30pre)

# =====================================================
# =====================================================
# Post 671
# =====================================================
# =====================================================
######################################################
# Engle-Granger second step: the error correction model
######################################################
# calculate the returns series: equivalent to first difference, so stationary
# and has general interest.

EPpost$rNomFX <- diff(EPpost$NomFX)
EPpost$rClose <- diff(EPpost$Close)
#saveRDS(EP, file = "Lab4.rds")
#
EG1 <- lm(EPpost$NomFX ~ EPpost$Close + c(1:length(EPpost$NomFX)))
summary(EG1)
#the short run regression with stationary variables
#general to specific methodology: starting from T^0.25 = 3 lags for each variable

EG2 <- lm(EPpost$rNomFX ~ EPpost$rClose + lag(EPpost$rClose, n=1) + lag(EPpost$rNomFX, n=1)+
             lag(EG1$residuals))
summary (EG2)
AIC(EG2)


######################################################
#Task 1: Finding the number of cointegrating vectors
######################################################

data <- cbind(EPpost$Close,EPpost$NomFX)
data_EG <- cbind(EPpost$Close,EPpost$NomFX)
# Select the lag order: Note that this is using level variables, not their returns
VARselect(data, lag.max=12, type = "both")
# Johansen test of cointegrating vectors
# We are using a linear trend specification in the cointegrating equation
# lags: K=2 (Should be 1)
# we are using the trace test
# specification is 'transitory' as the VECM coefficients will estimate the short run effects

# Case 4: restricted constant
Jtest4 <- ca.jo(data, type = "trace", ecdet = "const", K = 2, spec = "transitory")
summary (Jtest4)

# Case 3: unrestricted constant (preferable)
Jtest3 <- ca.jo(data, type = "trace", ecdet = "none", K =2, spec = "transitory")
summary (Jtest3)

# case 2: restricted trend
Jtest2 <- ca.jo(data, type = "trace", ecdet = "trend", K = 2, spec = "transitory")
summary (Jtest2)


######################################################
# Task 2: The estimation
######################################################
# one cointegrating relation, r=1
#######################

## VECM specification: LR equation normalization
vecmA <- cajorls(Jtest2, r = 1)
vecmA

res <- ts(vecmA$rlm$residuals)
plot(res,type="l")

###################################
## Note: there is no standard errors reported. Use VECM function with tsDyn

vecmA <- VECM(data, lag=1, r=1, include= c("const"), estim=c("ML"), LRinclude = c("trend"))
summary(vecmA)
vecm4 <- VECM(data, lag=1, r=1, include= c("none"), estim=c("ML"), LRinclude = c("const"))
summary(vecm4)
# present the ECM equations
coefs_all <- summary(vecmA)$coefMat
coefs_all
# check to confirm if the estimates are the same
beta <- cajorls(Jtest2, r = 1)$beta
all.equal(beta,
          coefB(vecmA), check.attributes = FALSE)

######################################################
# Vector Autogression
######################################################

######################################################
#Task 2a: Vector Autoregression (VAR)?
#(given the contradictory results regarding cointegration)
######################################################
data$rClose <- diff(data$Close)
data$rNomFX <- diff(data$NomFX)

# Diffferenced data for VAR
var_data <- cbind(data$rClose, data$rNomFX)
var_data <- na.omit(var_data)

# Select VAR lag order
VARselect(var_data, lag.max=12, type = "both")

# Estimate VAR
VAR_est <- VAR(y = var_data, p=1 ,lag.max=12, ic="AIC", type = c("const"))

# VAR regression output
coef(VAR_est$varresult$rClose)

coefplot(VAR_est$varresult$rClose)

######################################################
# Task 4: Interpretation
######################################################
################## Model (A)  #####

VARA <- vec2var(Jtest2, r=1)

#Prediction (forecast)
ttA <- predict(VARA, n.ahead = 12, ci = 0.95, dumvar = NULL)

#Forecast visualization
fanchart(ttA)

#The impulse response function (IRF) (slightly differs from Stata)
irfA <- irf(VARA, impulse="NomFX", response=c("Close"))
plot(irfA)

#The variance decomposition
fevd(VARA, n.ahead=12)

# ######################################################
# ## some optional graph code to decorate your irf!:
# irf1.prop<-ts(irfA$irf$NomFX)
# irf1.prop.l<-irfA$Lower$NomFX
# irf1.prop.u<-irfA$Upper$NomFX
# irf1.prop<-as.numeric(irf1.prop)
# 
# timet<-1:11
# df.irf<-data.frame(irf1.prop,irf1.prop.l,irf1.prop.u)
# df.irf<- df.irf
# df.irf<- data.frame(timet,df.irf)
# 
# 
# ggirf1.1<-ggplot(df.irf,aes(timet)) + xlab("Months") + ylab("% Change in interest rate")
# ggirf1.2<-ggirf1.1 +
#   geom_line(aes(y = irf1.prop, colour = "Close")) +
#   geom_line(aes(y = irf1.prop.l, colour = "95% CI")) +
#   geom_line(aes(y = irf1.prop.u, colour = "95% CI")) +
#   scale_colour_hue("Series:")+
#   theme(panel.background=element_rect(fill="white"), axis.line=element_line(color="black", size=1))+
#   ggtitle("Impulse variable NomFX, response variable Close")+
#   ggeasy::easy_center_title()
# ggirf1.2



######################################################
#diagnostic checks
######################################################

#any serial correlation? Type PT.adjusted for small sample
serial.test(VARA,lags.bg =1, type = "BG")

# Errors look normal?
normality.test(VARA, multivariate.only = TRUE)

# any arch effect?
arch.test(VARA,lags.multi=1, multivariate.only=TRUE)

#There is ARCH effect

data <- na.omit(data)
######################################################
#MGARCH
######################################################
###############################################################
# Task 1: Fitting GARCH(1,1) models for Oil price returns and S&P500 returns
###############################################################

garch_spec <- ugarchspec(variance.model=list(model="eGARCH", garchOrder= c(1,1)),
                         mean.model= list(armaOrder=c(0,0)))
garch_rClose <- ugarchfit(garch_spec,data=data$rClose)
garch_rNomFX <- ugarchfit(garch_spec,data=data$rNomFX)
names(garch_rClose@fit)
ht_rClose <- garch_rClose@fit$var
ht_rNomFX <- garch_rNomFX@fit$var
var_data <- cbind(ht_rClose, ht_rNomFX)


# Before proceeding to estimate VAR model of estimated volatility, we’ll need
# to check whether the predicted values are stationary or not (why?).

plot(ht_rClose, type="l")
plot(ht_rNomFX, type="l")

ADF_htClose <- ur.df(ht_rClose, type="trend", selectlags="AIC")
summary(ADF_htClose)
ADF_hNomFX <- ur.df(ht_rNomFX, type="trend", selectlags="AIC")
summary(ADF_hNomFX)

# #step2 of TS versus DS: regression on trend only (follow lab 2)
d <- na.omit(diff(ht_rNomFX, differences=1))
# L <- na.omit(lag(ht_rNomFX, n=1))
# LD <- na.omit(lag(d, n=1))
# tt <- cbind(d, L, LD)
# tt <- na.omit(tt)
# 
# ADF_hNomFX <- lm(tt$d ~ c(1:length(tt$d))+ttpre$Close.LD)
# summary(ADF_hNomFX)
# # trend not significant, need step 3
# 
# #step 3 of TS versus DS: ADF without trend
# ADF_hNomFX <- ur.df(EPpre$Close, type="drift", selectlags="AIC")
# summary(ADF_hNomFX)
# 
# #Step 4 of TS versus DS: DF without constant
# ADF_hNomFX <- ur.df(EPpre$Close, type="none")
# summary(ADF_hNomFX)

ADF_hNomFX <- ur.df(d, type="drift", selectlags="AIC")
summary(ADF_hNomFX)

VARselect(var_data, lag.max=30, type = "both")
VARvol <- VAR(y = var_data,lag.max=12, ic="AIC")


######################################################
# Task 4: Interpretation Garch Var Post
######################################################
################## Model (A)  #####

#Prediction (forecast)
ttA <- predict(VARvol, n.ahead = 12, ci = 0.95, dumvar = NULL)

#Forecast visualization
fanchart(ttA)

#The impulse response function (IRF) (slightly differs from Stata)
irfA <- irf(VARvol, impulse="ht_rNomFX", response=c("ht_rClose"))
plot(irfA)

irf <- irf(VARvol, impulse="ht_rClose", response=c("ht_rNomFX"))
plot(irf)


#The variance decomposition
fevd(VARvol, n.ahead=12)
causality(VARvol, cause="ht_rNomFX")$Granger
causality(VARvol, cause="ht_rClose")$Granger

###############################################################
# Task 2: Dynamic Conditional Correlations (DCC) GARCH model
###############################################################

rdata <- cbind(data$rClose, data$rNomFX)
uspec.n = multispec(replicate(2, ugarchspec(mean.model = list(armaOrder = c(1,0)))))
multf = multifit(uspec.n, rdata)
spec1 = dccspec(uspec = uspec.n, dccOrder = c(1, 1), distribution = 'mvnorm')
fit1 = dccfit(spec1, data = rdata, fit.control = list(eval.se = TRUE), fit = multf)
fit1

# Test of non-constant correlation based on Engle III and Sheppard (2001):
DCCtest(rdata, garchOrder = c(1,1), n.lags = 100, solver = "solnp",
        solver.control = list(), cluster = NULL, Z = NULL)
#forecast
dccforecast(fit1, n.ahead=10)
# Get the model based time varying covariance (arrays) and correlation matrices
cov1 = rcov(fit1) # extracts the covariance matrix within sample
cor1 = rcor(fit1) # extracts the correlation matrix within sample
dim(cor1) #dimension of dynamic conditional correlation
cor1[,,dim(cor1)[2]]
corr_rsp_rClose=cor1[2,1,450:583] # row 2, column 1, observations between 1100 to 1227
plot(corr_rsp_rClose, type="l")
cov1[,,dim(cov1)[2]]
cov_rsp_rClose=cov1[2,1,450:583] # row 2, column 1, observations between 1100 to 1227
plot(cov_rsp_rClose, type="l")

# # =====================================================
# # =====================================================
# # Pre 671
# # =====================================================
# # =====================================================


######################################################
# Vector Autogression Pre
######################################################

######################################################
#Task 2a: Vector Autoregression (VAR)? Pre
######################################################
predata <- cbind(EPpre$dtrdClose,EPpre$NomFX)
predata$rNomFX <- diff(predata$NomFX)
predata$rdtrd <- diff(predata$dtrdClose)

# Diffferenced data for VAR
var_predata <- cbind(predata$dtrdClose, predata$NomFX)
var_predata <- na.omit(var_predata)

# Select VAR lag order
VARselect(var_predata, lag.max=12, type = "both")

# Estimate VAR
VAR_est <- VAR(y = var_predata, p=2 ,lag.max=12, ic="AIC", type = c("const"))

# VAR regression output
coef(VAR_est$varresult$dtrdClose)

coefplot(VAR_est$varresult$dtrdClose)

######################################################
#Task 2a: ARDL Model Pre
######################################################

ARDLdata<-data.frame(EPpre)
ARDLdata<-ARDLdata[-3]
models <- auto_ardl(Close ~ NomFX, data = ARDLdata, max_order = 3)
models$top_orders
modelstrd<- auto_ardl(Close ~ NomFX + trend(Close), data = ARDLdata, max_order = 3)
modelstrd$top_orders

ardl_11 <- models$best_model
summary(ardl_11)

ardl_11_tr <- modelstrd$best_model
summary(ardl_11_tr)

ardl_11_tr <- ardl(Close ~ NomFX + trend(Close),
                   data = ARDLdata, order = c(1,1))
summary(ardl_11_tr)

ardl_11_tr$coefficients

######################################################
# Task 3b: Interpreting VAR output: Pre
######################################################
#Prediction (forecast)
ttpre<- predict(VAR_est, n.ahead = 12, ci = 0.95, dumvar = NULL)

#Forecast visualization
fanchart(ttpre)

#The impulse response function (IRF)
irf <- irf(VAR_est, impulse="NomFX", response=c("dtrdClose"))
plot(irf)


#The variance decomposition
vardec <- fevd(VAR_est, n.ahead=12)
vardec

#Granger Causality test
causality(VAR_est, cause="NomFX")$Granger
causality(VAR_est, cause="dtrdClose")$Granger

######################################################
# Task 3c: Diagnostic checks Pre
######################################################
#any serial correlation? Type PT.adjusted for small sample
serial.test(VAR_est,lags.bg = 1, type="BG")  # LM test

# Errors look normal?
normality.test(VAR_est, multivariate.only = TRUE)

# Parameter Stability test?
VAR_estStabil <- stability(VAR_est)
plot(VAR_estStabil)

arch.test(VAR_est,lags.multi=1, multivariate.only=TRUE)
#There is ARCH effect


######################################################
#MGARCH Pre
######################################################

predatagarch <- EPpre[,-3]
predatagarch$rClose <- diff(predatagarch$Close)
predatagarch$rNomFX <- diff(predatagarch$NomFX)
predatagarch <- na.omit(predatagarch)
###############################################################
# Task 1: Fitting GARCH(1,1) models for Oil price returns and S&P500 returns
###############################################################

garch_spec1 <- ugarchspec(variance.model=list(model="eGARCH", garchOrder= c(1,1)),
                         mean.model= list(armaOrder=c(1,0)))
garch_spec2 <- ugarchspec(variance.model=list(model="eGARCH", garchOrder= c(1,1)),
                         mean.model= list(armaOrder=c(2,2)))
garch_rClose <- ugarchfit(garch_spec1,data=predatagarch$rClose)
garch_rNomFX <- ugarchfit(garch_spec2,data=predatagarch$rNomFX)
names(garch_rClose@fit)
ht_rClose <- garch_rClose@fit$var
ht_rNomFX <- garch_rNomFX@fit$var
var_data <- cbind(ht_rClose, ht_rNomFX)


# Before proceeding to estimate VAR model of estimated volatility, we’ll need
# to check whether the predicted values are stationary or not (why?).

plot(ht_rClose, type="l")
plot(ht_rNomFX, type="l")


ADF_htClose <- ur.df(ht_rClose, type="trend", selectlags="AIC")
summary(ADF_htClose)
ADF_hNomFX <- ur.df(ht_rNomFX, type="trend", selectlags="AIC")
summary(ADF_hNomFX)

# Estimate the VAR model of the two volatility estimates, test using Granger
# causality test, interpret using impulse response function, and check for model
# diagnostics. Comment on the hedging capability 

VARselect(var_data, lag.max=30, type = "both")
VAR_est <- VAR(y = var_data,lag.max=12, ic="AIC")
#Granger Causality test
causality(VAR_est, cause="ht_rNomFX")$Granger
causality(VAR_est, cause="ht_rClose")$Granger
#The impulse response function (IRF)
irf1 <- irf(VAR_est, impulse="ht_rNomFX", response=c("ht_rClose"))
plot(irf1)
irfA <- irf(VAR_est, impulse="ht_rClose", response=c("ht_rNomFX"))
plot(irfA)


#Diagnostics
serial.test(VAR_est,lags.pt = 1, type = "BG")
normality.test(VAR_est, multivariate.only = TRUE)
VAR_estStabil <- stability(VAR_est)
plot(VAR_estStabil, type="l")
arch.test(VAR_est,lags.multi=1, multivariate.only=TRUE)

###############################################################
# Task 2: Dynamic Conditional Correlations (DCC) GARCH model
###############################################################
rdata<- cbind(predatagarch$rClose, predatagarch$rNomFX)
uspec.n = multispec(replicate(2, ugarchspec(mean.model = list(armaOrder = c(1,0)))))
multf = multifit(uspec.n, rdata)
spec1 = dccspec(uspec = uspec.n, dccOrder = c(1, 1), distribution = 'mvnorm')
fit1 = dccfit(spec1, data = rdata, fit.control = list(eval.se = TRUE), fit = multf)
fit1

# Test of non-constant correlation based on Engle III and Sheppard (2001):
DCCtest(rdata, garchOrder = c(1,1), n.lags = 50, solver = "solnp",
        solver.control = list(), cluster = NULL, Z = NULL)
#forecast
plot(dccforecast(fit1, n.ahead=10))

# Get the model based time varying covariance (arrays) and correlation matrices
cov1 = rcov(fit1) # extracts the covariance matrix within sample
cor1 = rcor(fit1) # extracts the correlation matrix within sample
dim(cor1) #dimension of dynamic conditional correlation
cor1[,,dim(cor1)[2]]
corr_rsp_rClose=cor1[2,1,500:669] # row 2, column 1, observations between 1100 to 1227
plot(corr_rsp_rClose, type="l")
cov1[,,dim(cov1)[2]]
cov_rsp_rClose=cov1[2,1,500:669]
plot(cov_rsp_rClose, type="l")



