#set working directory
setwd("")

# install.packages("aTSA")
# install.packages("tidyverse")
library("zoo")
library("xts")
library("ggplot2")
library("strucchange")
library("tidyverse")
library("urca")
library("lmtest")
library("forecast")

##############################################################
# EP processing
##############################################################
# load csv file, name it in R
OMXS30<-read.csv("_SE0000337842_2023-02-24.csv")
FXSwedenBIS<-read.csv("WS_EER_D_Export.csv")

# storge structure
class(OMXS30$Date)
class(FXSwedenBIS$Date)

#change to date
OMXS30$Date<-as.Date(OMXS30$Date, "%d/%m/%Y")
FXSwedenBIS$Date<-as.Date(FXSwedenBIS$Date,"%d/%m/%Y")
FXSwedenBIS <- na.omit(FXSwedenBIS) 

# Generate Subsample
OMXS30<-OMXS30[OMXS30$Date > "1993-01-04" & OMXS30$Date < "1997-12-30",]
FXSwedenBIS<-FXSwedenBIS[FXSwedenBIS$Date > "1993-01-04" & FXSwedenBIS$Date < "1997-12-30",]

#change to numeric
OMXS30$Closingprice <- as.numeric(as.character(OMXS30$Closingprice))
FXSwedenBIS$Narrow <- as.numeric(as.character(FXSwedenBIS$Narrow))

#xts: extensible time-series object.
OMXS30 <- xts(OMXS30$Closingprice, OMXS30$Date)
names(OMXS30)[1]<-"Close"

#xts: extensible time-series object.
FXSwedenBIS <- xts(FXSwedenBIS$Narrow, FXSwedenBIS$Date)
names(FXSwedenBIS)[1]<-"NomFX"

#############################################################
#merge the two EP frames
#############################################################
EP <- merge (OMXS30, FXSwedenBIS, join = "right", fill = NA)
EP <- na.omit(EP)

#############################################################
# time series properties of the EP
#############################################################

#check the time series plot to decide unit root test model
plot(index(EP),EP$Close, type="line")
plot(index(EP),EP$NomFX, type="line")

plot(log(EP$Close))
plot(log(EP$NomFX))



#plot with ggplot
p1<-ggplot(data=EP, aes(x=index(EP)))+
  geom_line(aes(y=log(OMXS30), color="OMSX30"), linetype="solid", color="black")
p1<-p1+theme(
             axis.text.x=element_text(color="black",size=12),
             axis.text.y=element_text(color="black", size=12),
             panel.background=element_rect(fill="white"),
             panel.grid.minor=element_blank(),
             panel.grid.major=element_blank(),
             axis.line=element_line(color="black", size=1),
             legend.position="none",
             plot.title=element_text(face="bold"))+ggeasy::easy_center_title()
p1<-p1+
  annotate("text", x=as.Date("1997-06-30"),y=6, label="Closing Price", fontface="bold",color="black")
p1<-p1+ xlab("Date")+ylab("Log Value")

# + labs(title="OMSX30")
p1
# ggsave("OMXS30.png")

p2<-ggplot(data=EP, aes(x=index(EP)))+
  geom_line(aes(y=log(NomFX),color="NomFX"), linetype="solid", color="blue")
p2<-p2+theme(
  axis.text.x=element_text(color="black",size=12),
  axis.text.y=element_text(color="black", size=12),
  panel.background=element_rect(fill="white"),
  panel.grid.minor=element_blank(),
  panel.grid.major=element_blank(),
  axis.line=element_line(color="black", size=1),
  legend.position="none",
  plot.title=element_text(face="bold"))+ggeasy::easy_center_title()
p2<-p2+annotate("text", x=as.Date("1997-06-30"),y=4.85, label="Nominal Exchange Rate", fontface="bold",color="blue")
p2<-p2+ xlab("Date")+ylab("Log Value")
# +labs(title="Nominal FX by CPI (BIS), Sweden")
p2
# ggsave("NomFX.png")

#two time series so specifying x axis in aes and then adding the two variables using geom_line
p3<-ggplot(data=EP, aes(x=index(EP)))+
  geom_line(aes(y=log(OMXS30), color="OMSX30"), linetype="solid", color="black")+
  geom_line(aes(y=log(NomFX),color="NomFX"), linetype="solid", color="blue")
p3<-p3+theme(
             axis.text.x=element_text(color="black",size=12),
             axis.text.y=element_text(color="black", size=12),
             panel.background=element_rect(fill="white"),
             panel.grid.minor=element_blank(),
             panel.grid.major=element_blank(),
             axis.line=element_line(color="black", size=1),
             legend.position="none",
             plot.title=element_text(face="bold"))+ggeasy::easy_center_title()

p3<-p3+annotate("text", x=as.Date("1997-06-30"),y=6, label="Closing Price", fontface="bold",color="black")
p3<-p3+annotate("text", x=as.Date("1997-06-30"),y=5, label="Nominal Exchange Rate", fontface="bold",color="blue")
p3<-p3+ xlab("Date")+ylab("Log Value")
# + labs(title="OMSX30 & Nominal FX by CPI (BIS), Sweden")
p3
# ggsave("OMXFX.png")

#scatter plot of any relation
ggplot(EP,aes(x=EP$Close,y=EP$NomFX))+geom_point() + xlab("Close Price")+ylab("Nominal Exchange Rate")
# ggsave("Scatter.png")


EP <- log(EP)
saveRDS(EP, file = "Coursework.rds")
# temp<-data.frame(EP)
# write.csv(temp, file = "coursework.csv")

#lets consider EP avoiding the known break-points of joining EU.
#we can check for further breakpoints roughly using the
#Zivot-Andrews unit root test. 

#Determine Breakpoints
#Autoregressive order 1 model (1 lag) on NomFX
dat <- tibble(NomFXlag0 = EP$NomFX,
              NomFXlag1 = lag(EP$NomFX)
)
drop_na()
qlr <- Fstats(NomFXlag0 ~ NomFXlag1, data = dat)
breakpoints(qlr)
sctest(qlr, type = "supF")


z_trend<-ur.za(EP$NomFX, model="both", lag=1)
summary(z_trend)

EPpre <- EP[0:670,]
EPpost <- EP[671:1254,]

# =====================================================
# =====================================================
# Post 671
# =====================================================
# =====================================================
######################################################
#Task 0: (b) Time series properties of the variables
######################################################

#plots to decide with or without trend for each variable's unit root test specification
plot(index(EPpost),EPpost$Close, type="l")
plot(index(EPpost),EPpost$NomFX, type="l")

# Unit root tests with trend
#step1 of TS versus DS

OMXS30post_adf <- ur.df(EPpost$Close, type="trend",selectlags="AIC")
summary(OMXS30post_adf)
# nonstationary
NomFXpost_adf <- ur.df(EPpost$NomFX, type="trend", selectlags="AIC")
summary(NomFXpost_adf)

#step2 of TS versus DS: regression on trend only (follow lab 2)
dOMXS30post <- na.omit(diff(EPpost$Close, differences=1))
dNomFXpost <- na.omit(diff(EPpost$NomFX, differences=1))
LOMXS30post <- na.omit(lag(EPpost$Close, n=1))
LNomFXpost <- na.omit(lag(EPpost$NomFX, n=1))
LDOMXS30post <- na.omit(lag(dOMXS30post, n=1))
LDNomFXpost <- na.omit(lag(dNomFXpost, n=1))
ttpost <- cbind(dOMXS30post, dNomFXpost, LOMXS30post, LNomFXpost, LDOMXS30post, LDNomFXpost)
ttpost <- na.omit(ttpost)

OMXS30post_adf <- lm(ttpost$Close ~ c(1:length(ttpost$Close))+ttpost$Close.2)
summary(OMXS30post_adf)
NomFXpost_adf <- lm(ttpost$NomFX ~ c(1:length(ttpost$NomFX))+ttpost$NomFX.2)
summary(NomFXpost_adf)


#step 3 of TS versus DS: ADF without trend
OMXS30post_adf <- ur.df(EPpost$Close, type="drift", selectlags="AIC")
summary(OMXS30post_adf)
NomFXpost_adf <- ur.df(EPpost$NomFX, type="drift", selectlags="AIC")
summary(NomFXpost_adf)

#Step 4 of TS versus DS: DF without constant
OMXS30post_adf <- ur.df(EPpost$Close, type="none")
summary(OMXS30post_adf)
NomFXpost_adf <- ur.df(EPpost$NomFX, type="none")
summary(NomFXpost_adf)

##################################################
#Task 0: (c) order of integration
##################################################
#Check if the variables are first difference stationary (I(1))
OMXS30post_adf <- ur.df(ttpost$Close, type="drift", selectlags="AIC")
summary(OMXS30post_adf)
NomFXpost_adf <- ur.df(ttpost$NomFX, type="drift", selectlags="AIC")
summary(NomFXpost_adf)
# All I(1)
Test<-na.omit(diff(EP$NomFX, differences = 2))
NomFXpost_adf <- ur.df(Test, type="drift", selectlags="AIC")
summary(NomFXpost_adf)
######################################################
#Task 1a: Engle-Granger first step
######################################################
##some useful graphs
ggplot(EPpost, aes(x=index(EPpost)))+
  geom_line(aes(y=EPpost$Close, color="Close"), linetype="solid", color="black")+
  geom_line(aes(y=EPpost$NomFX, color="NomFX"), linetype="solid", color="blue")



#scattposter plot of any relation
ggplot(EPpost,aes(x=EPpost$Close,y=EPpost$NomFX))+geom_point() + xlab("Close Price")+ylab("Nominal Exchange Rate")
# ggsave("ScatterPost.png")
ggplot(EPpre,aes(x=EPpre$Close,y=EPpre$NomFX))+geom_point() + xlab("Close Price")+ylab("Nominal Exchange Rate")
# ggsave("ScatterPre.png")

# OLS regression
EG1 <- lm(EPpost$NomFX ~ EPpost$Close + c(1:length(EPpost$NomFX)))
summary(EG1)
###########################################################
#durbin-Watson test of serial correlation in bi-variate case
# doesn't work with xts object
ttpost1 <- as.numeric(EPpost$NomFX)
ttpost2 <- as.numeric(EPpost$Close)
EGttpost <- lm(ttpost1 ~ ttpost2 + c(1:length(ttpost1)))

dwtest(EGttpost)
###########################################################

# plot residuals
EPpost$res <- EG1$residuals
plot(EPpost$res)

# Residual stationary?
ADF_res <- ur.df(EG1$residuals, type="none", selectlags="AIC")
summary(ADF_res)


##############################################################
# Fitting ARMA as the mean model 
##############################################################
# NomFX price returns
data <- cbind(EPpost$Close,EPpost$NomFX)
data$rClose <- diff(data$Close)
data$rNomFX <- diff(data$NomFX)
data <- na.omit(data)

acf(data$rNomFX, lag.max = "10")
pacf(data$rNomFX, lag.max = "10")
# ARMA (0,1)

NomFX01 <- arima(data$rNomFX, c(0,0,1), method= "ML")
NomFX01
NomFX00 <- auto.arima(data$rNomFX, stepwise = FALSE, method= "ML", ic=c("aic"))
NomFX00

lrtest (NomFX01, NomFX00)

# Close returns
acf(data$rClose, lag.max = "10")
pacf(data$rClose, lag.max = "10")

Close01 <- arima(data$rClose, c(0,0,1), method= "ML")
Close01
Close00 <- auto.arima(data$rClose, stepwise = FALSE, method= "ML", ic=c("aic"))
Close00

lrtest (Close00, Close01)

#############################################################

checkresiduals(arima(data$rClose, c(0,0,0), method= "ML"))

library ("aTSA")
##ARCH LM test: does the arma models have ARCH effect in the residuals?
arch.test(arima(data$rNomFX, c(0,0,0), method= "ML"))
arch.test(arima(data$rClose, c(0,0,0), method= "ML"))


##############################################################



# =====================================================
# =====================================================
# pre 671
# =====================================================
# =====================================================
######################################################
#Task 0: (b) Time series properties of the variables
######################################################

#plots to decide with or without trend for each variable's unit root test specification
plot(index(EPpre),EPpre$Close, type="l")
plot(index(EPpre),EPpre$NomFX, type="l")

# Unit root tests with trend
#step1 of TS versus DS

OMXS30pre_adf <- ur.df(EPpre$Close, type="trend",selectlags="AIC")
summary(OMXS30pre_adf)
# nonstationary
NomFXpre_adf <- ur.df(EPpre$NomFX, type="trend", selectlags="AIC")
summary(NomFXpre_adf)
# Stationary

#step2 of TS versus DS: regression on trend only (follow lab 2)
dOMXS30pre <- na.omit(diff(EPpre$Close, differences=1))
LOMXS30pre <- na.omit(lag(EPpre$Close, n=1))
LDOMXS30pre <- na.omit(lag(dOMXS30pre, n=1))
ttpre <- cbind(dOMXS30pre, LOMXS30pre, LDOMXS30pre)
ttpre <- na.omit(ttpre)

OMXS30pre_adf <- lm(ttpre$Close ~ c(1:length(ttpre$Close))+ttpre$Close.2)
summary(OMXS30pre_adf)

#step 3 of TS versus DS: ADF without trend
OMXS30pre_adf <- ur.df(EPpre$Close, type="drift", selectlags="AIC")
summary(OMXS30pre_adf)

#Step 4 of TS versus DS: DF without constant
OMXS30pre_adf <- ur.df(EPpre$Close, type="none")
summary(OMXS30pre_adf)

#========================================================
# Detrend the data sample for TS only (OMXS30pre)
#========================================================
dtrdOMXS30pre <- as.xts(residuals(lm(EPpre$Close ~ c(1:length(EPpre$Close)))))
plot(dtrdOMXS30pre)
names(dtrdOMXS30pre)[1]<-"dtrdClose"
EPpre<-cbind(EPpre,dtrdOMXS30pre)

#Generate the first difference of the series.
FD<-diff(dtrdOMXS30pre, differences=1) # first difference
names(FD)[1]<-"dtrdClose"
FD<-na.omit(FD)
plot(FD)

summary(ur.df(FD, type = "drift", selectlags = "AIC"))


###########################################################
#durbin-Watson test of serial correlation in bi-variate case
# doesn't work with xts object
ttpre1 <- as.numeric(EPpre$NomFX)
ttpre2 <- as.numeric(EPpre$Close)
EGttpost <- lm(ttpost1 ~ ttpre2 + c(1:length(ttpre1)))

dwtest(EGttpre)
###########################################################

##############################################################
# Fitting ARMA as the mean model 
##############################################################
# NomFX price returns
predata <- cbind(EPpre$Close,EPpre$NomFX)
predata$rClose <- diff(predata$Close)
predata$rNomFX <- diff(predata$NomFX)
predata <- na.omit(predata)

acf(predata$rNomFX, lag.max = "10")
pacf(predata$rNomFX, lag.max = "10")
# ARMA (2,1)

NomFX23 <- arima(predata$rNomFX, c(2,0,3), method= "ML")
NomFX23
NomFX22 <- auto.arima(predata$NomFX, stepwise = FALSE, method= "ML", ic=c("aic"))
NomFX22 <- arima(predata$rNomFX, c(2,0,2), method= "ML")

lrtest (NomFX22, NomFX23)

# Close returns
acf(predata$rClose, lag.max = "10")
pacf(predata$rClose, lag.max = "10")

Close11 <- arima(predata$rClose, c(1,0,1), method= "ML")
Close11
Close10 <- auto.arima(predata$rClose, stepwise = FALSE, method= "ML", ic=c("aic"))
Close10

lrtest (Close10, Close11)
# chosen model arma(1,0)
#############################################################

library ("aTSA")

##ARCH LM test: does the arma models have ARCH effect in the residuals?
arch.test(arima(predata$rNomFX, c(2,0,2), method= "ML"))
arch.test(arima(predata$rClose, c(1,0,0), method= "ML"))

# we reject the null and so there is arch effect and ARMA is not the best model
# for this data

##############################################################

