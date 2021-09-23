#Matthew Brigham
#STA 521 
#Final Project - GARCH
#Spring 2021

setwd("/Users/Matt/Library/Mobile Documents/com~apple~CloudDocs/Cleveland State/Courses/STA 521 - Time Series/521 Project/Datasets")


################################
###--------IMPORT DATA-------###
################################

dat <- read.csv("Daily TSLA.csv")
head(dat)
plot(dat$Close, type = "l", main = "Daily Stock Price for Tesla April 18 2016 - April 16 2021", ylab = "Closing Price (USD)")
dat.ts = ts(dat$Close, start = 1, end = nrow(dat))

#Calculate Daily Returns from Closing Stock Price
y = data.frame(NA)
colnames(y) = "Daily Returns"
for (i in 2:nrow(dat)) {
  y[i-1, 1] <- (dat$Close[i] - dat$Close[i-1])/dat$Close[i-1]
}
head(y)

#Create time series
y.ts_complete = ts(y, start = 1, end = nrow(y))
y.ts = ts(y, start = 1, end = nrow(y))

# #Divide data set into Training and Validation sets
# y.ts = y.ts_complete[1:880,]
# y.ts_test = y.ts_complete[881:nrow(y.ts_complete), ] #final 378 returns

################################
###-----DATA EXPLORATION-----###
################################
library(moments)

    #Descriptive Statistics
    summary(y.ts)
    kurtosis(y.ts) #leptokurtic 8.600 > 3 (normal=3) - means long, fat tails and tall/skinny peak
    skewness(y.ts) #0.302 approximately symmetric since between -0.5 and 0.5

    #Time Series Plot of Returns
    plot(y.ts, main = "Daily Returns for Tesla April 19 2016 - April 16 2021", xlab = "Time Index")

    #Time series plot of standardized returns
    s = sd(y.ts)
    mu = mean(y.ts)
    Z = (y.ts - mu)/s
    plot(Z)

    #Normality Check - shows not normal data
    shapiro.test(y.ts) 
    
    plot(density(y.ts*100), main = "Distribution of Returns TSLA", xlab = "Return %")
    curve(dnorm(x, mean=mean(y.ts*100), sd=sd(y.ts*100)), col="red", add = T)
    legend(7, 0.15, legend=c("Returns", "Normal Dist."),
           col=c("black", "red"), lty = 1:1, cex=0.8)
    
    #Correlograms - Stock Price - Not Stationary
    acf(dat.ts , type=c("correlation") , plot=T , main="ACF for Stock Price", lag.max = 50)$acf 
    
    #Correlograms - Returns
    par(mfrow=c(1,2))
    acf(y.ts , type=c("correlation") , plot=T , main="ACF for Returns", lag.max = 50)$acf #stationary, white noise
    pacf(y.ts, lag = 50, main = "PACF for Returns")

    #Correlograms - Squared Returns
    par(mfrow=c(1,2))
    acf(y.ts^2 , type=c("correlation") , plot=T , main="ACF for Squared Returns", lag.max = 50)$acf #not stationary
    pacf(y.ts^2, lag = 50, main = "PACF for Squared Returns")
    
    acf(diff(y.ts^2) , type=c("correlation") , plot=T , main="ACF for Squared Returns", lag.max = 50)$acf #not stationary
    pacf(diff(y.ts^2), lag = 50, main = "PACF for Squared Returns")
    
    #Correlograms for centered data - NOT HELPFUL CENTERING
    # acf(y.ts - mean(y.ts), type=c("correlation") , plot=T , main="ACF for Centered Returns", lag.max = 50)$acf #stationary, white noise
    # pacf(y.ts - mean(y.ts), lag = 50, main = "PACF for Centered Returns")
    # 
    # acf((y.ts - mean(y.ts))^2 , type=c("correlation") , plot=T , main="ACF for Squared Centered Returns", lag.max = 50)$acf #not stationary
    # pacf((y.ts - mean(y.ts))^2 , lag = 50, main = "PACF for Squared Centered Returns")

    #Correlograms for standardized returns - should have no correlation - assumption met
    acf(Z , type=c("correlation") , plot=T , main="ACF for Standardized Returns", lag.max = 50)$acf #stationary, white noise

    #Outlier detection
    s = sd(y.ts)
    length(y.ts[which(abs(y.ts) > 4*s)])  #number of outliers greater than 4 standard deviations
    y.outrm = y.ts[-which(abs(y.ts) > 3*s)] #TS with outliers removed
    
    
    
################################
###--------FIT MODELS--------###
################################    
    
library(rugarch)    
library(tseries)   
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  
#Fit GARCH(p,q) Models
    
#GARCH(1,1) - constant mean model

            #can specify different distribution.model, such as std (student), norm, sstd (skew student), etc
            #mean model assumes constant mean
            gspec11 <- ugarchspec(variance.model = list(model="sGARCH", garchOrder = c(1,1)), 
                                mean.model = list(armaOrder=c(0,0)), 
                                distribution.model="sstd")
        
            #Estimate GARCH model
            gfit11 <- ugarchfit(data = y.ts, spec = gspec11)
            plot(fitted(gfit11))
            plot(gfit11, which = 3, main = "Volatility of Model and Observed")
            plot(gfit11, which = "all")
            gfit11
            
            #Forecast volatility of returns
            forecast11 <- ugarchforecast(fitORspec = gfit11, n.ahead = 5)
            
            #Extract GARCH Info from Model
            gcoef11 <- coef(gfit11)           #coefficients
            guncvar11 <- uncvariance(gfit11)  #unconditional variance
            gmean11 <- fitted(gfit11)         #predicted mean
            gvol11 <- sigma(gfit11)          #predicted volatilities
            
            gcoef11
            
            #Plot Estimated Volatilities
            gfitvol11 = sigma(gfit11)
            plot(gfitvol11)

            #Forecasted volatilities - expect volatility to go down in coming days
            tail(gfitvol11, 1)
            sigma(forecast11)
            par(mfrow = c(1,2))
            plot(forecast11, which = 1)
            plot(forecast11, which = 3)
            
            ###########
            #Assessment
            ###########
            
            #Stationarity - stationarity met
            as.logical(gcoef11[3] + gcoef11[4] < 1)
            as.numeric(gcoef11[3] + gcoef11[4])
            
            #Significance tests - t-tests all parameters significant
            gfit11@fit$matcoef
            
            #Goodness of fit
            mspe11.mean <- mean(residuals(gfit11)^2) #mean squared prediction error of mean prediction
            mspe11.var <- mean((residuals(gfit11)^2 - sigma(gfit11)^2)^2) #variance prediction
            lik.gfit11 <- likelihood(gfit11) #likelihood
            
            mspe11.mean
            mspe11.var
            lik.gfit11
            
            #Information Criterion
            infocriteria(gfit11)[1] #AIC
            infocriteria(gfit11)[2] #BIC
            
            #Ljung-Box Test - residuals not correlated,  p>0.05, model adequate
            Box.test(residuals(gfit11), 22, type = "Ljung-Box")
            par(mfrow=c(2,1))
            acf(ts(as.numeric(residuals(gfit11))), plot=T , main="RSAC for GARCH(1,1)", lag.max = 50)$acf #no spikes, good
            pacf(ts(as.numeric(residuals(gfit11))), lag = 50, main="RSPAC ACF for GARCH(1,1)")
            
            # #Backtesting
            # groll11 = ugarchroll(gspec11, data = y.ts, n.start = 500, refit.window = "expanding", refit.every = 25)
            # g11pred = as.data.frame(groll11)
            # g11mspe.mean = mean((g11pred$Realized - g11pred$Mu)^2)
            # g11mspe.var = mean(((g11pred$Realized - g11pred$Mu)^2 - (g11pred$Sigma)^2)^2)
            # 
            # g11mspe.mean 
            # g11mspe.var
            # groll11@model$coef
            
            #Distribution Estimation Check
            library(fGarch)
            plot(density(y.ts*100), main = "Distribution of Returns TSLA", xlab = "Return %")
            curve(dsstd(x, mean=mean(y.ts*100), sd=sd(y.ts*100), nu = 3.501419 , xi = 1.042796), col="blue", add = T)
            curve(dnorm(x, mean=mean(y.ts*100), sd=sd(y.ts*100)), col="red", add = T)
            legend(6, 0.18, legend=c("Returns", "Normal Dist.", "SST Dist."),
                   col=c("black", "red", "blue"), lty = 1:1, cex=0.8)
            
            #Multi-step-ahead forecast (Out-of-Sample)
            
            gfit11rmse = matrix(nrow = 252, ncol = 1)
            
            for (i in 1:252) {
              y.outcast = y.ts_complete[i:(1005+i), ]
              gcast11i <- ugarchfit(data = y.outcast, spec = gspec11)
              gcast11_for <- ugarchforecast(fitORspec = gcast11i, n.ahead = 1)
              gfit11rmse[i, 1] = sigma(gcast11_for)[1]
            }
            
            gfit11rmse = as.data.frame(gfit11rmse)
            actual = as.data.frame(y[1006:1257, ])
            g11resid = actual[,1]-gfit11rmse[,1]
            rmse = sqrt(mean((g11resid)^2))
            rmse
            
            #Are GARCH residuals normal?
            res = as.numeric(residuals(gfit11))
            head(as.numeric(res))
            rownames(res) = as.numeric(1:nrow(residuals(gfit11)))
            head(res)
            plot(res, main = "Residual Plot of GARCH(1,1)", ylab = "Residual")
            shapiro.test(res)
                      
#GARCH(1,2) - some insignificant parameters
            
            gspec12 <- ugarchspec(variance.model = list(model="sGARCH", garchOrder = c(1,2)), 
                                mean.model = list(armaOrder=c(0,0)), 
                                distribution.model="sstd")
            gfit12 <- ugarchfit(data = y.ts, spec = gspec12)
            plot(fitted(gfit12))
            
            ###########
            #Assessment
            ###########
            
            #Significance tests - t-tests all parameters significant
            gfit12@fit$matcoef
            
            #Goodness of fit
            mspe12.mean <- mean(residuals(gfit12)^2) #mean squared prediction error of mean prediction
            mspe12.var <- mean((residuals(gfit12)^2 - sigma(gfit12)^2)^2) #variance prediction
            lik.gfit12 <- likelihood(gfit12) #likelihood
            
            mspe12.mean
            mspe12.var
            lik.gfit12
            
            #Information Criterion
            infocriteria(gfit12)[1] #AIC
            infocriteria(gfit12)[2] #BIC
            
            #Ljung-Box Test - residuals not correlated,  p>0.05, model adequate
            Box.test(residuals(gfit12), 22, type = "Ljung-Box")
            acf(residuals(gfit12), plot=T , main="RSAC for GARCH(1,2)", lag.max = 50)$acf #no spikes, good
            pacf(residuals(gfit12), lag = 50, main="RSPAC ACF for GARCH(1,2)") #no spikes, good
            
            # #Backtesting
            # groll12 = ugarchroll(gspec12, data = y.ts, n.start = 500, refit.window = "expanding", refit.every = 25)
            # g12pred = as.data.frame(groll12)
            # g12mspe.mean = mean((g12pred$Realized - g12pred$Mu)^2)
            # g12mspe.var = mean(((g12pred$Realized - g12pred$Mu)^2 - (g12pred$Sigma)^2)^2)
            # g12mspe.mean 
            # g12mspe.var
            
#GARCH(2,1) 
            
            gspec21 <- ugarchspec(variance.model = list(model="sGARCH", garchOrder = c(2,1)), 
                                  mean.model = list(armaOrder=c(0,0)), 
                                  distribution.model="sstd")
            gfit21 <- ugarchfit(data = y.ts, spec = gspec21)
            plot(fitted(gfit21))            
            
            ###########
            #Assessment
            ###########
            
            #Significance tests - t-tests all parameters significant
            gfit21@fit$matcoef
            
            #Goodness of fit
            mspe21.mean <- mean(residuals(gfit21)^2) #mean squared prediction error of mean prediction
            mspe21.var <- mean((residuals(gfit21)^2 - sigma(gfit21)^2)^2) #variance prediction
            lik.gfit21 <- likelihood(gfit21) #likelihood
            
            mspe21.mean
            mspe21.var
            lik.gfit21
            
            #Information Criterion
            infocriteria(gfit21)[1] #AIC
            infocriteria(gfit21)[2] #BIC
            
            #Ljung-Box Test - residuals not correlated,  p>0.05, model adequate
            Box.test(residuals(gfit21), 22, type = "Ljung-Box")
            acf(residuals(gfit21), plot=T , main="RSAC for GARCH(2,1)", lag.max = 50)$acf #no spikes, good
            pacf(residuals(gfit21), lag = 50, main="RSPAC ACF for GARCH(2,1)") #no spikes, good
            
            # #Backtesting
            # groll21 = ugarchroll(gspec21, data = y.ts, n.start = 500, refit.window = "expanding", refit.every = 25)
            # g21pred = as.data.frame(groll21)
            # g21mspe.mean = mean((g21pred$Realized - g21pred$Mu)^2)
            # g21mspe.var = mean(((g21pred$Realized - g21pred$Mu)^2 - (g21pred$Sigma)^2)^2)
            # g21mspe.mean 
            # g21mspe.var
            
#GARCH(2,2)             
            
            gspec22 <- ugarchspec(variance.model = list(model="sGARCH", garchOrder = c(2,2)), 
                                  mean.model = list(armaOrder=c(0,0)), 
                                  distribution.model="sstd")
            gfit22 <- ugarchfit(data = y.ts, spec = gspec22)
            plot(fitted(gfit22))            
            
            ###########
            #Assessment
            ###########
            
            #Significance tests - t-tests all parameters significant
            gfit22@fit$matcoef
            
            #Goodness of fit
            mspe22.mean <- mean(residuals(gfit22)^2) #mean squared prediction error of mean prediction
            mspe22.var <- mean((residuals(gfit22)^2 - sigma(gfit22)^2)^2) #variance prediction
            lik.gfit22 <- likelihood(gfit22) #likelihood
            
            mspe22.mean
            mspe22.var
            lik.gfit22
            
            #Information Criterion
            infocriteria(gfit22)[1] #AIC
            infocriteria(gfit22)[2] #BIC
            
            #Ljung-Box Test - residuals not correlated,  p>0.05, model adequate
            Box.test(residuals(gfit22), 22, type = "Ljung-Box")
            acf(residuals(gfit22), plot=T , main="RSAC for GARCH(2,2)", lag.max = 50)$acf #no spikes, good
            pacf(residuals(gfit22), lag = 50, main="RSPAC ACF for GARCH(2,2)") #no spikes, good
            
            # #Backtesting
            # groll22 = ugarchroll(gspec22, data = y.ts, n.start = 500, refit.window = "expanding", refit.every = 25)
            # g22pred = as.data.frame(groll22)
            # g22mspe.mean = mean((g22pred$Realized - g22pred$Mu)^2)
            # g22mspe.var = mean(((g22pred$Realized - g22pred$Mu)^2 - (g22pred$Sigma)^2)^2)
            # g22mspe.mean 
            # g22mspe.var

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~            
            
#Fitting an mcsGARCH model from GARCH(1,1)
            
            library(zoo)            
            library(xts)            
            
            df = as.data.frame(sigma(gfit11))    #dataframe of garch volatilities      
            mcsgfitsigma = as.xts(df)            #convert to xts 
            xtsy.ts = as.xts(y.ts)               #convert to xts
            
            mcsspecg11 = ugarchspec(mean.model = list(armaOrder = c(0,0)), 
                                    variance.model = list(model = 'mcsGARCH'),
                                    distribution.model = "sstd")
            mcsfitg11 = ugarchfit(data = xtsy.ts, spec = mcsspecg11, DailyVar = mcsgfitsigma^2)
            mcsfitg11
            
            plot(mcsfitg11, which = "all")
            
            #Plot Seasonal Daily Volatility
            plot(as.numeric(mcsfitg11@model$DiurnalVar^0.5), type = 'l', main = "Sigma[Diurnal]")
            
            #Ljung-Box Test
            Box.test(residuals(mcsfitg11), 22, type = "Ljung-Box")
            acf(residuals(mcsfitg11), plot=T , main="RSAC for mcsGARCH(1,1)", lag.max = 50)$acf #no spikes, good
            pacf(residuals(mcsfitg11), lag = 50, main="RSPAC ACF for mcsGARCH(1,1)")
            
            #Information Criterion
            infocriteria(mcsfitg11)[1] #AIC
            infocriteria(mcsfitg11)[2] #BIC
            
            lik.gfit11 <- likelihood(mcsfitg11) #likelihood
            lik.gfit11
            
            #Backtesting
            # mscgroll11 = ugarchroll(mcsfitg11, data = xtsy.ts, n.start = 500, refit.window = "expanding", refit.every = 25)
            # mscg11pred = as.data.frame(mscgroll11)
            # mscg11mspe.mean = mean((mscg11pred$Realized - mscg11pred$Mu)^2)
            # mscg11mspe.var = mean(((mscg11pred$Realized - mscg11pred$Mu)^2 - (mscg11pred$Sigma)^2)^2)
            # 
            # mscg11mspe.mean 
            # mscg11mspe.var
            # mscgroll11@model$coef
            
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
            
#Compute GJR-GARCH models
            
#GJR-GARCH(1,1)   
            
            gjrspec11 <- ugarchspec(variance.model = list(model="gjrGARCH", garchOrder = c(1,1)), 
                                  mean.model = list(armaOrder=c(0,0)), 
                                  distribution.model="sstd")
            
            #Estimate GJR-GARCH model
            gjrfit11 <- ugarchfit(data = y.ts, spec = gjrspec11)
            plot(fitted(gjrfit11))
            gjrfit11
            plot(gjrfit11, which = "all")
            ###########
            #Assessment
            ###########
            
            #Significance tests - t-tests all parameters significant
            gjrfit11@fit$matcoef
            
            #Goodness of fit
            gjrmspe11.mean <- mean(residuals(gjrfit11)^2) #mean squared prediction error of mean prediction
            gjrmspe11.var <- mean((residuals(gjrfit11)^2 - sigma(gjrfit11)^2)^2) #variance prediction
            lik.gjrfit11 <- likelihood(gjrfit11) #likelihood
            
            gjrmspe11.mean
            gjrmspe11.var
            lik.gjrfit11
            
            #Information Criterion
            infocriteria(gjrfit11)[1] #AIC
            infocriteria(gjrfit11)[2] #BIC
            
            #Ljung-Box Test 
            Box.test(residuals(gjrfit11), 22, type = "Ljung-Box")
            acf(residuals(gjrfit11), plot=T , main="RSAC for GJR-GARCH(1,1)", lag.max = 50)$acf 
            pacf(residuals(gjrfit11), lag = 50, main="RSPAC ACF for GJR-GARCH(1,1)")
            
            #Backtesting
            # gjrroll11 = ugarchroll(gjrspec11, data = y.ts, n.start = 500, refit.window = "expanding", refit.every = 25)
            # gjr11pred = as.data.frame(gjrroll11)
            # gjr11mspe.mean = mean((gjr11pred$Realized - gjr11pred$Mu)^2)
            # gjr11mspe.var = mean(((gjr11pred$Realized - gjr11pred$Mu)^2 - (gjr11pred$Sigma)^2)^2)
            # gjr11mspe.mean
            # gjr11mspe.var
            
#GJR-GARCH(1,2)   
            
            gjrspec12 <- ugarchspec(variance.model = list(model="gjrGARCH", garchOrder = c(1,2)), 
                                    mean.model = list(armaOrder=c(0,0)), 
                                    distribution.model="sstd")
            
            #Estimate GJR-GARCH model
            gjrfit12 <- ugarchfit(data = y.ts, spec = gjrspec12)
            plot(fitted(gjrfit12))
            gjrfit12
            
            ###########
            #Assessment
            ###########
            
            #Significance tests - t-tests all parameters significant
            gjrfit12@fit$matcoef
            
            #Goodness of fit
            gjrmspe12.mean <- mean(residuals(gjrfit12)^2) #mean squared prediction error of mean prediction
            gjrmspe12.var <- mean((residuals(gjrfit12)^2 - sigma(gjrfit12)^2)^2) #variance prediction
            lik.gjrfit12 <- likelihood(gjrfit12) #likelihood
            
            gjrmspe12.mean
            gjrmspe12.var
            lik.gjrfit12
            
            #Information Criterion
            infocriteria(gjrfit12)[1] #AIC
            infocriteria(gjrfit12)[2] #BIC
            
            #Ljung-Box Test 
            Box.test(residuals(gjrfit12), 22, type = "Ljung-Box")
            acf(residuals(gjrfit12), plot=T , main="RSAC for GJR-GARCH(1,2)", lag.max = 50)$acf 
            pacf(residuals(gjrfit12), lag = 50, main="RSPAC ACF for GJR-GARCH(1,2)")
            
            # #Backtesting
            # gjrroll12 = ugarchroll(gjrspec12, data = y.ts, n.start = 500, refit.window = "expanding", refit.every = 25)
            # gjr12pred = as.data.frame(gjrroll12)
            # gjr12mspe.mean = mean((gjr12pred$Realized - gjr12pred$Mu)^2)
            # gjr12mspe.var = mean(((gjr12pred$Realized - gjr12pred$Mu)^2 - (gjr12pred$Sigma)^2)^2)
            # gjr12mspe.mean
            # gjr12mspe.var
            
#GJR-GARCH(2,1)   
            
            gjrspec21 <- ugarchspec(variance.model = list(model="gjrGARCH", garchOrder = c(2,1)), 
                                    mean.model = list(armaOrder=c(0,0)), 
                                    distribution.model="sstd")
            
            #Estimate GJR-GARCH model
            gjrfit21 <- ugarchfit(data = y.ts, spec = gjrspec21)
            plot(fitted(gjrfit21))
            gjrfit21
            
            ###########
            #Assessment
            ###########
            
            #Significance tests - t-tests all parameters significant
            gjrfit21@fit$matcoef
            
            #Goodness of fit
            gjrmspe21.mean <- mean(residuals(gjrfit21)^2) #mean squared prediction error of mean prediction
            gjrmspe21.var <- mean((residuals(gjrfit21)^2 - sigma(gjrfit21)^2)^2) #variance prediction
            lik.gjrfit21 <- likelihood(gjrfit21) #likelihood
            
            gjrmspe21.mean
            gjrmspe21.var
            lik.gjrfit21
            
            #Information Criterion
            infocriteria(gjrfit21)[1] #AIC
            infocriteria(gjrfit21)[2] #BIC
            
            #Ljung-Box Test 
            Box.test(residuals(gjrfit21), 22, type = "Ljung-Box")
            acf(residuals(gjrfit21), plot=T , main="RSAC for GJR-GARCH(2,1)", lag.max = 50)$acf 
            pacf(residuals(gjrfit21), lag = 50, main="RSPAC ACF for GJR-GARCH(2,1)")            
            
            # #Backtesting
            # gjrroll21 = ugarchroll(gjrspec21, data = y.ts, n.start = 500, refit.window = "expanding", refit.every = 25)
            # gjr21pred = as.data.frame(gjrroll21)
            # gjr21mspe.mean = mean((gjr21pred$Realized - gjr21pred$Mu)^2)
            # gjr21mspe.var = mean(((gjr21pred$Realized - gjr21pred$Mu)^2 - (gjr21pred$Sigma)^2)^2)
            # gjr21mspe.mean
            # gjr21mspe.var
            
#GJR-GARCH(2,2)   
            
            gjrspec22 <- ugarchspec(variance.model = list(model="gjrGARCH", garchOrder = c(2,2)), 
                                    mean.model = list(armaOrder=c(0,0)), 
                                    distribution.model="sstd")
            
            #Estimate GJR-GARCH model
            gjrfit22 <- ugarchfit(data = y.ts, spec = gjrspec22)
            plot(fitted(gjrfit22))
            gjrfit22
            
            ###########
            #Assessment
            ###########
            
            #Significance tests - t-tests all parameters significant
            gjrfit22@fit$matcoef
            
            #Goodness of fit
            gjrmspe22.mean <- mean(residuals(gjrfit22)^2) #mean squared prediction error of mean prediction
            gjrmspe22.var <- mean((residuals(gjrfit22)^2 - sigma(gjrfit22)^2)^2) #variance prediction
            lik.gjrfit22 <- likelihood(gjrfit22) #likelihood
            
            gjrmspe22.mean
            gjrmspe22.var
            lik.gjrfit22
            
            #Information Criterion
            infocriteria(gjrfit22)[1] #AIC
            infocriteria(gjrfit22)[2] #BIC
            
            #Ljung-Box Test 
            Box.test(residuals(gjrfit22), 22, type = "Ljung-Box")
            acf(residuals(gjrfit22), plot=T , main="RSAC for GJR-GARCH(2,2)", lag.max = 50)$acf 
            pacf(residuals(gjrfit22), lag = 50, main="RSPAC ACF for GJR-GARCH(2,2)") 

            # #Backtesting
            # gjrroll22 = ugarchroll(gjrspec22, data = y.ts, n.start = 500, refit.window = "expanding", refit.every = 10)
            # gjr22pred = as.data.frame(gjrroll22)
            # gjr22mspe.mean = mean((gjr22pred$Realized - gjr22pred$Mu)^2)
            # gjr22mspe.var = mean(((gjr22pred$Realized - gjr22pred$Mu)^2 - (gjr22pred$Sigma)^2)^2)
            # gjr22mspe.mean
            # gjr22mspe.var
            
            #Are GJR-GARCH(2,2) residuals normal?
            resgjr22 = as.numeric(residuals(gjrfit22))
            plot(resgjr22)
            shapiro.test(resgjr22)
            
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

#EGARCH(1,1)   
            
            egspec11 <- ugarchspec(variance.model = list(model="eGARCH", garchOrder = c(1,1)), 
                                    mean.model = list(armaOrder=c(0,0)), 
                                    distribution.model="sstd")
            
            #Estimate EGARCH model
            egfit11 <- ugarchfit(data = y.ts, spec = egspec11)
            plot(fitted(egfit11))
            egfit11
            
            ###########
            #Assessment
            ###########
            
            #Significance tests - t-tests all parameters significant
            egfit11@fit$matcoef
            
            #Goodness of fit
            egmspe11.mean <- mean(residuals(egfit11)^2) #mean squared prediction error of mean prediction
            egmspe11.var <- mean((residuals(egfit11)^2 - sigma(egfit11)^2)^2) #variance prediction
            lik.egfit11 <- likelihood(egfit11) #likelihood
            
            egmspe11.mean
            egmspe11.var
            lik.egfit11
            
            #Information Criterion
            infocriteria(egfit11)[1] #AIC
            infocriteria(egfit11)[2] #BIC
            
            #Ljung-Box Test 
            Box.test(residuals(egfit11), 22, type = "Ljung-Box")
            acf(residuals(egfit11), plot=T , main="RSAC for EGARCH(1,1)", lag.max = 50)$acf 
            pacf(residuals(egfit11), lag = 50, main="RSPAC ACF for EGARCH(1,1)") 
            
            # #Backtesting
            # egroll11 = ugarchroll(egspec11, data = y.ts, n.start = 500, refit.window = "expanding", refit.every = 25)
            # eg11pred = as.data.frame(egroll11)
            # eg11mspe.mean = mean((eg11pred$Realized - eg11pred$Mu)^2)
            # eg11mspe.var = mean(((eg11pred$Realized - eg11pred$Mu)^2 - (eg11pred$Sigma)^2)^2)
            # eg11mspe.mean
            # eg11mspe.var
            
#EGARCH(1,2)   
            
            egspec12 <- ugarchspec(variance.model = list(model="eGARCH", garchOrder = c(1,2)), 
                                   mean.model = list(armaOrder=c(0,0)), 
                                   distribution.model="sstd")
            
            #Estimate EGARCH model
            egfit12 <- ugarchfit(data = y.ts, spec = egspec12)
            plot(fitted(egfit12))
            egfit12
            
            ###########
            #Assessment
            ###########
            
            #Significance tests - t-tests all parameters significant
            egfit12@fit$matcoef
            
            #Goodness of fit
            egmspe12.mean <- mean(residuals(egfit12)^2) #mean squared prediction error of mean prediction
            egmspe12.var <- mean((residuals(egfit12)^2 - sigma(egfit12)^2)^2) #variance prediction
            lik.egfit12 <- likelihood(egfit12) #likelihood
            
            egmspe12.mean
            egmspe12.var
            lik.egfit12
            
            #Information Criterion
            infocriteria(egfit12)[1] #AIC
            infocriteria(egfit12)[2] #BIC
            
            #Ljung-Box Test 
            Box.test(residuals(egfit12), 22, type = "Ljung-Box")
            acf(residuals(egfit12), plot=T , main="RSAC for EGARCH(1,2)", lag.max = 50)$acf 
            pacf(residuals(egfit12), lag = 50, main="RSPAC ACF for EGARCH(1,2)") 
            
            # #Backtesting
            # egroll12 = ugarchroll(egspec12, data = y.ts, n.start = 500, refit.window = "expanding", refit.every = 25)
            # eg12pred = as.data.frame(egroll12)
            # eg12mspe.mean = mean((eg12pred$Realized - eg12pred$Mu)^2)
            # eg12mspe.var = mean(((eg12pred$Realized - eg12pred$Mu)^2 - (eg12pred$Sigma)^2)^2)
            # eg12mspe.mean
            # eg12mspe.var
            
#EGARCH(2,1)   
            
            egspec21 <- ugarchspec(variance.model = list(model="eGARCH", garchOrder = c(2,1)), 
                                   mean.model = list(armaOrder=c(0,0)), 
                                   distribution.model="sstd")
            
            #Estimate EGARCH model
            egfit21 <- ugarchfit(data = y.ts, spec = egspec21)
            plot(fitted(egfit21))
            egfit21
            
            ###########
            #Assessment
            ###########
            
            #Significance tests - t-tests all parameters significant
            egfit21@fit$matcoef
            
            #Goodness of fit
            egmspe21.mean <- mean(residuals(egfit21)^2) #mean squared prediction error of mean prediction
            egmspe21.var <- mean((residuals(egfit21)^2 - sigma(egfit21)^2)^2) #variance prediction
            lik.egfit21 <- likelihood(egfit21) #likelihood
            
            egmspe21.mean
            egmspe21.var
            lik.egfit21
            
            #Information Criterion
            infocriteria(egfit21)[1] #AIC
            infocriteria(egfit21)[2] #BIC
            
            #Ljung-Box Test 
            Box.test(residuals(egfit21), 22, type = "Ljung-Box")
            acf(residuals(egfit21), plot=T , main="RSAC for EGARCH(2,1)", lag.max = 50)$acf 
            pacf(residuals(egfit21), lag = 50, main="RSPAC ACF for EGARCH(2,1)") 
            
            # #Backtesting
            # egroll21 = ugarchroll(egspec21, data = y.ts, n.start = 500, refit.window = "expanding", refit.every = 25)
            # eg21pred = as.data.frame(egroll21)
            # eg21mspe.mean = mean((eg21pred$Realized - eg21pred$Mu)^2)
            # eg21mspe.var = mean(((eg21pred$Realized - eg21pred$Mu)^2 - (eg21pred$Sigma)^2)^2)
            # eg21mspe.mean
            # eg21mspe.var
            
#EGARCH(2,2)   
            
            egspec22 <- ugarchspec(variance.model = list(model="eGARCH", garchOrder = c(2,2)), 
                                   mean.model = list(armaOrder=c(0,0)), 
                                   distribution.model="sstd")
            
            #Estimate EGARCH model
            egfit22 <- ugarchfit(data = y.ts, spec = egspec22)
            plot(fitted(egfit22))
            egfit22
            
            ###########
            #Assessment
            ###########
            
            #Significance tests - t-tests all parameters significant
            egfit22@fit$matcoef
            
            #Goodness of fit
            egmspe22.mean <- mean(residuals(egfit22)^2) #mean squared prediction error of mean prediction
            egmspe22.var <- mean((residuals(egfit22)^2 - sigma(egfit22)^2)^2) #variance prediction
            lik.egfit22 <- likelihood(egfit22) #likelihood
            
            egmspe22.mean
            egmspe22.var
            lik.egfit22
            
            #Information Criterion
            infocriteria(egfit22)[1] #AIC
            infocriteria(egfit22)[2] #BIC
            
            #Ljung-Box Test 
            Box.test(residuals(egfit22), 22, type = "Ljung-Box")
            acf(residuals(egfit22), plot=T , main="RSAC for EGARCH(2,2)", lag.max = 50)$acf 
            pacf(residuals(egfit22), lag = 50, main="RSPAC ACF for EGARCH(2,2)") 

            # #Backtesting
            # egroll22 = ugarchroll(egspec22, data = y.ts, n.start = 500, refit.window = "expanding", refit.every = 25)
            # eg22pred = as.data.frame(egroll22)
            # eg22mspe.mean = mean((eg22pred$Realized - eg22pred$Mu)^2)
            # eg22mspe.var = mean(((eg22pred$Realized - eg22pred$Mu)^2 - (eg22pred$Sigma)^2)^2)
            # eg22mspe.mean
            # eg22mspe.var
            