##################################  AMS 518 Code for Project Draft in R - Koushani Das

install.packages("quantmod")
install.packages("extRemes")
install.packages("PerformanceAnalytics")
install.packages("hqreg")
install.packages("reshape")
install.packages("rugarch")

library(extRemes)
library(PerformanceAnalytics)
library(quantmod)
library(quantreg)
library(ggplot2)
library(reshape)
library(rugarch)
library(dplyr)
library(readxl)
library(hqreg)
library(zoo)


p =0.99
start <- as.Date('2011-1-1')
end <- as.Date('2021-12-31')

# Read stock names from file Equity_List.csv
# The 30 stocks are randomly chosen from S&P500
# Portfolio is formed by equally weighted among 30 stocks
#tickers <- read.table("D:\\Equity_list.txt", header= F)


#print(tickers)
tickers <- c('IBM','AAPL','MSFT','AMZN','TSLA','GOOGL','UNH','JPM','JNJ','HD')
Data <- c()
for (i in 1:length(tickers)) {
  
  data <- getSymbols(tickers[i], src = "yahoo", from = start, to = end, auto.assign = FALSE)[,6]
  Data <- merge.xts(Data, data)
  
}
Data <- as.data.frame(Data)
print(Data)
#tickers <- list('MCD','HOG')
#stock_data <- data.frame(lapply(tickers$V1, function (x) Ad(getSymbols.yahoo(x, from = start, to = end, auto.assign = TRUE))))
stock_rets <- na.omit(ROC(Data, type = "discrete"))
port_ret <- apply(stock_rets, 1, mean)
print(port_ret)
port <- data.frame(date = as.Date(row.names(stock_rets)), port_ret, row.names = NULL)
print(port)
# Calculate 1-day VaR forecast with rolling window of 500 days
# Historical Simulation
HS <- function(mydata){
  res <- mydata[order(mydata)][5]
  return(res)
}

HS_VaR <- na.omit(lag(rollapply(port$port_ret, 500, HS)))

# Parametric VaR
# Normal distribution is used
cal_var <- function(x){
  return(-1*qnorm(p, mean(x), sd(x)))
}

para_VaR <- na.omit(lag(rollapply(port$port_ret, 500, cal_var)))
print(para_VaR)
#Garch 
# This block might take 5+ minutes to run
garch11.spec <- ugarchspec(mean.model = list(armaOrder=c(0,0)),
                           distribution.model = 'norm')

garchmodel <- function(x){
  garch11.fit <- ugarchfit(data = x,
                           spec = garch11.spec,
                           solver.control = list(tol = 1e-4, delta=1e-9))
  
  forecast <- ugarchforecast(garch11.fit)
  
  std <- forecast@forecast$sigmaFor[1]
  
  return(-1*qnorm(p, mean(x), std))
}

garch_var <- NULL
for(i in 500:length(port$port_ret)){
  garch_var <- rbind(garch_var, garchmodel(port$port_ret[1:i]))
}
garch_var <- na.omit(lag(garch_var))

# EVT VaR
EVT_model <- function(mydata){
  ret <- -1 * mydata
  thresh <- quantile(ret, p)
  
  # Log likelihood fucntion
  likelihood <- function(x){
    mu <- ret[ret > thresh] - thresh
    logL <- -1 * sum(log(((1+(x[2]*mu/x[1]))^(-1/x[2]-1))/x[1]))
    return(logL)
  }
  
  # Optimization process
  c <- c(0.1,0.1)   # Start point
  coef <- optim(c,likelihood)
  # Optimization results
  scale = coef$par[1]
  shape = coef$par[2]
  
  # Calculate VaR
  Nu <- length(ret[ret > thresh])
  N <- length(ret)
  EVT <- -1 * ( thresh + ( ((1-p)*N/Nu)^(-1*shape)-1 ) * scale / shape)
  
  return(EVT)
}

EVT_VaR <- na.omit(lag(rollapply(port$port_ret, 500, EVT_model)))
print(EVT_VaR)
# Combine the 3 VaR estimates together with portfolio returns
all_VaR <- data.frame(port$port_ret[501:length(port$port_ret)], HS_VaR, para_VaR, EVT_VaR)
names(all_VaR) <- c('return', 'HS', 'para', 'evt')
print(all_VaR)
# Estimate qr_VaR using quantile regression
# Unpenalized quantile regression
QRmodel <- function(mydata){
  mydata <- data.frame(mydata)
  col <- ncol(mydata)
  row <- nrow(mydata)
  X <- mydata[1:row-1, 1:4]
  y <- mydata[1:row-1, 4]
  fit <- rq(y~X$HS+X$para+X$evt, tau = 1-p, data = mydata)
  print(fit)
  b0 <- fit$coeff[1]
  print(b0)
  b1 <- fit$coeff[2]
  print(b1)
  b2 <- fit$coeff[3]
  b3 <- fit$coeff[4]
  b4 <- fit$coeff[5]
  res <- b0 + b1*mydata[row,1] + b2*mydata[row,2] + b3*mydata[row,3] + b4*mydata[row,4]
  return(res)
}

qr_VaR <- rollapply(all_VaR, 501, QRmodel, by.column = FALSE)
print(qr_VaR)

# Elastic Net Penalized quantile regression
QRmodel_Penalized <- function(mydata){
  col <- ncol(mydata)
  row <- nrow(mydata)
  X <- mydata[1:row-1, 1:4]
  y <- mydata[1:row-1, 4]
  
  fit <- hqreg_raw(X,y, tau = 1-p, method = 'quantile')
  print(fit)
  b0 <- fit$beta[1,100]
  print(b0)
  b1 <- fit$beta[2,100]
  b2 <- fit$beta[3,100]
  b3 <- fit$beta[4,100]
  b4 <- fit$beta[5,100]
  
  res <- b0 + b1*mydata[row,1] + b2*mydata[row,2] + b3*mydata[row,3] + b4*mydata[row,4]
  
  return(res)
}

qr_VaR_penalized <- rollapply(all_VaR, 501, QRmodel_Penalized, by.column = FALSE)
print(qr_VaR_penalized )
# Get all data together
VaR_data <- data.frame(port$date[501:length(port$date)], all_VaR)
print(VaR_data)
VaR_data <- data.frame(VaR_data[501:length(VaR_data$para), ], qr_VaR_penalized)

print(VaR_data)
names(VaR_data) <- c( 'Date', 'Return', 'HS', 'Para', 'EVT', 'QR_P')
print(names(VaR_data))
VaR_data <- mutate(VaR_data, Mean=(VaR_data$HS+VaR_data$Para+VaR_data$EVT)/4)
print(VaR_data)
# Backtesting
# We will test forecasts performance over 3 periods: 
# 1. Overall Period (jan 1, 2011 - Dec 29, 2021)
# 2. Fluctuant Period (Dec 24, 2011 - Dec 31, 2019)
# 3. Calm Period (Jan 1, 2019 - Dec 29, 2021)
# Eceptions are defined that the return is less than VaR forecasts
Exceptions <- data.frame(cbind(VaR_data$Return < VaR_data$Para,
                               VaR_data$Return < VaR_data$HS,
                               VaR_data$Return < VaR_data$EVT,
                               #VaR_data$Return < VaR_data$Garch,
                               VaR_data$Return < VaR_data$Mean,
                               #VaR_data$Return < VaR_data$QR,
                               VaR_data$Return < VaR_data$QR_P))
Exceptions_F <- data.frame(cbind(VaR_data$Return[1:1265] < VaR_data$Para[1:1265],
                                 VaR_data$Return[1:1265] < VaR_data$HS[1:1265],
                                 VaR_data$Return[1:1265] < VaR_data$EVT[1:1265],
                                # VaR_data$Return[1:1264] < VaR_data$Garch[1:1264],
                                 VaR_data$Return[1:1265] < VaR_data$Mean[1:1265],
                                 #VaR_data$Return[1:1264] < VaR_data$QR[1:1264],
                                 VaR_data$Return[1:1265] < VaR_data$QR_P[1:1265]))
Exceptions_C <- data.frame(cbind(VaR_data$Return[1265:2523] < VaR_data$Para[1265:2523],
                                 VaR_data$Return[1265:2523] < VaR_data$HS[1265:2523],
                                 VaR_data$Return[1265:2523] < VaR_data$EVT[1265:2523],
                                # VaR_data$Return[1265:2523] < VaR_data$Garch[1265:2523],
                                 VaR_data$Return[1265:2523] < VaR_data$Mean[1265:2523],
                                # VaR_data$Return[1265:2523] < VaR_data$QR[1265:2523],
                                 VaR_data$Return[1265:2523] < VaR_data$QR_P[1265:2523]))
colnames(Exceptions) <- c('Normal', 'HS', 'EVT', 'Mean', 'QR_P')
colnames(Exceptions_F) <- c('Normal', 'HS', 'EVT', 'Mean', 'QR_P')
colnames(Exceptions_C) <- c('Normal', 'HS', 'EVT', 'Mean', 'QR_P')
Backtest <- data.frame(apply(Exceptions, 2, sum))
print(Backtest)
Backtest_F <- data.frame(apply(Exceptions_F, 2, sum))
print(Backtest_F)
Backtest_C <- data.frame(apply(Exceptions_C, 2, sum))
print(Backtest_C)

# Unconditional Coverage
UC_test <- function(mydata, arg){
  T <- arg
  N <- mydata
  P <- N / T
  Q <- 1 - P
  q <- 1 - p
  LR.ration <- log(((Q^(T-N)*P^N)/(p^(T-N)*q^N))^2)
  return(LR.ration)
}

LR_UC <- apply(Backtest, 2, UC_test, arg = 2523)
UC <- (round((1 - pchisq(LR_UC, df = 1)), digits = 3))

LR_UC_F <- apply(Backtest_F, 2, UC_test, arg = 1264)
UC_F <- (round((1 - pchisq(LR_UC_F, df = 1)), digits = 3))

LR_UC_C <- apply(Backtest_C, 2, UC_test, arg = 2523-1264)
UC_C <- (round((1 - pchisq(LR_UC_C, df = 1)), digits = 3))

# Conditional Coverage
# Test independence
# Overall Period
temp <- apply(apply(apply(Exceptions, 2, as.numeric), 2, as.character), 2, paste, collapse = '')
print(temp)
T_11 <- NULL
for (i in temp) {
  z <- sapply(sapply(strsplit(i,'0'), strsplit, split = ''), length)
  T_11 <- rbind(T_11, sum(z)-length(z[z!=0]))
}
T_01 <- Backtest - T_11
T_10 <- T_01
T_00 <- 2523 - Backtest - T_10

CC_data <- data.frame(T_00, T_10, T_11, T_01)
names(CC_data) <- c('T_00', 'T_10', 'T_11', 'T_01')

CC_data <- CC_data %>%
  mutate(P = (T_01+T_11) / (T_00+T_01+T_10+T_11),
         P_01 = T_01 / (T_00 + T_01),
         P_11 = T_11 / (T_10 + T_11),
         LR.ration = -2 * log((1-P)^(T_00+T_10)*P^(T_01+T_11))
         + 2 * log((1-P_01)^T_00 * P_01^T_01 * (1-P_11)^T_10 * P_11^T_11))
LR_IND <- CC_data$LR.ration
IND <- (round((1 - pchisq(LR_IND, df = 1)), digits = 3))
print(IND)
#Fluctuant Period
temp <- apply(apply(apply(Exceptions_F, 2, as.numeric), 2, as.character), 2, paste, collapse = '')

T_11_F <- NULL
for (i in temp) {
  z <- sapply(sapply(strsplit(i,'0'), strsplit, split = ''), length)
  T_11_F <- rbind(T_11_F, sum(z)-length(z[z!=0]))
}
T_01_F <- Backtest_F - T_11_F
T_10_F <- T_01_F
T_00_F <- 1264 - Backtest_F - T_10_F

CC_data_F <- data.frame(T_00_F, T_10_F, T_11_F, T_01_F)
names(CC_data_F) <- c('T_00', 'T_10', 'T_11', 'T_01')

CC_data_F <- CC_data_F %>%
  mutate(P = (T_01+T_11) / (T_00+T_01+T_10+T_11),
         P_01 = T_01 / (T_00 + T_01),
         P_11 = T_11 / (T_10 + T_11),
         LR.ration = -2 * log((1-P)^(T_00+T_10)*P^(T_01+T_11))
         + 2 * log((1-P_01)^T_00 * P_01^T_01 * (1-P_11)^T_10 * P_11^T_11))
LR_IND_F <- CC_data_F$LR.ration
IND_F <- (round((1 - pchisq(LR_IND_F, df = 1)), digits = 3))
print(IND_F)
# Calm Period
temp <- apply(apply(apply(Exceptions_C, 2, as.numeric), 2, as.character), 2, paste, collapse = '')

T_11_C <- NULL
for (i in temp) {
  z <- sapply(sapply(strsplit(i,'0'), strsplit, split = ''), length)
  T_11_C <- rbind(T_11_C, sum(z)-length(z[z!=0]))
}
T_01_C <- Backtest_C - T_11_C
T_10_C <- T_01_C
T_00_C <- 1259 - Backtest_C - T_10_C

CC_data_C <- data.frame(T_00_C, T_10_C, T_11_C, T_01_C)
names(CC_data_C) <- c('T_00', 'T_10', 'T_11', 'T_01')

CC_data_C <- CC_data_C %>%
  mutate(P = (T_01+T_11) / (T_00+T_01+T_10+T_11),
         P_01 = T_01 / (T_00 + T_01),
         P_11 = T_11 / (T_10 + T_11),
         LR.ration = -2 * log((1-P)^(T_00+T_10)*P^(T_01+T_11))
         + 2 * log((1-P_01)^T_00 * P_01^T_01 * (1-P_11)^T_10 * P_11^T_11))
LR_IND_C <- CC_data_C$LR.ration
IND_C <- (round((1 - pchisq(LR_IND_C, df = 1)), digits = 3))

# LR_CC
LR_CC <- LR_UC + LR_IND
CC <- round(1 - pchisq(LR_CC, df = 2), digits = 3)

LR_CC_F <- LR_UC_F + LR_IND_F
CC_F <- round(1 - pchisq(LR_CC_F, df = 2), digits = 3)

LR_CC_C <- LR_UC_C + LR_IND_C
CC_C <- round(1 - pchisq(LR_CC_C, df = 2), digits = 3)

# Plot
f1 <- ggplot(VaR_data, aes(x = Date, y = Return, col = 'Portfolio Return')) + 
  geom_line() + 
  geom_line(aes(y = Para, col = 'Normal VaR')) + 
  geom_point(aes(y=Return, col='Exceptions'), 
             color = ifelse(VaR_data$Return < VaR_data$Para,'black','red'), 
             size = ifelse(VaR_data$Return < VaR_data$Para,1,-1)) + 
  labs(title = "Normal Model") + 
  theme(plot.title=element_text(hjust=0.6),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = 'top',
        axis.title.y = element_blank())
print(f1)

f2 <- ggplot(VaR_data, aes(x = Date, y = Return, col = 'Portfolio Return')) + 
  geom_line() + 
  geom_line(aes(y = HS, col = 'HS VaR')) + 
  geom_point(aes(y=Return, col='Exceptions'), 
             color = ifelse(VaR_data$Return < VaR_data$HS,'black','blue'), 
             size = ifelse(VaR_data$Return < VaR_data$HS,1,-1)) + 
  labs(title = "Historical Simulation") + 
  theme(plot.title=element_text(hjust=0.5),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = 'top',
        axis.title.y = element_blank())
print(f2)
f3 <- ggplot(VaR_data, aes(x = Date, y = Return, col = 'Portfolio Return')) + 
  geom_line() + 
  geom_line(aes(y = EVT, col = 'EVT VaR')) + 
  geom_point(aes(y=Return, col='Exceptions'), 
             color = ifelse(VaR_data$Return < VaR_data$EVT,'black','blue'), 
             size = ifelse(VaR_data$Return < VaR_data$EVT,1,-1)) + 
  labs(title = "Extreme Value Theory") + 
  theme(plot.title=element_text(hjust=0.5),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = 'top',
        axis.title.y = element_blank())
print(f3)
f4 <- ggplot(VaR_data, aes(x = Date, y = Return, col = 'Portfolio Return')) + 
  geom_line() + 
  geom_line(aes(y = Garch, col = 'GARCH VaR')) + 
  geom_point(aes(y=Return, col='Exceptions'), 
             color = ifelse(VaR_data$Return < VaR_data$Garch,'black','blue'), 
             size = ifelse(VaR_data$Return < VaR_data$Garch,1,-1)) + 
  labs(title = "GARCH(1, 1) Model") + 
  theme(plot.title=element_text(hjust=0.5),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = 'top',
        axis.title.y = element_blank())
print(f4)
f5 <- ggplot(VaR_data, aes(x = Date, y = Return, col = 'Portfolio Return')) + 
  geom_line() + 
  geom_line(aes(y = Mean, col = 'Simple Mean')) + 
  geom_point(aes(y=Return, col='Exceptions'), 
             color = ifelse(VaR_data$Return < VaR_data$Mean,'black','blue'), 
             size = ifelse(VaR_data$Return < VaR_data$Mean,1,-1)) + 
  labs(title = "Simple Mean") + 
  theme(plot.title=element_text(hjust=0.5),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = 'top',
        axis.title.y = element_blank())
print(f5)

f6 <- ggplot(VaR_data, aes(x = Date, y = Return, col = 'Portfolio Return')) + 
  geom_line() + 
  geom_line(aes(y = QR, col = 'QR Unpenalized')) + 
  geom_point(aes(y=Return, col='Exceptions'), 
             color = ifelse(VaR_data$Return < VaR_data$QR,'black','orange'), 
             size = ifelse(VaR_data$Return < VaR_data$QR,1,-1)) + 
  labs(title = "Unpenalized QR") + 
  theme(plot.title=element_text(hjust=0.5),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = 'top',
        axis.title.y = element_blank())

f7 <- ggplot(VaR_data, aes(x = Date, y = Return, col = 'Portfolio Return')) + 
  geom_line() + 
  geom_line(aes(y = QR_P, col = 'QR Penalized')) + 
  geom_point(aes(y=Return, col='Exceptions'), 
             color = ifelse(VaR_data$Return < VaR_data$QR_P,'black','orange'), 
             size = ifelse(VaR_data$Return < VaR_data$QR_P,1,-1)) + 
  labs(title = "Penalized QR") + 
  theme(plot.title=element_text(hjust=0.5),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = 'top',
        axis.title.y = element_blank())

print(f7)