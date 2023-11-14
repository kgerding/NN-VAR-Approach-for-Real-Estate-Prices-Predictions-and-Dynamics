#Author: Gerding, Kilian
#Date: 18.11.2023

# Contents
# Chapter 0: Set-up and Data Import
# Chapter 1: VAR
# Chapter 2: NN-VAR
# Chapter 3: Impulse Response

# Clear Environment
rm(list = ls(all = TRUE))

# graphics
par(mar=c(1,1,1,1))
par(mfrow = c(1,1))

# install package "easypackages" to be able to use the libraries function later on
install.packages("easypackages")
library(easypackages)

# ========== Chapter 0: Set-up and Data Import ==========

# Step 1: Install/Import Neccessary Packages

# define necessary packages
packages <- c("ggplot2",
              "dplyr",
              "tidyr",
              "vars",
              "AER",
              "lubridate",
              "tseries",
              "psych",
              "deSolve",
              "tsDyn",
              "TSPred",
              "imputeTS",
              "xts",
              "forecast",
              "data.table",
              "Metrics",
              "VAR.etp"
)

# eliminating unnecessary loading time by only installing packages not installed yet
lapply(packages, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)
})

# import defined packages
libraries(packages)

# modify path as needed, where the Swiss and US data reside
path <- 
setwd(path)

# Step 2: Read in data and data cleaning

# import data hpi, gdp, cpi, rate
df <- read.csv("swiss_quarterly_all.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

# transform date format
df$Date <- as.Date(df$Date, format = "%Y-%m-%d")

# switch to xts format
df <- as.xts(df[-1], order.by = df$Date)

# Step 3: Data wrangling and formatting

# HPI log transformation
df$hpi_log <- log(df$hpi)
df$hpi_log_chg <- df$hpi_log - stats::lag(df$hpi_log)

# de-trend and compute chages from the time series trend
lin.mod <- lm(df$hpi_log ~ time(df$hpi_log))
lin.trend <- lin.mod$fitted.values
df$linear_hpi <- xts(lin.trend, order.by = index(df))
df$hpi_lin_cycle <- df$hpi_log - df$linear_hpi
df$hpi_lin_cycle_df <- df$hpi_lin_cycle - stats::lag(df$hpi_lin_cycle)


# GDP log values
df$gdp_log <- log(df$gdp)
df$gdp_log_chg <- df$gdp_log - stats::lag(df$gdp_log)

# de-trend and compute chages from the time series trend
lin.mod <- lm(df$gdp_log ~ time(df$gdp_log))
lin.trend <- lin.mod$fitted.values
df$linear_gdp <- xts(lin.trend, order.by = index(df))
df$gdp_lin_cycle <- df$gdp_log - df$linear_gdp
df$gdp_lin_cycle_df <- df$gdp_lin_cycle - stats::lag(df$gdp_lin_cycle)

# get rate change
df$rate_chg <- df$rate - stats::lag(df$rate)

# get get CPI YoY changes
df$cpi_chg <- df$cpi - stats::lag(df$cpi)

# merge all CH/US data, ordered in the discussed order in the thesis

ch_data <- merge(df$gdp_lin_cycle_df,
                 df$cpi_chg,
                 df$rate,
                 df$hpi_lin_cycle_df,
                 all = FALSE)

colnames(ch_data) <- c("gdp", "cpi", "rate", "hpi")

us_data <- merge(df$gdp_lin_cycle_df,
                 df$cpi_chg,
                 df$rate_chg,
                 df$hpi_lin_cycle_df,
                 all = FALSE)

colnames(us_data) <- c("gdp", "cpi", "rate", "hpi")

# ========== Chapter 1: VAR  ==========

# amend if needed from CH to US Data
df <- na.omit(ch_data)

# split in training, validation and, testing

# t for training subset, set to 0.5, 0.7 and 0.8
t <- round(0.5*nrow(df))

# v for validation subset
v <- round((nrow(df)-t)/2)

# split data
df_train  <- df[1:t,]
df_val    <- df[(t+1):(t+v),]
df_train2 <- df[1:(t+v),]
df_test   <- df[(t+v+1):(t+2*v),]

time_index_all <- index(ch_data)
time_index_val <- index(ch_data[1:(t+v),])
time_index_test <- index(ch_data[1:(t+2*v),])

# ADF Test to test for stationary or non-stationary time series
for (i in 1:ncol(df_train)) {
  p <- adf.test(df_train[,i])
  print(p)
}

for (i in 1:ncol(df_train2)) {
  p <- adf.test(df_train2[,i])
  print(p)
}

# Specify VAR model
lagorder <- VARselect(df_train, lag.max = 4, type = "const")
lagorder

lagorder2 <- VARselect(df_train2, lag.max = 4, type = "const")
lagorder2

# plot acf of model variables
acf(df_train$cpi, main = "CPI")
acf(df_train$rate, main= "Rate")
acf(df_train$gdp, main = "GDP")
acf(df_train$hpi, main = "HPI")

# plot pacf of model variables
pacf(df_train$cpi, main = "CPI")
pacf(df_train$rate, main = "Rate")
pacf(df_train$gdp, main = "GDP")
pacf(df_train$hpi, main = "HPI")

# LM test

# estimate VAR model
model <- VAR(df_train, p =1, type = "const")
summary(model)

model2 <- VAR(df_train2, p =1, type = "const")
summary(model2)

# check roots
roots <- roots(model)
roots

roots2 <- roots(model2)
roots2

# check residuals autocorrelation
resid_auto <- serial.test(model)
resid_auto

resid_auto2 <- serial.test(model2)
resid_auto2

# check residuals ARCH effects
resid_arch <- arch.test(model)
resid_arch

resid_arch2 <- arch.test(model2)
resid_arch2

# check residuals distribution
resid_normal <- normality.test(model, multivariate.only = TRUE)
resid_normal

resid_normal2 <- normality.test(model2, multivariate.only = TRUE)
resid_normal2

# plot stability cumulative sum
stab <- stability(model, type = "OLS-CUSUM")
plot(stab)

stab2 <- stability(model2, type = "OLS-CUSUM")
plot(stab2)

# plot stability via fluctuations
stabf <- stability(model, type = "fluctuation")
plot(stabf)

stabf2 <- stability(model2, type = "fluctuation")
plot(stabf2)

# plot residuals, as well as histogram, acf and pacf of residuals and acf
# and pacf of squared residuals

lags <- 1 # change in accordance to chosen lag number

residuals_cpi <- xts(model$varresult$cpi$residuals, order.by = index(df_train[(lags+1):t-1,]))
plot(residuals_cpi)
hist(residuals_cpi)
acf(residuals_cpi)
pacf(residuals_cpi)

residuals_rate <- xts(model$varresult$rate$residuals, order.by = index(df_train[(lags+1):t-1,]))
plot(residuals_rate)
hist(residuals_rate)
acf(residuals_rate)
pacf(residuals_rate)

residuals_gdp <- xts(model$varresult$gdp$residuals, order.by = index(df_train[(lags+1):t-1,]))
plot(residuals_gdp)
hist(residuals_gdp)
acf(residuals_gdp)
pacf(residuals_gdp)

residuals_hpi <- xts(model$varresult$hpi$residuals, order.by = index(df_train[(lags+1):t-1,]))
plot(residuals_hpi)
hist(residuals_hpi)
acf(residuals_hpi)
pacf(residuals_hpi)

dev.off()
plot(model, names = "cpi")
dev.off()
plot(model, names = "rate")
dev.off()
plot(model, names = "gdp")
dev.off()
plot(model, names = "hpi")
dev.off()

# forecasting CH

# forecasting horizon is definied by the validation subsets of the data
ahead <- v

# predict n-ahead
forecast1 <- predict(model, n.ahead = ahead, ci = 0.90)
forecast2 <- predict(model2, n.ahead = ahead, ci = 0.90)

# extract validation forecasts
f_cpi <- forecast1$fcst$cpi[,1]
f_rate <- forecast1$fcst$rate[,1]
f_gdp <- forecast1$fcst$gdp[,1]
f_hpi <- forecast1$fcst$hpi[,1]
f_all <- cbind(f_cpi, f_rate,f_gdp, f_hpi)
f_all <- xts(f_all, order.by = index(df_val))

# extract validation forecasts lower bounds
f_cpi_l <- forecast1$fcst$cpi[,2]
f_rate_l <- forecast1$fcst$rate[,2]
f_gdp_l <- forecast1$fcst$gdp[,2]
f_hpi_l <- forecast1$fcst$hpi[,2]
f_all_l <- cbind(f_cpi_l, f_rate_l, f_gdp_l,  f_hpi_l)
f_all_l <- xts(f_all_l, order.by = index(df_val))

# extract validation forecasts upper bounds
f_cpi_u <- forecast1$fcst$cpi[,3]
f_rate_u <- forecast1$fcst$rate[,3]
f_gdp_u <- forecast1$fcst$gdp[,3]
f_hpi_u <- forecast1$fcst$hpi[,3]
f_all_u <- cbind(f_cpi_u, f_rate_u, f_gdp_u,  f_hpi_u)
f_all_u <- xts(f_all_u, order.by = index(df_val))

# extract testing forecasts
f_cpi <- forecast2$fcst$cpi[,1]
f_rate <- forecast2$fcst$rate[,1]
f_gdp <- forecast2$fcst$gdp[,1]
f_hpi <- forecast2$fcst$hpi[,1]
f_all2 <- cbind(f_cpi, f_rate, f_gdp, f_hpi)
f_all2 <- xts(f_all2, order.by = index(df_test))

# extract testing forecasts lower bounds
f_cpi_l <- forecast2$fcst$cpi[,2]
f_rate_l <- forecast2$fcst$rate[,2]
f_gdp_l <- forecast2$fcst$gdp[,2]
f_hpi_l <- forecast2$fcst$hpi[,2]
f_all_l2 <- cbind(f_cpi_l, f_rate_l, f_gdp_l, f_hpi_l)
f_all_l2 <- xts(f_all_l2, order.by = index(df_test))

# extract validation forecasts upper bounds
f_cpi_u <- forecast2$fcst$cpi[,3]
f_rate_u <- forecast2$fcst$rate[,3]
f_gdp_u <- forecast2$fcst$gdp[,3]
f_hpi_u <- forecast2$fcst$hpi[,3]
f_all_u2 <- cbind(f_cpi_u, f_rate_u, f_gdp_u, f_hpi_u)
f_all_u2 <- xts(f_all_u2, order.by = index(df_test))

# merge with actual data
comp <- merge.xts(df_train, f_all, f_all_l, f_all_u, f_all2, f_all_l2, f_all_u2, df_val, df_test, all = TRUE)
colnames(comp) <- c("GDP",
                    "CPI",
                    "Rate",
                    "HPI",
                    "Validation Forecast",
                    "Validation Forecast",
                    "Validation Forecast",
                    "Validation Forecast",
                    "Validation Lower Bound",
                    "Validation Lower Bound",
                    "Validation Lower Bound",
                    "Validation Lower Bound",
                    "Validation Upper Bound",
                    "Validation Upper Bound",
                    "Validation Upper Bound",
                    "Validation Upper Bound",
                    "Testing Forecast",
                    "Testing Forecast",
                    "Testing Forecast",
                    "Testing Forecast",
                    "Testing Lower Bound",
                    "Testing Lower Bound",
                    "Testing Lower Bound",
                    "Testing Lower Bound",
                    "Testing Upper Bound",
                    "Testing Upper Bound",
                    "Testing Upper Bound",
                    "Testing Upper Bound",
                    "1",
                    "2",
                    "3",
                    "4",
                    "5",
                    "6",
                    "7",
                    "8")

# for visual effectiveness
copy_comp <- comp
row <- t
copy_comp[row,c(5,6,7,8)] <- copy_comp[row,c(1,2,3,4)]
copy_comp[row,c(29,30,31,32)] <- copy_comp[row,c(1,2,3,4)]
copy_comp[row,c(9,10,11,12)] <- copy_comp[row+1,c(5,6,7,8)]
copy_comp[row,c(13,14,15,16)] <- copy_comp[row+1,c(5,6,7,8)]
copy_comp[row+ahead,c(17,18,19,20)] <- copy_comp[row+ahead,c(29,30,31,32)]
copy_comp[row+ahead,c(21,22,23,24)] <- copy_comp[row+ahead,c(17,18,19,20)]
copy_comp[row+ahead,c(25,26,27,28)] <- copy_comp[row+ahead,c(17,18,19,20)]
copy_comp[row+ahead,c(33,34,35,36)] <- copy_comp[row+ahead,c(29,30,31,32)]

# plot forecasts
plot(copy_comp[,c(1,5,17)], main = "GDP",
     cex = 0.5, grid.ticks.lwd = 0.15, yaxis.left = TRUE, yaxis.right = FALSE, legend.loc = "bottomleft")
lines(copy_comp[,c(29,33)], lwd=0.75, lty = 1, col='black')
lines(copy_comp[,c(9,13)], lwd=1, lty = 2, col='red')
lines(copy_comp[,c(21,25)], lwd=1, lty = 2, col='lightgreen')

plot(copy_comp[,c(2,6,18)], main = "CPI",
     cex = 0.5, grid.ticks.lwd = 0.15, yaxis.left = TRUE, yaxis.right = FALSE, legend.loc = "bottomleft")
lines(copy_comp[,c(30,34)], lwd=0.75, lty = 1, col='black')
lines(copy_comp[,c(10,14)], lwd=1, lty = 2, col='red')
lines(copy_comp[,c(22,26)], lwd=1, lty = 2, col='lightgreen')

plot(copy_comp[,c(3,7,19)], main = "Rate",
     cex = 0.5, grid.ticks.lwd = 0.15, yaxis.left = TRUE, yaxis.right = FALSE, legend.loc = "topright")
lines(copy_comp[,c(31,35)], lwd=0.75, lty = 1, col='black')
lines(copy_comp[,c(11,15)], lwd=1, lty = 2, col='red')
lines(copy_comp[,c(23,27)], lwd=1, lty = 2, col='lightgreen')

plot(copy_comp[,c(4,8,20)], main = "HPI",
     cex = 0.5, grid.ticks.lwd = 0.15, yaxis.left = TRUE, yaxis.right = FALSE, legend.loc = "topleft")
lines(copy_comp[,c(32,36)], lwd=0.75, lty = 1, col='black')
lines(copy_comp[,c(12,16)], lwd=1, lty = 2, col='red')
lines(copy_comp[,c(24,28)], lwd=1, lty = 2, col='lightgreen')

# prediction evaluation
f1 <- rbind(df_train, f_all)
f2 <- rbind(df_train, df_val, f_all2)
cum_comp <- merge(df, f1, f2)

starting <- as.xts(matrix(c(1,1,1,1,1,1,1,1,1,1,1,1), nrow = 1, ncol = 12),
                   as.Date("2000-06-01"))

colnames(starting) <- colnames(cum_comp)
cum_comp <- rbind(starting, cum_comp)
cum_comp[1:nrow(df_train),5:8] <- NA
cum_comp[1:(nrow(df_train)+nrow(df_val)),9:12] <- NA
colnames(cum_comp) <- c("CPI",
                        "Rate",
                        "GDP",
                        "HPI",
                        "Validation Forecast",
                        "Validation Forecast",
                        "Validation Forecast",
                        "Validation Forecast",
                        "Testing Forecast",
                        "Testing Forecast",
                        "Testing Forecast",
                        "TestingForecast")

# Squared Error and RMSE
v1 <- na.omit(cum_comp[,c(1,5)])
v2 <- na.omit(cum_comp[,c(2,6)])
v3 <- na.omit(cum_comp[,c(3,7)])
v4 <- na.omit(cum_comp[,c(4,8)])

v1_fe <- Metrics::rmse(v1$CPI, v1$`Validation Forecast`)
v2_fe <- Metrics::rmse(v2$Rate, v2$`Validation Forecast`)
v3_fe <- Metrics::rmse(v3$GDP, v3$`Validation Forecast`)
v4_fe <- Metrics::rmse(v4$HPI, v4$`Validation Forecast`)

v5 <- na.omit(cum_comp[,c(1,9)])
v6 <- na.omit(cum_comp[,c(2,10)])
v7 <- na.omit(cum_comp[,c(3,11)])
v8 <- na.omit(cum_comp[,c(4,12)])

v5_fe <- Metrics::rmse(v5$CPI, v5$`Testing Forecast`)
v6_fe <- Metrics::rmse(v6$Rate, v6$`Testing Forecast`)
v7_fe <- Metrics::rmse(v7$GDP, v7$`Testing Forecast`)
v8_fe <- Metrics::rmse(v8$HPI, v8$TestingForecast)

fe_all <- matrix(c(v1_fe,
                   v2_fe,
                   v3_fe,
                   v4_fe,
                   v5_fe,
                   v6_fe,
                   v7_fe,
                   v8_fe),
                 ncol = 4,
                 nrow = 2,
                 byrow = TRUE)
fe_all

# accumulated squared errors
sce <- merge(cum_comp, residuals_gdp, residuals_cpi, residuals_rate, residuals_hpi)
sce[,5] <- (sce[,5] - sce[,1])**2
sce[,6] <- (sce[,6] - sce[,2])**2
sce[,7] <- (sce[,7] - sce[,3])**2
sce[,8] <- (sce[,8] - sce[,4])**2

sce[,9] <- (sce[,9] - sce[,1])**2
sce[,10] <- (sce[,10] - sce[,2])**2
sce[,11] <- (sce[,11] - sce[,3])**2
sce[,12] <- (sce[,12] - sce[,4])**2

sce[,13] <- (sce[,13])**2
sce[,14] <- (sce[,14])**2
sce[,15] <- (sce[,15])**2
sce[,16] <- (sce[,16])**2

sce <- sce[-nrow(sce),-c(1,2,3,4)]
sce[is.na(sce)] <- 0
sce$gdp <- sce[,1] + sce[,5] + sce[,9]
sce$cpi <- sce[,2] + sce[,6] + sce[,10]
sce$rate <- sce[,3] + sce[,7] + sce[,11]
sce$hpi <- sce[,4] + sce[,8] + sce[,12]
sce <- sce[,c(13,14,15,16)]

starting <- as.xts(matrix(c(1,1,1,1), nrow = 1, ncol = 4),
                   as.Date("2000-06-01"))

colnames(starting) <- colnames(sce)
cum_sce <- rbind(starting, sce)
cum_sce <- cumsum(cum_sce)
cum_sce50 <- cum_sce

# cumulative squared error plots
plot(cum_sce[,1], main = "GDP Squared Errors")
lines(cum_sce[(t+1):(t+v),1], lwd=2, lty = 1, col='red')
lines(cum_sce[(t+v):nrow(cum_sce),1], lwd=2, lty = 1, col='lightgreen')

plot(cum_sce[,2], main = "CPI Squared Errors")
lines(cum_sce[(t+1):(t+v),2], lwd=2, lty = 1, col='red')
lines(cum_sce[(t+v):nrow(cum_sce),2], lwd=2, lty = 1, col='lightgreen')

plot(cum_sce[,3], main = "Rate Squared Errors")
lines(cum_sce[(t+1):(t+v),3], lwd=2, lty = 1, col='red')
lines(cum_sce[(t+v):nrow(cum_sce),3], lwd=2, lty = 1, col='lightgreen')

plot(cum_sce[,4], main = "HPI Squared Errors")
lines(cum_sce[(t+1):(t+v),4], lwd=2, lty = 1, col='red')
lines(cum_sce[(t+v):nrow(cum_sce),4], lwd=2, lty = 1, col='lightgreen')

# ========== Chapter 2: NN-VAR ==========

# python configuration
install.packages("tensorflow")

#path to python
library(reticulate)
path_to_python <- install_python()
virtualenv_create("r-reticulate", python = path_to_python)

# tensorflow
library(tensorflow)
install_tensorflow(envname = "r-reticulate")

# keras
install.packages("keras")
library(keras)
install_keras(envname = "r-reticulate")

# tf probability
install.packages("tfprobability")
library(tfprobability)
install_tfprobability(envname = "r-reticulate")

# amend if needed from CH to US data
df <- na.omit(ch_data)
tdf <- data.table(df)

# split df
t <- round(0.5*nrow(tdf))
v <- round((nrow(tdf)-t)/2)
tdf_train  <- tdf[1:t,]
tdf_val    <- tdf[(t+1):(t+v),]
tdf_train2 <- tdf[1:(t+v),]
tdf_test   <- tdf[(t+v+1):(t+2*v),]


# define object with variable names
names <- colnames(tdf)

# define the number of lags (identical to the VAR specification)
lags <- 1

# define the length of the data used for modelling
n <- nrow(tdf) - lags

# define the width of the data used for modelling
k <- ncol(tdf)

# Step 2: Standardize the data and reshape

# define scaling parameters
scaling_parameters_train <- list(
  means = tdf_train[ ,lapply(.SD, mean), .SDcols = names ],
  sd = tdf_train[ ,lapply(.SD, sd), .SDcols = names ]
)

scaling_parameters_train2 <- list(
  means = tdf_train[ ,lapply(.SD, mean), .SDcols = names ],
  sd = tdf_train[ ,lapply(.SD, sd), .SDcols = names ]
)

# standardize the data to mean = 0 and std deviation = 1
tdf_train <- tdf_train[ , (names):= lapply( .SD, function(i) {(i - mean( i )) / sd( i )}) ,.SDcols = names ]
tdf_train2 <- tdf_train2[ , (names):= lapply( .SD, function(i) {(i - mean( i )) / sd( i )}) ,.SDcols = names ]

tdf_val <- tdf_val[ , (names):= lapply( .SD, function(i) {(i - mean( i )) / sd( i )}) ,.SDcols = names ]
tdf_test <- tdf_test[ , (names):= lapply( .SD, function(i) {(i - mean( i )) / sd( i )}) ,.SDcols = names ]

# define independent and dependent variables
y_train  <- as.matrix( tdf_train[ ( lags + 1):.N, 1:k ])
x_train  <- embed( as.matrix(tdf_train), lags + 1 )[ , -c(1:k) ]

y_train2  <- as.matrix( tdf_train2[ ( lags + 1):.N, 1:k ])
x_train2  <- embed( as.matrix(tdf_train2), lags + 1 )[ , -c(1:k) ]

y_val  <- as.matrix( tdf_val[ ( lags + 1):.N, 1:k ])
x_val  <- embed( as.matrix(tdf_val), lags + 1 )[ , -c(1:k) ]

y_test  <- as.matrix( tdf_test[ ( lags + 1):.N, 1:k ])
x_test  <- embed( as.matrix(tdf_test), lags + 1 )[ , -c(1:k) ]

# create colnames such as all variable names are accompanied with their lag number
colnames(x_train) <- c( sapply( 1:lags , function(p) sprintf("%s_l%i", names, p)))
colnames(x_train2) <- c( sapply( 1:lags , function(p) sprintf("%s_l%i", names, p)))
colnames(x_val) <- c( sapply( 1:lags , function(p) sprintf("%s_l%i", names, p)))
colnames(x_test) <- c( sapply( 1:lags , function(p) sprintf("%s_l%i", names, p)))

# reshape independent and dependent variables in a multidimensional array 
x_train <- array_reshape(x_train, c(nrow(x_train), 1, ncol(x_train)))
y_train <- array_reshape(y_train, c(nrow(y_train), 1, ncol(y_train)))

x_train2 <- array_reshape(x_train2, c(nrow(x_train2),1,ncol(x_train2)))
y_train2 <- array_reshape(y_train2, c(nrow(y_train2),1,ncol(y_train2)))

x_val <- array_reshape(x_val, c(nrow(x_val),1,ncol(x_val)))
y_val <- array_reshape(y_val, c(nrow(y_val),1,ncol(y_val)))

x_test <- array_reshape(x_test, c(nrow(x_test),1,ncol(x_test)))
y_test <- array_reshape(y_test, c(nrow(y_test),1,ncol(y_test)))

# Step 3: Prepare the NNVAR model

# define the number of units
units <- 50

# define the number of layers
layers <- 2

# drop out threshold
drop_out <- 0.5

# define the learning rate
epsilon <- 0.001

# define the optimizer alogrithm
optim <- keras$optimizers$legacy$Adam(learning_rate=epsilon) #optimizer_adam(learning_rate=epsilon)

# define the dimensions
nnvar_dim <- dim(x_train)[2:3]

# define the LSTM Networks for each equation nested in a list

nnvar_list <- lapply(
  # for all variables
  1:k,
  function(k) {
    layers_list <- lapply(
      
      # define the pre-defined number of network layers per variable
      1:layers, function(l) {
        list(keras::layer_lstm(units = units,
                               return_sequences = ifelse(l<layers, TRUE, FALSE),
                               input_shape = nnvar_dim), 
             keras::layer_dropout(rate = drop_out))
      }
    )
    
    # call the defined function layers_list
    layers_list <- do.call(c, layers_list)
    
    # define the model 
    model <- keras::keras_model_sequential(layers_list) %>%
      keras::layer_dense(units = 2, activation = "linear") %>%
      
      # define mean and standard deviatino of distribution function
      tfprobability::layer_distribution_lambda(
        function(x) {
          tfprobability::tfd_normal(
            loc = x[, 1, drop = FALSE],
            scale = 1e-3 + tf$math$softplus(x[, 2, drop = FALSE])
          )
        }
      )
    
    # define loss function
    loglikelihood <- function(y, model) - (model %>% tfd_log_prob(y))
    
    # optimiye model with defined loss function and optimizer algorythm
    model %>% keras::compile( loss = loglikelihood,optimizer = optim)
  }
)

# Fit the LSTM networks models on the two datasets

nnvar_fitted <- lapply(
  1:k,
  function(k) {
    hist <- nnvar_list[[k]] %>%
      keras::fit( x = x_train, y = y_train[,,k], verbose = 1, epochs = 500)
    list( model = nnvar_list[[k]],
          history = hist
    )
  }
)

nnvar_fitted2 <- lapply(
  1:k,
  function(k) {
    hist <- nnvar_list[[k]] %>%
      keras::fit( x = x_train2, y = y_train2[,,k], verbose = 1, epochs = 500)
    list( model = nnvar_list[[k]],
          history = hist
    )
  }
)


# yield fitted values for both models

fitted_values <- lapply(
  1:length(nnvar_list),
  
  function(k) {
    
    mod <- nnvar_list[[k]]
    fitted <- mod(x_train)
    y_hat <- as.numeric(fitted %>% tfd_mean())
    sd <- as.numeric(fitted %>% tfd_stddev())
    
    # rescale to the original mean and standard deviation
    y_hat <- (y_hat + scaling_parameters_train$means[[k]]) * scaling_parameters_train$sd[[k]]
    sd <- (sd + scaling_parameters_train$means[[k]]) * scaling_parameters_train$sd[[k]]
    
    return(list(y_hat = unlist(y_hat), sd = unlist(sd)))
  }
)

fitted_values2 <- lapply(
  1:length(nnvar_list),
  
  function(k) {
    
    mod <- nnvar_list[[k]]
    fitted <- mod(x_train2)
    y_hat2 <- as.numeric(fitted %>% tfd_mean())
    sd2 <- as.numeric(fitted %>% tfd_stddev())
    
    # rescale to the original mean and standard deviation
    y_hat2 <- (y_hat2 + scaling_parameters_train2$means[[k]]) * scaling_parameters_train2$sd[[k]] 
    sd2 <- (sd2 + scaling_parameters_train2$means[[k]]) * scaling_parameters_train2$sd[[k]]
    
    return(list(y_hat2 = unlist(y_hat2), sd2 = unlist(sd2)))
  }
)


# transform the fitted values 
y_hat <- matrix(sapply(fitted_values, function(i) i$y_hat), ncol = k)
rownames(y_hat) <- NULL
colnames(y_hat) <- names

y_hat2 <- matrix(sapply(fitted_values2, function(i) i$y_hat2), ncol = k)
rownames(y_hat2) <- NULL
colnames(y_hat2) <- names

# compute residuals
data <- tdf_train

# rescale
data[,1] <- (data[,1] + scaling_parameters_train$means[[1]])  *  scaling_parameters_train$sd[[1]]
data[,2] <- (data[,2] + scaling_parameters_train$means[[2]]) *  scaling_parameters_train$sd[[2]]
data[,3] <- (data[,3] + scaling_parameters_train$means[[3]]) *  scaling_parameters_train$sd[[3]]
data[,4] <- (data[,4] + scaling_parameters_train$means[[4]]) *  scaling_parameters_train$sd[[4]]

res1 <- cbind(data[(lags+1):nrow(data),], y_hat)
res1 <- xts(res1, order.by = index(df_train[(lags+1):nrow(df_train)]))
res11 <- res1[,1:4] - res1[,5:8]

par(mfrow = c(2,1))
plot(res1[,c(1,5)],main = "GDP In-Sample")
plot(res11[,1], main = "GDP Residuals")
plot(res1[,c(2,6)], main = "CPI In-Sample")
plot(res11[,2], main = "CPI Residuals")
plot(res1[,c(3,7)], main = "Rate In-Sample")
plot(res11[,3], main = "Rate Residuals")
plot(res1[,c(4,8)], main = "HPI In-Sample")
plot(res11[,4], main = "HPI Residuals")

data <- tdf_train2

# rescale
data[,1] <- (data[,1] + scaling_parameters_train2$means[[1]])  *  scaling_parameters_train2$sd[[1]]
data[,2] <- (data[,2] + scaling_parameters_train2$means[[2]]) *  scaling_parameters_train2$sd[[2]]
data[,3] <- (data[,3] + scaling_parameters_train2$means[[3]]) *  scaling_parameters_train2$sd[[3]]
data[,4] <- (data[,4] + scaling_parameters_train2$means[[4]]) *  scaling_parameters_train2$sd[[4]]

res2 <- cbind(data[(lags+1):nrow(data),], y_hat2)
res2 <- xts(res2, order.by = index(df_train2[(lags+1):nrow(df_train2)]))
res22 <- res2[,1:4] - res2[,5:8]

par(mfrow = c(2,1))
plot(res2[,c(1,5)],main = "GDP In-Sample")
plot(res22[,1], main = "GDP Residuals")
plot(res2[,c(2,6)], main = "CPI In-Sample")
plot(res22[,2], main = "CPI Residuals")
plot(res2[,c(3,7)], main = "Rate In-Sample")
plot(res22[,3], main = "Rate Residuals")
plot(res2[,c(4,8)], main = "HPI In-Sample")
plot(res22[,4], main = "HPI Residuals")
par(mfrow = c(1,1))

# compute out of sample 1-step ahead recursive forecast

# define steps ahead and counter
ahead <- v
counter <- 1

# define base set
forecast <- data.table::copy(tdf_train[.N,])
uncertainty <- data.table::copy(tdf_train[.N,])
uncertainty[] <- 0

# change for loop
data <- tdf_train

# rescale
data[,1] <- (data[,1] + scaling_parameters_train$means[[1]])  *  scaling_parameters_train$sd[[1]]
data[,2] <- (data[,2] + scaling_parameters_train$means[[2]]) *  scaling_parameters_train$sd[[2]]
data[,3] <- (data[,3] + scaling_parameters_train$means[[3]]) *  scaling_parameters_train$sd[[3]]
data[,4] <- (data[,4] + scaling_parameters_train$means[[4]]) *  scaling_parameters_train$sd[[4]]

while(counter <= ahead) {
  
  # change data to explanatory data set
  explanatory = as.matrix(
    data[
      (.N-(lags-1)):.N, # take last p rows
      sapply(
        0:(lags-1),
        function(lag) {
          data.table::shift(.SD, lag)
        }
      )
    ][.N,]
  )
  
  explanatory <- array_reshape(explanatory, dim = c(1,1,ncol(explanatory)))
  
  # predict with fitted values per step
  fitted_values <- lapply(
    1:length(nnvar_list),
    
    function(k) {
      
      mod <- nnvar_list[[k]]
      fitted <- mod(explanatory)
      y_hat <- as.numeric(fitted %>% tfd_mean())
      sd <- as.numeric(fitted %>% tfd_stddev())
      
      # rescale to the original mean and standard deviation
      y_hat <- (y_hat + scaling_parameters_train$means[[k]]) * scaling_parameters_train$sd[[k]]
      sd <- (sd + scaling_parameters_train$means[[k]]) * scaling_parameters_train$sd[[k]]
      
      return(list(y_hat = unlist(y_hat), sd = unlist(sd)))
    }
  )
  
  
  # transform the fitted values 
  y_hat <- matrix(sapply(fitted_values, function(i) i$y_hat), ncol = k)
  rownames(y_hat) <- NULL
  colnames(y_hat) <- names
  
  # transform the fitted values' variance
  sd <- matrix(sapply(fitted_values, function(i) i$sd), ncol = k)
  rownames(sd) <- NULL
  colnames(sd) <- names
  
  # reformat predictions and std ("uncertainty")
  predictions <- data.table::melt(
    data.table(y_hat),
    measure.vars = names)
  
  uc <- data.table::melt(
    data.table(sd),
    measure.vars = names)
  
  # update variables for next iteration
  forecast_next_step <- data.table::dcast(predictions, .~variable)[,-1]
  forecast <- rbind(forecast, forecast_next_step)
  
  uncertainty_next_step <- data.table::dcast(uc, .~variable)[,-1]
  uncertainty <- rbind(uncertainty, uncertainty_next_step)
  
  data <- rbind(data, forecast_next_step)
  counter <- counter + 1
  
}

# define steps ahead and counter
ahead <- v
counter <- 1

# define base set
forecast <- data.table::copy(tdf_train2[.N,])
uncertainty <- data.table::copy(tdf_train2[.N,])
uncertainty[] <- 0

# change for loop
data2 <- tdf_train2

# rescale
data2[,1] <- (data2[,1] + scaling_parameters_train2$means[[1]])  *  scaling_parameters_train2$sd[[1]]
data2[,2] <- (data2[,2] + scaling_parameters_train2$means[[2]]) *  scaling_parameters_train2$sd[[2]]
data2[,3] <- (data2[,3] + scaling_parameters_train2$means[[3]]) *  scaling_parameters_train2$sd[[3]]
data2[,4] <- (data2[,4] + scaling_parameters_train2$means[[4]]) *  scaling_parameters_train2$sd[[4]]

while(counter <= ahead) {
  
  # change data to explanatory data set
  explanatory = as.matrix(
    data2[
      (.N-(lags-1)):.N, # take last p rows
      sapply(
        0:(lags-1),
        function(lag) {
          data.table::shift(.SD, lag)
        }
      )
    ][.N,]
  )
  
  explanatory <- array_reshape(explanatory, dim = c(1,1,ncol(explanatory)))
  
  # predict with fitted values per step
  fitted_values2 <- lapply(
    1:length(nnvar_list),
    
    function(k) {
      
      mod <- nnvar_list[[k]]
      fitted <- mod(explanatory)
      y_hat2 <- as.numeric(fitted %>% tfd_mean())
      sd2 <- as.numeric(fitted %>% tfd_stddev())
      
      # rescale to the original mean and standard deviation
      y_hat2 <- (y_hat2 + scaling_parameters_train2$means[[k]]) * scaling_parameters_train2$sd[[k]]
      sd2 <- (sd2 + scaling_parameters_train2$means[[k]]) * scaling_parameters_train2$sd[[k]]
      
      return(list(y_hat2 = unlist(y_hat2), sd2 = unlist(sd2)))
    }
  )
  
  
  # transform the fitted values 
  y_hat2 <- matrix(sapply(fitted_values2, function(i) i$y_hat2), ncol = k)
  rownames(y_hat2) <- NULL
  colnames(y_hat2) <- names
  
  # transform the fitted values' variance
  sd2 <- matrix(sapply(fitted_values2, function(i) i$sd2), ncol = k)
  rownames(sd2) <- NULL
  colnames(sd2) <- names
  
  # reformat predictions and std ("uncertainty")
  predictions2 <- data.table::melt(
    data.table(y_hat2),
    measure.vars = names)
  
  uncertainty2 <- data.table::melt(
    data.table(sd2),
    measure.vars = names)
  
  # update variabels for next iteration
  forecast_next_step <- data.table::dcast(predictions2, .~variable)[,-1]
  forecast <- rbind(forecast, forecast_next_step)
  
  uncertainty_next_step <- data.table::dcast(uncertainty2, .~variable)[,-1]
  uncertainty <- rbind(uncertainty, uncertainty_next_step)
  
  data2 <- rbind(data2, forecast_next_step)
  counter <- counter + 1
  
}

# evaluate forecasts
eval_forecast <- xts(data, order.by = index(df_train2))
eval_forecast2 <- xts(data2, order.by = index(df))

val_pred <- eval_forecast[(t+1):(t+v),]
test_pred <- eval_forecast2[(t+v+1):(nrow(df)),]

# combine
nn_pred <- merge(df_train, val_pred, test_pred, df_val, df_test, all = TRUE)

colnames(nn_pred) <- c("GDP",
                       "CPI",
                       "Rate",
                       "HPI",
                       "Validation Forecast",
                       "Validation Forecast",
                       "Validation Forecast",
                       "Validation Forecast",
                       "Testing Forecast",
                       "Testing Forecast",
                       "Testing Forecast",
                       "TestingForecast",
                       "1",
                       "2",
                       "3",
                       "4",
                       "5",
                       "6",
                       "7",
                       "8")

copy_nn_pred <- nn_pred
copy_nn_pred[row,c(5,6,7,8)] <- copy_nn_pred[row,c(1,2,3,4)]
copy_nn_pred[row+ahead,c(9,10,11,12)] <- copy_nn_pred[row+ahead,c(13,14,15,16)]
copy_nn_pred[row,c(13,14,15,16)] <- copy_nn_pred[row,c(1,2,3,4)]
copy_nn_pred[row+ahead,c(17,18,19,20)] <- copy_nn_pred[row+ahead,c(13,14,15,16)]

plot(copy_nn_pred[,c(1,5,9)], main = "GDP",
     cex = 0.5, grid.ticks.lwd = 0.15, yaxis.left = TRUE, yaxis.right = FALSE, legend.loc = "topleft")
lines(copy_nn_pred[,c(13,17)], lwd=0.75, lty = 1, col='black')

plot(copy_nn_pred[,c(2,6,10)], main = "CPI",
     cex = 0.5, grid.ticks.lwd = 0.15, yaxis.left = TRUE, yaxis.right = FALSE, legend.loc = "topleft")
lines(copy_nn_pred[,c(14,18)], lwd=0.75, lty = 1, col='black')

plot(copy_nn_pred[,c(3,7,11)], main = "Rate",
     cex = 0.5, grid.ticks.lwd = 0.15, yaxis.left = TRUE, yaxis.right = FALSE, legend.loc = "topleft")
lines(copy_nn_pred[,c(15,19)], lwd=0.75, lty = 1, col='black')

plot(copy_nn_pred[,c(4,8,12)], main = "HPI",
     cex = 0.5, grid.ticks.lwd = 0.15, yaxis.left = TRUE, yaxis.right = FALSE, legend.loc = "topleft")
lines(copy_nn_pred[,c(16,20)], lwd=0.75, lty = 1, col='black')


# visual prediction evaluation
nn_f1 <- rbind(df_train, val_pred)
nn_f2 <- rbind(df_train, df_val, test_pred)
nn_cum_comp <- merge(df, nn_f1, nn_f2)

starting <- as.xts(matrix(c(1,1,1,1,1,1,1,1,1,1,1,1), nrow = 1, ncol = 12),
                   as.Date("2000-06-01"))

colnames(starting) <- colnames(nn_cum_comp)
nn_cum_comp <- rbind(starting, nn_cum_comp)
nn_cum_comp[1:nrow(df_train),5:8] <- NA
nn_cum_comp[1:(nrow(df_train)+nrow(df_val)),9:12] <- NA
colnames(nn_cum_comp) <- c("GDP",
                           "CPI",
                           "Rate",
                           "HPI",
                           "Validation Forecast",
                           "Validation Forecast",
                           "Validation Forecast",
                           "Validation Forecast",
                           "Testing Forecast",
                           "Testing Forecast",
                           "Testing Forecast",
                           "TestingForecast")

# RMSE NNVAR
v1 <- na.omit(nn_cum_comp[,c(1,5)])
v2 <- na.omit(nn_cum_comp[,c(2,6)])
v3 <- na.omit(nn_cum_comp[,c(3,7)])
v4 <- na.omit(nn_cum_comp[,c(4,8)])

v1_fe <- Metrics::rmse(v1$GDP, v1$`Validation Forecast`)
v2_fe <- Metrics::rmse(v2$CPI, v2$`Validation Forecast`)
v3_fe <- Metrics::rmse(v3$Rate, v3$`Validation Forecast`)
v4_fe <- Metrics::rmse(v4$HPI, v4$`Validation Forecast`)

v5 <- na.omit(nn_cum_comp[,c(1,9)])
v6 <- na.omit(nn_cum_comp[,c(2,10)])
v7 <- na.omit(nn_cum_comp[,c(3,11)])
v8 <- na.omit(nn_cum_comp[,c(4,12)])

v5_fe <- Metrics::rmse(v5$GDP, v5$`Testing Forecast`)
v6_fe <- Metrics::rmse(v6$CPI, v6$`Testing Forecast`)
v7_fe <- Metrics::rmse(v7$Rate, v7$`Testing Forecast`)
v8_fe <- Metrics::rmse(v8$HPI, v8$TestingForecast)

nn_fe_all <- matrix(c(v1_fe,
                      v2_fe,
                      v3_fe,
                      v4_fe,
                      v5_fe,
                      v6_fe,
                      v7_fe,
                      v8_fe),
                    ncol = 4,
                    nrow = 2,
                    byrow = TRUE)

# compare RSME
nn_fe_all
fe_all

#compute residuals

# accumulated squared errors
nn_residuals <- xts(res11, order.by = index(df_train[(lags+1):(nrow(df_train))]))
nn_sce <- merge(nn_cum_comp, nn_residuals)
nn_sce[,5] <- (nn_sce[,5] - nn_sce[,1])**2
nn_sce[,6] <- (nn_sce[,6] - nn_sce[,2])**2
nn_sce[,7] <- (nn_sce[,7] - nn_sce[,3])**2
nn_sce[,8] <- (nn_sce[,8] - nn_sce[,4])**2

nn_sce[,9] <- (nn_sce[,9] - nn_sce[,1])**2
nn_sce[,10] <- (nn_sce[,10] - nn_sce[,2])**2
nn_sce[,11] <- (nn_sce[,11] - nn_sce[,3])**2
nn_sce[,12] <- (nn_sce[,12] - nn_sce[,4])**2

nn_sce[,13] <- (nn_sce[,13])**2
nn_sce[,14] <- (nn_sce[,14])**2
nn_sce[,15] <- (nn_sce[,15])**2
nn_sce[,16] <- (nn_sce[,16])**2

nn_sce <- nn_sce[-nrow(nn_sce),-c(1,2,3,4)]
nn_sce[is.na(nn_sce)] <- 0
nn_sce$gdp <- nn_sce[,1] + nn_sce[,5] + nn_sce[,9]
nn_sce$cpi <- nn_sce[,2] + nn_sce[,6] + nn_sce[,10]
nn_sce$rate <- nn_sce[,3] + nn_sce[,7] + nn_sce[,11]
nn_sce$hpi <- nn_sce[,4] + nn_sce[,8] + nn_sce[,12]
nn_sce <- nn_sce[,c(9,10,11,12)]

starting <- as.xts(matrix(c(1,1,1,1), nrow = 1, ncol = 4),
                   as.Date("2000-06-01"))

colnames(starting) <- colnames(nn_sce)
cum_nn_sce <- rbind(starting, nn_sce)
cum_nn_sce <- cumsum(cum_nn_sce)

# plot cumulative squared error
plot(cum_nn_sce[,1], main = "GDP Squared Errors")
lines(cum_nn_sce[(t+1):(t+v),1], lwd=2, lty = 1, col='red')
lines(cum_nn_sce[(t+v):nrow(cum_nn_sce),1], lwd=2, lty = 1, col='lightgreen')

plot(cum_nn_sce[,2], main = "CPI Squared Errors")
lines(cum_nn_sce[(t+1):(t+v),2], lwd=2, lty = 1, col='red')
lines(cum_nn_sce[(t+v):nrow(cum_nn_sce),2], lwd=2, lty = 1, col='lightgreen')

plot(cum_nn_sce[,3], main = "Rate Squared Errors")
lines(cum_nn_sce[(t+1):(t+v),3], lwd=2, lty = 1, col='red')
lines(cum_nn_sce[(t+v):nrow(cum_nn_sce),3], lwd=2, lty = 1, col='lightgreen')

plot(cum_nn_sce[,4], main = "HPI Squared Errors")
lines(cum_nn_sce[(t+1):(t+v),4], lwd=2, lty = 1, col='red')
lines(cum_nn_sce[(t+v):nrow(cum_nn_sce),4], lwd=2, lty = 1, col='lightgreen')

ss_all <- merge(cum_sce50, cum_nn_sce50)
ss_all[4,] <- 1
ss_all <- ss_all[-(1:3),]
colnames(ss_all) <- c("GDP VAR", "CPI VAR", "Rate VAR", "HPI VAR",
                      "GDP NNVAR", "CPI NNVAR", "Rate NNVAR", "HPI NNVAR"
)
plot(ss_all[,c(1,5)], main = "GDP Squared Errors", legend.loc ="topleft")
plot(ss_all[,c(2,6)], main = "CPI Squared Errors", legend.loc ="topleft")
plot(ss_all[,c(3,7)], main = "Rate Squared Errors", legend.loc ="topleft")
plot(ss_all[,c(4,8)], main = "HPI Squared Errors", legend.loc ="topleft")



# ========== Chapter 3: Impulse Response ==========

# classic IRF in VAR model
gdp_shock <- vars::irf(model, impulse = "gdp", response = c("cpi"),
                       n.ahead = 20, boot = TRUE, ortho = TRUE, runs = 1000, ci = 0.95)
plot(gdp_shock)

cpi_shock <- vars::irf(model2, impulse = "cpi", response = c("hpi"),
                       n.ahead = 20, boot = TRUE, ortho = TRUE, runs = 1000, ci = 0.95)
plot(cpi_shock)

rate_shock <- vars::irf(model2, impulse = "rate", response = c("hpi"),
                        n.ahead = 20, boot = TRUE, ortho = TRUE, runs = 1000, ci = 0.95)
plot(rate_shock)

# NNVAR IRF

# reshape
X_val_nnIRF <- x_val
X_val_nnIRF[,,1] <- scaling_parameters_train$means[[1]]
X_val_nnIRF[,,2] <- scaling_parameters_train$means[[2]]
X_val_nnIRF[,,3] <- scaling_parameters_train$means[[3]]
X_val_nnIRF[,,4] <- scaling_parameters_train$means[[4]]
#X_val_nnIRF[,,5] <- means_train[,2]
#X_val_nnIRF[,,6] <- means_train[,3]
#X_val_nnIRF[,,7] <- means_train[,1] 
#X_val_nnIRF[,,8] <- means_train[,2]
#X_val_nnIRF[,,9] <- means_train[,3]
#X_val_nnIRF[,,10] <- means_train[,1] 
#X_val_nnIRF[,,11] <- means_train[,2]
#X_val_nnIRF[,,12] <- means_train[,3]

#X_test_nnIRF <- keras::array_reshape(X_test,c(nrow(X_test),1,ncol(X_test)))

# shock interest rate
X_val_nnIRF[1,1,1] <- scaling_parameters_train$means[[1]] + 10
#X_val_nnIRF[2,1,6] <- means_train[,2] + 1


# create lists for impulse responses
tval_nnIRF_1 <- list()
tval_nnIRF_2 <- list()
tval_nnIRF_3 <- list()
tval_nnIRF_4 <- list()

runs <- 100

for (i in seq(1,runs,1)) {
  
  tval_nnIRF_1[[i]] <- nnvar_fitted[[1]]$model %>% predict(X_val_nnIRF)
  tval_nnIRF_2[[i]] <- nnvar_fitted[[2]]$model %>% predict(X_val_nnIRF)
  tval_nnIRF_3[[i]] <- nnvar_fitted[[3]]$model %>% predict(X_val_nnIRF)
  tval_nnIRF_4[[i]] <- nnvar_fitted[[4]]$model %>% predict(X_val_nnIRF)
  
  #test_pred[[1]] <- fitted_nnvar2[[1]]$model %>% predict(X_test)
  #test_pred[[2]] <- fitted_nnvar2[[2]]$model %>% predict(X_test)
  #test_pred[[3]] <- fitted_nnvar2[[3]]$model %>% predict(X_test)
  
}

tval_nnIRF_1 <- data.frame(t(sapply(tval_nnIRF_1,c)))
tval_nnIRF_2 <- data.frame(t(sapply(tval_nnIRF_2,c)))
tval_nnIRF_3 <- data.frame(t(sapply(tval_nnIRF_3,c)))
tval_nnIRF_4 <- data.frame(t(sapply(tval_nnIRF_4,c)))

tval_nnIRF_1m <- data.frame(colMeans(tval_nnIRF_1)) #* sds_train[[1]] - means_train[[1]]
tval_nnIRF_2m <- data.frame(colMeans(tval_nnIRF_2)) #* sds_train[[2]] - means_train[[2]]
tval_nnIRF_3m <- data.frame(colMeans(tval_nnIRF_3)) #* sds_train[[3]] - means_train[[3]]
tval_nnIRF_4m <- data.frame(colMeans(tval_nnIRF_4)) #* sds_train[[3]] - means_train[[3]]

tval_nnIRF_1m["5%"] <- sapply(tval_nnIRF_1, quantile, p= 0.05) #* sds_train[[1]] - means_train[[1]]
tval_nnIRF_2m["5%"] <- sapply(tval_nnIRF_2, quantile, p= 0.05) #* sds_train[[2]] - means_train[[2]]
tval_nnIRF_3m["5%"] <- sapply(tval_nnIRF_3, quantile, p= 0.05) #* sds_train[[3]] - means_train[[3]]
tval_nnIRF_4m["5%"] <- sapply(tval_nnIRF_4, quantile, p= 0.05) #* sds_train[[3]] - means_train[[3]]

tval_nnIRF_1m["95%"] <- sapply(tval_nnIRF_1, quantile, p= 0.95) #* sds_train[[1]] - means_train[[1]]
tval_nnIRF_2m["95%"] <- sapply(tval_nnIRF_2, quantile, p= 0.95) #* sds_train[[2]] - means_train[[2]]
tval_nnIRF_3m["95%"] <- sapply(tval_nnIRF_3, quantile, p= 0.95) #* sds_train[[3]] - means_train[[3]]
tval_nnIRF_4m["95%"] <- sapply(tval_nnIRF_4, quantile, p= 0.95) #* sds_train[[3]] - means_train[[3]]

tval_nnIRF_1m <- xts(tval_nnIRF_1m, order.by = index(df_val[(2:nrow(df_val)),]))
tval_nnIRF_2m <- xts(tval_nnIRF_2m, order.by = index(df_val[(2:nrow(df_val)),]))
tval_nnIRF_3m <- xts(tval_nnIRF_3m, order.by = index(df_val[(2:nrow(df_val)),]))
tval_nnIRF_4m <- xts(tval_nnIRF_4m, order.by = index(df_val[(2:nrow(df_val)),]))

plot(tval_nnIRF_1m)
plot(tval_nnIRF_2m)
plot(tval_nnIRF_3m)
plot(tval_nnIRF_4m)

###########
### END ###
###########