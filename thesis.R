library(fGarch)
library(astsa)
library(locfit)

## Here is a demo of how `garchFit` function works
## `fit` is in class `fGARCH`, the estimated parameters and log likelihood can be called in the following ways

x = rnorm(100)
fit = garchFit( ~ garch(1,1),data = x, cond.dist = "QMLE", trace = F)
unname(fit@fit$llh)
fit@fit$params$params[c("mu","alpha1","beta1")]

# Calculate Log-Likelihood

llh = function(sigmat, dat)
{
  return(-0.5 * sum(log(sigmat^2) + dat^2 / sigmat^2))
}

# Use Bootstrap to Obtain the Critical Values

## The theory is:
## 1. Assume we know the parameters, generate random sequence
## 3. Use the random sequence, we find the maximal difference of likelihood
## 4. Repeat step 1 ~ 2 multiple times, obtain the 95% percentile value as critical value

# This function returns the critical value with specific parameters

crit_value_1 = function(n, k, alpha, beta, omega)
# `n` is the total length of the inverval, `k` is the position in the interval;
# `alpha`, `beta`, `omega` are parameters.
{
  set.seed(12345)
  
  # Simulated Series
  diff = numeric(length = 200)
  for(i in seq(200))
  {
    simu = garchSim(spec = garchSpec
                    (model = list(omega = omega, alpha = alpha, beta = beta)),
                    n = n)
    temp_1 = simu[1:(n-k)]; temp_2 = simu[(n - k + 1):n]
    
    # Adapative and original model log-likelihoods
    ## New Fit
    lh_adapt_1 = garchFit(~ garch(1,1), data = temp_1, cond.dist = "QMLE", trace = F)@fit$llh
    lh_adapt_2 = garchFit(~ garch(1,1), data = temp_2, cond.dist = "QMLE", trace = F)@fit$llh
    lh_adapt = garchFit(~ garch(1,1), data = simu, cond.dist = "QMLE", trace = F)@fit$llh
    diff[i] = -2 * (lh_adapt_1 + lh_adapt_2 - lh_adapt)
  }
  return(sort(diff)[190])
}

# Calculating critical value among all possible parameters (95%)

gene_crit = function(n,params)
{
  alpha = params[1]; beta = params[2]; omega = params[3]
  pos = seq(n/10,9*n/10,n/10)
  tem = vector()
  for (i in seq(9)) {
    tem = c(tem,crit_value_1(n,pos[i],alpha,beta,omega))
    }
  max_crit = max(tem)
  return(max_crit)
}

## The problem of parametric bootstrap is that the process works too slow. 
## Therefore I came up an idea to simplify the process by making linear regression 
## to obtain the critical value under different circumstances. (parameters and length)
## Which turns out to be effective and approved by the advisor, saying that is the pattern used in industry.

params = matrix(data = c(0.1,0.8,0.1,0.2,0.7,0.1,0.3,0.6,0.1,0.7,0.1,0.1),
              byrow = T, nrow = 4, ncol = 3) # Some reasonable combinations of parameters

result = matrix(NA, nrow = 5, ncol = 4)
for(i in seq(2,5)){
  n = 2^(i-1) * 50
  for(j in seq(4)){
    result[i,j] = gene_crit(n = n, params = para[j,])
  }
}

# Screen plot
para = params[1,]
alpha = params[1]; beta = params[2]; omega = params[3]
pos = seq(10,90,10)
param_1 = vector(0)
for(i in seq(9)) {param_1 = c(param_1,crit_value_1(100,pos[i], alpha, beta, omega))}
## Same as param_2, param_3, param_4

par(mfrow = c(2,2))
plot(pos, param_1, pch = 20, main = "alpha = 0.1, beta = 0.8", ylab = "95% test statistics")
lines(ksmooth(pos, param_1, "normal",bandwidth = 20),col = 2)
plot(pos, param_2, pch = 20, main = "alpha = 0.2, beta = 0.7",ylab = "")
lines(ksmooth(pos, param_2, "normal",bandwidth = 20),col = 2)
plot(pos, param_3, pch = 20, main = "alpha = 0.3, beta = 0.6",ylab = "95% test statistics")
lines(ksmooth(pos, param_3, "normal",bandwidth = 20),col = 2)
plot(pos, param_4, pch = 20, main = "alpha = 0.7, beta = 0.1",ylab = "")
lines(ksmooth(pos, param_4, "normal",bandwidth = 20),col = 2)

row.names(result) = c("n = 50","n = 100", "n = 200", "n = 400", "n = 800")
colnames(result) = c("alpha = 0.1, beta = 0.8","alpha = 0.2, beta = 0.7",
                     "alpha = 0.3, beta = 0.6","alpha = 0.7, beta = 0.1")
result = t(result)
res = as.vector(result)
n = c(rep(50,4), rep(100,4), rep(200,4), rep(400,4), rep(800,4)) # length of the interval
a = rep(c(0.1,0.2,0.3,0.7), 5) # alpha
b = rep(c(0.8,0.7,0.6,0.1), 5) # beta
fit_lm = lm(res ~ n + a + b)
summary(fit_lm)

# Split the entire interval to multiple subintervals

max_K = function(length)
{
  K = floor(log(length / 10) / log(2))
  Ik = length
  for (i in seq(0, K))
  {
    Ik = c(Ik, ceiling(length - 10 * 2^i))
  }
  Ik = c(Ik, 1)
  return(Ik)
}

# Read in data and pre-processing
stock = read.csv("stock_data.csv", sep = ",")[c("Date","Open")]
stock$Date = as.Date(stock$Date, "%Y-%m-%d")

## Plot (stock price ~ date)
plot(Open ~ Date,stock,xaxt = "n",type = "l",ylab = "Index")
axis(1, stock$Date, labels = stock$Date[c(1,366,731,1097)],
     at = stock$Date[c(1,366,731,1097)], cex.axis = .7, las = 1)

## Plot (stock price returns ~ date)
price = log(stock[,2])
data = (embed(price,2)[,2] - embed(price,2)[,1])
plot(data ~ stock$Date[-1],xlab = "Date", xaxt = "n",type = "l", ylab = "Stock Index Return")
axis(1, stock$Date, labels = stock$Date[c(1,366,731,1097)],
     at = stock$Date[c(1,366,731,1097)], cex.axis = .7, las = 1)


# The following code shows three ways to fit the data:
# 1. Parametric model: The most common way to obtain the model, use entire sequence to fit the parameters and make forecasts
# 2. Stepwise parametric model: The reasonable way to obtain the model, make forecast only use data before the particular date
# 3. Adaptive parametric model: The new suggested algorithm, only use part of the former data to build forecast model

## Parametric GARCH model
fit_data = garchFit(~ garch(1,1),data = data, cond.dist = "QMLE", trace = F)
para_data = unname(fit_data@fit$params$params[c("omega","alpha1","beta1")])
omega_data = para_data[1]; alpha_data = para_data[2]; beta_data = para_data[3]
pred_data = fit_data@h.t

# Main function (Running the stepwise parametric and the new adaptive algorithm)
main = function(pos)
{
  data_main = data[1:pos]
  
  ## Fit the original model
  fit_orig = garchFit(~ garch(1,1), data = data_main, cond.dist = "QMLE", trace = F)
  para_orig = unname(fit_orig@fit$params$params[c("omega","alpha1","beta1")])
  omega_orig = para_orig[1]; alpha_orig = para_orig[2]; beta_orig = para_orig[3]
  
  ## Make predictions based on stepwise parametric method (named as `pred_orig`)
  pred_orig = omega_orig + alpha_orig * data[pos]^2 + beta_orig * tail(fit_orig@h.t,1)
  
  ## Adaptive parametric method
  step = max_K(pos)
  
  ### Divide the sequence into several subintervals
  for (i in seq(length(step)-1))
  {
    assign(paste("a",i,sep = ""), data_main[(step[i+1]+1):step[i]])
  }
  for(k in seq(3,length(step)-1))
  {
    test = get(paste("a",k,sep="")); len = length(test); diff = vector()
    breakpoint = -10.31 + 0.0021*len + 24.48*alpha_orig + 20.73*beta_orig # Use the fitted linear model to calculate critical value
    if (len >= 10) # In case the serie is too short to make reasonable estimation
    {for (j in seq(5,len-5))
      {
        main_lha_1 = garchFit(~ garch(1,1), data = test[1:j], cond.dist = "QMLE", trace = F)@fit$llh
        main_lha_2 = garchFit(~ garch(1,1), data = test[(j+1):len], cond.dist = "QMLE",trace = F)@fit$llh
        main_lho = garchFit(~ garch(1,1), data = test, cond.dist = "QMLE", trace = F)@fit$llh
        diff[j-4] = -2 * (main_lha_1+main_lha_2-main_lho)
      }
      max_diff = max(diff)
      if (max_diff > breakpoint) break}
  }
  
  ## If a breakpoint is detected, use the corresponding data, if not, use the entire data sequence
  data_new = data[step[k]:pos]
  fit_new = garchFit(~ garch(1,1), data = data_new, cond.dist = "QMLE", trace = F)
  para_new = unname(fit_new@fit$params$params[c("omega","alpha1","beta1")])
  omega_new = para_new[1]; alpha_new = para_new[2]; beta_new = para_new[3]
  
  ## Make predictions based on adaptive parametric method (named as `pred_new`)
  pred_new = omega_new + alpha_new * data[pos]^2 + beta_new * tail(fit_new@h.t,1)
  return(c(pred_orig,pred_new))
}

len = length(data)
main_pred = matrix(NA, ncol = 2, nrow = len - 50) # Taper the first 50 data as "training set"

## Predict at every position
for(i in 1:(len - 50))
{
  tryCatch({
    main_pred[i,] = main(i + 49)
  }, error = function(e){})
}
tempo2 = main_pred[1:(len - 50),]

result_raw = cbind(pred_data[51:len], main_pred)
colnames(result_raw) = c("data","orig","new")
write.csv(result_raw, file = "result.csv", row.names = F)

## Window smoothing (window size: 30)
size = 30
na = which(is.na(main_pred[,1]) == 1)
eff_len = len - 50 - length(na)
na_filter = c(1:(size/2 - 1), (eff_len - size/2 - 1):eff_len)

## Calculate prediction errors
error = function(pred)
{
  se = (pred[51:len] - data[51:len]^2)^2
  se = se[-na] # Remove NA values
  se = as.vector(filter(se, rep(1/30,30)))[-na_filter] # Window Smoothing
}
se_data = error(pred_data) # Parametric
se_new = error(main_pred[,2]) # Adaptive
se_orig = error(main_pred[,1]) # Stepwise

par(mfrow = c(3,1))
par(mar = c(3,4,2,1))

plot(pred_date, se_new/se_orig,xaxt = "n", ylab = "Ratio of errors",
        main = "Local GARCH to Stepwise GARCH", lty = 1)
plot(ksmooth(pred_date, se_new/se_data, "normal",bandwidth = 40), type="l", col = 2) # Kernel smoothing
axis(1, pred_date, labels = pred_date[c(1,366,731,1044)],
     at = pred_date[c(1,366,731,1044)], cex.axis = .7, las = 1)
abline(h = 1)
plot(pred_date, se_new/se_data,xaxt = "n", type = "l", ylab = "Ratio of errors",
     main = "Local GARCH to Stepwise GARCH", lty = 2)
axis(1, pred_date, labels = pred_date[c(1,366,731,1044)],
     at = pred_date[c(1,366,731,1044)], cex.axis = .7, las = 1)
abline(h = 1)
plot(pred_date, se_orig/se_data,xaxt = "n", type = "l", ylab = "Ratio of errors",
     main = "Local GARCH to Stepwise GARCH", lty = 3)
axis(1, pred_date, labels = pred_date[c(1,366,731,1044)],
     at = pred_date[c(1,366,731,1044)], cex.axis = .7, las = 1)
abline(h = 1)


# Nonparametric models

## Local linear regression
## 1. Use the whole sequence as known data (parametric)
library(locfit)
x = seq(length(data))
lo1 = locfit(data ~ lp(x, nn = 0.5, h = 0.5))
lo_data = unname(predict(lo1, data.frame(x = 51:len)))
se_lo1 = error(lo_data)

## 2. Use only former data (stepwise)
lo2_data = vector()
for (i in seq(1208))
{
  data_temp = data[1:(49+i)]; x_temp = seq(length(data_temp))
  lo2 = locfit(data_temp ~ lp(x_temp, nn = 0.5, h = 0.5))
  lo2_data[i] = unname(predict(lo2, data.frame(x_temp = (50+i))))
}
se_lo2 = error(lo2_data)

pred_date = (stock$Date + 50)[51:len]
pred_date = pred_date[-na]

par(mfrow = c(2,1))
par(mar = c(3,4,2,1))
plot(pred_date, se_lo1/se_data, xaxt = "n", #type = "l", 
     xlab = "Date", ylab = "Ratio of Errors", 
     main = "Local Linear to Parametric GARCH", cex.main = 1)
axis(1, pred_date, labels = pred_date[c(1,366,731,1044)],
     at = pred_date[c(1,366,731,1044)], cex.axis = .7, las = 1)
abline(h=1)
plot(pred_date, se_lo2/se_data, xaxt = "n", 
     type = "l", lty = 2, xlab = "Date", ylab = "Ratio of Errors")
axis(1, pred_date, labels = pred_date[c(1,366,731,1044)],
     at = pred_date[c(1,366,731,1044)], cex.axis = .7, las = 1)
abline(h=1)
