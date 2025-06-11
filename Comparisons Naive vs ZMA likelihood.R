library(rugarch)
library(xts)

year <- daily_rv_estimators_1year_with_log3

year$date <- as.Date(year$date)
year$return <- c(NA, diff(year$log_close))
year <- na.omit(year)

rets_xts     <- xts(year$return,     order.by = year$date)
rv_5min_xts <- xts(year$RV_5min,   order.by = year$date)
rv_zma_xts   <- xts(year$RV_ZMA,     order.by = year$date)

base.spec <- ugarchspec(
  variance.model     = list(model = "realGARCH", garchOrder = c(1,1)),
  mean.model         = list(armaOrder = c(1,0)),
  distribution.model = "norm"
)

roll_length <- 150
refit_every <- 10

roll.naive <- ugarchroll(
  spec            = base.spec,
  data            = rets_xts,
  n.ahead         = 1,
  forecast.length = roll_length,
  refit.every     = refit_every,
  refit.window    = "moving",
  calculate.VaR    = TRUE,
  VaR.alpha       = c(0.01),
  realizedVol     = rv_5min_xts
)

roll.zma <- ugarchroll(
  spec            = base.spec,
  data            = rets_xts,
  n.ahead         = 1,
  forecast.length = roll_length,
  refit.every     = refit_every,
  refit.window    = "moving",
  calculate.VaR    = TRUE,
  VaR.alpha       = c(0.01),
  realizedVol     = rv_zma_xts
)

fc_5min <- as.data.frame(roll.naive@forecast)
fc_zma   <- as.data.frame(roll.zma@forecast)

sigma_5min <- fc_5min[,"density.Sigma"]
r_5min     <- fc_5min[,"density.Realized"]

sigma_zma   <- fc_zma[,"density.Sigma"]
r_zma       <- fc_zma[,"density.Realized"]

logscore_5min <- dnorm(r_5min,
                        mean = 0,
                        sd   = sigma_5min,
                        log  = TRUE)

logscore_zma   <- dnorm(r_zma,
                        mean = 0,
                        sd   = sigma_zma,
                        log  = TRUE)

cumll_naive <- sum(logscore_5min, na.rm = TRUE)
cumll_zma   <- sum(logscore_zma,   na.rm = TRUE)

cat("Cumulative out-of-sample log-likelihood:\n",
    "  Naive RV →", round(cumll_naive,2), "\n",
    "  ZMA RV   →", round(cumll_zma,2),   "\n")

d_t <- logscore_5min - logscore_zma

mfit <- lm(d_t ~ 1)

n     <- length(d_t)
m     <- coef(mfit)[1]
se_iid<- sd(d_t) / sqrt(n)
Z_iid <- m / se_iid
p_iid <- 2 * (1 - pnorm(abs(Z_iid)))

cat("IID Vuong‐style Z =", round(Z_iid,3),
    " p-value =", round(p_iid,4), "\n")

library(sandwich)
library(lmtest)

nw_se <- sqrt(NeweyWest(mfit, prewhite = FALSE)[1,1])
Z_nw  <- m / nw_se
p_nw  <- 2 * (1 - pnorm(abs(Z_nw)))

cat("NW‐Vuong Z =", round(Z_nw,3),
    " p-value =", round(p_nw,4), "\n")
