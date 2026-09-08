################################################################################
## set up

set.seed(90743892)

library(ConformalIDR)
library(WeightedForecastVerification)
library(ggplot2)
source("scripts/ufuncs_sim.R")


################################################################################
## set up

N_ts <- 5000
x_ts <- runif(N_ts, 0, 10)

# isotonic example
y_ts <- rgamma(N_ts, shape = sqrt(x_ts), scale = pmin(pmax(x_ts, 1), 6))

# antitonic example
y_ts_at <- rgamma(N_ts, shape = sqrt(10 - x_ts), scale = pmin(pmax(10 - x_ts, 1), 6))

# less isotonic example
y_ts_le <- rnorm(N_ts, mean = 2*x_ts + 5*sin(x_ts), sd = 0.2*x_ts)

# plot example data
plot_example(x_ts[1:1000], y_ts[1:1000], filename = "plots/simstudy_data.png")
plot_example(x_ts[1:1000], y_ts_at[1:1000], filename = "plots/simstudy_data_at.png")
plot_example(x_ts[1:1000], y_ts_le[1:1000], ylims = c(0, 30), filename = "plots/simstudy_data_le.png")


################################################################################
## isotonic

##### results

res_100 <- get_results(N_tr = 100, x_ts, y_ts) # 30 secs
res_500 <- get_results(N_tr = 500, x_ts, y_ts) # 2 mins
res_1000 <- get_results(N_tr = 1000, x_ts, y_ts) # 4 mins
res_2000 <- get_results(N_tr = 2000, x_ts, y_ts) # 20 mins


##### evaluate

## probabilistic calibration of crisp CDF
plot_cal(res_100$pit, res_100$crps, filename = "plots/simstudy_100.png")
plot_cal(res_500$pit, res_500$crps, filename = "plots/simstudy_500.png")
plot_cal(res_1000$pit, res_1000$crps, filename = "plots/simstudy_1000.png")
plot_cal(res_2000$pit, res_2000$crps, filename = "plots/simstudy_2000.png")
plot_cal(res_2000$pit, res_2000$crps, vert = F, filename = "plots/simstudy_2000_ho.png")

## threshold calibration of crisp CDF
plot_cal(F_t = res_2000$F_t, obs = y_ts, type = "threshcal", filename = "plots/simstudy_2000_tc.png")

## unconditional coverage of prediction intervals
plot_cov(res_2000$cover, filename = "plots/simstudy_2000_cov.png")

## conditional coverage of prediction intervals
plot_cov(res_2000$cover, x_ts = x_ts, filename = "plots/simstudy_2000_ccov.png")

## interval score of prediction intervals
plot_is(res_2000$intsc, filename = "plots/simstudy_2000_is.png")

## unconditional average width of prediction intervals
plot_width(res_2000$widt, filename = "plots/simstudy_2000_wid.png")

## conditional average width of prediction intervals
plot_width(res_2000$widt, x_ts = x_ts, filename = "plots/simstudy_2000_cwid.png", type = "cond")

## average width of prediction intervals
plot_width(lapply(list(res_100, res_500, res_1000, res_2000), function(x) x$widt),
           filename = "plots/simstudy_2000_wid_samp.png", type = "samp")


################################################################################
## antitonic

##### results

res_100_at <- get_results(N_tr = 100, x_ts, y_ts_at, type = "anti")
res_500_at <- get_results(N_tr = 500, x_ts, y_ts_at, type = "anti")
res_1000_at <- get_results(N_tr = 1000, x_ts, y_ts_at, type = "anti")
res_2000_at <- get_results(N_tr = 2000, x_ts, y_ts_at, type = "anti")


##### evaluate

## probabilistic calibration of crisp CDF
plot_cal(res_100_at$pit, res_100_at$crps, filename = "plots/simstudy_100_at.png")
plot_cal(res_500_at$pit, res_500_at$crps, filename = "plots/simstudy_500_at.png")
plot_cal(res_1000_at$pit, res_1000_at$crps, filename = "plots/simstudy_1000_at.png")
plot_cal(res_2000_at$pit, res_2000_at$crps, filename = "plots/simstudy_2000_at.png")

## threshold calibration of crisp CDF
plot_cal(F_t = res_1000_at$F_t, obs = y_ts_at, type = "threshcal", filename = "plots/simstudy_1000_tc_at.png")

## unconditional coverage of prediction intervals
plot_cov(res_2000_at$cover, filename = "plots/simstudy_2000_cov_at.png")

## conditional coverage of prediction intervals
plot_cov(res_2000_at$cover, x_ts = x_ts, filename = "plots/simstudy_2000_ccov_at.png")

## interval score of prediction intervals
plot_is(res_2000_at$intsc, filename = "plots/simstudy_2000_is_at.png")

## unconditional average width of prediction intervals
plot_width(res_2000_at$widt, filename = "plots/simstudy_2000_wid_at.png")

## conditional average width of prediction intervals
plot_width(res_2000_at$widt, x_ts = x_ts, filename = "plots/simstudy_2000_cwid_at.png", type = "cond")

## average width of prediction intervals
plot_width(lapply(list(res_100_at, res_500_at, res_1000_at, res_2000_at), function(x) x$widt),
           filename = "plots/simstudy_2000_wid_samp_at.png", type = "samp")


################################################################################
## less isotonic

##### results

res_100_le <- get_results(N_tr = 100, x_ts, y_ts_le, type = "less")
res_500_le <- get_results(N_tr = 500, x_ts, y_ts_le, type = "less")
res_1000_le <- get_results(N_tr = 1000, x_ts, y_ts_le, type = "less")
res_2000_le <- get_results(N_tr = 2000, x_ts, y_ts_le, type = "less")


##### evaluate

## probabilistic calibration of crisp CDF
plot_cal(res_100_le$pit, res_100_le$crps, filename = "plots/simstudy_100_le.png")
plot_cal(res_500_le$pit, res_500_le$crps, filename = "plots/simstudy_500_le.png")
plot_cal(res_1000_le$pit, res_1000_le$crps, filename = "plots/simstudy_1000_le.png")
plot_cal(res_2000_le$pit, res_2000_le$crps, filename = "plots/simstudy_2000_le.png")
plot_cal(res_2000_le$pit, res_2000_le$crps, vert = F, filename = "plots/simstudy_2000_le_ho.png")

## threshold calibration of crisp CDF
plot_cal(F_t = res_2000_le$F_t, obs = y_ts_le, type = "threshcal", filename = "plots/simstudy_2000_tc_le.png")

## unconditional coverage of prediction intervals
plot_cov(res_2000_le$cover, filename = "plots/simstudy_2000_cov_le.png")

## conditional coverage of prediction intervals
plot_cov(res_2000_le$cover, x_ts = x_ts, filename = "plots/simstudy_2000_ccov_le.png")

## interval score of prediction intervals
plot_is(res_2000_le$intsc, filename = "plots/simstudy_2000_is_le.png")

## unconditional average width of prediction intervals
plot_width(res_2000_le$widt, filename = "plots/simstudy_2000_wid_le.png")

## conditional average width of prediction intervals
plot_width(res_2000_le$widt, x_ts = x_ts, filename = "plots/simstudy_2000_cwid_le.png", type = "cond")

## average width of prediction intervals
plot_width(lapply(list(res_100_le, res_500_le, res_1000_le, res_2000_le), function(x) x$widt),
           filename = "plots/simstudy_2000_wid_samp_le.png", type = "samp")


################################################################################
## compare thicknesses

df <- data.frame(thicc = c(res_100$thick$cidr, res_100_at$thick$cidr, res_100_le$thick$cidr),
                 mode = rep(c(" Isotonic", "Antitonic", " Less isotonic"), each = N_ts))
ggplot(df) + geom_boxplot(aes(x = mode, y = thicc)) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Thickness") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave("plots/simstudy_thick.png", width = 5, height = 3)


