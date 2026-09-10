################################################################################
##### EUPPBench temperature case study

set.seed(91)

library(ConformalIDR)
library(ncdf4)
library(ggplot2)
library(zoo)
library(scoringRules)
library(WeightedForecastVerification)
library(lubridate)

source("scripts/ufuncs_eupp.R")


################################################################################
## load data

# load processed data
load("scripts/data/dat_eupp.RData")

## plot predicted vs observed temperature
plot_pred(filename = "plots/EUMN_data.png")

## plot stations
plot_map(lons, lats, rowMeans(ts_obs), filename = "plots/EUMN_stations_alt.png")

## thresholds at which to evaluate threshold calibration
t_vec <- tr_obs |> quantile(c(0.1, 0.25, 0.5, 0.75, 0.9)) |> unname()

## levels at which to evaluate prediction intervals
a_vec <- seq(0.05, 0.95, 0.05)

## initialise lists to store verification data
verif_lists(ts_obs, t_vec, a_vec)

################################################################################
## prediction

# path to save the results
results_path <- "scripts/results/results_data_eupp.RData"

if (file.exists(results_path)) {
  # load results if they have already been generated
  load(results_path)
} else {
  ## fit models

  # number of bins for conformal binning
  k <- 20

  for (j in seq_along(stat_ids)) { # fit models separately for each station
    st <- stat_ids[j]
    for (i in seq_along(ts_times)) { # rolling training window
      print(paste0('Forecast at Station: ', st, ' (', j, ' from ', length(stat_ids), ') and Day: ', i))

      ### Get train data
      tr_ind <- roll_index(i, win_len = 45, tr_times, ts_times)
      train <- data.frame(obs = tr_obs[j, tr_ind], ens.mu = tr_fc_mn[j, tr_ind])
      test <- data.frame(obs = ts_obs[j, i], ens.mu = ts_fc_mn[j, i])


      ### LSPM
      lspm_preds <- conformal_lspm(x = train$ens.mu, y = train$obs, x_out = test$ens.mu, y_out = test$obs)
      pcal[['lspm']][j, i] <- lspm_preds$pit
      crpsc[['lspm']][j, i] <- lspm_preds$crps
      thick[['lspm']][j, i] <- lspm_preds$thick
      F_t[['lspm']][j, i, ] <- threshcal(lspm_preds, t_vec)
      cover[['lspm']][j, i, ] <- sapply(a_vec, function(a) coverage(lspm_preds, alpha = a, average = F))
      intsc[['lspm']][j, i, ] <- sapply(a_vec, function(a) int_score(lspm_preds, alpha = a))
      widt[['lspm']][j, i, ] <- sapply(a_vec, function(a) width(lspm_preds, alpha = a))


      ### CIDR
      cidr_preds <- conformal_idr(x = train$ens.mu, y = train$obs, x_out = test$ens.mu, y_out = test$obs)
      pcal[['cidr']][j, i] <- cidr_preds$pit
      crpsc[['cidr']][j, i] <- cidr_preds$crps
      thick[['cidr']][j, i] <- cidr_preds$thick
      F_t[['cidr']][j, i, ] <- threshcal(cidr_preds, t_vec)
      cover[['cidr']][j, i, ] <- sapply(a_vec, function(a) coverage(cidr_preds, alpha = a, average = F))
      intsc[['cidr']][j, i, ] <- sapply(a_vec, function(a) int_score(cidr_preds, alpha = a))
      widt[['cidr']][j, i, ] <- sapply(a_vec, function(a) width(cidr_preds, alpha = a))


      ### LB
      locb_preds <- conformal_bin(x = train$ens.mu, y = train$obs, x_out = test$ens.mu, y_out = test$obs, k = k)
      pcal[['locb']][j, i] <- sapply(locb_preds, function(x) x$pit)
      crpsc[['locb']][j, i] <- sapply(locb_preds, function(x) x$crps)
      thick[['locb']][j, i] <- sapply(locb_preds, function(x) x$thick)
      F_t[['locb']][j, i, ] <- sapply(locb_preds, function(x) threshcal(x, t_vec)) |> t()
      cover[['locb']][j, i, ] <- sapply(locb_preds, function(x) sapply(a_vec, function(a) coverage(x, alpha = a, average = F))) |> t()
      intsc[['locb']][j, i, ] <- sapply(locb_preds, function(x) sapply(a_vec, function(a) int_score(x, alpha = a))) |> t()
      widt[['locb']][j, i, ] <- sapply(locb_preds, function(x) sapply(a_vec, function(a) width(x, alpha = a))) |> t()

    }
  }
  # save data
  save(pcal, crpsc, F_t, thick, cover, intsc, widt, k, file = results_path)
  rm(i, j, st, tr_ind, train, test, lspm_preds, cidr_preds, locb_preds)
}


################################################################################
## results

## PIT histograms
plot_pit_hists(pcal, crpsc, filename = "plots/EUMN_pit_comp.png")

## PIT pp-plots
plot_pit_pp(pcal, crpsc, filename = "plots/EUMN_pitrd_comp.png")

## Threshold calibration diagrams
plot_tcal(F_t, ts_obs, t_vec, filename = "plots/EUMN_tcal_comp.png")

## Thickness
st <- sample(seq_along(stat_ids), 1)
th_all <- thick[['cidr']]
th_loc <- th_all[st, ]
plot_thick(th_loc, type = "traffic", obs = ts_obs[st, ], times = ts_times, filename = "plots/EUMN_thick_obs_ts.png")
plot_thick(th_loc, type = "traffic", obs = crpsc[['cidr']][st, ], times = ts_times, ylab = "CRPS", filename = "plots/EUMN_thick_crps_ts.png")
plot_thick(th_all, type = "hist", filename = "plots/EUMN_thick_cidr.png")
plot_thick(th_loc, type = "scatter", x = ts_fc_mn[st, ], filename = "plots/EUMN_thick_ens.png")

## Unconditional coverage
plot_cov_unc(cover, a_vec, filename = "plots/EUMN_cov_unc.png")

## Conditional coverage
plot_cov_con(cover, a_vec, x_ts = ts_fc_mn, filename = "plots/EUMN_cov_con.png")

## Average width
plot_width(widt, a_vec, filename = "plots/EUMN_wid.png")

## Interval score
plot_is(intsc, a_vec, filename = "plots/EUMN_intsc.png")

