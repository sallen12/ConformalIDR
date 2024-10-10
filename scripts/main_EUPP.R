################################################################################
## set up

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

load_data(na_prop = 0) # stations with any missing data are removed

## plot predicted vs observed temperature
plot_pred(filename = "plots/EUMN_data.png")

## plot stations
plot_map(lons, lats, rowMeans(ts_obs), filename = "plots/EUMN_stations.png")

## thresholds at which to evaluate threshold calibration
t_vec <- tr_obs |> quantile(c(0.1, 0.25, 0.5, 0.75, 0.9)) |> unname()

## initialise lists to store verification data
verif_lists(ts_obs, t_vec)


################################################################################
## prediction (seasonal)

## optimal number of bins for local binning
k <- local_binning_cv()

## fit models
for (j in seq_along(stat_ids)) {
  st <- stat_ids[j]
  for (i in 1:4) {
    s <- c("Wi", "Sp", "Su", "Au")[i]
    print(paste0('Forecast at Station: ', st, ' (', j, ' from ', length(stat_ids), ') and Season: ', s))

    ### Get train data
    seas_ind <- tr_seas == s
    train <- data.frame(obs = tr_obs[j, seas_ind], ens.mu = tr_fc_mn[j, seas_ind])
    seas_ind <- ts_seas == s
    test <- data.frame(obs = ts_obs[j, seas_ind], ens.mu = ts_fc_mn[j, seas_ind])

    ### LSPM
    lspm_preds <- conformal_lspm(x = train$ens.mu, y = train$obs, x_out = test$ens.mu, y_out = test$obs)
    pcal[['lspm']][j, seas_ind] <- lspm_preds$pit
    score[['lspm']][j, seas_ind] <- lspm_preds$crps
    thick[['lspm']][j, seas_ind] <- lspm_preds$thick
    F_t[['lspm']][j, seas_ind, ] <- threshcal(lspm_preds, t_vec)

    ### CIDR
    cidr_preds <- conformal_idr(x = train$ens.mu, y = train$obs, x_out = test$ens.mu, y_out = test$obs)
    pcal[['cidr']][j, seas_ind] <- cidr_preds$pit
    score[['cidr']][j, seas_ind] <- cidr_preds$crps
    thick[['cidr']][j, seas_ind] <- cidr_preds$thick
    F_t[['cidr']][j, seas_ind, ] <- threshcal(cidr_preds, t_vec)

    ### LB
    locb_preds <- conformal_bin(x = train$ens.mu, y = train$obs, x_out = test$ens.mu, y_out = test$obs, k[i, j])
    pcal[['locb']][j, seas_ind] <- sapply(locb_preds, function(x) x$pit)
    score[['locb']][j, seas_ind] <- sapply(locb_preds, function(x) x$crps)
    thick[['locb']][j, seas_ind] <- sapply(locb_preds, function(x) x$thick)
    F_t[['locb']][j, seas_ind, ] <- sapply(locb_preds, function(x) threshcal(x, t_vec)) |> t()
  }
}
rm(i, j, s, st, seas_ind, train, test, lspm_preds, cidr_preds, locb_preds, scores)


################################################################################
## prediction (rolling)

roll_index <- function(i, win_len, tr_times, ts_times) {
  ind <- (i - win_len):(i + win_len)
  roll_times <- ts_times[i] + days(c(-win_len, win_len))
  tr_ind <- logical(length(tr_times))
  for (k in -2:2) {
    roll_times_mod <- roll_times + years(k)
    c <- 0
    while (any(is.na(roll_times_mod))) {
      c <- c + 1
      roll_times_mod <- ts_times[i] + days(c(-win_len-c, win_len+c)) + years(k)
    }
    tr_ind <- tr_ind | (tr_times >= roll_times_mod[1]  & tr_times <= roll_times_mod[2])
  }
  tr_ind <- replicate(20, tr_ind) |> t() |> as.vector()
  return(tr_ind)
}

## fit models
win_len <- 45
for (j in 1) {#seq_along(stat_ids)) {
  st <- stat_ids[j]
  for (i in seq_along(ts_times)) {
    print(paste0('Forecast at Station: ', st, ' (', j, ' from ', length(stat_ids), ') and Day: ', i))

    ### Get train data
    #tr_ind <- lubridate::month(tr_times) %% 12
    #tr_ind <- tr_ind %in% ((lubridate::month(ts_times[i]) + -1:1) %% 12)
    #tr_ind <- replicate(20, tr_ind) |> t() |> as.vector()

    tr_ind <- roll_index(i, win_len, tr_times, ts_times)
    train <- data.frame(obs = tr_obs[j, tr_ind], ens.mu = tr_fc_mn[j, tr_ind])
    test <- data.frame(obs = ts_obs[j, i], ens.mu = ts_fc_mn[j, i])

    ### CIDR
    cidr_preds <- conformal_idr(x = train$ens.mu, y = train$obs, x_out = test$ens.mu)
    pcal[['cidr']][j, i] <- pit(cidr_preds, test$obs)
    score[['cidr']][j, i] <- crps(cidr_preds, test$obs)
    F_t[['cidr']][j, i, ] <- threshcal(cidr_preds, t_vec)
    thick[['cidr']][j, i] <- thickness(cidr_preds)

  }
}
rm(i, j, tr_ind, train, test, cidr_preds)



################################################################################
## results

## PIT histograms
plot_pit_hists(pcal, score, filename = "plots/EUMN_pitrd_comp.png")

## PIT pp-plots
plot_pit_pp(pcal, score, filename = "plots/EUMN_pit_comp.png")

## Threshold calibration diagrams
plot_tcal(F_t, ts_obs, t_vec, filename = "plots/EUMN_tcal_comp.png")

## Thickness
st <- 1
th_all <- thick[['cidr']]
th_loc <- th_all[st, ]
plot_thick(th_loc, type = "traffic", obs = ts_obs[st, ], times = ts_times, filename = "plots/EUMN_thick_obs_ts.png")
plot_thick(th_loc, type = "traffic", obs = score[['cidr']][st, ], times = ts_times, ylab = "CRPS", filename = "plots/EUMN_thick_crps_ts.png")
plot_thick(th_all, type = "hist", filename = "plots/EUMN_thick_cidr.png")
plot_thick(th_loc, type = "scatter", x = ts_fc_mn[st, ], filename = "plots/EUMN_thick_ens.png")

