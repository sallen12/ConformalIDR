################################################################################
## set up

set.seed(91)

library(ConformalIDR)
library(ncdf4)
library(ggplot2)
library(zoo)
library(scoringRules)
library(WeightedForecastVerification)

source("scripts/utility_funcs.R")


################################################################################
## load data

load_temp_data(na_prop = 0) # stations with any missing data are removed

## plot predicted vs observed temperature
plot_temp_pred(filename = "plots/EUMN_data.png")

## plot stations
plot_temp_map(lons, lats, rowMeans(ts_obs), filename = "plots/EUMN_stations.png")

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
    lspm_preds <- fit_lspm(y = train$obs, X = train$ens.mu, X_ts = test$ens.mu)
    scores <- eval_lspm(lspm_preds, test$obs, t_vec)
    pit[['lspm']][j, seas_ind] <- scores$pit
    score[['lspm']][j, seas_ind] <- scores$crps
    F_t[['lspm']][j, seas_ind, ] <- scores$F_t
    thick[['lspm']][j, seas_ind] <- scores$thick

    ### CIDR
    cidr_preds <- fit_cidr(y = train$obs, X = train$ens.mu, X_ts = test$ens.mu)
    scores <- eval_cidr(cidr_preds, test$obs, t_vec)
    pit[['cidr']][j, seas_ind] <- scores$pit
    score[['cidr']][j, seas_ind] <- scores$crps
    F_t[['cidr']][j, seas_ind, ] <- scores$F_t
    thick[['cidr']][j, seas_ind] <- scores$thick

    ### LB
    locb_preds <- fit_locb(y = train$obs, X = train$ens.mu, X_ts = test$ens.mu, k[i_s, j])
    scores <- eval_locb(locb_preds, test$obs, t_vec)
    pit[['locb']][j, seas_ind] <- scores$pit
    score[['locb']][j, seas_ind] <- scores$crps
    F_t[['locb']][j, seas_ind, ] <- scores$F_t
    thick[['locb']][j, seas_ind] <- scores$thick
  }
}
rm(i, j, s, st, seas_ind, train, test, lspm_preds, cidr_preds, locb_preds, scores)


################################################################################
## results

## PIT histograms
plot_pit_hists(pit, score, filename = "plots/EUMN_pitrd_comp.png")

## PIT pp-plots
plot_pit_pp(pit, score, filename = "plots/EUMN_pit_comp.png")

## Threshold calibration diagrams
plot_tcal(F_t, ts_obs, t_vec, filename = "plots/EUMN_tcal_comp.png")

## Thickness
st <- 1
th_all <- thick[['cidr']]
th_loc <- th_all[st, ]
plot_thick(th_loc, type = "traffic", obs = ts_obs[st, ], times = ts_times, filename = "plots/EUMN_thick_obs_ts.png")
plot_thick(th_all, type = "hist", filename = "plots/EUMN_thick_cidr.png")
plot_thick(th_loc, type = "scatter", x = ts_fc_mn[st, ], filename = "plots/EUMN_thick_ens.png")

