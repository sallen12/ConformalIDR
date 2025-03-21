################################################################################
## set up

set.seed(783)

library(ConformalIDR)
library(dplyr)
library(ggplot2)
library(zoo)
library(scoringRules)
library(WeightedForecastVerification)

source("scripts/ufuncs_icu.R")


################################################################################
## load data

load_data()

## initialise lists to store verification data
verif_lists(N_ts, n_t, t_vec)


################################################################################
## prediction

## optimal k for conformal binning
k <- local_binning_cv()

## fit models
for (i in seq_along(icu_vec)) {
  icu <- icu_vec[i]
  print(icu)

  ### Get train data
  train <- subset(data_tr, icuCode == icu)
  val <- subset(data_val, icuCode == icu)
  test <- subset(data_ts, icuCode == icu)
  icu_ind <- data_ts$icuCode == icu

  ### Get index from estimation set
  out <- get_index(val, train, test)
  train <- out$train
  test <- out$test

  ## LSPM
  start <- Sys.time()
  lspm_preds <- conformal_lspm(x = train$index, y = train$los, x_out = test$index, y_out = test$los)
  pcal[['lspm']][icu_ind] <- lspm_preds$pit
  score[['lspm']][icu_ind] <- lspm_preds$crps
  thick[['lspm']][icu_ind] <- lspm_preds$thick
  F_t[['lspm']][icu_ind, ] <- threshcal(lspm_preds, t_vec)
  end <- Sys.time()
  time_meth[which(icu_vec == icu), 1] <- difftime(end, start, units = "mins")

  ### CIDR
  start <- Sys.time()
  cidr_preds <- conformal_idr(x = train$index, y = train$los, x_out = test$index, y_out = test$los)
  pcal[['cidr']][icu_ind] <- cidr_preds$pit
  score[['cidr']][icu_ind] <- cidr_preds$crps
  thick[['cidr']][icu_ind] <- cidr_preds$thick
  F_t[['cidr']][icu_ind, ] <- threshcal(cidr_preds, t_vec)
  end <- Sys.time()
  time_meth[which(icu_vec == icu), 2] <- difftime(end, start, units = "mins")

  ### LB
  start <- Sys.time()
  locb_preds <- conformal_bin(x = train$index, y = train$los, x_out = test$index, y_out = test$los, x_est = val$index, y_est = val$los, k = k[i])
  pcal[['locb']][icu_ind] <- sapply(locb_preds, function(x) x$pit)
  score[['locb']][icu_ind] <- sapply(locb_preds, function(x) x$crps)
  thick[['locb']][icu_ind] <- sapply(locb_preds, function(x) x$thick)
  F_t[['locb']][icu_ind, ] <- sapply(locb_preds, function(x) threshcal(x, t_vec)) |> t()
  end <- Sys.time()
  time_meth[which(icu_vec == icu), 3] <- difftime(end, start, units = "mins")

  print(time_meth[which(icu_vec == icu), ])
}
preds <- list(cidr = cidr_preds, lspm = lspm_preds, locb = locb_preds) # save last predictions


################################################################################
## results

## save data

save(pcal, score, F_t, thick, k, preds, file = "scripts/results/results_data_icu.RData")
#load("scripts/results/results_data_icu.RData")

## PIT histograms
plot_pit_hists_icu(pcal$cidr, score$cidr, filename = "plots/ICU_pitrd_comp_ind.png")

## PIT pp-plots
plot_pit_pp(pcal, score, filename = "plots/ICU_pitrd_comp.png")

## PIT histograms
plot_pit_hists(pcal, score, filename = "plots/ICU_pit_comp.png")

## PIT histograms at ICUs
plot_pit_hists_icu(pcal[['cidr']], score[['cidr']], filename = "plots/ICU_pit_comp_ind.png")

## Threshold calibration diagrams
plot_tcal(F_t, data_ts$los, t_vec, filename = "plots/ICU_tcal_comp.png")

## Thickness
plot_thick(thick[['cidr']], filename = "plots/ICU_thick_cidr.png")
plot_thick(thick[['cidr']], type = "scatter", icu_tr = data_tr$icuCode, icu_ts = data_ts$icuCode, filename = "plots/ICU_thick_n.png")

## Example CDFs
plot_example(preds$cidr, preds$lspm, preds$locb, filename = "plots/ICU_dist_ex.png")

