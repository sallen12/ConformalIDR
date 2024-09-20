################################################################################
## set up

set.seed(783)

library(ConformalIDR)
library(dplyr)
library(ggplot2)
library(zoo)
library(scoringRules)
library(WeightedForecastVerification)

source("scripts/utility_funcs.R")


################################################################################
## load data

load_icu_data()

## initialise lists to store verification data
verif_icu_lists(N_ts, n_t, t_vec)


################################################################################
## prediction (seasonal)

## optimal k for Mondrian prediction
k <- local_binning_icu_cv()

## fit models
for (i in seq_along(icu_vec)) {
  icu <- icu_vec[i]
  print(icu)

  ### Get train data
  train <- subset(data_tr, icuCode == icu)
  val <- subset(data_val, icuCode == icu)
  test <- subset(data_ts, icuCode == icu)
  icu_ind <- data_ts$icuCode == icu

  ### Get index
  out <- get_index(val, train, test)
  train <- out$train
  test <- out$test

  ## LSPM
  start <- Sys.time()
  lspm_preds <- fit_lspm(train$los, train$index, test$index)
  scores <- eval_lspm(lspm_preds, test$los, t_vec)
  pit[['lspm']][icu_ind] <- scores$pit
  score[['lspm']][icu_ind] <- scores$crps
  thick[['lspm']][icu_ind] <- scores$thick
  F_t[['lspm']][icu_ind, ] <- scores$F_t
  end <- Sys.time()
  time_meth[which(icu_vec == icu), 1] <- difftime(end, start, units = "mins")

  ### CIDR
  start <- Sys.time()
  cidr_preds <- fit_cidr(train$los, train$index, test$index)
  scores <- eval_cidr(cidr_preds, test$los, t_vec)
  pit[['cidr']][icu_ind] <- scores$pit
  score[['cidr']][icu_ind] <- scores$crps
  thick[['cidr']][icu_ind] <- scores$thick
  F_t[['cidr']][icu_ind, ] <- scores$F_t
  end <- Sys.time()
  time_meth[which(icu_vec == icu), 2] <- difftime(end, start, units = "mins")

  ### LB
  start <- Sys.time()
  locb_preds <- fit_locb(train$los, train$index, test$index, k[i])
  scores <- eval_locb(locb_preds, test$los, t_vec)
  pit[['locb']][icu_ind] <- scores$pit
  score[['locb']][icu_ind] <- scores$crps
  thick[['locb']][icu_ind] <- scores$thick
  F_t[['locb']][icu_ind, ] <- scores$F_t
  end <- Sys.time()
  time_meth[which(icu_vec == icu), 3] <- difftime(end, start, units = "mins")

  print(time_meth[which(icu_vec == icu), ])
}


################################################################################
## results

## PIT histograms
plot_pit_hists(pit, score, filename = "plots/ICU_pitrd_comp.png")

## PIT pp-plots
plot_pit_pp(pit, score, filename = "plots/ICU_pit_comp.png")

## PIT histograms at ICUs
plot_pit_hists_icu(pit[['cidr']], score[['cidr']], filename = "plots/ICU_pit_comp_ind.png")

## Threshold calibration diagrams
plot_tcal(F_t, data_ts$los, t_vec, filename = "plots/ICU_tcal_comp.png")

## Thickness
plot_thick(thick[['cidr']], type = "hist", filename = "plots/ICU_thick_cidr.png")

## Example CDFs
plot_icu_example(cidr_preds, lspm_preds, filename = "plots/ICU_dist_ex.png")

