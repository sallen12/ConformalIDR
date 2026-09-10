################################################################################
## ICU length of stay case study

set.seed(783)

library(ConformalIDR)
library(dplyr)
library(ggplot2)
library(zoo)
library(scoringRules)
library(WeightedForecastVerification)
library(quantreg)
library(splines)

source("scripts/ufuncs_icu.R")


################################################################################
## load data

# load processed data
load("scripts/data/dat_icu.RData")

## initialise lists to store verification data
verif_lists(N_ts, n_t, t_vec)


################################################################################
## prediction

# path to save the results
results_path <- "scripts/results/results_data_icu.RData"

if (file.exists(results_path)) {
  # load results if they have already been generated
  load("scripts/results/results_data_icu.RData")
} else {
  ## optimal k for conformal binning
  k <- conformal_binning_cv()

  ## fit conformal predictive systems
  x_ts <- matrix(NA, N_ts, 3)
  colnames(x_ts) <- c("index", "obs", "icu")
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

    ### Store index and observation (for evaluation)
    x_ts[icu_ind, ] <- cbind(test$index, test$los, test$icuCode)

    ## LSPM
    start <- Sys.time()
    lspm_preds <- conformal_lspm(x = train$index, y = train$los, x_out = test$index, y_out = test$los)
    pcal[['lspm']][icu_ind] <- lspm_preds$pit
    score[['lspm']][icu_ind] <- lspm_preds$crps
    thick[['lspm']][icu_ind] <- lspm_preds$thick
    F_t[['lspm']][icu_ind, ] <- threshcal(lspm_preds, t_vec)
    cover[['lspm']][icu_ind, ] <- sapply(a_vec, function(a) coverage(lspm_preds, alpha = a, average = F))
    intsc[['lspm']][icu_ind, ] <- sapply(a_vec, function(a) int_score(lspm_preds, alpha = a))
    widt[['lspm']][icu_ind, ] <- sapply(a_vec, function(a) width(lspm_preds, alpha = a))
    end <- Sys.time()
    time_meth[which(icu_vec == icu), 1] <- difftime(end, start, units = "mins")

    ### CIDR
    start <- Sys.time()
    cidr_preds <- conformal_idr(x = train$index, y = train$los, x_out = test$index, y_out = test$los)
    pcal[['cidr']][icu_ind] <- cidr_preds$pit
    score[['cidr']][icu_ind] <- cidr_preds$crps
    thick[['cidr']][icu_ind] <- cidr_preds$thick
    F_t[['cidr']][icu_ind, ] <- threshcal(cidr_preds, t_vec)
    cover[['cidr']][icu_ind, ] <- sapply(a_vec, function(a) coverage(cidr_preds, alpha = a, average = F))
    intsc[['cidr']][icu_ind, ] <- sapply(a_vec, function(a) int_score(cidr_preds, alpha = a))
    widt[['cidr']][icu_ind, ] <- sapply(a_vec, function(a) width(cidr_preds, alpha = a))
    end <- Sys.time()
    time_meth[which(icu_vec == icu), 2] <- difftime(end, start, units = "mins")

    ### LB
    start <- Sys.time()
    locb_preds <- conformal_bin(x = train$index, y = train$los, x_out = test$index, y_out = test$los, x_est = val$index, y_est = val$los, k = k[i])
    pcal[['locb']][icu_ind] <- sapply(locb_preds, function(x) x$pit)
    score[['locb']][icu_ind] <- sapply(locb_preds, function(x) x$crps)
    thick[['locb']][icu_ind] <- sapply(locb_preds, function(x) x$thick)
    F_t[['locb']][icu_ind, ] <- sapply(locb_preds, function(x) threshcal(x, t_vec)) |> t()
    cover[['locb']][icu_ind, ] <- sapply(locb_preds, function(x) sapply(a_vec, function(a) coverage(x, alpha = a, average = F))) |> t()
    intsc[['locb']][icu_ind, ] <- sapply(locb_preds, function(x) sapply(a_vec, function(a) int_score(x, alpha = a))) |> t()
    widt[['locb']][icu_ind, ] <- sapply(locb_preds, function(x) sapply(a_vec, function(a) width(x, alpha = a))) |> t()
    end <- Sys.time()
    time_meth[which(icu_vec == icu), 3] <- difftime(end, start, units = "mins")

    print(time_meth[which(icu_vec == icu), ])
  }
  preds <- list(cidr = cidr_preds, lspm = lspm_preds, locb = locb_preds) # save last predictions


  ## fit conformalized quantile regression
  cover$cqr <- intsc$cqr <- widt$cqr <- matrix(NA, N_ts, n_a)
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
    val <- out$val
    n_cal <- nrow(train)

    # train on estimation data (val)
    for (j in seq_along(a_vec)) {
      a <- a_vec[j]
      print(a)

      # train on estimation data (val)
      fit_low <- rq(los ~ ns(index, df = 3), tau = a/2, data = val)
      fit_upp <- rq(los ~ ns(index, df = 3), tau = 1 - a/2, data = val)

      # get predictions on calibration data (train)
      cal_low <- predict(fit_low, newdata = data.frame(index = train$index))
      cal_upp <- predict(fit_upp, newdata = data.frame(index = train$index))

      cal_low <- pmax(cal_low, 0)
      cal_upp <- pmax(cal_upp, 0)

      # calculate scores
      scores <- pmax(cal_low - train$los, train$los - cal_upp)
      k_cqr <- ceiling((n_cal + 1) * (1 - a))
      q_hat <- sort(scores)[k_cqr]

      # get predictions on test data (test)
      pred_low <- predict(fit_low, newdata = data.frame(index = test$index))
      pred_upp <- predict(fit_upp, newdata = data.frame(index = test$index))

      pred_low <- pmax(pred_low, 0)
      pred_upp <- pmax(pred_upp, 0)

      # get prediction intervals on test data
      l <- pred_low - q_hat
      u <- pred_upp + q_hat
      l <- pmax(pmin(l, u), 0)
      u <- pmax(pmax(l, u), 0)

      # evaluate
      cover$cqr[icu_ind, j] <- as.numeric(test$los >= l & test$los <= u)
      widt$cqr[icu_ind, j] <- (u - l)
      intsc$cqr[icu_ind, j] <- scoringRules::ints_quantiles(test$los, l, u, 1 - a)
    }
  }
  ## save data
  save(pcal, score, F_t, thick, cover, intsc, widt, k, preds, x_ts, file = results_path)
}


################################################################################
## results

## PIT pp-plots for crisp CDF
plot_pit_pp(pcal, score, filename = "plots/ICU_pitrd_comp.png")

## PIT histograms for crisp CDF
plot_pit_hists(pcal, score, filename = "plots/ICU_pit_comp.png")

## PIT histograms for crisp CDF at ICUs
plot_pit_hists_icu(pcal[['cidr']], score[['cidr']], filename = "plots/ICU_pit_comp_ind.png")

## Threshold calibration diagrams for crisp CDF
plot_tcal(F_t, data_ts$los, t_vec, filename = "plots/ICU_tcal_comp.png")

## Thickness
plot_thick(thick[['cidr']], filename = "plots/ICU_thick_cidr.png")
plot_thick(thick[['cidr']], type = "scatter", icu_tr = data_tr$icuCode, icu_ts = data_ts$icuCode, filename = "plots/ICU_thick_n.png")

## Example crisp CDFs
plot_example(preds$cidr, preds$lspm, preds$locb, filename = "plots/ICU_dist_ex.png")

## Unconditional coverage of prediction intervals
plot_cov_unc(cover, a_vec, filename = "plots/ICU_cov_unc.png")

## Conditional coverage of prediction intervals
plot_cov_con(cover, x_ts = x_ts[, 1], icu = x_ts[, 3], a_vec, filename = "plots/ICU_cov_con.png")

## Average width of prediction intervals
plot_width(widt, a_vec, filename = "plots/ICU_wid.png")

## Interval score of prediction intervals
plot_is(intsc, a_vec, filename = "plots/ICU_intsc.png")



