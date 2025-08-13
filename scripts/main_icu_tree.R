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

# wrapper to find the optimal hyperparameters for conformal binning
local_binning_cv_tree <- function(cp_vec = c(1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1)) {

  score_mat <- matrix(NA, nrow = n_icu, ncol = length(cp_vec))
  for (j in seq_along(icu_vec)) {
    print(icu_vec[j])

    ### Get train data
    dat_icu <- data_val %>% subset(icuCode == icu_vec[j])
    train <- dat_icu %>% sample_frac(size = 0.5)
    val <- dat_icu %>% anti_join(train, by = "id")

    for (i in seq_along(cp_vec)) {
      cp <- cp_vec[i]
      scores <- conformal_bin(x = train$index, y = train$los,
                              x_out = val$index, y_out = val$los,
                              x_est = dat_icu$index, y_est = dat_icu$los,
                              binning = "tree", cp = cp)
      score_mat[j, i] <- sapply(scores, function(x) x$crps) |> mean()
    }
  }
  cp <- cp_vec[apply(score_mat, 1, which.min)]

  return(cp)
}

# wrapper to plot pit pp-plots
plot_pit_pp_tree <- function(pit, score, filename = NULL) {
  lspm_plot <- pit_reldiag(pit[['lspm']],
                           title = paste("LSPM: CRPS =", round(mean(score[['lspm']], na.rm = T), 3)))
  cidr_plot <- pit_reldiag(pit[['cidr']],
                           title = paste("CIDR: CRPS =", round(mean(score[['cidr']], na.rm = T), 3)))
  locb_plot <- pit_reldiag(pit[['lspm']],
                           title = paste("CB (kmeans): CRPS =", round(mean(score[['locb']], na.rm = T), 3)))
  tree_plot <- pit_reldiag(pit[['tree']],
                           title = paste("CB (tree): CRPS =", round(mean(score[['tree']], na.rm = T), 3)))
  pit_plot <- gridExtra::grid.arrange(lspm_plot, cidr_plot, locb_plot, tree_plot, nrow = 1)
  if (!is.null(filename)) {
    ggsave(plot = pit_plot, filename, width = 1.7*7.5, height = 1.3*2.5)
  } else {
    return(pit_plot)
  }
}

# wrapper to plot threshold calibration
plot_tcal_tree <- function(F_t, ts_obs, t_vec, filename = NULL) {
  lspm_plot <- tc_reldiag(F_t[['lspm']], ts_obs, t_vec, xlab = "F(x)", ylab = "P(Y \u2264 x | F(x))", title = "LSPM")
  lspm_plot <- lspm_plot + guides(col = guide_legend(ncol = 2)) # split LSPM legend over two columns
  cidr_plot <- tc_reldiag(F_t[['cidr']], ts_obs, t_vec, xlab = "F(x)", ylab = "P(Y \u2264 x | F(x))", title = "CIDR")
  locb_plot <- tc_reldiag(F_t[['locb']], ts_obs, t_vec, xlab = "F(x)", ylab = "P(Y \u2264 x | F(x))", title = "CB (kmeans)")
  tree_plot <- tc_reldiag(F_t[['tree']], ts_obs, t_vec, xlab = "F(x)", ylab = "P(Y \u2264 x | F(x))", title = "CB (tree)")
  tc_plot <- gridExtra::grid.arrange(lspm_plot, cidr_plot, locb_plot, tree_plot, nrow = 1)
  if (!is.null(filename)) {
    ggsave(plot = tc_plot, filename, width = 1.7*7.5, height = 1.3*2.5)
  } else {
    return(tc_plot)
  }
}


################################################################################
## load data

load_data()

load("scripts/results/results_data_icu.RData")

pcal$tree <- score$tree <- thick$tree <- NA
F_t[['tree']] <- array(NA, dim(F_t[['lspm']]))


################################################################################
## prediction

## optimal k for conformal binning
cp <- local_binning_cv_tree()

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

  ### LB tree
  start <- Sys.time()
  locb_preds <- conformal_bin(x = train$index, y = train$los, x_out = test$index,
                              y_out = test$los, x_est = val$index, y_est = val$los,
                              method = "tree", cp = cp[i])
  pcal[['tree']][icu_ind] <- sapply(locb_preds, function(x) x$pit)
  score[['tree']][icu_ind] <- sapply(locb_preds, function(x) x$crps)
  thick[['tree']][icu_ind] <- sapply(locb_preds, function(x) x$thick)
  F_t[['tree']][icu_ind, ] <- sapply(locb_preds, function(x) threshcal(x, t_vec)) |> t()
  end <- Sys.time()

  print(end - start)
}
preds[['tree']] <- locb_preds # save last predictions


################################################################################
## results

## save data

save(pcal, score, F_t, thick, k, cp, preds, file = "scripts/results/results_data_icu_tree.RData")
#load("scripts/results/results_data_icu_tree.RData")

## PIT pp-plots
plot_pit_pp_tree(pcal, score, filename = "plots/ICU_pitrd_tree.png")

## Threshold calibration diagrams
plot_tcal_tree(F_t, data_ts$los, t_vec, filename = "plots/ICU_tcal_tree.png")


