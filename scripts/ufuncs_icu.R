################################################################################
##### utility functions for ICU length of stay case study

# load data
load_data <- function() {
  load("C:/Users/sa20i493/Documents/Data/icu_data/mdsi.rda")

  icu_vec <<- data$icuCode %>% unique()
  n_icu <<- icu_vec %>% length()

  ## split into estimation (val), calibration (tr), and test (ts) data
  data <- data %>% arrange(outDate)
  data_tr <- data %>% slice_head(n = round(nrow(data)*0.75))
  data_ts <- data %>% anti_join(data_tr, by = "id")
  data_val <- data_tr %>% sample_frac(size = 2/3)
  data_tr <- data_tr %>% anti_join(data_val, by = "id")

  t_vec <<- data_ts$los %>% quantile(c(0.1, 0.25, 0.5, 0.75, 0.9)) %>% unname()
  n_t <<- t_vec %>% length()
  N_ts <<- data_ts %>% nrow()

  time_meth <<- matrix(NA, nrow = n_icu, ncol = 3)
  colnames(time_meth) <<- c("lspm", "cidr", "mond")

  ## add index
  data_val$index <- NA
  data_tr$index <- NA
  data_ts$index <- NA
  for (icu in icu_vec) {
    print(icu)

    ### Get train data
    train <- subset(data_tr, icuCode == icu)
    val <- subset(data_val, icuCode == icu)
    test <- subset(data_ts, icuCode == icu)

    ### Get index
    out <- get_index(val, train, test)

    icu_ind <- data_val$icuCode == icu
    data_val[icu_ind, ] <- out$val

    icu_ind <- data_tr$icuCode == icu
    data_tr[icu_ind, ] <- out$train

    icu_ind <- data_ts$icuCode == icu
    data_ts[icu_ind, ] <- out$test
  }

  data_val <<- data_val
  data_tr <<- data_tr
  data_ts <<- data_ts
}

# function to get the index from the covariates
get_index <- function(val, train, test) {
  out <- tryCatch({
    fit <- mgcv::gam(log(los) ~ s(age, k = 3) + sex + planned + readmission + from + diag1 +
                       s(nems1, k = 3) + s(severity, k = 3) + interv,
                     data = val)
    val$index <- predict(fit, val) |> as.vector() |> exp()
    train$index <- predict(fit, train) |> as.vector() |> exp()
    test$index <- predict(fit, test) |> as.vector() |> exp()
    out <- list(val = val, train = train, test = test)
  },
  error = function(cond) {
    fit <- mgcv::gam(log(los) ~ s(age, k = 3) + sex + planned + readmission +
                       s(nems1, k = 3) + s(severity, k = 3),
                     data = val)
    val$index <- predict(fit, val) |> as.vector() |> exp()
    train$index <- predict(fit, train) |> as.vector() |> exp()
    test$index <- predict(fit, test) |> as.vector() |> exp()
    return(list(val = val, train = train, test = test))
  })
  return(out)
}

# initialise lists to store verification data
verif_lists <- function(N_ts, n_t, t_vec) {

  pcal <<- data.frame(replicate(3, numeric(N_ts)))
  colnames(pcal) <<- c("lspm", "cidr", "locb")
  score <<- pcal
  thick <<- pcal

  F_t <<- list(lspm = matrix(NA, N_ts, n_t),
               cidr = matrix(NA, N_ts, n_t),
               locb = matrix(NA, N_ts, n_t))
}

# perform cross validation to find the optimal number of bins at each station
conformal_binning_cv <- function(k_vec = c(1, seq(10, 100, 10))) {

  score_mat <- matrix(NA, nrow = n_icu, ncol = length(k_vec))
  for (j in seq_along(icu_vec)) {
    print(icu_vec[j])

    ### Get train data
    dat_icu <- data_val %>% subset(icuCode == icu_vec[j])
    train <- dat_icu %>% sample_frac(size = 0.5)
    val <- dat_icu %>% anti_join(train, by = "id")

    for (i in seq_along(k_vec)) {
      k <- k_vec[i]
      scores <- conformal_bin(x = train$index, y = train$los, x_out = val$index, y_out = val$los, x_est = dat_icu$index, y_est = dat_icu$los, k = k)
      score_mat[j, i] <- sapply(scores, function(x) x$crps) |> mean()
    }
  }
  k <- k_vec[apply(score_mat, 1, which.min)]

  return(k)
}

# wrapper to plot pit histograms
plot_pit_hists <- function(pit, score, filename = NULL) {
  lspm_plot <- pit_hist(pit[['lspm']], ranks = F, ymax = 0.4, xlab = NULL, xticks = F,
                        title = paste("LSPM: CRPS =", round(mean(score[['lspm']], na.rm = T), 3)))
  cidr_plot <- pit_hist(pit[['cidr']], ranks = F, ymax = 0.4, xlab = NULL, xticks = F,
                        title = paste("CIDR: CRPS =", round(mean(score[['cidr']], na.rm = T), 3)))
  locb_plot <- pit_hist(pit[['locb']], ranks = F, ymax = 0.4, xlab = NULL, xticks = F,
                        title = paste("CB: CRPS =", round(mean(score[['locb']], na.rm = T), 3)))
  pit_plot <- gridExtra::grid.arrange(lspm_plot, cidr_plot, locb_plot, nrow = 1)
  if (!is.null(filename)) {
    ggsave(plot = pit_plot, filename, width = 7.5, height = 2.5)
  } else {
    return(pit_plot)
  }
}

# plot PIT histograms at individual ICUs
plot_pit_hists_icu <- function(pit, score, ids = c(44, 65, 76, 77), filename = NULL) {
  code_vec <- paste0("ICU", ids)
  plot_list <- vector("list", 4)
  for (i in 1:4) {
    ind <- data_ts$icuCode == code_vec[i]
    plot_list[[i]] <- pit_hist(pit[ind], ranks = F, bins = 20, ymax = 0.08, xlab = NULL, xticks = F,
                               title = paste0(code_vec[i], ": CRPS = ", round(mean(score[ind], na.rm = T), 3)))
  }
  plot_list[['nrow']] <- 1
  pit_plot <- do.call(gridExtra::grid.arrange, plot_list)

  if (!is.null(filename)) {
    ggsave(plot = pit_plot, filename, width = 10, height = 2.5)
  }
}

# wrapper to plot pit pp-plots
plot_pit_pp <- function(pit, score, filename = NULL) {
  lspm_plot <- pit_reldiag(pit[['lspm']],
                           title = paste("LSPM: CRPS =", round(mean(score[['lspm']], na.rm = T), 3)))
  cidr_plot <- pit_reldiag(pit[['cidr']],
                           title = paste("CIDR: CRPS =", round(mean(score[['cidr']], na.rm = T), 3)))
  locb_plot <- pit_reldiag(pit[['lspm']],
                           title = paste("CB: CRPS =", round(mean(score[['locb']], na.rm = T), 3)))
  pit_plot <- gridExtra::grid.arrange(lspm_plot, cidr_plot, locb_plot, nrow = 1)
  if (!is.null(filename)) {
    ggsave(plot = pit_plot, filename, width = 1.3*7.5, height = 1.3*2.5)
  } else {
    return(pit_plot)
  }
}

# wrapper to plot tail calibration
plot_tcal <- function(F_t, ts_obs, t_vec, filename = NULL) {
  lspm_plot <- tc_reldiag(F_t[['lspm']], ts_obs, t_vec, xlab = "F(x)", ylab = "P(Y \u2264 x | F(x))", title = "LSPM")
  lspm_plot <- lspm_plot + guides(col = guide_legend(ncol = 2)) # split LSPM legend over two columns
  cidr_plot <- tc_reldiag(F_t[['cidr']], ts_obs, t_vec, xlab = "F(x)", ylab = "P(Y \u2264 x | F(x))", title = "CIDR")
  locb_plot <- tc_reldiag(F_t[['locb']], ts_obs, t_vec, xlab = "F(x)", ylab = "P(Y \u2264 x | F(x))", title = "CB")
  tc_plot <- gridExtra::grid.arrange(lspm_plot, cidr_plot, locb_plot, nrow = 1)
  if (!is.null(filename)) {
    ggsave(plot = tc_plot, filename, width = 1.3*7.5, height = 1.3*2.5)
  } else {
    return(tc_plot)
  }
}

# wrapper to plot thickness of conformal IDR bands
plot_thick <- function(thick, type = "hist", icu_tr = NULL, icu_ts = NULL, ylab = "Temperature", filename = NULL) {

  if (is.matrix(thick)) thick <- as.vector(thick)

  if (type == "hist") {
    ## histogram of thicknesses
    plot_obj <- ggplot(data.frame(x = thick)) +
      geom_histogram(aes(x = x, y = after_stat(count) / sum(after_stat(count))),
                     boundary = 0, binwidth = 0.025, fill = "lightgrey", col = "darkgrey") +
      scale_x_continuous(name = "Thickness", limits = c(0, 1), expand = c(0, 0)) +
      scale_y_continuous(name = "Relative frequency", expand = expansion(c(0, 0.15))) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            plot.margin = margin(c(5.5, 10.5, 5.5, 5.5)))
    if (!is.null(filename)) {
      ggsave(filename, width = 3.7, height = 2.5)
    } else {
      return(plot_obj)
    }

  } else if (type == "scatter") {
    ## scatter plot vs sample size

    if (!identical(levels(icu_tr), levels(icu_ts))) stop("levels of 'icu_tr' and 'icu_ts' should be the same")
    df <- sapply(levels(icu_tr), function(icu) {
      n <- sum(icu_tr == icu)
      th <- mean(thick[icu_ts == icu])
      c(n = n, th = th)
    }) |> t() |> as.data.frame()

    plot_obj <- ggplot(df) +
      geom_point(aes(x = n, y = th)) +
      scale_x_continuous(name = "Sample size") +
      scale_y_continuous(name = "Average thickness") +
      theme_bw() +
      theme(panel.grid = element_blank())

    if (!is.null(filename)) {
      ggsave(filename, width = 3.7, height = 2.5)
    } else {
      return(plot_obj)
    }

  } else {
    stop("'type' must be one of 'traffic', 'hist', and 'scatter'")
  }

}

# plot examples of predictive cdf's
plot_example <- function(cidr_preds, lspm_preds, locb_preds, filename = NULL) {

  ind <- 359 #25
  ind <- sample(1:ncol(cidr_preds$points), 1) #440/1202/1708/2335/4577/4997
  x <- seq(0, 10, 0.01)
  lspm_cdf <- stepfun(lspm_preds$points[-1, ind], lspm_preds$cdf_crisp[, ind])(x)
  cidr_cdf <- stepfun(cidr_preds$points[-1, ind], cidr_preds$cdf_crisp[, ind])(x)
  locb_cdf <- stepfun(locb_preds[[ind]]$points[-1], locb_preds[[ind]]$cdf_crisp)(x)
  df <- data.frame(x = x, y = c(lspm_cdf, cidr_cdf, locb_cdf), mth = rep(c("LSPM", "CIDR", "CB"), each = length(x)))
  plot_obj <- ggplot(df) + geom_step(aes(x = x, y = y, col = mth)) +
    scale_x_continuous(name = "LOS (hours)", limits = c(0, 10)) +
    scale_y_continuous(name = "Predictive CDF") +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = "bottom")
  plot_obj

  if (!is.null(filename)) {
    ggsave(plot = plot_obj, filename, width = 3.2, height = 2.7)
  } else {
    return(plot_obj)
  }

}

