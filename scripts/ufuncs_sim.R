################################################################################
## simulation study helper functions


# wrapper to fit and evaluate the different conformal predictive systems and CQR
get_results <- function(N_tr = 500, x_ts, y_ts, type = c("iso", "anti", "less"), k = 10) {
  start <- Sys.time()
  type <- match.arg(type)

  x_tr <- runif(N_tr, 0, 10)
  if (type == "anti") {
    y_tr <- rgamma(N_tr, shape = sqrt(10 - x_tr), scale = pmin(pmax(10 - x_tr, 1), 6))
  } else if (type == "less") {
    y_tr <- rnorm(N_tr, mean = 2*x_tr + 5*sin(x_tr), sd = 0.2*x_tr)
  } else {
    y_tr <- rgamma(N_tr, shape = sqrt(x_tr), scale = pmin(pmax(x_tr, 1), 6))
  }

  N_ts <- length(y_ts)
  pcal <- data.frame(replicate(3, numeric(N_ts)))
  colnames(pcal) <- c("lspm", "cidr", "locb")
  thick <- crps <- pcal

  t_vec <- unname(quantile(y_ts, c(0.1, 0.25, 0.5, 0.75, 0.9)))
  n_t <- length(t_vec)
  F_t <- list(lspm = matrix(NA, N_ts, n_t),
              cidr = matrix(NA, N_ts, n_t),
              locb = matrix(NA, N_ts, n_t))

  a_vec <- seq(0.05, 0.95, 0.05)
  n_a <- length(a_vec)
  cover <- intsc <- widt <- list(lspm = matrix(NA, N_ts, n_a),
                                 cidr = matrix(NA, N_ts, n_a),
                                 locb = matrix(NA, N_ts, n_a))

  ### LSPM
  lspm_preds <- conformal_lspm(x = x_tr, y = y_tr, x_out = x_ts, y_out = y_ts)
  pcal[['lspm']] <- lspm_preds$pit
  crps[['lspm']] <- lspm_preds$crps
  thick[['lspm']] <- lspm_preds$thick
  F_t[['lspm']] <- threshcal(lspm_preds, t_vec)
  cover[['lspm']] <- sapply(a_vec, function(a) coverage(lspm_preds, alpha = a, average = F))
  intsc[['lspm']] <- sapply(a_vec, function(a) int_score(lspm_preds, alpha = a))
  widt[['lspm']] <- sapply(a_vec, function(a) width(lspm_preds, alpha = a))

  ### CIDR
  cidr_preds <- conformal_idr(x = x_tr, y = y_tr, x_out = x_ts, y_out = y_ts)
  pcal[['cidr']] <- cidr_preds$pit
  crps[['cidr']] <- cidr_preds$crps
  thick[['cidr']] <- cidr_preds$thick
  F_t[['cidr']] <- threshcal(cidr_preds, t_vec)
  cover[['cidr']] <- sapply(a_vec, function(a) coverage(cidr_preds, alpha = a, average = F))
  intsc[['cidr']] <- sapply(a_vec, function(a) int_score(cidr_preds, alpha = a))
  widt[['cidr']] <- sapply(a_vec, function(a) width(cidr_preds, alpha = a))

  ### CB
  locb_preds <- conformal_bin(x = x_tr, y = y_tr, x_out = x_ts, y_out = y_ts, k = k)
  pcal[['locb']] <- sapply(locb_preds, function(x) x$pit)
  crps[['locb']] <- sapply(locb_preds, function(x) x$crps)
  thick[['locb']] <- sapply(locb_preds, function(x) x$thick)
  F_t[['locb']] <- sapply(locb_preds, function(x) threshcal(x, t_vec)) |> t()
  cover[['locb']] <- sapply(locb_preds, function(x) sapply(a_vec, function(a) coverage(x, alpha = a, average = F))) |> t()
  intsc[['locb']] <- sapply(locb_preds, function(x) sapply(a_vec, function(a) int_score(x, alpha = a))) |> t()
  widt[['locb']] <- sapply(locb_preds, function(x) sapply(a_vec, function(a) width(x, alpha = a))) |> t()

  ### CQR
  cqr_preds <- get_cqr_results(N_tr, x_ts, y_ts, type)
  cover$cqr <- cqr_preds$cover
  intsc$cqr <- cqr_preds$intsc
  widt$cqr <- cqr_preds$widt

  end <- Sys.time()
  print(end - start)

  return(list(pit = pcal, crps = crps, F_t = F_t, thick = thick, cover = cover, intsc = intsc, widt = widt))
}

# function to fit and evaluate conformalized quantile regression (CQR)
get_cqr_results <- function(N_tr = 500, x_ts, y_ts, type = c("iso", "anti", "less"), N_val = N_tr) {
    start <- Sys.time()
    type <- match.arg(type)

    x_tr <- runif(N_tr, 0, 10)
    x_val <- runif(N_val, 0, 10)
    if (type == "anti") {
      y_tr <- rgamma(N_tr, shape = sqrt(10 - x_tr), scale = pmin(pmax(10 - x_tr, 1), 6))
      y_val <- rgamma(N_val, shape = sqrt(10 - x_val), scale = pmin(pmax(10 - x_val, 1), 6))
    } else if (type == "less") {
      y_tr <- rnorm(N_tr, mean = 2*x_tr + 5*sin(x_tr), sd = 0.2*x_tr)
      y_val <- rnorm(N_val, mean = 2*x_val + 5*sin(x_val), sd = 0.2*x_val)
    } else {
      y_tr <- rgamma(N_tr, shape = sqrt(x_tr), scale = pmin(pmax(x_tr, 1), 6))
      y_val <- rgamma(N_val, shape = sqrt(x_val), scale = pmin(pmax(x_val, 1), 6))
    }

    N_ts <- length(y_ts)

    a_vec <- seq(0.05, 0.95, 0.05)
    n_a <- length(a_vec)
    cover <- intsc <- widt <- matrix(NA, N_ts, n_a)

    ## CQR
    for (j in seq_along(a_vec)) {
      a <- a_vec[j]
      print(a)

      # train on estimation data (val)
      val_df <- data.frame(x = x_val, y = y_val)
      fit_low <- rq(y ~ ns(x, df = 3), tau = a/2, data = val_df)
      fit_upp <- rq(y ~ ns(x, df = 3), tau = 1 - a/2, data = val_df)

      # get predictions on calibration data (train)
      cal_low <- predict(fit_low, newdata = data.frame(x = x_tr))
      cal_upp <- predict(fit_upp, newdata = data.frame(x = x_tr))

      cal_low <- pmax(cal_low, 0)
      cal_upp <- pmax(cal_upp, 0)

      # calculate scores
      scores <- pmax(cal_low - y_tr, y_tr - cal_upp)
      k_cqr <- ceiling((N_tr + 1) * (1 - a))
      q_hat <- sort(scores)[k_cqr]

      # get predictions on test data (test)
      pred_low <- predict(fit_low, newdata = data.frame(x = x_ts))
      pred_upp <- predict(fit_upp, newdata = data.frame(x = x_ts))

      pred_low <- pmax(pred_low, 0)
      pred_upp <- pmax(pred_upp, 0)

      # get prediction intervals on test data
      l <- pred_low - q_hat
      u <- pred_upp + q_hat
      l <- pmax(pmin(l, u), 0)
      u <- pmax(pmax(l, u), 0)

      # evaluate
      cover[, j] <- as.numeric(y_ts >= l & y_ts <= u)
      widt[, j] <- (u - l)
      intsc[, j] <- scoringRules::ints_quantiles(y_ts, l, u, 1 - a)
    }
    return(list(cover = cover, intsc = intsc, widt = widt))
}

# wrapper to plot and save PIT histograms, pp-plots, and threshold calibration plots
plot_cal <- function(pcal, score, F_t = NULL, obs = NULL, type = "pitpp", vert = T, filename = NULL) {
  ncol <- if (vert) {1} else {3}
  if (type == "pithist") {
    ## PIT histograms
    lspm_plot <- pit_hist(pcal[['lspm']], ranks = F, ymax = 0.4, xlab = NULL, xticks = F,
                          title = paste("LSPM: CRPS =", round(mean(score[['lspm']], na.rm = T), 3)))
    cidr_plot <- pit_hist(pcal[['cidr']], ranks = F, ymax = 0.4, xlab = NULL, xticks = F,
                          title = paste("CIDR: CRPS =", round(mean(score[['cidr']], na.rm = T), 3)))
    locb_plot <- pit_hist(pcal[['locb']], ranks = F, ymax = 0.4, xlab = NULL, xticks = F,
                          title = paste("CB: CRPS =", round(mean(score[['locb']], na.rm = T), 3)))
  } else if (type == "pitpp") {
    ## PIT pp-plots
    lspm_plot <- pit_reldiag(pcal[['lspm']], title = paste("LSPM: CRPS =", round(mean(score[['lspm']], na.rm = T), 3)))
    cidr_plot <- pit_reldiag(pcal[['cidr']], title = paste("CIDR: CRPS =", round(mean(score[['cidr']], na.rm = T), 3)))
    locb_plot <- pit_reldiag(pcal[['locb']], title = paste("CB: CRPS =", round(mean(score[['locb']], na.rm = T), 3)))
  } else if (type == "threshcal") {
    ## threshold calibration
    t_vec <- quantile(obs, c(0.1, 0.25, 0.5, 0.75, 0.9))
    lspm_plot <- tc_reldiag(F_t[['lspm']], obs, t_vec, xlab = "F(x)", ylab = "P(Y \u2264 x | F(x))", title = "LSPM")
    cidr_plot <- tc_reldiag(F_t[['cidr']], obs, t_vec, xlab = "F(x)", ylab = "P(Y \u2264 x | F(x))", title = "CIDR")
    locb_plot <- tc_reldiag(F_t[['locb']], obs, t_vec, xlab = "F(x)", ylab = "P(Y \u2264 x | F(x))", title = "CB")
  } else {
    stop("argument 'type' must be one of 'pithist', 'pitpp', and 'threshcal'")
  }

  cal_plot <- gridExtra::grid.arrange(lspm_plot, cidr_plot, locb_plot, ncol = ncol)

  if (!is.null(filename)) {
    if (type == "threshcal" || vert == F) {
      ggsave(plot = cal_plot, filename, width = 10.5, height = 3.5, dpi = 300)
    } else {
      ggsave(plot = cal_plot, filename, width = 2.5, height = 7.5, dpi = 300)
    }
  }
}

# wrapper to plot and save interval score plots
plot_is <- function(is, alpha = NULL, filename = NULL) {
  if (is.null(alpha)) alpha <- seq(0.05, 0.95, 0.05)
  scores <- sapply(is, colMeans)
  df <- data.frame(s = as.vector(scores), a = 1 - alpha, mth = rep(c("LSPM", "CIDR", "CB", "CQR"), each = length(alpha)))
  is_plot <- ggplot(df) + geom_point(aes(x = a, y = s, col = mth), size = 2) +
    geom_line(aes(x = a, y = s, col = mth), linewidth = 1) +
    scale_x_continuous(name = expression(paste("Nominal level (", 1 - alpha, ")")), limits = c(0, 1)) +
    scale_y_continuous(name = "Interval Score") +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.justification = c(0, 1),
          legend.position = c(0.01, 0.99))
  if (!is.null(filename)) {
    ggsave(plot = is_plot, filename, width = 5, height = 3, dpi = 300)
  }
  return(is_plot)
}

# wrapper to plot and save unconditional and conditional coverage plots
plot_cov <- function(cov, alpha = NULL, x_ts = NULL, n_bins = 10, filename = NULL) {
  if (is.null(alpha)) alpha <- seq(0.05, 0.95, 0.05)
  cov_mat <- sapply(cov, colMeans)

  if (is.null(x_ts)) {
    df <- data.frame(s = as.vector(cov_mat), a = 1 - alpha, mth = rep(c("LSPM", "CIDR", "CB", "CQR"), each = length(alpha)))
    cov_plot <- ggplot(df) +
      geom_abline(aes(intercept = 0, slope = 1), lty = "dotted") +
      geom_point(aes(x = a, y = s, col = mth), size = 2) +
      geom_line(aes(x = a, y = s, col = mth), linewidth = 1) +
      scale_x_continuous(name = expression(paste("Nominal level (", 1 - alpha, ")")), limits = c(0, 1)) +
      scale_y_continuous(name = "Empirical coverage", limits = c(0, 1)) +
      theme_bw() +
      theme(legend.title = element_blank(),
            legend.justification = c(0, 1),
            legend.position = c(0.01, 0.99))
  } else {
    breaks <- seq(0, 10, length.out = n_bins + 1)
    cond_cal_05 <- sapply(1:n_bins, function(i) {
      ind <- (x_ts >= breaks[i] & x_ts < breaks[i + 1])
      sapply(cov, function(x) mean(x[ind, which(abs(alpha - 0.5) < 1e-10)]))
    })
    cond_cal_09 <- sapply(1:n_bins, function(i) {
      ind <- (x_ts >= breaks[i] & x_ts < breaks[i + 1])
      sapply(cov, function(x) mean(x[ind, which(abs(alpha - 0.1) < 1e-10)]))
    })
    df <- data.frame(x = breaks[1:n_bins] + diff(breaks)/2,
                     c1 = as.vector(t(cond_cal_05)),
                     c2 = as.vector(t(cond_cal_09)),
                     mth = rep(c("LSPM", "CIDR", "CB", "CQR"), each = n_bins))
    cov_plot <- ggplot(df) +
      geom_hline(aes(yintercept = 0.5), lty = "dotted") +
      geom_hline(aes(yintercept = 0.9), lty = "dotted") +
      geom_point(aes(x = x, y = c1, col = mth), size = 2) +
      geom_line(aes(x = x, y = c1, col = mth), linewidth = 1, lty = "dashed") +
      geom_point(aes(x = x, y = c2, col = mth), size = 2) +
      geom_line(aes(x = x, y = c2, col = mth), linewidth = 1) +
      scale_x_continuous(name = "X", breaks = breaks, limits = c(0, 10)) +
      scale_y_continuous(name = "Empirical coverage", limits = c(0, 1)) +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            legend.title = element_blank(),
            legend.justification = c(0, 0),
            legend.position = c(0.01, 0.01)) +
      guides(colour = guide_legend(nrow = 1))
  }

  if (!is.null(filename)) {
    ggsave(plot = cov_plot, filename, width = 5, height = 3, dpi = 300)
  }
  return(cov_plot)
}

# wrapper to plot unconditional and conditional interval width plots
plot_width <- function(width, alpha = NULL, x_ts = NULL, n_bins = 10, type = c("ave", "cond", "samp"), filename = NULL) {
  type <- match.arg(type)
  if (is.null(alpha)) alpha <- seq(0.05, 0.95, 0.05)

  if (type == "ave") {
    av_width <- sapply(width, colMeans)
    df <- data.frame(s = as.vector(av_width), a = 1 - alpha, mth = rep(c("LSPM", "CIDR", "CB", "CQR"), each = length(alpha)))
    wid_plot <- ggplot(df) +
      geom_point(aes(x = a, y = s, col = mth), size = 2) +
      geom_line(aes(x = a, y = s, col = mth), linewidth = 1) +
      scale_x_continuous(name = expression(paste("Nominal level (", 1 - alpha, ")")), limits = c(0, 1)) +
      scale_y_continuous(name = "Average width") +
      theme_bw() +
      theme(legend.title = element_blank(),
            legend.justification = c(0, 1),
            legend.position = c(0.01, 0.99))
  } else if (type == "cond") {
    breaks <- seq(0, 10, length.out = n_bins + 1)
    cond_wid_05 <- sapply(1:n_bins, function(i) {
      ind <- (x_ts >= breaks[i] & x_ts < breaks[i + 1])
      sapply(width, function(x) mean(x[ind, which(abs(alpha - 0.5) < 1e-10)]))
    })
    cond_wid_09 <- sapply(1:n_bins, function(i) {
      ind <- (x_ts >= breaks[i] & x_ts < breaks[i + 1])
      sapply(width, function(x) mean(x[ind, which(abs(alpha - 0.1) < 1e-10)]))
    })
    df <- data.frame(x = breaks[1:n_bins] + diff(breaks)/2,
                     w1 = as.vector(t(cond_wid_05)),
                     w2 = as.vector(t(cond_wid_09)),
                     mth = rep(c("LSPM", "CIDR", "CB", "CQR"), each = n_bins))
    wid_plot <- ggplot(df) +
      geom_point(aes(x = x, y = w1, col = mth), size = 2) +
      geom_line(aes(x = x, y = w1, col = mth), linewidth = 1, lty = "dashed") +
      geom_point(aes(x = x, y = w2, col = mth), size = 2) +
      geom_line(aes(x = x, y = w2, col = mth), linewidth = 1) +
      scale_x_continuous(name = "X", breaks = breaks, limits = c(0, 10)) +
      scale_y_continuous(name = "Average width") +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            legend.title = element_blank(),
            legend.justification = c(0, 1),
            legend.position = c(0.01, 0.99))
  } else {
    a_ind <- which(abs(alpha - 0.1) < 1e-10)
    df <- lapply(seq_along(width), function(i) {
      id <- c(" 100", " 500", "1000", "2000")[i]
      mat_list <- width[[i]]
      inner_df <- lapply(seq_along(mat_list), function(j) {
        data.frame(value = mat_list[[j]][, a_ind],
                   id = id,
                   mth = c("LSPM", "CIDR", "CB", "CQR")[j])
      })
      do.call(rbind, inner_df)
    }) |> do.call(what = rbind)

    wid_plot <- ggplot(df, aes(x = id, y = value, fill = mth)) +
      geom_boxplot(width = 0.7) +
      scale_x_discrete(name = "Sample size") +
      scale_y_continuous(name = "Width") +
      theme_bw() +
      theme(legend.position = "bottom",
            legend.title = element_blank(),
            panel.grid.minor = element_blank())
  }

  if (!is.null(filename)) {
    ggsave(plot = wid_plot, filename, width = 5, height = 3, dpi = 300)
  }
  return(wid_plot)
}

# function to plot example data
plot_example <- function(x, y, filename = NULL, ylims = c(0, 80)) {
  df <- data.frame(x = x, y = y)
  plot_obj <- ggplot(df) + geom_point(aes(x = x, y = y), size = 0.5) +
    scale_x_continuous(name = "X", limits = c(0, 10)) +
    scale_y_continuous(name = "Y", limits = ylims) +
    theme_bw() +
    theme(panel.grid = element_blank())
  if (!is.null(filename)) {
    ggsave(filename, plot_obj, width = 3, height = 3)
  }
  return(plot_obj)
}

# wrapper to train and predict from conformalised quantile regression model
cqr_pred <- function(x, y, alpha = NULL) {
  if (is.null(alpha)) alpha <- seq(0.05, 0.95, 0.05)

  for (j in seq_along(a_vec)) {
    a <- a_vec[j]
    print(a)

    # train on estimation data (val)
    fit_low <- rq(y ~ ns(x, df = 3), tau = a/2, data = val)
    fit_upp <- rq(y ~ ns(x, df = 3), tau = 1 - a/2, data = val)

    # get predictions on calibration data (train)
    #cal_low <- predict(qrf_mod, newdata = data.frame(index = train$index), what = a/2)
    #cal_upp <- predict(qrf_mod, newdata = data.frame(index = train$index), what = 1 - a/2)
    cal_low <- predict(fit_low, newdata = data.frame(index = train$index))
    cal_upp <- predict(fit_upp, newdata = data.frame(index = train$index))

    cal_low <- pmax(cal_low, 0)
    cal_upp <- pmax(cal_upp, 0)

    # calculate scores
    scores <- pmax(cal_low - train$los, train$los - cal_upp)
    k_cqr <- ceiling((n_cal + 1) * (1 - a))
    q_hat <- sort(scores)[k_cqr]

    # get predictions on test data (test)
    #pred_low <- predict(qrf_mod, newdata = as.matrix(x = test$index), what = a/2)
    #pred_upp <- predict(qrf_mod, newdata = as.matrix(x = test$index), what = 1 - a/2)
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

