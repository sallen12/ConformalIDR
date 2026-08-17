################################################################################
## simulation study

set.seed(90743892)

library(ConformalIDR)
library(WeightedForecastVerification)
library(ggplot2)


# wrapper to fit and evaluate the different prediction methods
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
  cov <- is <- width <- list(lspm = matrix(NA, N_ts, n_a),
                        cidr = matrix(NA, N_ts, n_a),
                        locb = matrix(NA, N_ts, n_a))

  ### LSPM
  lspm_preds <- conformal_lspm(x = x_tr, y = y_tr, x_out = x_ts, y_out = y_ts)
  pcal[['lspm']] <- lspm_preds$pit
  crps[['lspm']] <- lspm_preds$crps
  thick[['lspm']] <- lspm_preds$thick
  F_t[['lspm']] <- threshcal(lspm_preds, t_vec)
  cov[['lspm']] <- sapply(a_vec, function(a) coverage(lspm_preds, alpha = a, average = F))
  is[['lspm']] <- sapply(a_vec, function(a) int_score(lspm_preds, alpha = a))
  width[['lspm']] <- sapply(a_vec, function(a) width(lspm_preds, alpha = a))

  ### CIDR
  cidr_preds <- conformal_idr(x = x_tr, y = y_tr, x_out = x_ts, y_out = y_ts)
  pcal[['cidr']] <- cidr_preds$pit
  crps[['cidr']] <- cidr_preds$crps
  thick[['cidr']] <- cidr_preds$thick
  F_t[['cidr']] <- threshcal(cidr_preds, t_vec)
  cov[['cidr']] <- sapply(a_vec, function(a) coverage(cidr_preds, alpha = a, average = F))
  is[['cidr']] <- sapply(a_vec, function(a) int_score(cidr_preds, alpha = a))
  width[['cidr']] <- sapply(a_vec, function(a) width(cidr_preds, alpha = a))

  ### CB
  locb_preds <- conformal_bin(x = x_tr, y = y_tr, x_out = x_ts, y_out = y_ts, k = k)
  pcal[['locb']] <- sapply(locb_preds, function(x) x$pit)
  crps[['locb']] <- sapply(locb_preds, function(x) x$crps)
  thick[['locb']] <- sapply(locb_preds, function(x) x$thick)
  F_t[['locb']] <- sapply(locb_preds, function(x) threshcal(x, t_vec)) |> t()
  cov[['locb']] <- sapply(locb_preds, function(x) sapply(a_vec, function(a) coverage(x, alpha = a, average = F))) |> t()
  is[['locb']] <- sapply(locb_preds, function(x) sapply(a_vec, function(a) int_score(x, alpha = a))) |> t()
  width[['locb']] <- sapply(locb_preds, function(x) sapply(a_vec, function(a) width(x, alpha = a))) |> t()

  end <- Sys.time()
  print(end - start)

  return(list(pit = pcal, crps = crps, F_t = F_t, thick = thick, cov = cov, is = is, width = width))
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
  df <- data.frame(s = as.vector(scores), a = 1 - alpha, mth = rep(c("LSPM", "CIDR", "CB"), each = length(alpha)))
  is_plot <- ggplot(df) + geom_point(aes(x = a, y = s, col = mth), size = 2) +
    geom_line(aes(x = a, y = s, col = mth), linewidth = 1) +
    scale_x_continuous(name = expression(1 - alpha), limits = c(0, 1)) +
    scale_y_continuous(name = "Interval Score") +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.justification = c(0, 1),
          legend.position = c(0.01, 0.99))
  if (!is.null(filename)) {
    ggsave(plot = is_plot, filename, width = 5, height = 3, dpi = 300)
  }
}

# wrapper to plot and save unconditional and conditional coverage plots
plot_cov <- function(cov, alpha = NULL, x_ts = NULL, n_bins = 10, filename = NULL) {
  if (is.null(alpha)) alpha <- seq(0.05, 0.95, 0.05)
  cov_mat <- sapply(cov, colMeans)

  if (is.null(x_ts)) {
    df <- data.frame(s = as.vector(cov_mat), a = 1 - alpha, mth = rep(c("LSPM", "CIDR", "CB"), each = length(alpha)))
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
                     mth = rep(c("LSPM", "CIDR", "CB"), each = n_bins))
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

}

# wrapper to plot unconditional and conditional interval width plots
plot_width <- function(width, alpha = NULL, x_ts = NULL, n_bins = 10, type = c("ave", "cond", "samp"), filename = NULL) {
  type <- match.arg(type)
  if (is.null(alpha)) alpha <- seq(0.05, 0.95, 0.05)

  if (type == "ave") {
    av_width <- sapply(width, colMeans)
    df <- data.frame(s = as.vector(av_width), a = 1 - alpha, mth = rep(c("LSPM", "CIDR", "CB"), each = length(alpha)))
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
                     mth = rep(c("LSPM", "CIDR", "CB"), each = n_bins))
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
                         mth = c("LSPM", "CIDR", "CB")[j])
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
    ggsave(plot = wid_plot, filename, width = 5, height = 3.3, dpi = 300)
  }

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

res_100 <- get_results(N_tr = 100, x_ts, y_ts) # 10 secs
res_500 <- get_results(N_tr = 500, x_ts, y_ts) # 40 secs
res_1000 <- get_results(N_tr = 1000, x_ts, y_ts) # 2 mins
res_2000 <- get_results(N_tr = 2000, x_ts, y_ts) # 10 mins


##### evaluate

## probabilistic calibration
plot_cal(res_100$pit, res_100$crps, filename = "plots/simstudy_100.png")
plot_cal(res_500$pit, res_500$crps, filename = "plots/simstudy_500.png")
plot_cal(res_1000$pit, res_1000$crps, filename = "plots/simstudy_1000.png")
plot_cal(res_2000$pit, res_2000$crps, filename = "plots/simstudy_2000.png")

plot_cal(res_2000$pit, res_2000$crps, vert = F, filename = "plots/simstudy_2000_ho.png")

## threshold calibration
plot_cal(F_t = res_2000$F_t, obs = y_ts, type = "threshcal", filename = "plots/simstudy_2000_tc.png")

## unconditional coverage of prediction intervals
plot_cov(res_2000$cov, filename = "plots/simstudy_2000_cov.png")

## conditional coverage of prediction intervals
plot_cov(res_2000$cov, x_ts = x_ts, filename = "plots/simstudy_2000_ccov.png")

## interval score of prediction intervals
plot_is(res_2000$is, filename = "plots/simstudy_2000_is.png")

## unconditional average width of prediction intervals
plot_width(res_2000$width, filename = "plots/simstudy_2000_wid.png")

## conditional average width of prediction intervals
plot_width(res_2000$width, x_ts = x_ts, filename = "plots/simstudy_2000_cwid.png", type = "cond")

## average width of prediction intervals
plot_width(lapply(list(res_100, res_500, res_1000, res_2000), function(x) x$width),
           filename = "plots/simstudy_2000_wid_samp.png", type = "samp")


################################################################################
## antitonic

##### results

res_100_at <- get_results(N_tr = 100, x_ts, y_ts_at, type = "anti")
res_500_at <- get_results(N_tr = 500, x_ts, y_ts_at, type = "anti")
res_1000_at <- get_results(N_tr = 1000, x_ts, y_ts_at, type = "anti")
res_2000_at <- get_results(N_tr = 2000, x_ts, y_ts_at, type = "anti")


##### evaluate

## probabilistic calibration
plot_cal(res_100_at$pit, res_100_at$crps, filename = "plots/simstudy_100_at.png")
plot_cal(res_500_at$pit, res_500_at$crps, filename = "plots/simstudy_500_at.png")
plot_cal(res_1000_at$pit, res_1000_at$crps, filename = "plots/simstudy_1000_at.png")
plot_cal(res_2000_at$pit, res_2000_at$crps, filename = "plots/simstudy_2000_at.png")

## threshold calibration
plot_cal(F_t = res_1000_at$F_t, obs = y_ts_at, type = "threshcal", filename = "plots/simstudy_1000_tc_at.png")

## unconditional coverage of prediction intervals
plot_cov(res_2000_at$cov, filename = "plots/simstudy_2000_cov_at.png")

## conditional coverage of prediction intervals
plot_cov(res_2000_at$cov, x_ts = x_ts, filename = "plots/simstudy_2000_ccov_at.png")

## interval score of prediction intervals
plot_is(res_2000_at$is, filename = "plots/simstudy_2000_is_at.png")

## unconditional average width of prediction intervals
plot_width(res_2000_at$width, filename = "plots/simstudy_2000_wid_at.png")

## conditional average width of prediction intervals
plot_width(res_2000_at$width, x_ts = x_ts, filename = "plots/simstudy_2000_cwid_at.png", type = "cond")

## average width of prediction intervals
plot_width(lapply(list(res_100_at, res_500_at, res_1000_at, res_2000_at), function(x) x$width),
           filename = "plots/simstudy_2000_wid_samp_at.png", type = "samp")


################################################################################
## less isotonic

##### results

res_100_le <- get_results(N_tr = 100, x_ts, y_ts_le, type = "less")
res_500_le <- get_results(N_tr = 500, x_ts, y_ts_le, type = "less")
res_1000_le <- get_results(N_tr = 1000, x_ts, y_ts_le, type = "less")
res_2000_le <- get_results(N_tr = 2000, x_ts, y_ts_le, type = "less")


##### evaluate

## probabilistic calibration
plot_cal(res_100_le$pit, res_100_le$crps, filename = "plots/simstudy_100_le.png")
plot_cal(res_500_le$pit, res_500_le$crps, filename = "plots/simstudy_500_le.png")
plot_cal(res_1000_le$pit, res_1000_le$crps, filename = "plots/simstudy_1000_le.png")
plot_cal(res_2000_le$pit, res_2000_le$crps, filename = "plots/simstudy_2000_le.png")

plot_cal(res_2000_le$pit, res_2000_le$crps, vert = F, filename = "plots/simstudy_2000_le_ho.png")

## threshold calibration
plot_cal(F_t = res_2000_le$F_t, obs = y_ts_le, type = "threshcal", filename = "plots/simstudy_2000_tc_le.png")


## unconditional coverage of prediction intervals
plot_cov(res_2000_le$cov, filename = "plots/simstudy_2000_cov_le.png")

## conditional coverage of prediction intervals
plot_cov(res_2000_le$cov, x_ts = x_ts, filename = "plots/simstudy_2000_ccov_le.png")

## interval score of prediction intervals
plot_is(res_2000_le$is, filename = "plots/simstudy_2000_is_le.png")

## unconditional average width of prediction intervals
plot_width(res_2000_le$width, filename = "plots/simstudy_2000_wid_le.png")

## conditional average width of prediction intervals
plot_width(res_2000_le$width, x_ts = x_ts, filename = "plots/simstudy_2000_cwid_le.png", type = "cond")

## average width of prediction intervals
plot_width(lapply(list(res_100_le, res_500_le, res_1000_le, res_2000_le), function(x) x$width),
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


