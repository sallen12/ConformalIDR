################################################################################
## simulation study

set.seed(90743892)

library(ConformalIDR)
library(WeightedForecastVerification)
library(ggplot2)


# wrapper to fit and evaluate the different prediction methods
get_results <- function(N_tr = 500, x_ts, y_ts, anti = FALSE, k = 10) {
  start <- Sys.time()

  x_tr <- runif(N_tr, 0, 10)
  y_tr <- rgamma(N_tr, shape = sqrt(x_tr), scale = pmin(pmax(x_tr, 1), 6))
  if (anti) x_tr <- 10 - x_tr

  N_ts <- length(y_ts)
  pcal <- data.frame(replicate(3, numeric(N_ts)))
  colnames(pcal) <- c("lspm", "cidr", "locb")
  thick <- score <- pcal

  t_vec <- unname(quantile(y_ts, c(0.1, 0.25, 0.5, 0.75, 0.9)))
  n_t <- length(t_vec)
  F_t <- list(lspm = matrix(NA, N_ts, n_t),
              cidr = matrix(NA, N_ts, n_t),
              locb = matrix(NA, N_ts, n_t))

  ### LSPM
  lspm_preds <- conformal_lspm(x = x_tr, y = y_tr, x_out = x_ts, y_out = y_ts)
  pcal[['lspm']] <- lspm_preds$pit
  score[['lspm']] <- lspm_preds$crps
  thick[['lspm']] <- lspm_preds$thick
  F_t[['lspm']] <- threshcal(lspm_preds, t_vec)

  ### CIDR
  cidr_preds <- conformal_idr(x = x_tr, y = y_tr, x_out = x_ts, y_out = y_ts)
  pcal[['cidr']] <- cidr_preds$pit
  score[['cidr']] <- cidr_preds$crps
  thick[['cidr']] <- cidr_preds$thick
  F_t[['cidr']] <- threshcal(cidr_preds, t_vec)

  ### LB
  locb_preds <- conformal_bin(x = x_tr, y = y_tr, x_out = x_ts, y_out = y_ts, k = k)
  pcal[['locb']] <- sapply(locb_preds, function(x) x$pit)
  score[['locb']] <- sapply(locb_preds, function(x) x$crps)
  thick[['locb']] <- sapply(locb_preds, function(x) x$thick)
  F_t[['locb']] <- sapply(locb_preds, function(x) threshcal(x, t_vec)) |> t()

  end <- Sys.time()
  print(end - start)

  return(list(pit = pcal, crps = score, F_t = F_t, thick = thick))
}

# wrapper to plot and save PIT histograms, pp-plots, and threshold calibration plots
plot_cal <- function(pcal, score, F_t = NULL, obs = NULL, type = "pitpp", filename = NULL) {
  if (type == "pithist") {
    ## PIT histograms
    lspm_plot <- pit_hist(pcal[['lspm']], ranks = F, ymax = 0.4, xlab = NULL, xticks = F,
                          title = paste("LSPM: CRPS =", round(mean(score[['lspm']], na.rm = T), 3)))
    cidr_plot <- pit_hist(pcal[['cidr']], ranks = F, ymax = 0.4, xlab = NULL, xticks = F,
                          title = paste("CIDR: CRPS =", round(mean(score[['cidr']], na.rm = T), 3)))
    locb_plot <- pit_hist(pcal[['locb']], ranks = F, ymax = 0.4, xlab = NULL, xticks = F,
                          title = paste("LB: CRPS =", round(mean(score[['locb']], na.rm = T), 3)))
    cal_plot <- gridExtra::grid.arrange(lspm_plot, cidr_plot, locb_plot, ncol = 1)
  } else if (type == "pitpp") {
    ## PIT pp-plots
    lspm_plot <- pit_reldiag(pcal[['lspm']], title = paste("LSPM: CRPS =", round(mean(score[['lspm']], na.rm = T), 3)))
    cidr_plot <- pit_reldiag(pcal[['cidr']], title = paste("CIDR: CRPS =", round(mean(score[['cidr']], na.rm = T), 3)))
    locb_plot <- pit_reldiag(pcal[['locb']], title = paste("LB: CRPS =", round(mean(score[['locb']], na.rm = T), 3)))
    cal_plot <- gridExtra::grid.arrange(lspm_plot, cidr_plot, locb_plot, ncol = 1)
  } else if (type == "threshcal") {
    ## threshold calibration
    t_vec <- quantile(obs, c(0.1, 0.25, 0.5, 0.75, 0.9))
    lspm_plot <- tc_reldiag(F_t[['lspm']], obs, t_vec, title = "LSPM")
    cidr_plot <- tc_reldiag(F_t[['cidr']], obs, t_vec, title = "CIDR")
    locb_plot <- tc_reldiag(F_t[['locb']], obs, t_vec, title = "LB")
    cal_plot <- gridExtra::grid.arrange(lspm_plot, cidr_plot, locb_plot, nrow = 1)
  } else {
    stop("argument 'type' must be one of 'pithist', 'pitpp', and 'threshcal'")
  }

  if (!is.null(filename)) {
    if (type == "threshcal") {
      ggsave(plot = cal_plot, filename, width = 10.5, height = 3.5, dpi = 300)
    } else {
      ggsave(plot = cal_plot, filename, width = 2.5, height = 7.5, dpi = 300)
    }

  }
}

# function to plot example data
plot_example <- function(x, y, filename = NULL) {
  df <- data.frame(x = x, y = y)
  plot_obj <- ggplot(df) + geom_point(aes(x = x, y = y), size = 0.5) +
    scale_x_continuous(name = "X", limits = c(0, 10)) +
    scale_y_continuous(name = "Y", limits = c(0, 80)) +
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

# plot example data
plot_example(x_ts[1:1000], y_ts[1:1000], filename = "plots/simstudy_data.png")
plot_example(x_ts[1:1000], y_ts_at[1:1000])


################################################################################
## isotonic

##### results

res_100 <- get_results(N_tr = 100, x_ts, y_ts) # 4 secs
res_500 <- get_results(N_tr = 500, x_ts, y_ts) # 20 secs
res_1000 <- get_results(N_tr = 1000, x_ts, y_ts) # 1 min
res_2000 <- get_results(N_tr = 2000, x_ts, y_ts) # 5 mins


##### evaluate

## probabilistic calibration
plot_cal(res_100$pit, res_100$crps, filename = "plots/simstudy_100.png")
plot_cal(res_500$pit, res_500$crps, filename = "plots/simstudy_500.png")
plot_cal(res_1000$pit, res_1000$crps, filename = "plots/simstudy_1000.png")
plot_cal(res_2000$pit, res_2000$crps, filename = "plots/simstudy_2000.png")

## threshold calibration
plot_cal(F_t = res_2000$F_t, obs = y_ts, type = "threshcal", filename = "plots/simstudy_2000_tc.png")


################################################################################
## antitonic

##### results

res_100_at <- get_results(N_tr = 100, x_ts, y_ts_at, anti = T) # 4 secs
res_500_at <- get_results(N_tr = 500, x_ts, y_ts_at, anti = T) # 20 secs
res_1000_at <- get_results(N_tr = 1000, x_ts, y_ts_at, anti = T) # 1 min
res_2000_at <- get_results(N_tr = 2000, x_ts, y_ts_at, anti = T) # 5 mins


##### evaluate

## probabilistic calibration
plot_cal(res_100_at$pit, res_100_at$crps, filename = "plots/simstudy_100_at.png")
plot_cal(res_500_at$pit, res_500_at$crps, filename = "plots/simstudy_500_at.png")
plot_cal(res_1000_at$pit, res_1000_at$crps, filename = "plots/simstudy_1000_at.png")
plot_cal(res_2000_at$pit, res_2000_at$crps, filename = "plots/simstudy_2000_at.png")

## threshold calibration
plot_cal(F_t = res_2000_at$F_t, obs = y_ts_at, type = "threshcal", filename = "plots/simstudy_2000_tc_at.png")


##### compare thicknesses
df <- data.frame(thicc = c(res_100$thick$cidr, res_100_at$thick$cidr),
                 mode = rep(c("Isotonic", "Antitonic"), each = N_ts))
ggplot(df) + geom_boxplot(aes(x = mode, y = thicc)) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Thickness") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave("plots/simstudy_thick.png", width = 3, height = 3)
