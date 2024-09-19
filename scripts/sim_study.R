################################################################################
## simulation study

set.seed(90743892)

library(scoringRules)
library(WeightedForecastVerification)
library(ggplot2)

source("scripts/utility_funcs.R")

# wrapper to fit and evaluate the different prediction methods
get_results <- function(N_tr = 500, x_ts, y_ts, anti = FALSE, k = 10) {
  start <- Sys.time()

  x_tr <- runif(N_tr, 0, 10)
  y_tr <- rgamma(N_tr, shape = sqrt(x_tr), scale = pmin(pmax(x_tr, 1), 6))
  if (anti) x_tr <- 10 - x_tr

  N_ts <- length(y_ts)
  pit <- data.frame(replicate(3, numeric(N_ts)))
  colnames(pit) <- c("lspm", "cidr", "mond")
  thicc <- crps <- pit

  t_vec <- unname(quantile(y_ts, c(0.1, 0.25, 0.5, 0.75, 0.9)))
  n_t <- length(t_vec)
  F_t <- list(lspm = matrix(NA, N_ts, n_t),
              cidr = matrix(NA, N_ts, n_t),
              mond = matrix(NA, N_ts, n_t))

  ### LSPM
  lspm_preds <- fit_lspm(y = y_tr, X = x_tr, X_ts = x_ts)
  scores <- eval_lspm(lspm_preds, y_ts, t_vec)
  pit[['lspm']] <- scores$pit
  crps[['lspm']] <- scores$crps
  F_t[['lspm']] <- scores$F_t
  thicc[['lspm']] <- scores$thick

  ### CIDR
  cidr_preds <- fit_cidr(y = y_tr, X = x_tr, X_ts = x_ts)
  scores <- eval_cidr(cidr_preds, y_ts, t_vec)
  pit[['cidr']] <- scores$pit
  crps[['cidr']] <- scores$crps
  F_t[['cidr']] <- scores$F_t
  thicc[['cidr']] <- scores$thick

  ### Mondrian
  mond_preds <- fit_mond(y = y_tr, X = x_tr, X_ts = x_ts, k)
  scores <- eval_mond(mond_preds, y_ts, t_vec)
  pit[['mond']] <- scores$pit
  crps[['mond']] <- scores$crps
  F_t[['mond']] <- scores$F_t
  thicc[['mond']] <- scores$thick

  end <- Sys.time()
  print(end - start)

  return(list(pit = pit, crps = crps, F_t = F_t, thick = thicc))
}

# wrapper to plot and save PIT histograms, pp-plots, and threshold calibration plots
plot_cal <- function(pit, crps, F_t = NULL, obs = NULL, type = "pitpp", filename = NULL) {
  if (type == "pithist") {
    ## PIT histograms
    id_plot <- pit_hist(pit[['id']], ranks = F, ymax = 0.4, xlab = NULL, xticks = F,
                        title = paste("Ideal: CRPS =", round(mean(crps[['id']], na.rm = T), 3)))
    lspm_plot <- pit_hist(pit[['lspm']], ranks = F, ymax = 0.4, xlab = NULL, xticks = F,
                          title = paste("LSPM: CRPS =", round(mean(crps[['lspm']], na.rm = T), 3)))
    cidr_plot <- pit_hist(pit[['cidr']], ranks = F, ymax = 0.4, xlab = NULL, xticks = F,
                          title = paste("CIDR: CRPS =", round(mean(crps[['cidr']], na.rm = T), 3)))
    mond_plot <- pit_hist(pit[['mond']], ranks = F, ymax = 0.4, xlab = NULL, xticks = F,
                          title = paste("LB: CRPS =", round(mean(crps[['mond']], na.rm = T), 3)))
    cal_plot <- gridExtra::grid.arrange(lspm_plot, cidr_plot, mond_plot, ncol = 1)
  } else if (type == "pitpp") {
    ## PIT pp-plots
    id_plot <- pit_reldiag(pit[['id']], title = paste("Ideal: CRPS =", round(mean(crps[['id']], na.rm = T), 3)))
    lspm_plot <- pit_reldiag(pit[['lspm']], title = paste("LSPM: CRPS =", round(mean(crps[['lspm']], na.rm = T), 3)))
    cidr_plot <- pit_reldiag(pit[['cidr']], title = paste("CIDR: CRPS =", round(mean(crps[['cidr']], na.rm = T), 3)))
    mond_plot <- pit_reldiag(pit[['mond']], title = paste("LB: CRPS =", round(mean(crps[['mond']], na.rm = T), 3)))
    cal_plot <- gridExtra::grid.arrange(lspm_plot, cidr_plot, mond_plot, ncol = 1)
  } else if (type == "threshcal") {
    ## threshold calibration
    t_vec <- quantile(obs, c(0.1, 0.25, 0.5, 0.75, 0.9))
    lspm_plot <- threshreldiag(F_t[['lspm']], obs, t_vec, title = "LSPM")
    cidr_plot <- threshreldiag(F_t[['cidr']], obs, t_vec, title = "CIDR")
    mond_plot <- threshreldiag(F_t[['mond']], obs, t_vec, title = "LB")
    cal_plot <- gridExtra::grid.arrange(lspm_plot, cidr_plot, mond_plot, nrow = 1)
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


################################################################################
## set up

N_ts <- 5000
x_ts <- runif(N_ts, 0, 10)

# isotonic example
y_ts <- rgamma(N_ts, shape = sqrt(x_ts), scale = pmin(pmax(x_ts, 1), 6))

# antitonic example
y_ts_at <- rgamma(N_ts, shape = sqrt(10 - x_ts), scale = pmin(pmax(10 - x_ts, 1), 6))

df <- data.frame(x = x_ts[1:1000], y = y_ts[1:1000])
ggplot(df) + geom_point(aes(x = x, y = y), size = 0.5) +
  scale_x_continuous(name = "X", limits = c(0, 10)) +
  scale_y_continuous(name = "Y", limits = c(0, 80)) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave("plots/simstudy_data.png", width = 3, height = 3)


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
