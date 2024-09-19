################################################################################
## set up

set.seed(783)

library(ncdf4)
library(ggplot2)
library(zoo)
library(scoringRules)
library(WeightedForecastVerification)
library(isodistrreg)
library(Rcpp)
library(mgcv)

source("online_idr.R")
sourceCpp("online_idr_computation_weights.cpp")

source("utility_funcs.R")

get_index <- function(val, train, test) {
  out <- tryCatch({
    fit <- gam(log(los + 1) ~ s(age, k = 3) + sex + planned + readmission + from + diag1 +
                 s(nems1, k = 3) + s(severity, k = 3) + interv,
               data = val)
    train$index <- exp(as.vector(predict(fit, train)) - 1)
    test$index <- exp(as.vector(predict(fit, test)) - 1)
    out <- list(train = train, test = test)
  },
  error = function(cond) {
    fit <- gam(log(los + 1) ~ s(age, k = 3) + sex + planned + readmission +
                 s(nems1, k = 3) + s(severity, k = 3),
               data = val)
    train$index <- exp(as.vector(predict(fit, train)) - 1)
    test$index <- exp(as.vector(predict(fit, test)) - 1)
    return(list(train = train, test = test))
  })
  return(out)
}


################################################################################
## load data

load("C:/Users/sa20i493/Documents/Data/icu_data/mdsi.rda")

icu_vec <- unique(data$icuCode)
n_icu <- length(icu_vec)


## filter for efficiency
data <- data[seq(1, nrow(data), 2), ]


## split into train, validation, and test data
data_ts <- tail(data, floor(nrow(data)*0.2))
data_tr <- head(data, ceiling(nrow(data)*0.8))
val_ind <- sample(c(rep(TRUE, floor(nrow(data_tr)*0.5)),
                    rep(FALSE, ceiling(nrow(data_tr)*0.5))))
data_val <- data_tr[val_ind, ]
data_tr <- data_tr[!val_ind, ]

rm(data)


################################################################################
## forecasting

t_vec <- unname(quantile(data_ts$los, c(0.1, 0.25, 0.5, 0.75, 0.9)))
n_t <- length(t_vec)
N_ts <- nrow(data_ts)

pit <- crps <- thick <- data.frame(replicate(3, numeric(N_ts)))
colnames(pit) <- colnames(crps) <- colnames(thick) <- c("lspm", "cidr", "mond")

F_t <- list(lspm = matrix(NA, N_ts, n_t),
            cidr = matrix(NA, N_ts, n_t),
            mond = matrix(NA, N_ts, n_t))


## optimal k for Mondrian prediction
k_vec <- c(1, seq(10, 100, 10))
mondcrps <- matrix(NA, nrow = n_icu, ncol = length(k_vec))
for (j in seq_along(icu_vec)) {
  print(icu_vec[j])

  ### Get train data
  train <- subset(data_tr, icuCode == icu_vec[j])
  val <- subset(data_val, icuCode == icu_vec[j])

  ### Get index
  out <- get_index(val, train, val)
  train <- out$train
  val <- out$test

  for (k in seq_along(k_vec)) {
    ### Mondrian
    mond_preds <- fit_mond(train$los, train$index, val$index, k)
    scores <- eval_mond(mond_preds, val$los, t_vec)
    mondcrps[j, k] <- mean(scores$crps)
  }
}
plot(k_vec, colMeans(mondcrps))
#k <- k_vec[apply(mondcrps, 1, which.min)]
k <- 20


################################################################################
## fit models

time_meth <- matrix(NA, nrow = n_icu, ncol = 3)
colnames(time_meth) <- c("lspm", "cidr", "mond")
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
  lspm_preds <- fit_lspm(train$los, as.matrix(train$index), as.matrix(test$index))
  scores <- eval_lspm(lspm_preds, test$los, t_vec)
  pit[['lspm']][icu_ind] <- scores$pit
  crps[['lspm']][icu_ind] <- scores$crps
  thick[['lspm']][icu_ind] <- scores$thick
  F_t[['lspm']][icu_ind, ] <- scores$F_t
  end <- Sys.time()
  time_meth[which(icu_vec == icu), 1] <- end - start

  ### CIDR
  start <- Sys.time()
  cidr_preds <- fit_cidr(train$los, train$index, test$index)
  scores <- eval_cidr(cidr_preds, test$los, t_vec)
  pit[['cidr']][icu_ind] <- scores$pit
  crps[['cidr']][icu_ind] <- scores$crps
  thick[['cidr']][icu_ind] <- scores$thick
  F_t[['cidr']][icu_ind, ] <- scores$F_t
  end <- Sys.time()
  time_meth[which(icu_vec == icu), 2] <- end - start

  ### Mondrian
  start <- Sys.time()
  mond_preds <- fit_mond(train$los, train$index, test$index, k[i])
  scores <- eval_mond(mond_preds, test$los, t_vec)
  pit[['mond']][icu_ind] <- scores$pit
  crps[['mond']][icu_ind] <- scores$crps
  thick[['mond']][icu_ind] <- scores$thick
  F_t[['mond']][icu_ind, ] <- scores$F_t
  end <- Sys.time()
  time_meth[which(icu_vec == icu), 3] <- end - start

  print(time_meth[which(icu_vec == icu), ])
}


################################################################################
## results

## PIT histograms
lspm_plot <- pit_hist(pit[['lspm']], ranks = F, ymax = 0.4, xlab = NULL, xticks = F,
                      title = paste("LSPM: CRPS =", round(mean(crps[['lspm']], na.rm = T), 3)))
cidr_plot <- pit_hist(pit[['cidr']], ranks = F, ymax = 0.4, xlab = NULL, xticks = F,
                      title = paste("CIDR: CRPS =", round(mean(crps[['cidr']], na.rm = T), 3)))
mond_plot <- pit_hist(pit[['mond']], ranks = F, ymax = 0.4, xlab = NULL, xticks = F,
                      title = paste("LB: CRPS =", round(mean(crps[['mond']], na.rm = T), 3)))
pit_plot <- gridExtra::grid.arrange(lspm_plot, cidr_plot, mond_plot, nrow = 1)
#ggsave(plot = pit_plot, "ICU_pit_comp.png", width = 7.5, height = 2.5)


## PIT histograms at ICUs
code_vec <- paste0("ICU", c(44, 65, 76, 77))
plot_list <- vector("list", 4)
for (i in 1:4) {
  ind <- data_ts$icuCode == code_vec[i]
  plot_list[[i]] <- pit_hist(pit[['cidr']][ind], ranks = F, bins = 20, ymax = 0.08, xlab = NULL, xticks = F,
                             title = paste0(code_vec[i], ": CRPS = ", round(mean(crps[['cidr']][ind], na.rm = T), 3)))
}
plot_list[['nrow']] <- 1
pit_plot <- do.call(gridExtra::grid.arrange, plot_list)
#ggsave(plot = pit_plot, "ICU_pit_comp_ind.png", width = 10, height = 2.5)


## PIT reliability diagrams
lspm_plot <- pit_reldiag(pit[['lspm']],
                      title = paste("LSPM: CRPS =", round(mean(crps[['lspm']], na.rm = T), 3)))
cidr_plot <- pit_reldiag(pit[['cidr']],
                         title = paste("CIDR: CRPS =", round(mean(crps[['cidr']], na.rm = T), 3)))
mond_plot <- pit_reldiag(pit[['lspm']],
                         title = paste("LB: CRPS =", round(mean(crps[['mond']], na.rm = T), 3)))
pit_plot <- gridExtra::grid.arrange(lspm_plot, cidr_plot, mond_plot, nrow = 1)
#ggsave(plot = pit_plot, "ICU_pitrd_comp.png", width = 7.5, height = 2.5)


## Threshold calibration
lspm_plot <- threshreldiag(F_t[['lspm']], data_ts$los, t = t_vec, title = "LSPM")
cidr_plot <- threshreldiag(F_t[['cidr']], data_ts$los, t = t_vec, title = "CIDR")
mond_plot <- threshreldiag(F_t[['mond']], data_ts$los, t = t_vec, title = "LB")
pit_plot <- gridExtra::grid.arrange(lspm_plot, cidr_plot, mond_plot, nrow = 1)
#ggsave(plot = pit_plot, "ICU_tcal_comp.png", width = 10.5, height = 3.5)


## CRPS
crps_mat <- t(sapply(icu_vec, function(icu) {
  ind <- data_ts$icuCode == icu
  c(mean(crps[['lspm']][ind]), mean(crps[['cidr']][ind]), mean(crps[['mond']][ind]))
}))
rownames(crps_mat) <- icu_vec
colnames(crps_mat) <- c("LSPM", "CIDR", "LB")
crps_mat[order(icu_vec), ]


## Thickness of CIDR
ggplot(data.frame(x = thick[['cidr']])) +
  geom_histogram(aes(x = x, y = ..count.. / sum(..count..)), boundary = 0, binwidth = 0.025, fill = "lightgrey", col = "darkgrey") +
  scale_x_continuous(name = "Thickness", limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(name = "Relative frequency", expand = expansion(c(0, 0.15))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.margin = margin(c(5.5, 10.5, 5.5, 5.5)))
#ggsave("ICU_thick_cidr.png", width = 3.7, height = 2.5)


## Thickness vs covariates
plot(data_ts$los, thick[['cidr']])

df <- data.frame(th = sapply(icu_vec, function(i) mean(thick[['cidr']][data_ts$icuCode == i])),
                 n = sapply(icu_vec, function(i) sum(data_tr$icuCode == i)))
ggplot(df) + geom_point(aes(x = n, y = th)) +
  scale_x_continuous(name = "Sample size") +
  scale_y_continuous(name = "Average thickness") +
  theme_bw() +
  theme(panel.grid = element_blank())
#ggsave("ICU_thick_n.png", width = 3.7, height = 2.5)


## plot example cdf's

ind <- sample(1:2528, 1) #440/1202/1708
x <- cidr_preds$points
y <- ecdf(lspm_preds[ind, ][lspm_preds[ind, ] > 0])(x)
df <- data.frame(x = x, y = c(y, cidr_preds$cdf_oos[, ind]), mth = rep(c("LSPM", "CIDR"), each = length(x)))
ggplot(df) + geom_step(aes(x = x, y = y, col = mth)) +
  scale_x_continuous(name = "LOS", limits = c(0, 10)) +
  scale_y_continuous(name = "Predictive CDF") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")
#ggsave("ICU_dist_ex.png", width = 3.2, height = 2.7)
