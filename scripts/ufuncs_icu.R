################################################################################
##### utility functions for ICU length of stay case study

# load data
load_data <- function() {
  load("C:/Users/sa20i493/Documents/Data/icu_data/mdsi.rda")

  icu_vec <<- data$icuCode %>% unique()
  n_icu <<- icu_vec %>% length()

  ## split into train, estimation, and test data
  data <- data %>% arrange(outDate)
  data_tr <- data %>% slice_head(n = nrow(data)*0.8)
  data_ts <- data %>% slice_tail(n = nrow(data)*0.2)
  data_val <- data_tr %>% sample_frac(size = 0.5)
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
    fit <- mgcv::gam(log(los + 1) ~ s(age, k = 3) + sex + planned + readmission +
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
local_binning_cv <- function(k_vec = c(1, seq(10, 100, 10))) {

  score_mat <- matrix(NA, nrow = n_icu, ncol = length(k_vec))
  for (j in seq_along(icu_vec)) {
    print(icu_vec[j])

    ### Get train data
    train <- subset(data_tr, icuCode == icu_vec[j])
    val <- subset(data_val, icuCode == icu_vec[j])

    for (k in seq_along(k_vec)) {
      locb_preds <- conformal_bin(train$los, train$index, val$index, k)
      scores <- eval_locb(locb_preds, val$los, t_vec)
      score_mat[j, k] <- mean(scores$crps)
    }
  }
  plot(k_vec, colMeans(score_mat))
  k <- k_vec[apply(score_mat, 1, which.min)]

  return(k)
}

# plot PIT histograms at individual ICUs
plot_pit_hists <- function(pit, score, ids = c(44, 65, 76, 77), filename = NULL) {
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

# plot examples of predictive cdf's
plot_example <- function(cidr_preds, lspm_preds, filename = NULL) {

  ind <- sample(1:2528, 1) #440/1202/1708
  x <- cidr_preds$points
  y <- ecdf(lspm_preds[ind, ][lspm_preds[ind, ] > 0])(x)
  df <- data.frame(x = x, y = c(y, cidr_preds$cdf_oos[, ind]), mth = rep(c("LSPM", "CIDR"), each = length(x)))
  plot_obj <- ggplot(df) + geom_step(aes(x = x, y = y, col = mth)) +
    scale_x_continuous(name = "LOS", limits = c(0, 10)) +
    scale_y_continuous(name = "Predictive CDF") +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = "bottom")

  if (!is.null(filename)) {
    ggsave(plot_obj, filename, width = 3.2, height = 2.7)
  } else {
    return(plot_obj)
  }

}

