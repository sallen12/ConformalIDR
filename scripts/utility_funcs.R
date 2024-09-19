################################################################################
##### EUPPBench temperature data

# load data
load_temp_data <- function(na_prop = 10) {

  path <- "C:/Users/sa20i493/Documents/Data/EUMetNet/t2m_station_"

  ### train data

  ## fcst
  fcst_file <- nc_open(paste0(path, "refo_fc.ncdf4"))
  train_stid_fc <- ncvar_get(fcst_file, varid = "station_id")
  train_lon_fc <- ncvar_get(fcst_file, varid = "station_longitude")
  train_lat_fc <- ncvar_get(fcst_file, varid = "station_latitude")
  train_year_fc <- ncvar_get(fcst_file, varid = "year")
  train_time_fc <- as.POSIXct(ncvar_get(fcst_file, varid = "time")*(24*60*60), origin = '2017-01-02')
  train_lt_fc <- ncvar_get(fcst_file, varid = "step")
  train_ens_fc <- ncvar_get(fcst_file, varid = "number")
  tr_fc <<- ncvar_get(fcst_file, "t2m")
  tr_fc <<- tr_fc - 273 # convert to celcius
  tr_fc <<- aperm(tr_fc, c(5, 1, 2, 4, 3))
  nc_close(fcst_file)

  ## obs
  obs_file <- nc_open(paste0(path, "refo_obs.ncdf4"))
  train_stid_obs <- ncvar_get(obs_file, varid = "station_id")
  train_year_obs <- ncvar_get(obs_file, varid = "year")
  train_time_obs <- as.POSIXct(ncvar_get(obs_file, varid = "time")*(24*60*60), origin = '2017-01-02')
  train_lt_obs <- ncvar_get(obs_file, varid = "step")
  tr_obs <<- ncvar_get(obs_file, "t2m")
  tr_obs <<- tr_obs - 273 # convert to celcius
  nc_close(obs_file)


  ### test data

  ## fcst
  fcst_file <- nc_open(paste0(path, "1718_fc.ncdf4"))
  test_stid_fc <- ncvar_get(fcst_file, varid = "station_id")
  test_lon_fc <- ncvar_get(fcst_file, varid = "station_longitude")
  test_lat_fc <- ncvar_get(fcst_file, varid = "station_latitude")
  test_time_fc <- as.POSIXct(ncvar_get(fcst_file, varid = "time"), origin = '1970-01-01')
  test_lt_fc <- ncvar_get(fcst_file, varid = "step")
  test_ens_fc <- ncvar_get(fcst_file, varid = "number")
  ts_fc <<- ncvar_get(fcst_file, "t2m")
  ts_fc <<- ts_fc - 273 # convert to celcius
  ts_fc <<- aperm(ts_fc, c(4, 1, 2, 3))
  nc_close(fcst_file)

  ## obs
  obs_file <- nc_open(paste0(path, "1718_obs.ncdf4"))
  test_stid_obs <- ncvar_get(obs_file, varid = "station_id")
  test_time_obs <- as.POSIXct(ncvar_get(obs_file, varid = "time"), origin = '1970-01-01')
  test_lt_obs <- ncvar_get(obs_file, varid = "step")
  ts_obs <<- ncvar_get(obs_file, "t2m")
  ts_obs <<- ts_obs - 273 # convert to celcius
  nc_close(obs_file)


  ### checks
  if (identical(test_lt_fc, test_lt_obs) &
      identical(train_lt_fc, train_lt_obs) &
      identical(test_lt_fc, train_lt_fc)) {
    lead_times <<- test_lt_obs
  } else {
    stop("Lead times in forecast and observation data do not match")
  }
  if (identical(test_stid_fc, test_stid_obs) &
      identical(train_stid_fc, train_stid_obs) &
      identical(test_stid_fc, train_stid_obs)) {
    stat_ids <<- test_stid_obs
  } else {
    stop("Station IDs in forecast and observation data do not match")
  }
  if (identical(test_lon_fc, train_lon_fc)) {
    lons <<- test_lon_fc
  } else {
    stop("Station longitudes in train and test data do not match")
  }
  if (identical(test_lat_fc, train_lat_fc)) {
    lats <<- test_lat_fc
  } else {
    stop("Station latitudes in train and test data do not match")
  }
  if (identical(test_time_fc, test_time_obs)) {
    ts_times <<- test_time_obs
  } else {
    stop("Forecast reference times in forecast and observation test data do not match")
  }
  if (identical(train_time_fc, train_time_obs)) {
    tr_times <<- train_time_obs
  } else {
    stop("Forecast reference times in forecast and observation training data do not match")
  }
  if (identical(train_year_fc, train_year_obs)) {
    tr_years <<- train_year_obs
  } else {
    stop("Years in forecast and observation training data do not match")
  }
  n_ens <<- length(test_ens_fc)
  tr_n_ens <<- length(train_ens_fc)

  ### remove stations with missing data
  i <- 1
  while (i <= length(stat_ids)) {
    id <- stat_ids[i]
    train_na <- sapply(1:length(lead_times), function(lt) mean(is.na(tr_obs[i, lt, , ])))
    test_na <- sapply(1:length(lead_times), function(lt) mean(is.na(ts_obs[i, lt, ])))
    if (any(100*test_na > na_prop) | any(100*train_na > na_prop)) {
      ts_obs <<- ts_obs[-i, , ]
      ts_fc <<- ts_fc[-i, , , ]
      tr_obs <<- tr_obs[-i, , , ]
      tr_fc <<- tr_fc[-i, , , , ]
      stat_ids <<- stat_ids[-i]
      lons <<- lons[-i]
      lats <<- lats[-i]
      print(paste("Station", id, "has been removed due to a high proportion of missing values"))
    } else {
      i <- i + 1
    }
  }
  n_loc <<- length(lats)

  ## plot example forecast trajectory
  plot_obj <- plot_temp_example()


  ## restrict attention to 24h forecasts
  lead_time <- 24
  ts_fc <<- ts_fc[, which(lead_times == lead_time), , ]
  ts_obs <<- ts_obs[, which(lead_times == lead_time), ]
  tr_fc <<- tr_fc[, which(lead_times == lead_time), , , ]
  tr_fc <<- array(tr_fc, c(n_loc, length(tr_years)*length(tr_times), tr_n_ens))
  tr_obs <<- tr_obs[, which(lead_times == lead_time), , ]
  tr_obs <<- array(tr_obs, c(n_loc, length(tr_years)*length(tr_times)))


  ## get ensemble mean forecast
  ts_fc_mn <<- apply(ts_fc, c(1, 2), mean)
  tr_fc_mn <<- apply(tr_fc, c(1, 2), mean)


  ## get month and season data
  tr_month <<- lubridate::month(tr_times)
  ts_month <<- lubridate::month(ts_times)

  tr_seas <<- tr_month
  tr_seas[tr_month %in% c(12, 1, 2)] <<- "Wi"
  tr_seas[tr_month %in% c(3, 4, 5)] <<- "Sp"
  tr_seas[tr_month %in% c(6, 7, 8)] <<- "Su"
  tr_seas[tr_month %in% c(9, 10, 11)] <<- "Au"

  ts_seas <<- ts_month
  ts_seas[ts_month %in% c(12, 1, 2)] <<- "Wi"
  ts_seas[ts_month %in% c(3, 4, 5)] <<- "Sp"
  ts_seas[ts_month %in% c(6, 7, 8)] <<- "Su"
  ts_seas[ts_month %in% c(9, 10, 11)] <<- "Au"

  tr_month <<- rep(tr_month, each = length(tr_years))
  tr_seas <<- rep(tr_seas, each = length(tr_years))

  return(plot_obj)
}

# plot example temperature ensemble forecast
plot_temp_example <- function() {
  s <- sample(seq_along(stat_ids), 1)
  t <- sample(seq_along(ts_times), 1)
  df <- data.frame(lt = lead_times,
                   y = c(as.vector(ts_fc[s, , t, ]), ts_obs[s, , t]),
                   m = as.factor(rep(0:n_ens, each = length(lead_times))))

  ggplot(df) + geom_line(aes(x = lt, y = y, col = m)) +
    scale_x_continuous(name = "Lead time (hours)", expand = c(0, 0)) +
    scale_y_continuous(name = "Temperature (C)") +
    scale_color_manual(values = c(rep("grey", n_ens), "black")) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "none")
}

# plot predicted vs observed temperature
plot_temp_pred <- function(filename = NULL) {
  plot_obj <- ggplot(data.frame(x = ts_fc_mn[1, ], y = ts_obs[1, ])) +
    geom_point(aes(x = x, y = y), size = 0.5) +
    scale_x_continuous(name = "Ensemble mean temperature") +
    scale_y_continuous(name = "Observed temperature") +
    theme_bw() +
    theme(panel.grid = element_blank())
  if (!is.null(filename)) {
    ggsave(filename, plot_obj, width = 3, height = 3)
  } else {
    return(plot_obj)
  }
}

# plot map of stations
plot_temp_map <- function(lons, lats, z, filename = NULL, title = NULL){

  world <- map_data("world")
  ind <- (world$long >= 0.95*min(lons) & world$long < 1.05*max(lons)) &
    (world$lat >= 0.95*min(lats) & world$lat <= 1.05*max(lats))
  world <- world[ind, ]

  df <- data.frame(Lat = lats, Lon = lons, z = z)
  plot_obj <- ggplot() + borders("world") +
    geom_point(data = df, aes(x = Lon, y = Lat, fill = z), shape = 21, size = 2) +
    coord_fixed(ylim = range(lats), xlim = range(lons)) +
    scale_fill_fermenter(breaks = seq(5, 12, 1), name = "", palette = "Reds", direction = 1) +
    theme_void() + theme(legend.title = element_blank(), legend.position = "bottom",
                         panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
                         legend.key.width = unit(0.3, "in")) +
    ggtitle(title)


  if (!is.null(filename)) {
    ggsave(filename, plot_obj, width = 2, height = 2.5)
  }

  return(plot_obj)

}

# initialise lists to store verification data
verif_lists <- function(ts_obs, t_vec) {

  pit <<- list(lspm = array(NA, dim(ts_obs)),
               cidr = array(NA, dim(ts_obs)),
               locb = array(NA, dim(ts_obs)))
  score <<- pit
  thick <<- pit

  n_t <<- length(t_vec)
  F_t <<- list(lspm = array(NA, c(dim(ts_obs), n_t)),
               cidr = array(NA, c(dim(ts_obs), n_t)),
               locb = array(NA, c(dim(ts_obs), n_t)))
}

# perform cross validation to find the optimal number of bins at each station
local_binning_cv <- function(k_vec = c(1, seq(5, 100, 5))) {

  score_mat <- array(NA, c(n_loc, length(k_vec), 4))
  n <- ncol(tr_obs)
  val_ind <- sample(c(rep(T, 0.25*n), rep(F, 0.75*n)))
  for (j in seq_along(stat_ids)) {
    st <- stat_ids[j]
    for (i in 1:4) {
      s <- c("Wi", "Sp", "Su", "Au")[i]
      print(paste0('Forecast at Station: ', st, ' (', j, ' from ', length(stat_ids), ') and Season: ', s))

      ### Get train data
      seas_ind <- (tr_seas == s) & !val_ind
      train <- data.frame(obs = tr_obs[j, seas_ind], ens.mu = tr_fc_mn[j, seas_ind])
      seas_ind <- (tr_seas == s) & val_ind
      test <- data.frame(obs = tr_obs[j, seas_ind], ens.mu = tr_fc_mn[j, seas_ind])

      for (k in seq_along(k_vec)) {
        mond_preds <- fit_mond(y = train$obs, X = train$ens.mu, X_ts = test$ens.mu, k_vec[k])
        scores <- eval_mond(mond_preds, test$obs, t_vec)
        score_mat[j, k, i] <- mean(scores$crps)
      }
    }
  }
  plot(k_vec, apply(score_mat, 2, mean))
  k <- sapply(1:n_loc, function(i) sapply(1:4, function(j) k_vec[which.min(score_mat[i, , j])]))

  return(k)
}

# wrapper to plot pit histograms
plot_pit_hists <- function(pit, score, filename = NULL) {
  lspm_plot <- pit_hist(pit[['lspm']], ranks = F, ymax = 0.4, xlab = NULL, xticks = F,
                        title = paste("LSPM: score =", round(mean(score[['lspm']], na.rm = T), 3)))
  cidr_plot <- pit_hist(pit[['cidr']], ranks = F, ymax = 0.4, xlab = NULL, xticks = F,
                        title = paste("CIDR: score =", round(mean(score[['cidr']], na.rm = T), 3)))
  locb_plot <- pit_hist(pit[['locb']], ranks = F, ymax = 0.4, xlab = NULL, xticks = F,
                        title = paste("LB: score =", round(mean(score[['locb']], na.rm = T), 3)))
  pit_plot <- gridExtra::grid.arrange(lspm_plot, cidr_plot, locb_plot, nrow = 1)
  if (!is.null(filename)) {
    ggsave(plot = pit_plot, filename, width = 7.5, height = 2.5)
  } else {
    return(pit_plot)
  }
}

# wrapper to plot pit pp-plots
plot_pit_pp <- function(pit, score, filename = NULL) {
  lspm_plot <- pit_reldiag(pit[['lspm']],
                           title = paste("LSPM: score =", round(mean(score[['lspm']], na.rm = T), 3)))
  cidr_plot <- pit_reldiag(pit[['cidr']],
                           title = paste("CIDR: score =", round(mean(score[['cidr']], na.rm = T), 3)))
  locb_plot <- pit_reldiag(pit[['lspm']],
                           title = paste("LB: score =", round(mean(score[['locb']], na.rm = T), 3)))
  pit_plot <- gridExtra::grid.arrange(lspm_plot, cidr_plot, locb_plot, nrow = 1)
  if (!is.null(filename)) {
    ggsave(plot = pit_plot, filename, width = 7.5, height = 2.5)
  } else {
    return(pit_plot)
  }
}

# wrapper to plot pit pp-plots
plot_tcal <- function(F_t, ts_obs, t_vec, filename = NULL) {
  lspm_plot <- threshreldiag(F_t[['lspm']], ts_obs, t_vec, title = "LSPM")
  cidr_plot <- threshreldiag(F_t[['cidr']], ts_obs, t_vec, title = "CIDR")
  locb_plot <- threshreldiag(F_t[['locb']], ts_obs, t_vec, title = "LB")
  tc_plot <- gridExtra::grid.arrange(lspm_plot, cidr_plot, locb_plot, nrow = 1)
  if (!is.null(filename)) {
    ggsave(plot = tc_plot, , width = 1.3*7.5, height = 1.3*2.5)
  } else {
    return(tc_plot)
  }
}

# wrapper to plot thickness of conformal IDR bands
plot_thick(thick, type = "traffic", filename = NULL) {

  if (type == "traffic") {
    ## traffic light plot of thickness and temperature over time
    df <- data.frame(time = ts_times, y = thick[['cidr']][1, ], s = score[['cidr']][1, ], o = ts_obs[1, ])
    df$group <- c("Low", "Medium", "High")[1 + (df$y > 0.25) + (df$y > 0.5)]
    df$group <- factor(df$group, levels = c("Low", "Medium", "High"))
    plot_obj <- ggplot(na.omit(df)) +
      geom_point(aes(x = time, y = o, fill = group), shape = 21) +
      scale_x_datetime(name = NULL) +
      scale_y_continuous(name = "Temperature") +
      scale_fill_manual(name = "Thickness", values = c("green4", "orange3", "red3")) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            legend.position = "bottom")
    if (!is.null(filename)) {
      ggsave(filename, plot_obj, width = 3.7, height = 3)
    } else {
      return(plot_obj)
    }

  } else if (type == "hist") {
    ## histogram of thicknesses
    ggplot(data.frame(x = as.vector(thick[['cidr']]))) +
      geom_histogram(aes(x = x, y = ..count.. / sum(..count..)), boundary = 0, binwidth = 0.025, fill = "lightgrey", col = "darkgrey") +
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
    ## scatter plot vs covariates
    df <- data.frame(x = ts_fc_mn[1, ], th = as.vector(thick[['cidr']][1, ]))
    ggplot(df) +
      geom_point(aes(x = x, y = th)) +
      scale_x_continuous(name = "Ensemble mean") +
      scale_y_continuous(name = "Thickness") +
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


fit_cidr <- function(y, X, X_ts) {
  fit <- conformal_idr(x = X, y = y, x_out = X_ts, y_out = NA, online = FALSE)
  return(fit)
}

fit_locb <- function(y, X, X_ts, k = 1) {
  out <- kmeans(X, k)
  tr_cl <- out$cluster
  ts_cl <- clue::cl_predict(out, as.matrix(X_ts))
  return(list(tr = data.frame(y = y, cl = tr_cl), ts = as.vector(ts_cl)))
}


eval_lspm <- function(preds, obs, t_vec) {
  rank_y <- rowSums(obs > preds) + 1
  n <- ncol(preds) + 1
  pit <- (rank_y/n) +  runif(length(obs))*((rank_y + 1)/n - rank_y/n)
  crps <- crps_sample(as.numeric(obs), preds) # actual CRPS is infinite
  if (length(obs) > 1) {
    F_t <- sapply(t_vec, function(t) rowMeans(preds <= t))
  } else {
    F_t <- sapply(t_vec, function(t) mean(preds <= t))
  }
  thicc <- 1/sum(!is.na(obs))
  return(list(pit = pit, crps = crps, F_t = F_t, thick = thicc))
}

eval_cidr <- function(preds, obs, t_vec) {
  n <- length(obs)
  pit0 <- function(x, y, z) stepfun(x = x, y = c(0, y))(z)
  pit <- sapply(1:length(obs), function(i) pit0(preds$points, preds$cdf_oos[, i], obs[i]))

  ens <- preds$points[-c(1, length(preds$points))]
  if (length(obs) > 1) {
    w <- pmax(apply(preds$cdf_oos, 2, diff)[-1, ], 0)
    crps <- sapply(1:length(obs), function(i)
      crps_sample(obs[i], ens, w = w[, i]))
  } else {
    w <- pmax(diff(preds$cdf_oos)[-1], 0)
    crps <- crps_sample(obs, ens, w = w)
  }

  F_t <- sapply(t_vec, function(t) {
    sapply(1:length(obs), function(i) {
      pit0(preds$points, preds$cdf_oos[, i], t)
    })
  })

  thicc <- apply(abs(preds$cdf_lcnf - preds$cdf_ucnf), 2, max)

  return(list(pit = pit, crps = crps, F_t = F_t, thick = thicc))
}

eval_locb <- function(preds, obs, t_vec) {
  k <- sort(unique(preds[['tr']]$cl))
  F_x <- lapply(k, function(i) ecdf(preds[['tr']]$y[preds[['tr']]$cl == i]))
  out <- lapply(1:length(obs), function(i) {
    cl_y <- preds[['ts']][i]
    ens <- preds[['tr']]$y[preds[['tr']]$cl == cl_y]
    pit <- F_x[[cl_y]](obs[i])
    crps <- crps_sample(obs[i], ens)
    F_t <- F_x[[cl_y]](t_vec)
    thicc <- 1/length(ens)
    return(list(pit = pit, crps = crps, F_t = F_t, thick = thicc))
  })
  pit <- sapply(out, function(z) z$pit)
  crps <- sapply(out, function(z) z$crps)
  F_t <- t(sapply(out, function(z) z$F_t))
  thicc <- sapply(out, function(z) z$thick)
  return(list(pit = pit, crps = crps, F_t = F_t, thick = thicc))
}


threshreldiag <- function(x, y, t_vec, title = NULL, xlab = "x", ylab = "x_rc", pointSize = NULL, textSize = NULL, spaceLegend = NULL){

  y <- sapply(t_vec, function(t) as.numeric(y <= t))

  na_ind <- is.na(x)
  x <- matrix(x[!na_ind], ncol = length(t_vec))
  y <- matrix(y[!na_ind], ncol = length(t_vec))

  x_rc <- sapply(seq_along(t_vec), function(i) isoreg(x[, i], y[, i])$yf) # values correspond to ORDERED forecast values!
  x <- apply(x, 2, sort)

  df <- data.frame(x = as.vector(x), x_rc = as.vector(x_rc), t = rep(round(t_vec, 2), each = nrow(x)))
  plt <- ggplot(df) + geom_abline(aes(intercept = 0, slope = 1), lty = "dotted") +
    geom_line(aes(x = x, y = x_rc, col = as.factor(t))) +
    scale_x_continuous(name = xlab, limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(name = ylab, limits = c(0, 1), expand = c(0, 0)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.justification = c(0, 1),
          legend.position = c(0.01, 0.99),
          plot.margin = margin(c(5.5, 10.5, 5.5, 5.5))) +
    ggtitle(title)

  if (!is.null(pointSize) && !is.null(textSize) && !is.null(spaceLegend)) {
    plt <- plt + guides(shape = guide_legend(override.aes = list(size = pointSize)),
                        color = guide_legend(override.aes = list(size = pointSize))) +
      theme(legend.title = element_text(size = textSize),
            legend.text  = element_text(size = textSize),
            legend.key.size = unit(spaceLegend, "lines"))
  }

  return(plt)
}
