################################################################################
##### EUPPBench temperature case study

# load data
load_data <- function(na_prop = 10) {

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
  plot_obj <- plot_example()


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
plot_example <- function() {
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
plot_pred <- function(filename = NULL) {
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
plot_map <- function(lons, lats, z, filename = NULL){
  if (is.matrix(z)) z <- as.vector(z)

    ## elevation data
  dem <- geodata::elevation_global(res = 0.5, path = tempdir())

  xmin <- 2.5  # Western edge
  xmax <- 11 # Eastern edge
  ymin <- 45 # Southern edge
  ymax <- 54 # Northern edge

  region_extent <- terra::ext(xmin, xmax, ymin, ymax)
  cropped_dem <- terra::crop(dem, region_extent)
  dem_df <- terra::as.data.frame(cropped_dem, xy = TRUE)

  colnames(dem_df) <- c("lon", "lat", "elevation")

  ## boundary data
  eur_bords <- rnaturalearth::ne_countries(scale = "large", continent = "europe", returnclass = "sf")

  ## plot data
  df <- data.frame(lat = lats, lon = lons, z = z)

  # Plot the elevation data
  plot_obj <- ggplot() +
    geom_raster(data = dem_df, aes(x = lon, y = lat, fill = elevation)) +
    geom_sf(data = eur_bords, fill = NA, color = "black") +
    geom_point(data = df, aes(lon, lat, color = z)) +
    geom_point(data = df, aes(lon, lat), color = "black", shape = 21) +
    scale_x_continuous(name = "Longitude", expand = c(0, 0), limits = c(xmin, xmax)) +
    scale_y_continuous(name = "Latitude", expand = c(0, 0), limits = c(ymin, ymax)) +
    scale_fill_gradient(name = "Elevation (m)", low = "white", high = "grey40") +
    scale_color_gradientn(name = "Temp.", limits = c(4, 14), colors = colorRampPalette(c("white", "darkred"))(6)) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA)) +
    guides(fill = "none")

  if (!is.null(filename)) {
    ggsave(filename, plot_obj, width = 3, height = 4)
  } else {
    return(plot_obj)
  }
}

# initialise lists to store verification data
verif_lists <- function(ts_obs, t_vec) {

  pcal <<- list(lspm = array(NA, dim(ts_obs)),
                cidr = array(NA, dim(ts_obs)),
                locb = array(NA, dim(ts_obs)))
  score <<- pcal
  thick <<- pcal

  n_t <<- length(t_vec)
  F_t <<- list(lspm = array(NA, c(dim(ts_obs), n_t)),
               cidr = array(NA, c(dim(ts_obs), n_t)),
               locb = array(NA, c(dim(ts_obs), n_t)))
}

# perform cross validation to find the optimal number of bins at each station
conformal_binning_cv <- function(k_vec = c(1, seq(5, 50, 5))) {

  score_mat <- array(NA, c(n_loc, length(k_vec), 4))
  n <- ncol(tr_obs)
  val_ind <- sample(c(rep(T, 0.25*n), rep(F, 0.75*n)))
  for (j in seq_along(stat_ids)) {
    st <- stat_ids[j]
    for (i in 1:4) {
      s <- c("Wi", "Sp", "Su", "Au")[i]
      print(paste0('Forecast at Station: ', st, ' (', j, ' from ', length(stat_ids), ') and Season: ', s))

      ### Get train data
      ind <- (tr_seas == s) & !val_ind
      train <- data.frame(obs = tr_obs[j, ind], ens.mu = tr_fc_mn[j, ind])
      ind <- (tr_seas == s) & val_ind
      val <- data.frame(obs = tr_obs[j, ind], ens.mu = tr_fc_mn[j, ind])

      for (l in seq_along(k_vec)) {
        k <- k_vec[l]
        scores <- conformal_bin(x = train$ens.mu, y = train$obs, x_out = val$ens.mu, y_out = val$obs, k = k)
        score_mat[j, l, i] <- sapply(scores, function(x) x$crps) |> mean()
      }

    }
  }
  k <- sapply(1:n_loc, function(i) sapply(1:4, function(j) k_vec[which.min(score_mat[i, , j])]))

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
    ggsave(plot = pit_plot, filename, width = 1.3*7.5, height = 1.3*2.5)
  } else {
    return(pit_plot)
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
plot_thick <- function(thick, type = "traffic", obs = NULL, times = NULL, x = NULL, ylab = "Temperature", filename = NULL) {

  if (is.matrix(thick)) thick <- as.vector(thick)

  if (type == "traffic") {
    ## traffic light plot of thickness and temperature over time
    df <- data.frame(time = times, y = thick, o = obs)
    df$group <- c("Low", "Medium", "High")[1 + (df$y > 0.25) + (df$y > 0.5)]
    df$group <- factor(df$group, levels = c("Low", "Medium", "High"))
    plot_obj <- ggplot(na.omit(df)) +
      geom_point(aes(x = time, y = o, fill = group), shape = 21) +
      scale_x_datetime(name = NULL) +
      scale_y_continuous(name = ylab) +
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
    ## scatter plot vs covariates
    df <- data.frame(x = x, th = thick)
    plot_obj <- ggplot(df) +
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

