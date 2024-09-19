################################################################################
## set up

set.seed(91)

library(ncdf4)
library(ggplot2)
library(zoo)
library(scoringRules)
library(WeightedForecastVerification)
library(isodistrreg)
library(Rcpp)

source("online_idr.R")
sourceCpp("online_idr_computation_weights.cpp")

source("utility_funcs.R")


################################################################################
## load data

path <- "C:/Users/sa20i493/Documents/Data/EUMetNet/t2m_station_"
#path <- "C:/Users/sa20i493/Documents/Data/EUMetNet/tp6_station_"
load_data <- function(path, var = "t2m", na_prop = 10) {

  ### train data

  ## fcst
  fcst_file <- nc_open(paste0(path, "refo_fc.ncdf4"))
  train_stid_fc <- ncvar_get(fcst_file, varid = "station_id")
  train_lon_fc <- ncvar_get(fcst_file, varid = "station_longitude")
  train_lat_fc <- ncvar_get(fcst_file, varid = "station_latitude")
  train_year_fc <- ncvar_get(fcst_file, varid = "year")
  train_time_fc <- as.POSIXct(ncvar_get(fcst_file, varid = "time")*(24*60*60), origin = '2017-01-02')
  #train_time_fc <- ncvar_get(fcst_file, varid = "time")
  train_lt_fc <- ncvar_get(fcst_file, varid = "step")
  train_ens_fc <- ncvar_get(fcst_file, varid = "number")
  tr_fc <<- ncvar_get(fcst_file, var)
  if (var == "t2m") tr_fc <<- tr_fc - 273 # convert to celcius
  if (var == "tp6") tr_fc <<- pmax(tr_fc, 0)*1000 # convert from m to mm
  tr_fc <<- aperm(tr_fc, c(5, 1, 2, 4, 3))
  nc_close(fcst_file)

  ## obs
  obs_file <- nc_open(paste0(path, "refo_obs.ncdf4"))
  train_stid_obs <- ncvar_get(obs_file, varid = "station_id")
  train_year_obs <- ncvar_get(obs_file, varid = "year")
  train_time_obs <- as.POSIXct(ncvar_get(obs_file, varid = "time")*(24*60*60), origin = '2017-01-02')
  #train_time_obs <- ncvar_get(obs_file, varid = "time")
  train_lt_obs <- ncvar_get(obs_file, varid = "step")
  tr_obs <<- ncvar_get(obs_file, var)
  if (var == "t2m") tr_obs <<- tr_obs - 273 # convert to celcius
  if (var == "tp6") tr_obs <<- pmax(tr_obs, 0)*1000 # convert from m to mm
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
  ts_fc <<- ncvar_get(fcst_file, var)
  if (var == "t2m") ts_fc <<- ts_fc - 273 # convert to celcius
  if (var == "tp6") ts_fc <<- pmax(ts_fc, 0)*1000 # convert from m to mm
  ts_fc <<- aperm(ts_fc, c(4, 1, 2, 3))
  nc_close(fcst_file)

  ## obs
  obs_file <- nc_open(paste0(path, "1718_obs.ncdf4"))
  test_stid_obs <- ncvar_get(obs_file, varid = "station_id")
  test_time_obs <- as.POSIXct(ncvar_get(obs_file, varid = "time"), origin = '1970-01-01')
  test_lt_obs <- ncvar_get(obs_file, varid = "step")
  ts_obs <<- ncvar_get(obs_file, var)
  if (var == "t2m") ts_obs <<- ts_obs - 273 # convert to celcius
  if (var == "tp6") ts_obs <<- pmax(ts_obs, 0)*1000 # convert from m to mm
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

}
load_data(path, na_prop = 0)

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
plot_example()


## get month and season data
tr_month <- lubridate::month(tr_times)
ts_month <- lubridate::month(ts_times)

tr_seas <- tr_month
tr_seas[tr_month %in% c(12, 1, 2)] <- "Wi"
tr_seas[tr_month %in% c(3, 4, 5)] <- "Sp"
tr_seas[tr_month %in% c(6, 7, 8)] <- "Su"
tr_seas[tr_month %in% c(9, 10, 11)] <- "Au"

ts_seas <- ts_month
ts_seas[ts_month %in% c(12, 1, 2)] <- "Wi"
ts_seas[ts_month %in% c(3, 4, 5)] <- "Sp"
ts_seas[ts_month %in% c(6, 7, 8)] <- "Su"
ts_seas[ts_month %in% c(9, 10, 11)] <- "Au"

table(tr_seas)*length(tr_years) # roughly 1000 training data points in each season
table(tr_month)*length(tr_years) # roughly 350 training data points in each month


## restrict attention to 24h forecasts
n_loc <- length(lats)
lead_time <- 24
stat_ind <- sample(1:length(stat_ids), n_loc)
stat_ids <- stat_ids[stat_ind]
ts_fc <- ts_fc[stat_ind, which(lead_times == lead_time), , ]
ts_obs <- ts_obs[stat_ind, which(lead_times == lead_time), ]
tr_fc <- tr_fc[stat_ind, which(lead_times == lead_time), , , ]
tr_fc <- array(tr_fc, c(n_loc, length(tr_years)*length(tr_times), tr_n_ens))
tr_obs <- tr_obs[stat_ind, which(lead_times == lead_time), , ]
tr_obs <- array(tr_obs, c(n_loc, length(tr_years)*length(tr_times)))
tr_month <- rep(tr_month, each = length(tr_years))
tr_seas <- rep(tr_seas, each = length(tr_years))
rm(lead_times, path, stat_ind)


## get ensemble mean forecast
ts_fc_mn <- apply(ts_fc, c(1, 2), mean)
tr_fc_mn <- apply(tr_fc, c(1, 2), mean)


## plot predicted vs observed
ggplot(data.frame(x = ts_fc_mn[1, ], y = ts_obs[1, ])) +
  geom_point(aes(x = x, y = y), size = 0.5) +
  scale_x_continuous(name = "Ensemble mean temperature") +
  scale_y_continuous(name = "Observed temperature") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave("EUMN_data.png", width = 3, height = 3)



## plot stations
plot_map <- function(lons, lats, z, filename = NULL, title = NULL){

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
plot_map(lons, lats, rowMeans(ts_obs), filename = "EUMN_stations.png")


################################################################################
## prediction

pit <- crps <- list(ifs = array(NA, dim(ts_obs)),
                    emos = array(NA, dim(ts_obs)),
                    lspm = array(NA, dim(ts_obs)),
                    idr = array(NA, dim(ts_obs)),
                    cidr = array(NA, dim(ts_obs)),
                    mond = array(NA, dim(ts_obs)))

thick <- list(cidr = array(NA, dim(ts_obs)))


t_vec <- unname(quantile(tr_obs, c(0.1, 0.25, 0.5, 0.75, 0.9)))
n_t <- length(t_vec)
F_t <- list(ifs = array(NA, c(dim(ts_obs), n_t)),
            emos = array(NA, c(dim(ts_obs), n_t)),
            lspm = array(NA, c(dim(ts_obs), n_t)),
            idr = array(NA, c(dim(ts_obs), n_t)),
            cidr = array(NA, c(dim(ts_obs), n_t)),
            mond = array(NA, c(dim(ts_obs), n_t)))



################################################################################
## rolling window

## optimal k for Mondrian prediction
win_len <- 50
k_vec <- c(1, seq(5, win_len, 5))
mondcrps <- matrix(NA, nrow = n_loc, ncol = length(k_vec))
for (j in seq_along(stat_ids)) {
  st <- stat_ids[j]
  mondcrps_temp <- matrix(NA, nrow = length(times), ncol = length(k_vec))
  for (i in seq_along(times)) {
    print(paste0('Forecast at Station: ', st, ' (', j, ' from ', length(stat_ids), ') and Day: ', i))

    if (i > win_len) {
      ### Get train data
      tr_ind <- (i - win_len):(i - 1)
      train <- data.frame(obs = ts_obs[j, tr_ind], ens.mu = ts_fc_mn[j, tr_ind])
      test <- data.frame(obs = ts_obs[j, i], ens.mu = ts_fc_mn[j, i])

      for (k in seq_along(k_vec)) {
        ### Mondrian
        mond_preds <- fit_mond(y = train$obs, X = train$ens.mu, X_ts = test$ens.mu, k_vec[k])
        scores <- eval_mond(mond_preds, test$obs)
        mondcrps_temp[i, k] <- scores$crps
      }
    }
  }
  mondcrps[j, ] <- colMeans(mondcrps_temp, na.rm=T)
}
plot(k_vec, colMeans(mondcrps))
k <- 10

## fit models
win_len <- 50
for (j in seq_along(stat_ids)) {
  st <- stat_ids[j]
  for (i in seq_along(times)) {
    print(paste0('Forecast at Station: ', st, ' (', j, ' from ', length(stat_ids), ') and Day: ', i))

    if (i > win_len) {
      ### Get train data
      tr_ind <- (i - win_len):(i - 1)
      train <- data.frame(obs = ts_obs[j, tr_ind], ens.mu = ts_fc_mn[j, tr_ind])
      test <- data.frame(obs = tS_obs[j, i], ens.mu = ts_fc_mn[j, i])

      ### IFS
      scores <- eval_ifs(ts_fc[j, i, ], ts_obs[j, i])
      pit[['ifs']][j, i] <- scores$pit
      crps[['ifs']][j, i] <- scores$crps
      F_t[['ifs']][j, i, ] <- scores$F_t

      ### EMOS
      emos_preds <- fit_emos(y = train$obs, X = train$ens.mu, X_ts = test$ens.mu)
      scores <- eval_emos(emos_preds, as.numeric(test$obs))
      pit[['emos']][j, i] <- scores$pit
      crps[['emos']][j, i] <- scores$crps
      F_t[['emos']][j, i, ] <- scores$F_t

      ### LSPM
      lspm_preds <- fit_lspm(y = train$obs, X = as.matrix(train$ens.mu), X_ts = as.matrix(test$ens.mu))
      scores <- eval_lspm(lspm_preds, test$obs)
      pit[['lspm']][j, i] <- scores$pit
      crps[['lspm']][j, i] <- scores$crps
      F_t[['lspm']][j, i, ] <- scores$F_t

      ### IDR
      idr_preds <- fit_idr(y = train$obs, X = train$ens.mu, X_ts = test$ens.mu)
      scores <- eval_idr(idr_preds, test$obs)
      pit[['idr']][j, i] <- scores$pit
      crps[['idr']][j, i] <- scores$crps
      F_t[['idr']][j, i, ] <- scores$F_t

      ### CIDR
      cidr_preds <- fit_cidr(y = train$obs, X = train$ens.mu, X_ts = test$ens.mu)
      scores <- eval_cidr(cidr_preds, test$obs)
      pit[['cidr']][j, i] <- scores$pit
      crps[['cidr']][j, i] <- scores$crps
      F_t[['cidr']][j, i, ] <- scores$F_t
      thick[['cidr']][j, i] <- scores$thick

      ### Mondrian
      mond_preds <- fit_mond(y = train$obs, X = train$ens.mu, X_ts = test$ens.mu, k)
      scores <- eval_mond(mond_preds, test$obs)
      pit[['mond']][j, i] <- scores$pit
      crps[['mond']][j, i] <- scores$crps
      F_t[['mond']][j, i, ] <- scores$F_t
    }

  }
}


################################################################################
## monthly window

## optimal k for Mondrian prediction
k_vec <- c(1:4, seq(5, 50, 5))
mondcrps <- array(NA, c(n_loc, length(k_vec), 12))
val_ind <- sample(c(rep(T, 0.25*ncol(tr_obs)), rep(F, 0.75*ncol(tr_obs))))
for (j in seq_along(stat_ids)) {
  st <- stat_ids[j]
  for (m in 1:12) {
    print(paste0('Forecast at Station: ', st, ' (', j, ' from ', length(stat_ids), ') and Month: ', m))

    ### Get train data
    mon_ind <- (tr_month == m) & !val_ind
    train <- data.frame(obs = tr_obs[j, mon_ind], ens.mu = tr_fc_mn[j, mon_ind])
    mon_ind <- (tr_month == m) & val_ind
    test <- data.frame(obs = tr_obs[j, mon_ind], ens.mu = tr_fc_mn[j, mon_ind])

    for (k in seq_along(k_vec)) {
      ### Mondrian
      mond_preds <- fit_mond(y = train$obs, X = train$ens.mu, X_ts = test$ens.mu, k_vec[k])
      scores <- eval_mond(mond_preds, test$obs)
      mondcrps[j, k, m] <- mean(scores$crps)
    }
  }
}
plot(k_vec, apply(mondcrps, 2, mean))
k <- sapply(1:35, function(i) sapply(1:12, function(m) k_vec[which.min(mondcrps[i, , m])]))

## fit models
for (j in seq_along(stat_ids)) {
  st <- stat_ids[j]
  for (m in 1:12) {
    print(paste0('Forecast at Station: ', st, ' (', j, ' from ', length(stat_ids), ') and Month: ', m))

    ### Get train data
    mon_ind <- tr_month == m
    train <- data.frame(obs = tr_obs[j, mon_ind], ens.mu = tr_fc_mn[j, mon_ind])
    mon_ind <- ts_month == m
    test <- data.frame(obs = ts_obs[j, mon_ind], ens.mu = ts_fc_mn[j, mon_ind])

    ### IFS
    scores <- eval_ifs(ts_fc[j, mon_ind, ], ts_obs[j, mon_ind])
    pit[['ifs']][j, mon_ind] <- scores$pit
    crps[['ifs']][j, mon_ind] <- scores$crps
    F_t[['ifs']][j, mon_ind, ] <- scores$F_t

    ### EMOS
    emos_preds <- fit_emos(y = train$obs, X = train$ens.mu, X_ts = test$ens.mu)
    scores <- eval_emos(emos_preds, as.numeric(test$obs))
    pit[['emos']][j, mon_ind] <- scores$pit
    crps[['emos']][j, mon_ind] <- scores$crps
    F_t[['emos']][j, mon_ind, ] <- scores$F_t

    ### LSPM
    lspm_preds <- fit_lspm(y = train$obs, X = as.matrix(train$ens.mu), X_ts = as.matrix(test$ens.mu))
    scores <- eval_lspm(lspm_preds, test$obs)
    pit[['lspm']][j, mon_ind] <- scores$pit
    crps[['lspm']][j, mon_ind] <- scores$crps
    F_t[['lspm']][j, mon_ind, ] <- scores$F_t

    ### IDR
    idr_preds <- fit_idr(y = train$obs, X = train$ens.mu, X_ts = test$ens.mu)
    scores <- eval_idr(idr_preds, test$obs)
    pit[['idr']][j, mon_ind] <- scores$pit
    crps[['idr']][j, mon_ind] <- scores$crps
    F_t[['idr']][j, mon_ind, ] <- scores$F_t

    ### CIDR
    cidr_preds <- fit_cidr(y = train$obs, X = train$ens.mu, X_ts = test$ens.mu)
    scores <- eval_cidr(cidr_preds, test$obs)
    pit[['cidr']][j, mon_ind] <- scores$pit
    crps[['cidr']][j, mon_ind] <- scores$crps
    F_t[['cidr']][j, mon_ind, ] <- scores$F_t
    thick[['cidr']][j, mon_ind] <- scores$thick

    ### Mondrian
    mond_preds <- fit_mond(y = train$obs, X = train$ens.mu, X_ts = test$ens.mu, k)
    scores <- eval_mond(mond_preds, test$obs)
    pit[['mond']][j, mon_ind] <- scores$pit
    crps[['mond']][j, mon_ind] <- scores$crps
    F_t[['mond']][j, mon_ind, ] <- scores$F_t
  }

}


################################################################################
## seasonal window

## optimal k for Mondrian prediction
k_vec <- c(1, seq(5, 100, 5))
mondcrps <- array(NA, c(n_loc, length(k_vec), 4))
val_ind <- sample(c(rep(T, 0.25*ncol(tr_obs)), rep(F, 0.75*ncol(tr_obs))))
for (j in seq_along(stat_ids)) {
  st <- stat_ids[j]
  for (i_s in 1:4) {
    s <- c("Wi", "Sp", "Su", "Au")[i_s]
    print(paste0('Forecast at Station: ', st, ' (', j, ' from ', length(stat_ids), ') and Season: ', s))

    ### Get train data
    seas_ind <- (tr_seas == s) & !val_ind
    train <- data.frame(obs = tr_obs[j, seas_ind], ens.mu = tr_fc_mn[j, seas_ind])
    seas_ind <- (tr_seas == s) & val_ind
    test <- data.frame(obs = tr_obs[j, seas_ind], ens.mu = tr_fc_mn[j, seas_ind])

    for (k in seq_along(k_vec)) {
      ### Mondrian
      mond_preds <- fit_mond(y = train$obs, X = train$ens.mu, X_ts = test$ens.mu, k_vec[k])
      scores <- eval_mond(mond_preds, test$obs, t_vec)
      mondcrps[j, k, i_s] <- mean(scores$crps)
    }
  }
}
plot(k_vec, apply(mondcrps, 2, mean))
k <- sapply(1:35, function(i) sapply(1:4, function(j) k_vec[which.min(mondcrps[i, , j])]))

## fit models
for (j in seq_along(stat_ids)) {
  st <- stat_ids[j]
  for (i_s in 1:4) {
    s <- c("Wi", "Sp", "Su", "Au")[i_s]
    print(paste0('Forecast at Station: ', st, ' (', j, ' from ', length(stat_ids), ') and Season: ', s))

    ### Get train data
    seas_ind <- tr_seas == s
    train <- data.frame(obs = tr_obs[j, seas_ind], ens.mu = tr_fc_mn[j, seas_ind])
    seas_ind <- ts_seas == s
    test <- data.frame(obs = ts_obs[j, seas_ind], ens.mu = ts_fc_mn[j, seas_ind])

    ### IFS
    scores <- eval_ifs(ts_fc[j, seas_ind, ], ts_obs[j, seas_ind], t_vec)
    pit[['ifs']][j, seas_ind] <- scores$pit
    crps[['ifs']][j, seas_ind] <- scores$crps
    F_t[['ifs']][j, seas_ind, ] <- scores$F_t

    ### EMOS
    emos_preds <- fit_emos(y = train$obs, X = train$ens.mu, X_ts = test$ens.mu)
    scores <- eval_emos(emos_preds, as.numeric(test$obs), t_vec)
    pit[['emos']][j, seas_ind] <- scores$pit
    crps[['emos']][j, seas_ind] <- scores$crps
    F_t[['emos']][j, seas_ind, ] <- scores$F_t

    ### LSPM
    lspm_preds <- fit_lspm(y = train$obs, X = as.matrix(train$ens.mu), X_ts = as.matrix(test$ens.mu))
    scores <- eval_lspm(lspm_preds, test$obs, t_vec)
    pit[['lspm']][j, seas_ind] <- scores$pit
    crps[['lspm']][j, seas_ind] <- scores$crps
    F_t[['lspm']][j, seas_ind, ] <- scores$F_t

    ### IDR
    idr_preds <- fit_idr(y = train$obs, X = train$ens.mu, X_ts = test$ens.mu)
    scores <- eval_idr(idr_preds, test$obs, t_vec)
    pit[['idr']][j, seas_ind] <- scores$pit
    crps[['idr']][j, seas_ind] <- scores$crps
    F_t[['idr']][j, seas_ind, ] <- scores$F_t

    ### CIDR
    cidr_preds <- fit_cidr(y = train$obs, X = train$ens.mu, X_ts = test$ens.mu)
    scores <- eval_cidr(cidr_preds, test$obs, t_vec)
    pit[['cidr']][j, seas_ind] <- scores$pit
    crps[['cidr']][j, seas_ind] <- scores$crps
    F_t[['cidr']][j, seas_ind, ] <- scores$F_t
    thick[['cidr']][j, seas_ind] <- scores$thick

    ### Mondrian
    mond_preds <- fit_mond(y = train$obs, X = train$ens.mu, X_ts = test$ens.mu, k[i_s, j])
    scores <- eval_mond(mond_preds, test$obs, t_vec)
    pit[['mond']][j, seas_ind] <- scores$pit
    crps[['mond']][j, seas_ind] <- scores$crps
    F_t[['mond']][j, seas_ind, ] <- scores$F_t
  }

}


################################################################################
## fixed window

## optimal k for Mondrian prediction
k_vec <- c(1, seq(10, 200, 10))
mondcrps <- matrix(NA, nrow = n_loc, ncol = length(k_vec))
val_ind <- sample(c(rep(T, 0.25*ncol(tr_obs)), rep(F, 0.75*ncol(tr_obs))))
for (j in seq_along(stat_ids)) {
  st <- stat_ids[j]
  print(paste0('Forecast at Station: ', st, ' (', j, ' from ', length(stat_ids), ')'))

  ### Get train data
  train <- data.frame(obs = tr_obs[j, !val_ind], ens.mu = tr_fc_mn[j, !val_ind])
  test <- data.frame(obs = tr_obs[j, val_ind], ens.mu = tr_fc_mn[j, val_ind])

  for (k in seq_along(k_vec)) {
    ### Mondrian
    mond_preds <- fit_mond(y = train$obs, X = train$ens.mu, X_ts = test$ens.mu, k_vec[k])
    scores <- eval_mond(mond_preds, test$obs)
    mondcrps[j, k] <- mean(scores$crps)
  }
}
plot(k_vec, colMeans(mondcrps))
k <- sapply(1:35, function(i) k_vec[which.min(mondcrps[i, ])])

## fit models
for (j in seq_along(stat_ids)) {
  st <- stat_ids[j]
  print(paste0('Forecast at Station: ', st, ' (', j, ' from ', length(stat_ids), ')'))

  ### Get train data
  train <- data.frame(obs = tr_obs[j, ], ens.mu = tr_fc_mn[j, ])
  test <- data.frame(obs = ts_obs[j, ], ens.mu = ts_fc_mn[j, ])

  ### IFS
  scores <- eval_ifs(ts_fc[j, , ], ts_obs[j, ])
  pit[['ifs']][j, ] <- scores$pit
  crps[['ifs']][j, ] <- scores$crps
  F_t[['ifs']][j, , ] <- scores$F_t

  ### EMOS
  emos_preds <- fit_emos(y = train$obs, X = train$ens.mu, X_ts = test$ens.mu)
  scores <- eval_emos(emos_preds, as.numeric(test$obs))
  pit[['emos']][j, ] <- scores$pit
  crps[['emos']][j, ] <- scores$crps
  F_t[['emos']][j, , ] <- scores$F_t

  ### LSPM
  lspm_preds <- fit_lspm(y = train$obs, X = as.matrix(train$ens.mu), X_ts = as.matrix(test$ens.mu))
  scores <- eval_lspm(lspm_preds, test$obs)
  pit[['lspm']][j, ] <- scores$pit
  crps[['lspm']][j, ] <- scores$crps
  F_t[['lspm']][j, , ] <- scores$F_t

  ### IDR
  idr_preds <- fit_idr(y = train$obs, X = train$ens.mu, X_ts = test$ens.mu)
  scores <- eval_idr(idr_preds, test$obs)
  pit[['idr']][j, ] <- scores$pit
  crps[['idr']][j, ] <- scores$crps
  F_t[['idr']][j, , ] <- scores$F_t

  ### CIDR
  cidr_preds <- fit_cidr(y = train$obs, X = train$ens.mu, X_ts = test$ens.mu)
  scores <- eval_cidr(cidr_preds, test$obs)
  pit[['cidr']][j, ] <- scores$pit
  crps[['cidr']][j, ] <- scores$crps
  F_t[['cidr']][j, , ] <- scores$F_t
  thick[['cidr']][j, ] <- scores$thick

  ### Mondrian
  mond_preds <- fit_mond(y = train$obs, X = train$ens.mu, X_ts = test$ens.mu, k[j])
  scores <- eval_mond(mond_preds, test$obs)
  pit[['mond']][j, ] <- scores$pit
  crps[['mond']][j, ] <- scores$crps
  F_t[['mond']][j, , ] <- scores$F_t
}


################################################################################
## results

## PIT histograms
ifs_plot <- pit_hist(pit[['ifs']], ranks = F, ymax = 0.4, xlab = NULL, xticks = F,
                     title = paste("IFS: CRPS =", round(mean(crps[['ifs']], na.rm = T), 3)))
emos_plot <- pit_hist(pit[['emos']], ranks = F, ymax = 0.4, xlab = NULL, xticks = F,
                      title = paste("EMOS: CRPS =", round(mean(crps[['emos']], na.rm = T), 3)))
lspm_plot <- pit_hist(pit[['lspm']], ranks = F, ymax = 0.4, xlab = NULL, xticks = F,
                      title = paste("LSPM: CRPS =", round(mean(crps[['lspm']], na.rm = T), 3)))
idr_plot <- pit_hist(pit[['idr']], ranks = F, ymax = 0.4, xlab = NULL, xticks = F,
                      title = paste("IDR: CRPS =", round(mean(crps[['idr']], na.rm = T), 3)))
cidr_plot <- pit_hist(pit[['cidr']], ranks = F, ymax = 0.4, xlab = NULL, xticks = F,
                      title = paste("CIDR: CRPS =", round(mean(crps[['cidr']], na.rm = T), 3)))
mond_plot <- pit_hist(pit[['mond']], ranks = F, ymax = 0.4, xlab = NULL, xticks = F,
                      title = paste("LB: CRPS =", round(mean(crps[['mond']], na.rm = T), 3)))
pit_plot <- gridExtra::grid.arrange(ifs_plot, emos_plot, lspm_plot,
                                    idr_plot, cidr_plot, mond_plot, nrow = 2)
#ggsave(plot = pit_plot, "EUMN_pit_comp.png", width = 8, height = 5, dpi = 300)
#pit_plot <- gridExtra::grid.arrange(lspm_plot, cidr_plot, mond_plot, nrow = 1)
#ggsave(plot = pit_plot, "EUMN_pit_comp.png", width = 7.5, height = 2.5, dpi = 300)


## PIT reliability diagrams
lspm_plot <- pit_reldiag(pit[['lspm']],
                         title = paste("LSPM: CRPS =", round(mean(crps[['lspm']], na.rm = T), 3)))
cidr_plot <- pit_reldiag(pit[['cidr']],
                         title = paste("CIDR: CRPS =", round(mean(crps[['cidr']], na.rm = T), 3)))
mond_plot <- pit_reldiag(pit[['lspm']],
                         title = paste("LB: CRPS =", round(mean(crps[['mond']], na.rm = T), 3)))
pit_plot <- gridExtra::grid.arrange(lspm_plot, cidr_plot, mond_plot, nrow = 1)
#ggsave(plot = pit_plot, "EUMN_pitrd_comp.png", width = 7.5, height = 2.5)


## Threshold calibration diagrams
ifs_plot <- threshreldiag(F_t[['ifs']], ts_obs, t_vec, title = "IFS")
emos_plot <- threshreldiag(F_t[['emos']], ts_obs, t_vec, title = "EMOS")
lspm_plot <- threshreldiag(F_t[['lspm']], ts_obs, t_vec, title = "LSPM")
idr_plot <- threshreldiag(F_t[['idr']], ts_obs, t_vec, title = "IDR")
cidr_plot <- threshreldiag(F_t[['cidr']], ts_obs, t_vec, title = "CIDR")
mond_plot <- threshreldiag(F_t[['mond']], ts_obs, t_vec, title = "LB")
tc_plot <- gridExtra::grid.arrange(ifs_plot, emos_plot, lspm_plot,
                                   idr_plot, cidr_plot, mond_plot, nrow = 2)
#tc_plot <- gridExtra::grid.arrange(lspm_plot, cidr_plot, mond_plot, nrow = 1)
#ggsave(plot = tc_plot, "EUMN_tcal_comp.png", width = 1.3*7.5, height = 1.3*2.5)


## Thickness
df <- data.frame(time = ts_times, y = colMeans(thick[['cidr']]))
plot_1 <- ggplot(na.omit(df)) + geom_point(aes(x = time, y = y)) +
  scale_x_datetime(name = NULL) +
  scale_y_continuous(name = "Thickness", limits = c(0, 1), expand = c(0, 0.01)) +
  theme_bw() +
  theme(panel.grid = element_blank())

df <- data.frame(y = as.vector(crps[['cidr']]), x = as.vector(thick[['cidr']]))
plot_2 <- ggplot(na.omit(df)) + geom_point(aes(x = x, y = y)) +
  scale_x_continuous(name = "Thickness", limits = c(0, 1), expand = c(0, 0.01)) +
  scale_y_continuous(name = "CRPS", expand = c(0, 0.1)) +
  theme_bw() +
  theme(panel.grid = element_blank())

df <- data.frame(y = as.vector(ts_obs), x = as.vector(thick[['cidr']]))
plot_3 <- ggplot(na.omit(df)) + geom_point(aes(x = x, y = y)) +
  scale_x_continuous(name = "Thickness", limits = c(0, 1), expand = c(0, 0.01)) +
  scale_y_continuous(name = "Temperature", expand = c(0, 1)) +
  theme_bw() +
  theme(panel.grid = element_blank())

df <- data.frame(x = as.vector(ts_obs), y = as.vector(crps[['cidr']]))
plot_4 <- ggplot(na.omit(df)) + geom_point(aes(x = x, y = y)) +
  scale_x_continuous(name = "Temperature", expand = c(0, 1)) +
  scale_y_continuous(name = "CRPS", expand = c(0, 0.1)) +
  theme_bw() +
  theme(panel.grid = element_blank())

df <- data.frame(time = ts_times, y = thick[['cidr']][1, ], s = crps[['cidr']][1, ], o = ts_obs[1, ])
#df$group <- c("Low", "Medium", "High")[ceiling(3*df$y)]
df$group <- c("Low", "Medium", "High")[1 + (df$y > 0.25) + (df$y > 0.5)]
df$group <- factor(df$group, levels = c("Low", "Medium", "High"))
plot_5 <- ggplot(na.omit(df)) + geom_point(aes(x = time, y = o, fill = group), shape = 21) +
  scale_x_datetime(name = NULL) +
  scale_y_continuous(name = "Temperature") +
  scale_fill_manual(name = "Thickness", values = c("green4", "orange3", "red3")) +
  #scale_fill_brewer(name = "Thick.") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")
ggsave("EUMN_thick_obs_ts.png", width = 3.7, height = 3)

plot_6 <- ggplot(na.omit(df)) + geom_point(aes(x = time, y = s, fill = group), shape = 21) +
  scale_x_datetime(name = NULL) +
  scale_y_continuous(name = "CRPS") +
  scale_fill_manual(name = "Thickness", values = c("green4", "orange3", "red3")) +
  #scale_fill_brewer(name = "Thick.") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")
ggsave("EUMN_thick_crps_ts.png", width = 3.7, height = 3)

## Thickness of CIDR
ggplot(data.frame(x = as.vector(thick[['cidr']]))) +
  geom_histogram(aes(x = x, y = ..count.. / sum(..count..)), boundary = 0, binwidth = 0.025, fill = "lightgrey", col = "darkgrey") +
  scale_x_continuous(name = "Thickness", limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(name = "Relative frequency", expand = expansion(c(0, 0.15))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.margin = margin(c(5.5, 10.5, 5.5, 5.5)))
#ggsave("EUMN_thick_cidr.png", width = 3.7, height = 2.5)


## Thickness vs obs
df <- data.frame(x = as.vector(ts_obs), th = as.vector(thick[['cidr']]))
ggplot(df) + geom_point(aes(x = x, y = th)) +
  scale_x_continuous(name = "Observation") +
  scale_y_continuous(name = "Thickness") +
  theme_bw() +
  theme(panel.grid = element_blank())
#ggsave("EUMN_thick_obs.png", width = 3.7, height = 2.5)


## Thickness vs covariates
df <- data.frame(x = as.vector(ts_fc_mn), th = as.vector(thick[['cidr']]))
ggplot(df) +
  #geom_histogram(aes(x = x, y = 10 * ..count.. / sum(..count..)), fill = "lightgrey", col = "darkgrey", alpha = 0.1) +
  geom_point(aes(x = x, y = th)) +
  scale_x_continuous(name = "Ensemble mean") +
  scale_y_continuous(name = "Thickness") +
  theme_bw() +
  theme(panel.grid = element_blank())
#ggsave("EUMN_thick_ens.png", width = 3.7, height = 2.5)




