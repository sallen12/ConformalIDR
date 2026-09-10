################################################################################
## script to load, process and save original EUPP data

set.seed(286791)

library(ncdf4)
library(ggplot2)
library(zoo)
library(lubridate)

na_prop <- 0 # stations with any missing data are removed

path <- "Data/EUMetNet/t2m_station_"

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
tr_fc <- ncvar_get(fcst_file, "t2m")
tr_fc <- tr_fc - 273 # convert to celcius
tr_fc <- aperm(tr_fc, c(5, 1, 2, 4, 3))
nc_close(fcst_file)

## obs
obs_file <- nc_open(paste0(path, "refo_obs.ncdf4"))
train_stid_obs <- ncvar_get(obs_file, varid = "station_id")
train_year_obs <- ncvar_get(obs_file, varid = "year")
train_time_obs <- as.POSIXct(ncvar_get(obs_file, varid = "time")*(24*60*60), origin = '2017-01-02')
train_lt_obs <- ncvar_get(obs_file, varid = "step")
tr_obs <- ncvar_get(obs_file, "t2m")
tr_obs <- tr_obs - 273 # convert to celcius
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
ts_fc <- ncvar_get(fcst_file, "t2m")
ts_fc <- ts_fc - 273 # convert to celcius
ts_fc <- aperm(ts_fc, c(4, 1, 2, 3))
nc_close(fcst_file)

## obs
obs_file <- nc_open(paste0(path, "1718_obs.ncdf4"))
test_stid_obs <- ncvar_get(obs_file, varid = "station_id")
test_time_obs <- as.POSIXct(ncvar_get(obs_file, varid = "time"), origin = '1970-01-01')
test_lt_obs <- ncvar_get(obs_file, varid = "step")
ts_obs <- ncvar_get(obs_file, "t2m")
ts_obs <- ts_obs - 273 # convert to celcius
nc_close(obs_file)


### checks
if (identical(test_lt_fc, test_lt_obs) &
    identical(train_lt_fc, train_lt_obs) &
    identical(test_lt_fc, train_lt_fc)) {
  lead_times <- test_lt_obs
} else {
  stop("Lead times in forecast and observation data do not match")
}
if (identical(test_stid_fc, test_stid_obs) &
    identical(train_stid_fc, train_stid_obs) &
    identical(test_stid_fc, train_stid_obs)) {
  stat_ids <- test_stid_obs
} else {
  stop("Station IDs in forecast and observation data do not match")
}
if (identical(test_lon_fc, train_lon_fc)) {
  lons <- test_lon_fc
} else {
  stop("Station longitudes in train and test data do not match")
}
if (identical(test_lat_fc, train_lat_fc)) {
  lats <- test_lat_fc
} else {
  stop("Station latitudes in train and test data do not match")
}
if (identical(test_time_fc, test_time_obs)) {
  ts_times <- test_time_obs
} else {
  stop("Forecast reference times in forecast and observation test data do not match")
}
if (identical(train_time_fc, train_time_obs)) {
  tr_times <- train_time_obs
} else {
  stop("Forecast reference times in forecast and observation training data do not match")
}
if (identical(train_year_fc, train_year_obs)) {
  tr_years <- train_year_obs
} else {
  stop("Years in forecast and observation training data do not match")
}
n_ens <- length(test_ens_fc)
tr_n_ens <- length(train_ens_fc)

### remove stations with missing data
i <- 1
while (i <= length(stat_ids)) {
  id <- stat_ids[i]
  train_na <- sapply(1:length(lead_times), function(lt) mean(is.na(tr_obs[i, lt, , ])))
  test_na <- sapply(1:length(lead_times), function(lt) mean(is.na(ts_obs[i, lt, ])))
  if (any(100*test_na > na_prop) | any(100*train_na > na_prop)) {
    ts_obs <- ts_obs[-i, , ]
    ts_fc <- ts_fc[-i, , , ]
    tr_obs <- tr_obs[-i, , , ]
    tr_fc <- tr_fc[-i, , , , ]
    stat_ids <- stat_ids[-i]
    lons <- lons[-i]
    lats <- lats[-i]
    print(paste("Station", id, "has been removed due to a high proportion of missing values"))
  } else {
    i <- i + 1
  }
}
n_loc <- length(lats)

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

## restrict attention to 24h forecasts
lead_time <- 24
ts_fc <- ts_fc[, which(lead_times == lead_time), , ]
ts_obs <- ts_obs[, which(lead_times == lead_time), ]
tr_fc <- tr_fc[, which(lead_times == lead_time), , , ]
tr_fc <- array(tr_fc, c(n_loc, length(tr_years)*length(tr_times), tr_n_ens))
tr_obs <- tr_obs[, which(lead_times == lead_time), , ]
tr_obs <- array(tr_obs, c(n_loc, length(tr_years)*length(tr_times)))

## get ensemble mean forecast
ts_fc_mn <- apply(ts_fc, c(1, 2), mean)
tr_fc_mn <- apply(tr_fc, c(1, 2), mean)


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

tr_month <- rep(tr_month, each = length(tr_years))
tr_seas <- rep(tr_seas, each = length(tr_years))

save(ts_fc, ts_obs, tr_fc, tr_obs, tr_fc_mn, ts_fc_mn, lats, lons, stat_ids, tr_times, ts_times, file = "scripts/data/dat_eupp.RData")
