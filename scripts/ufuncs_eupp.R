################################################################################
##### utility functions for EUPPBench temperature case study

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
verif_lists <- function(ts_obs, t_vec, a_vec) {

  pcal <<- list(lspm = array(NA, dim(ts_obs)),
                cidr = array(NA, dim(ts_obs)),
                locb = array(NA, dim(ts_obs)))
  crpsc <<- pcal
  thick <<- pcal

  n_t <<- length(t_vec)
  F_t <<- list(lspm = array(NA, c(dim(ts_obs), n_t)),
               cidr = array(NA, c(dim(ts_obs), n_t)),
               locb = array(NA, c(dim(ts_obs), n_t)))

  n_a <<- length(a_vec)
  cover <<- list(lspm = array(NA, c(dim(ts_obs), n_a)),
                 cidr = array(NA, c(dim(ts_obs), n_a)),
                 locb = array(NA, c(dim(ts_obs), n_a)))
  intsc <<- cover
  widt <<- cover
}

# function to get indices of all dates within 90 days of target date
roll_index <- function(i, win_len, tr_times, ts_times) {
  ind <- (i - win_len):(i + win_len)
  roll_times <- ts_times[i] + days(c(-win_len, win_len))
  tr_ind <- logical(length(tr_times))
  for (k in -2:2) {
    roll_times_mod <- roll_times + years(k)
    c <- 0
    while (any(is.na(roll_times_mod))) {
      c <- c + 1
      roll_times_mod <- ts_times[i] + days(c(-win_len-c, win_len+c)) + years(k)
    }
    tr_ind <- tr_ind | (tr_times >= roll_times_mod[1]  & tr_times <= roll_times_mod[2])
  }
  tr_ind <- replicate(20, tr_ind) |> t() |> as.vector()
  return(tr_ind)
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
  }
  return(tc_plot)
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

# wrapper to plot average interval score
plot_is <- function(is, alpha = NULL, filename = NULL) {
  av_score <- sapply(intsc, apply, 3, mean)
  df <- data.frame(s = as.vector(av_score), a = 1 - alpha, mth = rep(c("LSPM", "CIDR", "CB"), each = length(alpha)))
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
  } else {
    return(is_plot)
  }
}

# wrapper to plot unconditional coverage
plot_cov_unc <- function(cov, alpha = NULL, filename = NULL) {
  if (is.null(alpha)) alpha <- seq(0.05, 0.95, 0.05)
  cov_mat <- sapply(cover, apply, 3, mean)
  df <- data.frame(s = as.vector(cov_mat),
                   a = 1 - alpha,
                   mth = rep(c("LSPM", "CIDR", "CB"), each = length(alpha)))
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

  if (!is.null(filename)) {
    ggsave(plot = cov_plot, filename, width = 5, height = 3, dpi = 300)
  } else {
    return(cov_plot)
  }
}

# wrapper to plot conditional coverage
plot_cov_con <- function(cov, alpha = NULL, x_ts = NULL, n_bins = 10, filename = NULL) {
  if (is.null(alpha)) alpha <- seq(0.05, 0.95, 0.05)
  n_loc <- nrow(x_ts)
  breaks <- seq(0, 1, length.out = n_bins + 1)
  cond_cal_05 <- sapply(1:n_bins, function(i) {
    cov_mat <- sapply(1:n_loc, function(j) {
      q <- quantile(ts_fc_mn[j, ], c(breaks[i], breaks[i + 1])) |> unname()
      if (i == n_bins) {
        ind <- (x_ts[j, ] >= q[1] & x_ts[j, ] <= q[2])
      } else {
        ind <- (x_ts[j, ] >= q[1] & x_ts[j, ] < q[2])
      }
      sapply(cov, function(x) mean(x[j, ind, which(abs(alpha - 0.5) < 1e-10)]))
    })
    rowMeans(cov_mat)
  })
  cond_cal_09 <- sapply(1:n_bins, function(i) {
    cov_mat <- sapply(1:n_loc, function(j) {
      q <- quantile(tr_obs[j, ], c(breaks[i], breaks[i + 1])) |> unname()
      if (i == n_bins) {
        ind <- (x_ts[j, ] >= q[1] & x_ts[j, ] <= q[2])
      } else {
        ind <- (x_ts[j, ] >= q[1] & x_ts[j, ] < q[2])
      }
      sapply(cov, function(x) mean(x[j, ind, which(abs(alpha - 0.1) < 1e-10)]))
    })
    rowMeans(cov_mat)
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
    scale_x_continuous(name = "Ensemble mean (as quantile)", breaks = breaks, limits = c(0, 1)) +
    scale_y_continuous(name = "Empirical coverage", limits = c(0, 1)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.justification = c(0, 0),
          legend.position = c(0.01, 0.01)) +
    guides(colour = guide_legend(nrow = 1))

  if (!is.null(filename)) {
    ggsave(plot = cov_plot, filename, width = 5, height = 3, dpi = 300)
  } else {
    return(cov_plot)
  }
}

# wrapper to plot average prediction interval length
plot_width <- function(width, alpha = NULL, filename = NULL) {
  if (is.null(alpha)) alpha <- seq(0.05, 0.95, 0.05)

  av_width <- sapply(width, apply, 3, mean)
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


  if (!is.null(filename)) {
    ggsave(plot = wid_plot, filename, width = 5, height = 3, dpi = 300)
  } else {
    return(wid_plot)
  }

}
