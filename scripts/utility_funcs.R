fit_emos <- function(y, X, X_ts) {
  fit <- crch::crch(y ~ X)
  mu <- predict(fit, newdata = data.frame(X = X_ts))
  sig <- predict(fit, newdata = data.frame(X = X_ts), type = 'scale')
  pred <- cbind(as.numeric(mu), as.numeric(sig))
  names(pred) <- c("mu", "sig")
  return(pred)
}

fit_lspm <- function(y, X, X_ts) {
  fit <- fit_lspm(y, as.matrix(X), as.matrix(X_ts))
  return(fit)
}

fit_idr <- function(y, X, X_ts) {
  fit <- idr(y = y, X = data.frame(X = X))
  pred <- predict(fit, data = data.frame(X = X_ts))
  return(pred)
}

fit_cidr <- function(y, X, X_ts) {
  fit <- conformal_idr(x = X, y = y, x_out = X_ts, y_out = NA, online = FALSE)
  return(fit)
}

fit_mond <- function(y, X, X_ts, k = 1) {
  out <- kmeans(X, k)
  tr_cl <- out$cluster
  ts_cl <- clue::cl_predict(out, as.matrix(X_ts))
  return(list(tr = data.frame(y = y, cl = tr_cl), ts = as.vector(ts_cl)))
}


eval_ifs <- function(ens, obs, t_vec) {
  if (is.vector(ens)) {
    M <- length(ens)
    pit <- rank(c(obs, ens), ties.method = "random")[1]
    F_t <- sapply(t_vec, function(t) mean(ens <= t))
  } else {
    M <- ncol(ens)
    pit <- apply(cbind(obs, ens), 1, rank, ties.method = "random")[1, ]
    F_t <- sapply(t_vec, function(t) rowMeans(ens <= t))
  }
  pit <- runif(length(obs), (pit - 1)/(M + 1), pit/(M + 1))
  crps <- crps_sample(obs, ens)
  return(list(pit = pit, crps = crps, F_t = F_t))
}

eval_emos <- function(preds, obs, t_vec) {
  mu <- preds[, 1]
  sig <- preds[, 2]
  pit <- pnorm(obs, mu, sig)
  crps <- crps_norm(as.numeric(obs), mu, sig)
  F_t <- sapply(t_vec, function(t) pnorm(t, mu, sig))
  return(list(pit = pit, crps = crps, F_t = F_t))
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

eval_idr <- function(preds, obs, t_vec) {
  pit <- pit(preds, obs)
  crps <- crps(preds, obs)
  F_t <- sapply(t_vec, function(t) {
    sapply(1:length(preds), function(i) {
      r <- sum(t >= preds[[i]]$points) + 1
      if (r == 1) {
        return(0)
      } else {
        return(preds[[i]]$cdf[r - 1])
      }
    })
  })
  return(list(pit = pit, crps = crps, F_t = F_t))
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

eval_mond <- function(preds, obs, t_vec) {
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

# evaluate ideal forecaster in gamma simulation study
eval_ideal <- function(preds, obs, anti = FALSE) {
  if (anti) {
    pit <- pgamma(obs, shape = sqrt(10 - preds), scale = pmin(pmax(10 - preds, 1), 6))
    crps <- crps_gamma(obs, shape = sqrt(10 - preds), scale = pmin(pmax(10 - preds, 1), 6))
  } else {
    pit <- pgamma(obs, shape = sqrt(preds), scale = pmin(pmax(preds, 1), 6))
    crps <- crps_gamma(obs, shape = sqrt(preds), scale = pmin(pmax(preds, 1), 6))
  }
  return(list(pit = pit, crps = crps))
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
