
library(ConformalIDR)
library(ggplot2)

set.seed(2651015)

n <- 512
x <- rnorm(n)
y <- rnorm(n, x)

x_n1 <- seq(-3, 3, 0.75)

fit <- conformal_idr(x = x, y = y, x_out = x_n1, y_out = NA, online = FALSE)

#fit$points[1] <- -5
#fit$points[n + 2] <- 5
#fit$cdf_lcnf[1, ] <- 0
#fit$cdf_ucnf[n + 2, ] <- 1

df <- data.frame(x = fit$points,
                 y = c(as.vector(fit$cdf_oos), as.vector(fit$cdf_ucnf), as.vector(fit$cdf_lcnf)),
                 b = rep(c("Crisp CDF", "Conformal bounds", "Conformal bound2"), each = 9*length(fit$points)),
                 xn = rep(x_n1, each = length(fit$points)))
ggplot(df) + geom_step(aes(x = x, y = y, col = b)) + facet_wrap(vars(xn)) +
  scale_x_continuous(name = "Threshold", limits = c(-5, 5)) +
  scale_y_continuous(name = "CDF") +
  scale_color_manual(values = c("Conformal bounds" = "#F8766D", "Conformal bound2" = "#F8766D", "Crisp CDF" = "#619CFF"),
                     breaks = c("Conformal bounds", "Crisp CDF")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")
ggsave("plots/simstudy_norm.png", width = 8, height = 6)

