
set.seed(2651015)

library(ConformalIDR)
library(ggplot2)

n <- 512
x <- rnorm(n)
y <- rnorm(n, x)

x_n1 <- seq(-3, 3, 0.75)

fit <- conformal_idr(x = x, y = y, x_out = x_n1)
N <- nrow(fit$points)

df <- data.frame(x = as.vector(fit$points),
                 y = c(as.vector(fit$cdf_crisp), as.vector(fit$cdf_upper), as.vector(fit$cdf_lower)),
                 b = rep(c("Crisp CDF", "Conformal bounds", "Conformal bound2"), each = 9*N),
                 xn = rep(x_n1, each = N))
ggplot(df) + geom_step(aes(x = x, y = y, col = b)) + facet_wrap(vars(xn)) +
  scale_x_continuous(name = "Threshold", limits = c(-5, 5)) +
  scale_y_continuous(name = "CDF") +
  scale_color_manual(values = c("Conformal bounds" = "#F8766D",
                                "Conformal bound2" = "#F8766D",
                                "Crisp CDF" = "#619CFF"),
                     breaks = c("Conformal bounds", "Crisp CDF")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

ggsave("plots/simstudy_norm.png", width = 8, height = 6)
