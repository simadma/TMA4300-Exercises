library(tidyverse)  # Collection of R packages designed for data science

library(INLA)       # Full Bayesian Analysis of Latent Gaussian Models using Integrated
# Nested Laplace Approximations

data(coal, package = "boot")

dexp_restricted <- function(x, rate, min, max) {
  if (rate > 0) {
    return(
      if_else(
        condition = min <= x & x <= max,
        true      = rate*exp(-(rate*(x - min) + log1p(-exp(-rate*(max - min))))),
        false     = 0
      )
    )
  } else if (rate < 0) {
    return(
      if_else(
        condition = min <= x & x <= max,
        true      = -rate*exp(-(rate*(x - max) + log1p(-exp(rate*(max - min))))),
        false     = 0
      )
    )
  } else {
    return(dunif(x, min, max))
  }
}

rexp_restricted <- function(n, rate, min, max) {
  u <- runif(n)
  if (rate > 0) {
    return(min - (1/rate)*(log(1 - u) + log1p(u/(1 - u)*exp(-rate*(max - min)))))
  } else if (rate < 0) {
    return(max - (1/rate)*(log(u) + log1p((1 - u)/u*exp(rate*(max - min)))))
  } else {
    return(min + u*(max - min))
  }
}

MH_step_single_site <- function(s, lambda, n) {
  
  alpha <- min(1, exp(log_ratio))
  if (runif(1) < alpha) {
    return(t_0 + t_star*(t_2 - t_0))
  }
  s
}

single_site_McMC <- function(data, t_1, lambda_0, lambda_1, beta, n = 1, n_iter = 5000) {
  x <- data$date
  y <- length(x)  # Number of disasters
  t_0 <- x[1]     # Start date
  t_2 <- x[y]     # End date
  theta <- matrix(nrow = n_iter, ncol = 4,
             dimnames = list(NULL, c("t_1", "lambda_0", "lambda_1", "beta"))
           )
  theta[1, ] <- c(t_1, lambda_0, lambda_1, beta)
  if (n_iter > 1) {
    for (i in 2:n_iter) {
      # Metropolis Hastings step
      # t_star <- rnorm(1, mean = t_1, sd = n)
      # if (t_0 < t_star & t_star < t_2) {
      #   alpha <- min(1, exp(-(lambda_0 - lambda_1)*(t_star - t_1)))
      #   if (runif(1) < alpha) {
      #     t_1 <- t_star
      #   }
      # }
      t_1 <- rexp_restricted(1, lambda_0 - lambda_1, t_0, t_2)
      
      y_0 <- sum(x < t_1)  # Number of disasters before t_1
      y_1 <- y - y_0       # Number of disasters after t_1
      
      # Gibbs step, sample lambda_0
      lambda_0 <- rgamma(1,
                    shape = y_0 + 2,
                    rate  = t_1 - t_0 + 1/beta
                  )
      # Gibbs step, sample lambda_1
      lambda_1 <- rgamma(1,
                    shape = y_1 + 2,
                    rate  = t_2 - t_1 + 1/beta
                  )
      # Gibbs step, sample beta
      beta <- 1 / rgamma(1,
                    shape = 4,
                    rate  = lambda_0 + lambda_1 + 1
                  )
      theta[i, ] <- c(t_1, lambda_0, lambda_1, beta)  # Update
    }
  }
  as_tibble(theta)
}


set.seed(731)
theta <- single_site_McMC(coal, 1890, 10, 10, 10, n = 1, n_iter = 5000)

theta$t_1

## Plot
plot(theta$t_1, type="l")
hist(theta$t_1[(nrow(theta)-1000):nrow(theta)], probability = TRUE, breaks = 50)

# 
# mean(theta$t_1[(nrow(theta)-1000):nrow(theta)])
# plot(theta$beta, type = "l")
# 
# dbeta_four <- function(x, shape1, shape2, start = 0, end = 1) {
#   y <- (x - start)/(end - start)
#   dbeta(y, shape1, shape2)/(end - start)
# }
# 
# rbeta_four <- function(n, shape1, shape2, start = 0, end = 1) {
#   start + rbeta(n, shape1, shape2)*(end - start)
# }
# 
# a <- coal$date[1]
# b <- coal$date[nrow(coal)]
# x <- seq(a, b, l=1000)
# 
# 
# 
# lambda <- 10
# tun_param <- 0.5
# alph <- exp(-tun_param*lambda)
# bet <- exp(tun_param*lambda)
# ye <- dexp_restricted(x, lambda, a, b)
# yb <- dbeta_four(x, alph, bet, a, b)
# 
# plot(x, ye, col = "red", lwd=2, type="l", ylim=c(0, 0.05))
# lines(x, yb, col = "blue", lwd=2, type="l")
# 
# 
# yb <- dbeta_four(x, 1, 1, a, b)
# plot(x, yb, col = "blue", lwd=2, type="l", ylim=c(0, max(yb[yb!=Inf])))
# lines(x, dbeta_four(x, exp(lambda), exp(-lambda), a, b), col = "red", lwd=2, type="l")

set.seed(5)
rgamma(1, 2, 10)
xrgamma(1, 2, rate=10)
