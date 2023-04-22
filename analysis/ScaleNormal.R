library(usethis)
library(rlang)
library(Scale2) # nolint: missing_package_linter.
library(tidyverse)

create_normal_data <- function() {

  ## SMALL
  set.seed(150)
  sigma <- 5
  mu <- 5
  n <- 10
  x_data <- rnorm(n, mean = mu, sd = sigma)

  dist_data <- MeanNormalData$new(x_data, mu_true = mu, sigma_true = sigma)

  small_normal_n10_m5_s5 <- list(
    n = n, mu = mu, sigma = sigma, x_data = x_data
  )

  use_data(small_normal_n10_m5_s5, overwrite = FALSE)

  ## LARGE
  set.seed(150)
  sigma <- 1
  mu <- -2
  n <- 100
  x_data <- rnorm(n, mean = mu, sd = sigma)

  large_normal_n100_mneg2_s1 <- list(
    n = n, mu = mu, sigma = sigma, x_data = x_data
  )

  use_data(large_normal_n100_mneg2_s1, overwrite = FALSE)

  ## Binomial NORMAL DATA
  set.seed(150)
  beta_true <- c(5, -2)
  n <- 20
  d <- 2
  x_1 <- rep(1, n)
  x_2 <- seq(0, 10, length.out = n)
  y <- beta_true[1]*x_1 + beta_true[2] * x_2 + as.matrix(rnorm(n, mean = 0, sd = 1))
  x <- matrix(c(x_1, x_2), ncol = 2)

  lm1 <- lm(y ~ x - 1)
  summary(lm1)
  bivariate_normal_n20_b5neg2 <- list(
    n = n, d = d, x = x, y = y, beta_true = beta_true
  )
  use_data(bivariate_normal_n20_b5neg2, overwrite = TRUE)

}

test_SCALE_mean_normal <- function() {

  dist_data <- small_normal_n10_m5_s5

  t_inc <- 0.01
  theta <- sqrt(t_inc)
  iterations <- 1000
  kill_time <- iterations * t_inc
  num_particles <- 1000
  rescale <- TRUE
  subsample <- TRUE

  dist_data <- MeanNormalData$new(data_x, mu_true = mu, sigma_true = sigma)

  filename <- glue("2318NormalN(mu={mu},s={sigma}){t_inc=} {iterations=} {kill_time=} {num_particles=} {rescale=} {subsample=}.rdata", .transformer = vv_transformer)
  print(paste(filename))
  # tryCatch(
  #   error = function(cnd) {
  #     print(cnd)
  #     print(filename)
  #   },
  SCALE_info <- SCALE(num_particles = num_particles, d = 1, theta = theta, num_meshes = iterations, kill_time = kill_time, data = dist_data,
                      ess_thresh = 0.5, parallel = FALSE, resample_every = Inf, rescale = rescale, subsample = subsample)

  # )

  save(SCALE_info, dist_data, file = filename)
  print(paste("Successfuly Saved", filename))

  list2env(SCALE_info, envir=global_env())

  # debug_hist <- SCALE_info$debug_hist





  idx_sample <- seq_along(debug_hist)
  id_trace_plot <- id_trace_path(debug_hist, idx = idx_sample, show_lines = FALSE, show_points = TRUE, sample_particles = TRUE, show_resample = TRUE)
  plot(id_trace_plot)

  idx_sample <- floor(seq(from = length(debug_hist) / 10, to = length(debug_hist)))
  density_plot <- norm_density_estimate(debug_hist, idx = idx_sample, dist_data, unscale=TRUE)
  trace_plot <- trace_path(debug_hist)
  ess_plot <- plot_ess(debug_hist, show_ess_thresh = TRUE)
  idx_sample <- floor(seq(from = 1, to = length(debug_hist), length.out = 10))
  density_path <- normal_density_path(debug_hist, idx = idx_sample, dist_data, unscale=TRUE)
  weight_path <- normal_weight_point_path(debug_hist, idx = idx_sample, dist_data, unscale=TRUE)

  #library(plotly)
  # ggplotly(density_path)
  # plot(density_path)
  scale_plot <- plot_grid(trace_plot, ess_plot, weight_path, density_plot)
  plot(scale_plot)
  save(SCALE_info, weight_path, id_trace_plot, dist_data, density_path, trace_plot, density_plot, scale_plot, ess_plot, file = filename)

  print(paste("Successfuly Saved with Plot", filename))

  return()

  ggplotly(density_path)
  ggplotly(density_plot)
  ggplotly(weight_path)
  plot(trace_plot)
  plot(ess_plot)

  # idx <- floor(seq(from = 2, to = length(debug_hist), length.out = 1000))
  get_mean_sd <- function(debug_hist, idx, unscale=TRUE) {
    x_unscale_vectoriser <- function(x) {
      map_dbl(x, dist_data$x_unscale)
    }
    debug_trbl <- debug_hist[idx] %>%
      map(as_tibble) %>%
      imap(., ~ mutate(.x, path_curr = x_unscale_vectoriser(path_curr), norm_weight = norm_weight / length(idx), mesh_idx = idx[.y], .keep='used')) %>%
      map(as_tibble) %>%
      reduce(add_row) %>% mutate(mesh_idx = as_factor(mesh_idx))


    sample_density <- density(debug_trbl$path_curr, weights=debug_trbl$norm_weight)
    # print(paste("Mean:", mean(sample_density$y), "SD", sd(sample_density$y)))

    print(str(trap_mean_sd(sample_density$x, sample_density$y)))
    # plot(sample_density)
  }
  get_mean_sd(debug_hist, idx)




}

bivariate_SCALE_normal <- function() {
  list2env(bivariate_normal_n20_b5neg2, rlang::current_env())
  # n = n, d = d, x = x, y = y, beta_true = beta_true
  bi_normal <- Normal$new(y, x)

  t_inc <- 0.001
  theta <- sqrt(t_inc)
  iterations <- 100
  kill_time <- iterations * t_inc
  num_particles <- 1000
  rescale <- TRUE
  subsample <- FALSE

  set.seed(150)
  SCALE_info <- SCALE(num_particles = num_particles, d = 2, theta = theta, num_meshes = iterations, kill_time = kill_time, data = bi_normal,
                      ess_thresh = 0.5, resample_every = Inf, rescale = rescale, subsample = subsample, print_updates = TRUE)

  idx <- floor(seq(1, iterations, length.out = 10))
  density_path(SCALE_info, idx = idx, unscale = TRUE)
  density_total(SCALE_info, idx = 10:iterations, unscale = TRUE)
  plot_ess(SCALE_info, show_ess_thresh = TRUE)
  GLM_trace_path(SCALE_info)
}

SCALE_mean_normal <- function() {

  data <- small_normal_n10_m5_s5

  t_inc <- 0.01
  theta <- sqrt(t_inc)
  iterations <- 100
  kill_time <- iterations * t_inc
  num_particles <- 1000
  rescale <- TRUE
  subsample <- TRUE

  dist_data <- MeanNormalData$new(data$x_data, mu_true = data$mu, sigma_true = data$sigma)

  # filename <- glue("NormalN(mu={mu},s={sigma}){t_inc=} {iterations=} {kill_time=} {num_particles=} {rescale=} {subsample=}.rdata", .transformer = vv_transformer)
  filename <- "small_normal_n10_m5_s5.rdata"
  set.seed(150)
  SCALE_info <- SCALE(num_particles = num_particles, d = 1, theta = theta, num_meshes = iterations, kill_time = kill_time, data = dist_data,
                      ess_thresh = 0.5, parallel = FALSE, resample_every = Inf, rescale = rescale, subsample = subsample, print_updates = TRUE)

  save(SCALE_info, dist_data, file = filename)
  print(paste("Successfuly Saved", filename))
  list2env(SCALE_info, envir=rlang::global_env())

  GLM_trace_path(SCALE_info, unscale = TRUE)
  plot_ess(SCALE_info, show_ess_thresh = TRUE)
  
  iterations <- SCALE_info$parameters$num_meshes
  idx <- floor(seq(1, iterations, length.out = 10))
  density_path(SCALE_info, idx = idx, unscale = TRUE)
  idx <- floor(seq(iterations * 0.1, iterations))
  density_total(SCALE_info, idx =idx, unscale = TRUE)
}

SCALE_binomial_normal <- function() {
  list2env(small_logistic_example, rlang::current_env())
  model_glm <- glm(y ~ x - 1, family = binomial)

  binomial_data <- Binomial$new(y = y, x = x)

  t_inc <- 0.1
  theta <- sqrt(t_inc)
  iterations <- 100
  kill_time <- iterations * t_inc
  num_particles <- 100
  rescale <- TRUE
  subsample <- FALSE

  set.seed(150)
  filename <- "SCALE_small_logistic_example"
  SCALE_info <- SCALE(num_particles = num_particles, d = 2, theta = theta, num_meshes = iterations, kill_time = kill_time, data = binomial_data,
                      ess_thresh = 0.5, resample_every = Inf, rescale = rescale, subsample = subsample, print_updates = TRUE)
  
  save(SCALE_info, binomial_data, file = filename)
  print(paste("Successfuly Saved", filename))
  list2env(SCALE_info, envir=rlang::global_env())

  GLM_trace_path(SCALE_info, unscale = TRUE)
  plot_ess(SCALE_info, show_ess_thresh = TRUE)
  
  iterations <- SCALE_info$parameters$num_meshes
  idx <- floor(seq(1, iterations, length.out = 10))
  density_path(SCALE_info, idx = idx, unscale = TRUE)
  idx <- floor(seq(iterations * 0.1, iterations))
  density_total(SCALE_info, idx =idx, unscale = TRUE)

  idx_sample <- floor(seq(from = length(debug_hist) / 10, to = length(debug_hist)))
  
}

bivariate_SCALE_binomial <- function() {
  set.seed(1)
  dsz             <<- 10                       # Size of data set
  dimen           <<- 2                       # Dimensionality (>1)
  examp.design <<- matrix(c(rep(1,dsz),(-1)^gtools::odd(1:(dsz*(dimen-1)))/c(1:dsz)),nrow=dsz,ncol=2,byrow=FALSE)
  examp.data <<- c(rep(1,2),rep(0,8))

  
  small_logistic_example <- list(
    n = 10, d = 2, x = examp.design, y = examp.data
  )

  model_glm <- glm(examp.data ~ examp.design - 1, family = binomial)
  use_data(small_logistic_example, overwrite = FALSE)


}