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

  set.seed(150)
  beta_true <- c(5, -2)
  n <- 20
  d <- 2
  x_1 <- rep(1, n)
  x_2 <- rnorm(n, mean=0, sd=1)
  y <- beta_true[1]*x_1 + beta_true[2] * x_2
  x <- matrix(c(x_1, x_2), ncol = 2)

  lm1 <- lm(y ~ x - 1)
  bivariate_normal_n20_b5neg2 <- list(
    n = n, d = d, x = x, y = y, beta_true = beta_true
  )
  use_data(bivariate_normal_n20_b5neg2, overwrite = FALSE)

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

  # filename <- glue("2318NormalN(mu={mu},s={sigma}){t_inc=} {iterations=} {kill_time=} {num_particles=} {rescale=} {subsample=}.rdata", .transformer = vv_transformer)
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

