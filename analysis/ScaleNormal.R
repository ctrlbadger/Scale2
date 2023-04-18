library(usethis)
library(rlang)


create_normal_data <- function() {

  ## SMALL
  set.seed(150)
  sigma <- 5
  mu <- 5
  n <- 10
  x_data <- rnorm(n, mean = mu, sd = sigma)

  dist_data <- MeanNormalData$new(x_data, mu_true = mu, sigma_true = sigma)

  small_normal_n10_m5_s5 <- list(
    n = n, mu = mu, sigma = sigma, x_data = x_data, dist_data
  )

  use_data(small_normal_n10_m5_s5, overwrite = TRUE)

  ## SMALL
  set.seed(150)
  sigma <- 1
  mu <- -2
  n <- 100
  x_data <- rnorm(n, mean = mu, sd = sigma)

  large_normal_n100_mneg2_s1 <- list(
    n = n, mu = mu, sigma = sigma, x_data = x_data
  )

  use_data(large_normal_n100_mneg2_s1, overwrite = TRUE)
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

normal_weight_point_path <- function(debug_hist, idx, dist_data, unscale = FALSE) {
  lim_max <- 10
  lim_length <- 1000
  mu_x <- seq(-lim_max, lim_max, length.out = lim_length)


  if (unscale) {
    x_unscale_vectoriser <- function(x) {
      map_dbl(x, dist_data$x_unscale)
      #x * dist_data$sigma_est + dist_data$x_hat
    }
  }else {
    x_unscale_vectoriser <- function(x) x
  }
  debug_trbl <- debug_hist[idx] %>%
    map(as_tibble) %>%
    imap(., ~ mutate(.x, path_curr = x_unscale_vectoriser(path_curr), norm_weight = norm_weight / length(idx), mesh_idx = idx[.y], .keep='used')) %>%
    map(as_tibble) %>%
    reduce(add_row) %>% mutate(mesh_idx = as_factor(mesh_idx))

  print(head(debug_trbl))

  pi_actual <- map_dbl(mu_x, dist_data$pi_actual)
  pi_observed <- map_dbl(mu_x, dist_data$pi_observed)
  pi_vals <- tibble(mu_x, pi_actual, pi_observed) %>%
    pivot_longer(cols = -1, names_to = "pi_method", values_to="pi") %>% mutate(pi_method = as.factor(pi_method))

  ggplot() +
    geom_line(data = pi_vals, aes(x=mu_x, y=pi, linetype=pi_method), inherit.aes = FALSE) +
    geom_point(data = debug_trbl, aes(x = path_curr, y=0, size = norm_weight, colour=mesh_idx), alpha = 0.2, inherit.aes = FALSE) +
    geom_density(data = debug_trbl, aes(x = path_curr, weight = norm_weight, fill=mesh_idx), alpha = 0.2, inherit.aes = FALSE) +
    xlim(-5, 10) +
    guides(size = "none")
  # geom_boxplot(data = debug_trbl, aes(x = path_curr, weight = norm_weight, colour=mesh_idx), inherit.aes = FALSE)
}

id_trace_path <- function(debug_hist, idx = seq_along(debug_hist), show_lines = TRUE, show_points = TRUE, sample_particles = TRUE, show_resample = TRUE) {
  id <- map(debug_hist[idx], "id")
  resample_every <- parameters$resample_every

  debug_trbl <- debug_hist[idx]

  sample_idx<- 1:100
  if (sample_particles)  {
    debug_trbl <- debug_trbl %>% map(., ~ as_tibble(.x)[sample_idx, ])
    id <- id %>% map(., ~ (.x)[sample_idx])
  } else {
    debug_trbl <- debug_trbl %>% map(., ~ as_tibble(.x))
  }

  debug_trbl  <- debug_trbl %>%
    imap(., ~ mutate(.x, path_curr, .keep='used', mesh_idx = idx[.y], location = seq_along(path_curr), norm_weight = norm_weight, resample = resample)) %>%
    map(as_tibble) %>%
    reduce(add_row) %>% mutate(id = as_factor(unlist(id)), mesh_idx = mesh_idx, location = as_factor(location))

  resample_points <- idx[map_lgl(debug_hist, "resample")]


  mean_trace <- debug_trbl %>%
    group_by(mesh_idx) %>% summarise(avg_path = mean(path_curr))

  temp_plot <- ggplot()
  if (show_lines) {
    temp_plot <- temp_plot +
      geom_line(data=debug_trbl, aes(x=mesh_idx, y=path_curr, colour = location), linewidth = 0.95, alpha  = 0.5)
  }
  if (show_points) {
    temp_plot <- temp_plot +
      geom_point(data=debug_trbl, aes(x = mesh_idx, y=path_curr, size=(10*norm_weight)^2, colour = location), alpha  = 0.5)
  }
  if (show_resample) {
    temp_plot <- temp_plot +
      geom_vline(xintercept = resample_points, linetype = 'dotted', color = "red" )
  }
  temp_plot +
    geom_line(data=mean_trace, aes(x = mesh_idx, y=avg_path), linewidth = 1, linetype='longdash') +
    scale_fill_continuous(type = "viridis") +
    guides(color = "none", size="none")
}

plot_ess <- function(debug_hist, show_ess_thresh = FALSE) {
  if (show_ess_thresh) {
    temp_ess_thresh <- NULL
    temp_resample_every <- NULL
    temp_ess_thresh <- parameters$ess_thresh
    # try(temp_ess_thresh <- as.list(parameters[[1]])$ess_thresh)
    # try(temp_resample_every <- as.list(parameters[[1]])$resample_every)
    temp_resample_every <- parameters$resample_every
    if (is.numeric(temp_ess_thresh)) ess_thresh <- temp_ess_thresh
    if (is.numeric(temp_resample_every)) ess_thresh <- temp_ess_thresh
  }
  # print(ess_thresh)
  ess_tbl <- tibble(mesh_idx = seq_along(debug_hist),
                    ess = map_dbl(debug_hist, "ess"),
                    resample =  as.factor(map_lgl(debug_hist, "resample")))

  ess_plot <- ggplot(data=ess_tbl, aes(x = mesh_idx, y = ess, colour = resample)) +
    geom_point()

  if (show_ess_thresh && (ess_thresh <= 1)) {
    ess_plot <- ess_plot + geom_hline(yintercept = ess_thresh, linetype="longdash", linewidth=1)
  }

  ess_plot
}
