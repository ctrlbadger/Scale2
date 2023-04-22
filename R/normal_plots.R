

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


normal_density_estimate <- function(debug_hist, idx = seq_along(debug_hist), dist_data, unscale = TRUE) {
  if (unscale) {
    x_unscale_vectoriser <- function(x) {
      map_dbl(x, dist_data$x_unscale)
      #x * dist_data$sigma_est + dist_data$x_hat
    }
  }
  else {
    x_unscale_vectoriser <- function(x) x
  }


  lim_max <- 10
  lim_length <- 1000
  mu_x <- seq(-lim_max, lim_max, length.out = lim_length)

  debug_trbl <- debug_hist[idx] %>%
    map(as_tibble) %>%
    imap(., ~ mutate(.x, path_curr = x_unscale_vectoriser(path_curr), norm_weight = norm_weight / length(idx), .keep='used', mesh_idx = idx[.y])) %>%
    map(as_tibble) %>%
    reduce(add_row)

  pi_actual <- map_dbl(mu_x, dist_data$pi_actual)
  pi_observed <- map_dbl(mu_x, dist_data$pi_observed)
  pi_vals <- tibble(mu_x, pi_actual, pi_observed) %>%
    pivot_longer(cols = -1, names_to = "pi_method", values_to="pi") %>% mutate(pi_method = as.factor(pi_method))


  ggplot() +
    geom_line(data = pi_vals, aes(x=mu_x, y=pi, linetype=pi_method), colour="red") +
    geom_density(data = debug_trbl, aes(x = path_curr, weight = norm_weight), alpha = 0.2) +
    ggplot2::xlim(-5, 10)
}

normal_density_path <- function(debug_hist, idx, dist_data, unscale = FALSE) {
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


  pi_actual <- map_dbl(mu_x, dist_data$pi_actual)
  pi_observed <- map_dbl(mu_x, dist_data$pi_observed)
  pi_vals <- tibble(mu_x, pi_actual, pi_observed) %>%
    pivot_longer(cols = -1, names_to = "pi_method", values_to="pi") %>% mutate(pi_method = as.factor(pi_method))

  ggplot() +
    geom_line(data = pi_vals, aes(x=mu_x, y=pi, linetype=pi_method)) +
    geom_density(data = debug_trbl, aes(x = path_curr, weight = norm_weight, colour=mesh_idx, fill=mesh_idx), alpha = 0.2) +
    xlim(-5, 10)
}
