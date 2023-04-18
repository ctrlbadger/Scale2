cauchy.scale <- function(x, cauchy.x_0 = 0, cauchy.gamma=1) (x - cauchy.x_0)/ cauchy.gamma

cauchy.pi <- function(x, cauchy.x_0=0, cauchy.gamma=1) {
  (3.1415 * cauchy.gamma * (1 + cauchy.scale(x, cauchy.x_0, cauchy.gamma)^2))^-1
}

cauchy_density_path <- function(path_hist, weight_hist, debug_hist, idx = seq_along(debug_hist), cauchy.x_0, cauchy.gamma) {
  text.res <- 10000
  x_lap <- seq(from = -15, to = 15, length.out = text.res)
  y_lap <- cauchy.pi(x_lap, cauchy.x_0, cauchy.gamma)

  cauchy_trbl <- tibble(x_lap, y_lap)

  debug_trbl <- debug_hist[idx] %>%
    map(as_tibble) %>%
    reduce(add_row) %>%
    mutate(mesh_idx = factor(mesh_idx))

  ggplot() +
    geom_density(data = debug_trbl, aes(x = path_curr, weight = norm_weight, color = mesh_idx, fill = mesh_idx), alpha = 0.2) +
    geom_line(data = cauchy_trbl, aes(x = x_lap, y = y_lap), linewidth = 1)


}

cauchy_density_estimate <- function(debug_hist, idx = seq_along(debug_hist), cauchy.x_0, cauchy.gamma) {
  text.res <- 10000
  x_lap <- seq(from = -15, to = 15, length.out = text.res)
  y_lap <- cauchy.pi(x_lap, cauchy.x_0, cauchy.gamma)

  cauchy_trbl <- tibble(x_lap, y_lap)

  debug_trbl <- debug_hist[idx] %>%
    map(as_tibble) %>%
    imap(., ~ mutate(.x, path_curr, norm_weight = norm_weight / length(idx), .keep='used', mesh_idx = idx[.y])) %>%
    map(as_tibble) %>%
    reduce(add_row)



  ggplot() +
    geom_density(data = debug_trbl, aes(x = path_curr, weight = norm_weight), alpha = 0.2) +
    geom_line(data = cauchy_trbl, aes(x = x_lap, y = y_lap), linewidth = 1)
}





