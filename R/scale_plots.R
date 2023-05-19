
#' Trace Path with ID's of Chain
#'
#' @param debug_hist debug_hist of chain
#' @param idx Indexes of chain to show
#' @param show_lines Logical. Show Lines
#' @param show_points Logical. Show Points
#' @param sample_particles Logical. Show ID's
#' @param show_resample Logical. Show Resample Points
#'
#' @return
#' @export
id_trace_path <- function(debug_hist, idx = seq_along(debug_hist), show_lines = TRUE, show_points = TRUE, sample_particles = TRUE, show_resample = TRUE) {

  id <- map(debug_hist[idx], "id")
  resample_every <- parameters$resample_every

  debug_trbl <- debug_hist[idx]

  sample_idx <- 1:100
  if (sample_particles)  {
    debug_trbl <- debug_trbl %>% map(., ~ as_tibble(.x)[sample_idx, ])
    id <- id %>% purrr::map(., ~ (.x)[sample_idx])
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

#' ESS Plot of ScaLE
#'
#' @param SCALE_info ScaLE return object
#' @param show_ess_thresh Logical. Show ESS threshold
#'
#' @return
#' @export
plot_ess <- function(SCALE_info, show_ess_thresh = FALSE) {
  debug_hist <- SCALE_info$debug_hist
  if (show_ess_thresh) {
    ess_thresh <- SCALE_info$parameters$ess_thresh
    resample_every <- SCALE_info$parameters$resample_every
  }
  # print(ess_thresh)
  ess_tbl <- tibble(mesh_idx = seq_along(debug_hist),
                    ess = map_dbl(debug_hist, "ess"),
                    resample =  as.factor(purrr::map_lgl(debug_hist, "resample")))

  ess_plot <- ggplot(data=ess_tbl, aes(x = mesh_idx, y = ess, colour = resample)) +
    ggplot2::geom_point()

  if (show_ess_thresh && (ess_thresh <= 1)) {
    ess_plot <- ess_plot + ggplot2::geom_hline(yintercept = ess_thresh, linetype="longdash", linewidth=1)
  }

  ess_plot
}

#' Trace Path of ScaLE
#'
#' @param debug_hist debug_hist list from SCALE return
#' @param idx Index vector of chains to show
#'
#' @return
#' @export
#'
#' @examples
trace_path <- function(debug_hist, idx = seq_along(debug_hist)) {
  debug_trbl <- debug_hist[idx] %>%
    map(as_tibble) %>%
    imap(., ~ mutate(.x, path_curr, .keep='used', mesh_idx = idx[.y])) %>%
    map(as_tibble) %>%
    reduce(add_row)


  mean_trace <- debug_trbl %>%
    group_by(mesh_idx) %>% summarise(avg_path = mean(path_curr))
  ybins <- abs(diff(range(debug_trbl$path_curr))) / 10
  ggplot(debug_trbl, aes(x=mesh_idx, y=path_curr)) + geom_bin_2d(binwidth=c(1, ybins)) +
    geom_line(data=mean_trace, aes(x = mesh_idx, y=avg_path), linewidth = 1, color='yellow', alpha=0.5) +
    scale_fill_continuous(type = "viridis")

}

#' Trace Path for GLM ScaLE Models
#'
#' @param SCALE_info SCALE return object
#' @param idx Index of mesh points to show
#' @param unscale Logical. Unscale Particles?
#'
#' @return
#' @export
GLM_trace_path <- function(SCALE_info, idx = seq_along(SCALE_info$debug_hist), unscale = TRUE) {
  list2env(SCALE_info, rlang::current_env())
  debug_trbl <- debug_hist[idx] %>%
    map(as_tibble) %>%
    imap(., ~ mutate(.x, path_curr, .keep='used', mesh_idx = idx[.y])) %>%
    map(as_tibble) %>%
    reduce(add_row)



  n <- data$n
  d <- data$d
  num_p <- parameters$num_particles
  iters <- length(idx)

  path_str <- ifelse(unscale, "rescaled_path_curr", "path_curr")
  debug_trbl <- map2(debug_hist[idx], idx, ~
                       tibble(beta = unlist((pluck(.x, path_str))),
                              beta_dim = factor(rep(1:d, num_p)),
                              mesh_idx = factor(.y),
                              norm_weight = rep(pluck(.x , 'norm_weight'), each = d),
                              mesh_time = .x$mesh_time)) %>%
    reduce(add_row) %>%
    group_by(mesh_idx, beta_dim) %>% mutate(beta_average = mean(beta)) %>% dplyr::ungroup()

  ggplot(debug_trbl) +
    ggplot2::geom_bin_2d(aes(x = mesh_time, y=beta)) +
    ggplot2::geom_line(aes(x = mesh_time, y=beta_average), linewidth = 1, color='yellow', alpha=0.5) +
    ggplot2::facet_grid(rows = ggplot2::vars(beta_dim), scales = "free", labeller = "label_both") +
    ggplot2::scale_fill_continuous(type = "viridis")
}
