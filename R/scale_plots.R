
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
