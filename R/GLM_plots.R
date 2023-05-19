#' Kernel Density Estimates Plot across mesh points
#'
#' @param SCALE_info Return object from SCALE
#' @param idx Mesh indexes to display
#' @param unscale Logical. Unscale particles
#' @param ... Additional Parameters for ggplot
#'
#' @return
#' @export
density_path <- function(SCALE_info, idx = seq_along(debug_hist), unscale = TRUE, ...) {
    # list(debug_hist = debug_hist, parameters = parameters, data = data)

    list2env(SCALE_info, rlang::current_env())

    n <- data$n
    d <- data$d
    num_p <- parameters$num_particles
    iters <- length(idx)

    path_str <- ifelse(unscale, "rescaled_path_curr", "path_curr")
    debug_trbl <- map2(debug_hist[idx], idx, ~
        tibble(beta = unlist((pluck(.x, path_str))),
                beta_dim = factor(rep(1:d, num_particles)),
                mesh_idx = factor(.y),
                norm_weight = rep(pluck(.x , 'norm_weight'), each = d))) %>%
                reduce(add_row)

    ggplot(debug_trbl, aes(x=beta, weight = norm_weight, fill = mesh_idx), ...) +
        geom_density(alpha = 0.25) +
        ggplot2::facet_wrap(ggplot2::vars(beta_dim), scales = "free", labeller = "label_both")
}

#' Total Kernel Density Estimates Plot for ScaLE
#'
#' @param SCALE_info Return object from SCALE
#' @param idx Mesh indexes to include in total
#' @param unscale Logical. Unscale particles
#' @param ... Additional Parameters for ggplot
#'
#' @return
#' @export
density_total <- function(SCALE_info, idx = seq_along(debug_hist), unscale = TRUE, ...) {
    # list(debug_hist = debug_hist, parameters = parameters, data = data)

    list2env(SCALE_info, rlang::current_env())

    n <- data$n
    d <- data$d
    num_p <- parameters$num_particles
    iters <- length(idx)

    path_str <- ifelse(unscale, "rescaled_path_curr", "path_curr")

    debug_trbl <- map2(debug_hist[idx], idx, ~
        tibble(beta = unlist((pluck(.x, path_str))),
                beta_dim = factor(rep(1:d, num_particles)),
                mesh_idx = factor(.y),
                norm_weight = rep(pluck(.x , 'norm_weight') / iters, each = d))) %>%
                reduce(add_row)

    ggplot(debug_trbl, aes(x=beta, weight = norm_weight), ...) +
        geom_density(alpha = 0.25) +
        ggplot2::facet_wrap(ggplot2::vars(beta_dim), scales = "free", labeller = "label_both")
}

