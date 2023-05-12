

# Particle Constructor
new_particle <- function(d, id = NULL) {
  path_vec <- vector("double", length = d)
  structure(
    list(
      d = d, # Dimension of distribution
      path_curr = path_vec, # Current BM Path Position
      path_u = path_vec, # BM upper hyper cube boundary
      path_l = path_vec, # BM lower hyper cube boundary
      path_tau = path_vec, # R^d BM hypercube hitting time location
      minI = path_vec, # {-1, 1}^d => for each dimension: ifelse ( path_tau > path_curr, minI <- 1, minI <- 1 )  )
      theta = path_vec, # BM Hypercube boundary size > 0
      time_tau = path_vec, # {R+}^d BM hypercube hitting time
      bdry = path_vec, # BM hypercube opposite of hitting time location bdry <- path_curr + theta * minI
      tau_order = 1:d, # Index of hitting times sorted in ascending order
      log_weight = 0,  # Currently accrued log weight of particle
      time_curr = 0,   # Current BM time of particle
      time_next = 0,   # Next BM event time of particle
      time_next_is_tau = FALSE, # Boolean for whether at next event time the hypercube needs to be updated
      iter_count = 0L, # Poisson Iteration Count
      phi_l = 0,       # Lower bound of killing rate function phi
      phi_u = 0        # Upper bound of killing rate function phi
    ),
    id = id,             # Integer ID of particle
    class = "particle"
  )
}

# Define particle constructor Helper
particle <- function(d, path_curr, theta, log_weight = NULL, id = NULL) {
  # Get new particle class
  p <- new_particle(d, id)

  # Initialise current path and set upper and lower bounds of hypercube
  p$path_curr <- path_curr

  p$theta <- rep(theta[1], d)


  # Initialise log weight
  if (!is.null(log_weight))  p$log_weight <- log_weight


  p <- center_hypercube(p)



  p
}

## Initialises each particle at time 0 with even weights
init_particles <- function(num_particles, d, theta, data) {
  # Initialise a list of particles with a normal distributed path
  # path_curr <- map(seq_along(d), rnorm())
  # purrr::map(1:num_particles, ~ particle(d, path_curr = rep(runif(1, min=-1, max = 1), d), theta, log_weight = - log(num_particles), id = as.integer(.x)))
  map(1:num_particles, ~ particle(d, path_curr = rnorm(d, mean = 0, sd = 1), theta, log_weight = - log(num_particles), id = as.integer(.x)))
}


# Get new particle event as minimum of BM kingman local bound poisson proccess time, hitting time or time_target
# Changes time_next & time_next_is_tau
get_next_event <- function(p, time_target) {
  # kingman poisson process update event
  kingman_event <- rexp(1, rate = p$phi_u - p$phi_l) + p$time_curr
  # Next tau hitting time
  next_tau_event <- p$time_tau[p$tau_order[1]]

  # Next event is minimum of all 3 events
  p$time_next <- min(kingman_event, time_target, next_tau_event)

  # Is next event a hitting time
  # If true we update the hypercube at next step
  p$time_next_is_tau <- p$time_next == next_tau_event
  p
}

# Bubble sorter to resort tau_order according to which dimension has the first hitting time tau
reorder_tau <- function(tau_order, time_tau, d = length(tau_order)) { # nolint
  ## Bubble sorter to resort tau_order according to which dimension has the first hitting time tau
  if (d == 1) return(1)

  for (i in 2:d) {
    if (time_tau[tau_order[i - 1]] < time_tau[tau_order[i]]) {
      return(tau_order) # return sorted
    } else {
      temp <- tau_order[i - 1]
      tau_order[i - 1] <- tau_order[i]
      tau_order[i] <- temp
    }
  }
  tau_order
}

# Updates hypercube for current brownian motion path
# Updates single dimension of the hypercube
update_hypercube <- function(p) {
  # Get dimension of hitting time boundary reached
  dim_idx <- p$tau_order[1]

  # Get next brownian motion hitting time
  bm_pass_info <- scale::bm.pass(s = p$time_curr, x = p$path_curr[dim_idx], p$theta[dim_idx])
  p$time_tau[dim_idx] <- bm_pass_info$tau
  p$minI[dim_idx] <- bm_pass_info$minI
  p$path_tau[dim_idx] <- bm_pass_info$y

  # Reorder hitting times
  p$tau_order <- reorder_tau(p$tau_order, p$time_tau, p$d)

  # Update lower and upper hitting time boundaries
  p$path_l[dim_idx] <- p$path_curr[dim_idx] - p$theta[dim_idx]
  p$path_u[dim_idx] <- p$path_curr[dim_idx] + p$theta[dim_idx]

  p$bdry <- (p$minI * (p$path_u - p$path_l) + p$path_l + p$path_u) / 2

  p
}

# Centers All dimensions of hypercube around path curr
# Used in initialising
center_hypercube <- function(p) {
  # For each dimension, find hitting time and location and initialise them
  # This function calls to the original scale reference paper
  bm_pass_info <- purrr::map2(p$path_curr, p$theta, ~ bm.pass(s = p$time_curr, x = .x, theta = .y)) |>
    purrr::list_transpose(simplify = TRUE)

  p$path_l <- p$path_curr - p$theta
  p$path_u <- p$path_curr + p$theta

  p$time_tau <- bm_pass_info$tau
  p$path_tau <- bm_pass_info$y
  p$minI <- bm_pass_info$minI

  # Sorted index of hitting times tau
  p$tau_order <- order(p$time_tau)

  # Find opposite boundary to hitting boundary
  p$bdry <- (p$minI * (p$path_u - p$path_l) + p$path_l + p$path_u) / 2

  p
}

# Updates trajectory of particle based conditionally on the hypercube
update_trajectory <- function(p) {
  path_proposal <- numeric(p$d)

  dimensions_to_update <- 1:p$d

  # if (p$time_next_is_tau) {
  #   idx <- p$tau_order[1]
  #   path_proposal[idx] <- p$path_tau[idx]
  #
  #   # Everything but tau
  #   dimensions_to_update <- dimensions_to_update[-idx]
  # }


  for (i in seq_along(dimensions_to_update)) {

    path_proposal[i] <- scale::scale.midpt(
      q = p$time_next,
      s = p$time_curr,
      tau = p$time_tau[i],
      x = p$path_curr[i],
      y = p$path_tau[i],
      minI = p$minI[i],
      bdry = p$bdry[i])$w

  }
  return(path_proposal)
}



# Helper function for returning info about the algorithm
particle_snapshot <- function(particles, ...) {
  purrr::list_transpose(particles) |>  append(values = list(...))
}



#' ScaLE (Scaleable Langevin Exact Algorithm)
#'
#' @param num_particles Integer number of particles
#' @param d Integer dimension of distribution
#' @param theta Brownian motion hypercube size. Positive single value
#' @param num_meshes Integer. Equivalent to iteration count
#' @param kill_time Stop time of the algorithm. Positive numeric
#' @param data An R6 [Data()] class specifying the distribution to be sampled from.
#' @param ess_thresh Effective sample size threshold. numeric from [0, 1]. If empirical ess is lower than threshold, particles will be resampled
#' @param parallel Whether to parallelise the importance sampling step at each iteration. Does not currently have better performance
#' @param resample_every A second way to resample particles at every nth iteration. Numeric integer less than `num_meshes`
#' @param rescale Logical. If TRUE, BM path will be centered at 0 but all calls to phi and boundedness will be scaled by preconditioning matrix and centered, should help with mixing. If FALSE no such scaling will take place.
#' @param subsample Logical. If TRUE, [get_phi_info()] will sample 2 data points for each phi update. If FALSE, use full data set for each phi_bound.
#' @param print_updates Logical. Display updates
#'
#' @return A [list()] containing all input `parameters` and a second list `mesh_hist` of length `num_meshes`.
#' Each `mesh_hist` for any given iteration includes:
#' - `particles` list containing all `particle` classes at
#' - `iter_counter` what iteration count
#' - `incr_log_weight` a vector of incremental log weight for each particle
#' - `norm_weight` A vector of normalised weights
#' - `ess` The empirical effective sample size of the current particles
#' - `resample` Logical, whether the particles were resampled at current iteration
#' - `mesh_idx` Index of current mesh time
#' - `mesh_time` Mesh time of Brownian Motion path particles
#' - `id` A vector of id's of all particles
#' @export
#'
#' @examples
SCALE <- function(num_particles, d, theta, num_meshes, kill_time, data, ess_thresh = 0, parallel = FALSE, resample_every = 10, rescale = FALSE, subsample = FALSE, print_updates = FALSE) {
  parameters <- as.list(environment())
  force(parameters)

  # Create times at which the BM chain is sampled at
  mesh_times <- seq(from = kill_time / num_meshes, to = kill_time, length.out = num_meshes)

  # Initialise Particles
  particles <- init_particles(num_particles, d, theta, data)

  mesh_hist <- list()


  if (print_updates) print(paste("Num Particles:", num_particles, "Mesh Number:", num_meshes, "Kill Time:", kill_time, "Threshold:", ess_thresh))


  for (mesh_idx in seq_along(mesh_times)) {
    propogate_info <- purrr::map(particles, ~ propogate_to_time_target(.x, mesh_times[mesh_idx], data, rescale, subsample), .progress=print_updates) |>
      purrr::list_transpose()

    particles <- propogate_info$p
    incr_log_weight <- propogate_info$incr_log_weight
    iter_counter <- propogate_info$iter_counter


    # Perform Resampling step if needed and calculate log_weight
    if ((mesh_idx %% resample_every) == 0) resample_overide <- TRUE
    else resample_overide <- FALSE

    subsample_info <- resample_particles(num_particles, particles, incr_log_weight, ess_thresh, resample_overide, mesh_idx)
    particles <- subsample_info$particles
    ess <- subsample_info$ess
    norm_weight <- subsample_info$norm_weight
    id <- subsample_info$id
    resample <- subsample_info$resample
    print_msg <- subsample_info$print_msg


    # Create snapshot instance
    if (rescale) {
      for (i in seq_along(particles)) {
        particles[[i]]$rescaled_path_curr <- data$x_unscale(particles[[i]]$path_curr)
      }
    }
    mesh_snapshot <- particle_snapshot(particles,
                                        iter_counter = iter_counter,
                                        incr_log_weight = incr_log_weight,
                                        norm_weight = norm_weight, ess = ess, resample = resample,
                                        mesh_idx = mesh_idx, mesh_time = mesh_times[mesh_idx], id = id)
    mesh_hist[[length(mesh_hist) + 1]] <- mesh_snapshot
    if (print_updates) print(print_msg)
  }
  invisible(NULL)
  return(list(debug_hist = mesh_hist, parameters = parameters, data=data))
}
