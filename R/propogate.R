

# For a given particle propogate forward the BM path until we have reached time_target
# This algorithm is the importance sampler step between particle resampling steps
# Returns updated particle, count of iteration time and the log weight between initial and final time of the particle
propogate_to_time_target <- function(p, time_target, data, rescale, subsample) {
  ## Propogates a particle to specified time by updating BM path
  iter_counter <- 0L
  incr_log_weight <- 0

  # Get initial bounds on phi and phi function
  phi_info <- get_phi_info(data, p$path_l, p$path_u, rescale, subsample)
  phi_pi <- phi_info$phi_estimator_hyp
  p$phi_u <- phi_info$phi_u
  p$phi_l <- phi_info$phi_l

  if ((p$phi_u - p$phi_l) < 0) stop("Somethings wrong")

  p <- center_hypercube(p)
  # Get next event
  p <- get_next_event(p, time_target)
  # Propogate particle until time_next is time_target
  while (p$time_next < time_target) {

    # Propogate particle to next event time and calculate weight and new hypercube if necessary
    p$path_curr <- update_trajectory(p)

    incr_log_weight <- incr_log_weight - p$phi_l * (p$time_next - p$time_curr)

    p$time_curr <- p$time_next

    # If particle has reached hitting time tau then update hypercube
    if (p$time_next_is_tau) {
      p <- update_hypercube(p)

      # Get initial bounds on phi and phi function
      phi_info <- get_phi_info(data, p$path_l, p$path_u, rescale, subsample)
      phi_pi <- phi_info$phi_estimator_hyp
      p$phi_u <- phi_info$phi_u
      p$phi_l <- phi_info$phi_l
    } else {
      if ((p$phi_u - phi_pi(p$path_curr)) <= 0) {
        print("Somethings Wrong")
        print(str(p))
      }
      incr_log_weight <- incr_log_weight +
        log(p$phi_u - phi_pi(p$path_curr)) - log(p$phi_u - p$phi_l)
    }


    # Get next event time for particle
    p <- get_next_event(p, time_target)

    iter_counter <- iter_counter + 1L
    p$iter_count <- p$iter_count + 1L

  }
  if (p$time_next != time_target) {
    print("Somethings wrong")
  }

  # Update weights to time_target
  incr_log_weight <- incr_log_weight - p$phi_l * (p$time_next - p$time_curr)

  # Update BM path to final time_target
  p$path_curr <- update_trajectory(p)
  p$time_curr <- p$time_next

  return(list(p = p, iter_counter = iter_counter, incr_log_weight = incr_log_weight))
}
