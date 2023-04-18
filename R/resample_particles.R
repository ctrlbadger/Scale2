

resample_particles <- function(num_particles, particles, incr_log_weight, ess_thresh, resample_overide, mesh_idx) {
  # Get previous log weights for each particle
  prev_log_weight <- map_dbl(particles, "log_weight")


  # Get unormalised weights from previous weight and recent propogation
  unorm_log_weight <- prev_log_weight + incr_log_weight
  # norm_weight <- unorm_weight / sum(unorm_weight)

  max_log_weight <- max(unorm_log_weight) # Find the maximum log weight
  log_sum_weight <- max_log_weight + c(logSumExp(unorm_log_weight - max_log_weight)) # Calculate the log of the sum of the exponentiated weights
  norm_log_weight <- c(unorm_log_weight - log_sum_weight)
  norm_weight <- c(exp(norm_log_weight))

  ess <- 1 / (num_particles * sum(norm_weight ^ 2))

  # Resample Particles if Effective Sample < User Threshold
  resample <- (ess < ess_thresh) || resample_overide
  # print(paste(ess, "resample:", resample || resample_overide))
  # print(paste("Weight Range", min(norm_weight), max(norm_weight)))
  # print(paste("Weight Sums", sum(norm_weight)))

  #print(paste(resample || resample_overide, round(ess, 5), round(min(norm_weight), 7), round(max(norm_weight), 5)))

  debug_msg <- glue("\n{mesh_idx}: ESS: {round(ess, 5)} Resampling: {resample || resample_overide}\n",
                    "Norm Weights: {round(min(norm_weight), 7)} {round(max(norm_weight), 5)} \n",
                    "Log Incr Weights: {round(min(incr_log_weight), 7)} {round(max(incr_log_weight), 5)} \n")

  if (resample || resample_overide) {

    ## Perform Sampling Step

    # MULTINOMIAL RESAMPLING
    sample_idx <- sample.int(num_particles, replace = TRUE, prob = norm_weight)

    # Stratified Resampling
    # sample_idx <- strat.resamp(norm_weight, num_particles)$p.idx

    ## Loop

    # Residual Resampling
    # sample_idx <- resid.resamp(norm_weight, num_particles)$p.idx


    particles <- particles[sample_idx]




    norm_weight <- rep(1 / num_particles, num_particles)

    # This function is slow (let's rewrite it)
    # walk(1:num_particles, ~ (pluck(particles, .x, "log_weight") <<- -log(num_particles))) # Old
    lw_temp <- -log(num_particles)
    for (i in seq_along(particles)) {
      particles[[i]]$log_weight <- lw_temp
    }
    # walk(particles, ~ list_assign(.x, "log_weight" = lw_temp)) # New

    current_paths <- map_dbl(particles, "path_curr")

    debug_msg <- debug_msg %>%
      glue("Unique: {length(unique(sample_idx)) / num_particles}\n",
           "Mean {mean(current_paths * norm_weight)} SD {sd(current_paths)}")

  } else {
    # Update log weights of particles
    for (i in seq_along(particles)) {
      particles[[i]]$log_weight <- norm_log_weight[i]
    }
    print("")
    current_paths <- map_dbl(particles, "path_curr")
    debug_msg <- debug_msg %>%
      glue("Mean {mean(current_paths * norm_weight)} SD {sd(current_paths)}")
  }
  id <- map_int(1:num_particles, ~ pluck(particles, .x, attr_getter("id")))
  print("")
  print(debug_msg)
  list(particles = particles, ess = ess, norm_weight = norm_weight, resample = resample, id=id)
}
