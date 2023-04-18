test_that("Cauchy Centered", {
  # Create data
  set.seed(150)
  mu <- 0
  gamma <- 1

  cauchy_data <- CauchyData$new(mu, gamma)
  
  t_inc <- 1
  iterations <- 2
  kill_time <- iterations * t_inc
  num_particles <- 10
  rescale <- FALSE
  subsample <- FALSE

  SCALE_info <- SCALE(num_particles = num_particles, d = 1, theta = 5, num_meshes = iterations, kill_time = kill_time, data = cauchy_data, 
                        ess_thresh = 0, parallel = FALSE, resample_every = Inf, rescale = rescale, subsample = subsample)
  
  expect_snapshot(SCALE_info)
})

test_that("Cauchy Scaled", {
  set.seed(150)
    # Create data
  mu <- -2
  gamma <- 5

  cauchy_data <- CauchyData$new(mu, gamma)
  
  t_inc <- 1
  iterations <- 2
  kill_time <- iterations * t_inc
  num_particles <- 10
  rescale <- TRUE
  subsample <- FALSE

  SCALE_info <- SCALE(num_particles = num_particles, d = 1, theta = 5, num_meshes = iterations, kill_time = kill_time, data = cauchy_data, 
                        ess_thresh = 0, parallel = FALSE, resample_every = Inf, rescale = rescale, subsample = subsample)

  expect_snapshot(SCALE_info)
}) 

test_that("Cauchy Resample", {
    # Create data
  set.seed(150)
  mu <- 0
  gamma <- 1

  cauchy_data <- CauchyData$new(mu, gamma)
  
  t_inc <- 1
  iterations <- 2
  kill_time <- iterations * t_inc
  num_particles <- 10
  rescale <- FALSE
  subsample <- FALSE

  SCALE_info <- SCALE(num_particles = num_particles, d = 1, theta = 5, num_meshes = iterations, kill_time = kill_time, data = cauchy_data, 
                        ess_thresh = 1, parallel = FALSE, resample_every = Inf, rescale = rescale, subsample = subsample)

  expect_snapshot(SCALE_info)

})

test_that("Normal Rescale Sampled", {
  list2env(large_normal_n100_mneg2_s1, rlang::current_env())

  t_inc <- 0.001
  theta <- sqrt(t_inc)
  iterations <- 2
  kill_time <- iterations * t_inc
  num_particles <- 10
  rescale <- TRUE
  subsample <- TRUE
  
  set.seed(150)
  dist_data <- MeanNormalData$new(x_data, mu_true = mu, sigma_true = sigma)
    

  SCALE_info <- SCALE(num_particles = num_particles, d = 1, theta = theta, num_meshes = iterations, kill_time = kill_time, data = dist_data, 
                      ess_thresh = 00, resample_every = Inf, rescale = rescale, subsample = subsample)
  expect_snapshot(SCALE_info)
}) 

test_that("Normal Rescale Unsampled", {
  list2env(large_normal_n100_mneg2_s1, rlang::current_env())

  t_inc <- 0.001
  theta <- sqrt(t_inc)
  iterations <- 2
  kill_time <- iterations * t_inc
  num_particles <- 10
  rescale <- TRUE
  subsample <- FALSE
  
  set.seed(150)
  dist_data <- MeanNormalData$new(x_data, mu_true = mu, sigma_true = sigma)
    

  SCALE_info <- SCALE(num_particles = num_particles, d = 1, theta = theta, num_meshes = iterations, kill_time = kill_time, data = dist_data, 
                      ess_thresh = 00, resample_every = Inf, rescale = rescale, subsample = subsample)
  expect_snapshot(SCALE_info)
}) 
