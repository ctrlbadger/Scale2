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