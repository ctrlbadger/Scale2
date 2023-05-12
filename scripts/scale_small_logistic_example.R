library(Scale2)

list2env(small_logistic_example, rlang::current_env())
model_glm <- glm(y ~ x - 1, family = binomial)

binomial_data <- Binomial$new(y = y, x = x)


num_particles <- 1024
iterations <- 1000
kill_time <- 50
t_inc <- kill_time / iterations
theta <- sqrt(t_inc)
rescale <- TRUE
subsample <- FALSE
ess_thresh <- 0.5

set.seed(150)
filename <- "SCALE_small_logistic_example_2.rda"
SCALE_info <- SCALE(num_particles = num_particles, d = 2, theta = theta, num_meshes = iterations, kill_time = kill_time, data = binomial_data, ess_thresh = ess_thresh, resample_every = Inf, rescale = rescale, subsample = subsample, print_updates = TRUE)
save(SCALE_info, binomial_data, file = filename)
