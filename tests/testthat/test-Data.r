test_that("Small Normal Data:", {
  list2env(small_normal_n10_m5_s5, rlang::current_env())

  dist_data <- MeanNormalData$new(x_data, mu_true = mu, sigma_true = sigma)

  expect_snapshot(
    dist_data
  )


  lim_max <- 10
  lim_length <- 1000
  mu_x <- seq(-lim_max, lim_max, length.out = lim_length)


  # Testing Phi Functions
  phi_test_plot <- function() {
    mu_x <- seq(-lim_max, lim_max, length.out = lim_length)
    phi_actual <- map_dbl(mu_x, ~ dist_data$phi(dist_data$x_unscale(.x)))
    phi_estimator <- map_dbl(mu_x, ~ dist_data$phi_estimator(dist_data$x_unscale(.x), 3, 1))
    phi_lb <- map_dbl(mu_x, ~ dist_data$phi_estimator_bounds(c(dist_data$lambda %*% as.matrix(.x)))$phi_l)
    phi_ub <- map_dbl(mu_x, ~ dist_data$phi_estimator_bounds(c(dist_data$lambda %*% as.matrix(.x)))$phi_u)
    tibble(mu_x, phi_actual, phi_estimator, phi_lb, phi_ub) %>%
      pivot_longer(cols=-1, names_to="phi_method", values_to = "phi") %>%
      ggplot(aes(x=mu_x, y=phi, colour=phi_method)) + geom_line() # + ylim(c(-1000, 1000))
  }

  alpha_test_plot <- function() {
    alpha_actual <- map_dbl(mu_x, ~ dist_data$alpha(dist_data$x_unscale(.x)))
    alpha_estimator <- map_dbl(mu_x, ~ dist_data$alpha_subsample(dist_data$x_unscale(.x), 1))
    alpha_bound <- map_dbl(mu_x, ~ ((dist_data$n + 1) * dist_data$hessian_bound * abs(mu - .x)))
    tibble(mu_x, alpha_actual, alpha_estimator, alpha_bound) %>%
      pivot_longer(cols=-1, names_to="alpha", values_to = "x") %>%
      ggplot(aes(x=mu_x, y = x, colour = alpha)) + geom_line()
  }


  pi_test_plot <- function() {
    pi_actual <- map_dbl(mu_x, ~ dist_data$pi_actual(.x))
    pi_observed <- map_dbl(mu_x, ~ dist_data$pi_observed(.x))
    tibble(mu_x, pi_actual, pi_observed) %>%
      pivot_longer(cols=-1, names_to="pi_method", values_to = "pi") %>%
      ggplot(aes(x=mu_x, y=pi, colour=pi_method)) + geom_line()

  }

  plot(alpha_test_plot())

  plot(phi_test_plot())

  plot(pi_test_plot())
})


test_that("Large Normal Data:", {
  list2env(large_normal_n100_mneg2_s1, rlang::current_env())

  dist_data <- MeanNormalData$new(x_data, mu_true = mu, sigma_true = sigma)

  expect_snapshot(
    dist_data
  )
})


test_that("Centered Cauchy Example", {
  list2env(cauchy_centered, rlang::current_env())

  dist_data <- CauchyData$new(mu, gamma)
  expect_snapshot(
    dist_data
  )
})

test_that("Deviated Cauchy Example", {
  list2env(cauchy_centered, rlang::current_env())

  dist_data <- CauchyData$new(mu, gamma)
  expect_snapshot(
    dist_data
  )
})
