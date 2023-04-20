test_that("Normal GLM", {
  list2env(small_normal_n10_m5_s5, rlang::current_env())
  y <- x_data
  x <- rep(1, n)



  plot_log_observed_i <- function() {
    beta_test <- seq(from=-5, to=5, length.out=100)
    log_pi_observed_1 <- map_dbl(beta_test, ~ glm_data$log_pi_observed_i(.x, 1))
    log_pi_observed <- map_dbl(beta_test, ~ glm_data$log_pi_observed(.x))
    tbl <- tibble(beta_test, log_pi_observed_1, log_pi_observed) %>%
      pivot_longer(cols=-1, names_to="log_pi_method", values_to="log_y")
    ggplot(tbl, aes(x=beta_test, y=log_y, color = log_pi_method)) +
        geom_line()
  }

  list2env(small_normal_n10_m5_s5, rlang::current_env())

  glm_data <- Normal$new(y, x)
  glm_data$dispersion <- var(y)
  glm_data$v(1)
  glm_data$link(1)
  glm_data$logistic(1)

  glm_data$grad_ll(2, 1) |> print()
  glm_data$v(1)

  glm_data$data_y[1]
  glm_data$data_x[1,]

  # eta_i <- glm_data$eta(beta=1, x=1)
  # glm_data$grad_ll(0, 1)

  # GLM_grad_ll_1 <- map_dbl(beta_test, ~ glm_data$grad_ll(.x, 1))
  # norma_grad_ll_1 <- map_dbl(beta_test, ~ glm_data$grad_ll(.x, 1))
  tol = 10^-4
  list2env(small_normal_n10_m5_s5, rlang::current_env())
  normal_data <- MeanNormalData$new(x_data, mu_true = mu, sigma_true = sigma)
  glm_data$grad_ll(1, 1)
  expect_equal(glm_data$grad_ll(1, 1), -0.3392704, tolerance = tol)
  expect_equal(glm_data$total_grad_ll(1), 2.407614, tolerance = tol)
  normal_data$total_lap_ll(1)
  normal_data$lap_ll(1, 1)
  normal_data$total_lap_ll(1)
  expect_equal(glm_data$lap_ll(1, 1), -0.08152505, tolerance = tol)
  expect_equal(glm_data$total_lap_ll(1), -0.8152505, tolerance = tol)

  expect_equal(glm_data$get_dispersion_param(), 12.26617, tolerance = tol)

  glm_data$get_hessian_bound()

  glm_data$get_preconditioning_matrix()

  phi_test_plot <- function(dist_data) {
    lim_max <- 10
    lim_length <- 1000
    mu_x <- seq(-lim_max, lim_max, length.out = lim_length)
    phi_actual <- map_dbl(mu_x, ~ dist_data$phi(dist_data$x_unscale(.x)))
    phi_estimator <- map_dbl(mu_x, ~ dist_data$phi_estimator(dist_data$x_unscale(.x), 3, 1))
    phi_lb <- map_dbl(mu_x, ~ dist_data$phi_estimator_bounds(c(dist_data$lambda %*% as.matrix(.x)))$phi_l)
    phi_ub <- map_dbl(mu_x, ~ dist_data$phi_estimator_bounds(c(dist_data$lambda %*% as.matrix(.x)))$phi_u)
    tibble(mu_x, phi_actual, phi_estimator, phi_lb, phi_ub) %>%
      pivot_longer(cols=-1, names_to="phi_method", values_to = "phi") %>%
      ggplot(aes(x=mu_x, y=phi, colour=phi_method)) + geom_line()
  }

  phi_test_plot(glm_data)
})


test_that("Bivariate Normal GLM", {
  print("Bivariate Example")
  list2env(bivariate_normal_n20_b5neg2, rlang::current_env())
  # n = n, d = d, x = x, y = y, beta_true = beta_true
  bi_normal <- Normal$new(y, x)
  # expect_snapshot(bi_normal)

  beta <- c(1, 6)
  bi_normal$log_pi_observed_i(beta, 1)

  expect_equal(bi_normal$eta(c(3, 2), c(5, 10)), c(15, 20))
  bi_normal$x_hat
  bi_normal$hessian_bound
  bi_normal$lambda
  bi_normal$inv_lambda
  bi_normal$grad_ll(beta, 1)
  bi_normal$lap_ll(beta, 1)
  bi_normal$total_lap_ll(beta)
  bi_normal$total_grad_ll(beta)
  bi_normal$phi(beta)
  bi_normal$phi_estimator(beta, 1, 2)
  bi_normal$phi_estimator_bounds(beta)


}
)

phi_plot_3d <- function() {
  library(plotly)
  list2env(bivariate_normal_n20_b5neg2, rlang::current_env())
  # n = n, d = d, x = x, y = y, beta_true = beta_true
  dist_data <- Normal$new(y, x)

  lim_max <- 15
  lim_length <- 100

  mu <- seq(-lim_max, lim_max, length.out = lim_length)
  sd <- seq(-lim_max, lim_max, length.out = lim_length)

  dist_data$phi_estimator(dist_data$x_hat, 1, 1)

  z_phi <- matrix(0, nrow=lim_length, ncol=lim_length)
  for (x in seq_along(mu)) {
    for (y in seq_along(sd)) {
      z_phi[x, y] <- dist_data$phi(dist_data$x_unscale(c(mu[x], sd[y])))
    }
  }

  z_ub <- matrix(0, nrow=lim_length, ncol=lim_length)
  z_lb <- matrix(0, nrow=lim_length, ncol=lim_length)
  for (x in seq_along(mu)) {
    for (y in seq_along(sd)) {
      phi_info <- dist_data$phi_estimator_bounds(abs(dist_data$lambda %*% c(mu[x], sd[y])))
      z_ub[x, y] <- phi_info$phi_u
      z_lb[x, y] <- phi_info$phi_l
    }
  }

  fig <- plot_ly(
    x = ~mu,
    y = ~sd,
  ) %>%
    add_surface(z = ~z_phi) %>%
    add_surface(z = ~z_ub) %>%
    add_surface(z = ~z_lb)

  fig
}
