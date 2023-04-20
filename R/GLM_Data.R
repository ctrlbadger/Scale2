
GLM <- R6::R6Class("GLM", inherit = Data,
    public = list(
    dispersion = NULL,
    offset = NULL,
    prior_weights = NULL,
    data_y = NULL,
    data_x = NULL,
    log_pi_observed = function(beta) {
      sum(map_dbl(1:self$n, ~ self$log_pi_observed_i(beta, .x)))
    },
    unit_deviance = function(mu, y) {

    },
    v = function(mu) { # Variance Function
    },
    log_pi_observed_i = function(beta, i) {
      y <- self$data_y[i]
      x <- self$data_x[i, ]

      eta_i <- self$eta(beta, y)
      mu <- self$logistic(eta_i)

      -0.5*self$dispersion * self$unit_deviance(mu, y)
    },
    link = function(mu) { # g(mu) = eta

    },
    grad_link = function(mu) {
    },
    grad2_link = function(mu) {
    },
    logistic = function(eta) { # mu = g^-1(eta)
    },
    eta = function(beta, x) {
      (beta * x)
    },
    pi_actual = function(x) {

    },
    pi_observed = function(x) {
    },
    grad_ll = function(beta, i) {
      if (i == 0) return(0)
      y <- self$data_y[i]
      x <- self$data_x[i, ]
      eta_i <- self$eta(beta, x)

      mu <- self$logistic(eta_i)
      (y - mu) * (self$prior_weights[i] * x) / (self$dispersion * self$v(mu) * self$grad_link(mu))
    },
    lap_ll = function(beta, i) {
      if (i == 0) return(0)
      y <- self$data_y[i]
      x <- self$data_x[i, ]
      eta_i <- self$eta(beta, x)
      mu <- self$logistic(eta_i)
      f1 <- -self$prior_weights[i]  / (self$dispersion * self$v(mu) * self$grad_link(mu)^2) * x * x
      # f2 <- (y - mu) * self$prior_weights[i] * x^2 / (self$dispersion * self$v(mu) * self$grad2_link(mu))
      f2 <- 0
      f1 + f2
    },
    fisher_information = function(i) { # Fisher information
      if (i == 0) return(0)
      y <- self$data_y[i]
      x <- self$data_x[i, ]
      eta_i <- self$eta(self$x_hat, x)
      mu <- self$logistic(eta_i)

      self$prior_weights[i] * x^2 / (self$dispersion * self$v(mu) * self$grad_link(mu)^2)
    },
    get_hessian_bound = function() { # Find hessian bound
        get_radial_bound <- function(i) { # Function for finding maximum eigenvalue for the ith hessian log likelihood
            (eigen(diag(self$fisher_information(i), self$d))$values %>% abs() %>% max())
        }

        # Bound based on direct computation of n hessians at each datum
        hessian_bound <- map_dbl(1:self$n, ~ get_radial_bound(.x)) %>% max()

        hessian_bound
    },
    get_dispersion_param = function() {
      temp <- 0
      for (i in 1:self$n) {
        y <- self$data_y[i]
        x <- self$data_x[i, ]
        eta_i <- self$eta(self$x_hat, x)
        mu <- self$logistic(eta_i)
        temp <- temp + self$prior_weights[i] * (y - mu)^2 / self$v(mu)
      }

      temp / (self$n - self$d)
    },
    get_preconditioning_matrix = function() {
      fi <- Reduce('+', map(1:self$n, ~ diag(self$fisher_information(.x), self$d))) %>% sqrt()
      self$inv_lambda <- fi
      self$lambda <- solve(fi)

      list(lambda = self$lambda, inv_lambda = self$inv_lambda)
    },
    initialize = function(n, d) {
        if (is.null(self$dispersion)) self$dispersion <- self$get_dispersion_param()
        if (is.null(self$offset)) self$offset <- rep(0, n)
        if (is.null(self$prior_weights)) self$prior_weights = rep(1, n)
        if (is.null(self$hessian_bound)) self$hessian_bound <- self$get_hessian_bound()
        if (is.null(self$lambda)) self$get_preconditioning_matrix()
        super$initialize(n, d)
    }
    )
)


Normal <- R6::R6Class("Binomial", inherit = GLM,
    public = list(
    dispersion = 1,
    offset = 0,
    weights = NULL,
    unit_deviance = function(mu, y) {
      (y - mu)^2
    },
    v = function(mu) { # Variance Function
      1
    },
    link = function(mu) { # g(mu) = eta
      mu
    },
    grad_link = function(mu) {
      1
    },
    grad2_link = function(mu) {
      0
    },
    logistic = function(eta) { # mu = g^-1(eta)
      eta
    },
    pi_actual = function(x) {

    },
    pi_observed = function(x) {
    },
    initialize = function(y, x) {
      self$data_y <- y
      self$data_x <- as.matrix(x)



      self$n <- dim(self$data_x)[1]
      self$d <- dim(self$data_x)[2]

      if (self$d == 1) self$dispersion <- var(y)

      beta_hat <- coef(lm(y ~ x - 1))[1:self$d]
      names(beta_hat) <- NULL
      self$x_hat <- beta_hat


      super$initialize(self$n, self$d)
      invisible(self)
    }
  )
)

