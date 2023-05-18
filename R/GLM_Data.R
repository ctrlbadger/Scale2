GLM <- R6::R6Class("GLM", inherit = Data,
    public = list(
    dispersion = NULL,
    offset = NULL,
    prior_weights = NULL,
    data_y = NULL,
    data_x = NULL,
    constant_marginal = NULL,
    log_pi_observed = function(beta) {
      sum(map_dbl(1:self$n, ~ self$log_pi_observed_i(beta, .x)))
    },
    log_pi_observed_i = function(beta, i) {
      y <- self$data_y[i]
      x <- self$data_x[i, ]

      eta_i <- self$eta(beta, y)
      mu <- self$logistic(eta_i)

      -0.5*self$dispersion * self$unit_deviance(mu, y)
    },
    eta = function(beta, x) {
      sum(beta * x)
    },
    grad_ll = function(beta, i) {
      if (i == 0) return(rep(0, self$d))
      y <- self$data_y[i]
      x <- self$data_x[i, ]
      eta_i <- self$eta(beta, x)

      mu <- self$logistic(eta_i)
      (y - mu) * (self$prior_weights[i] * x) / (self$dispersion * self$v(mu) * self$grad_link(mu))
    },
    lap_ll = function(beta, i) {
      if (i == 0) return(rep(0, self$d))
      y <- self$data_y[i]
      x <- self$data_x[i, ]
      eta_i <- self$eta(beta, x)
      mu <- self$logistic(eta_i)
      f1 <- -self$prior_weights[i]  / (self$dispersion * self$v(mu) * self$grad_link(mu)^2) * x * x
      f2 <- (y - mu) * self$prior_weights[i] * x^2 / (self$dispersion * self$v(mu) * self$grad2_link(mu))
      # f2 <- 0
      f1 + f2
    },
    fisher_information = function(i) { # Fisher information
      if (i == 0) return(rep(0, self$d))
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
        if (is.null(self$prior_weights)) self$prior_weights <- rep(1, n)
        if (is.null(self$hessian_bound)) self$hessian_bound <- self$get_hessian_bound()
        if (is.null(self$lambda)) self$get_preconditioning_matrix()
        super$initialize(n, d)
    }
    ),
    private = list(
      link = function(mu) { # g(mu) = eta

      },
      grad_link = function(mu) {
      },
      grad2_link = function(mu) {
      },
      logistic = function(eta) { # mu = g^-1(eta)
      },
      unit_deviance = function(mu, y) {

      },
      v = function(mu) { # Variance Function
      },
      pi_actual = function(x) {

      },
      pi_observed = function(x) {
      }
    )
)

#' Normal GLM Model R6 Class
#' @description A normal linear model class for specifying the distribution to the scale algorithm
#' @examples
#' list2env(bivariate_normal_n20_b5neg2, rlang::current_env())
#' bi_normal <- Normal$new(y, x)
#' @export
Normal <- R6::R6Class("Normal", inherit = GLM,
    public = list(
    #' @description The unit deviance of a Normal model
    #'
    #' @param mu Mean value mu
    #' @param y Dependent variable y
    unit_deviance = function(mu, y) {
      (y - mu)^2
    },
    #' @description Variance Function
    #'
    #' @param mu Mean value mu
    v = function(mu) { # Variance Function
      1
    },
    #' @description Link Function g(mu) = eta
    #'
    #' @param mu Mean value mu
    link = function(mu) { # g(mu) = eta
      mu
    },
    #' @description Derivative of Link Function g(mu) = eta
    #'
    #' @param mu Mean value mu
    grad_link = function(mu) {
      1
    },
    #' @description Second derivative of Link Function
    #'
    #' @param mu Mean value mu
    grad2_link = function(mu) {
      0
    },
    #' @description Logistic Function - Inverse of Link Function
    #'
    #' @param eta Eta. Inverse of value mu
    logistic = function(eta) { # mu = g^-1(eta)
      eta
    },
    #' @description Laplacian of log likelihood at datum i
    #'
    #' @param beta - Unkown parameter $\beta$
    #' @param i - Index of data matrix
    lap_ll = function(beta, i) {
      if (i == 0) return(rep(0, self$d))
      y <- self$data_y[i]
      x <- self$data_x[i, ]
      eta_i <- self$eta(beta, x)
      mu <- self$logistic(eta_i)
      f1 <- -self$prior_weights[i]  / (self$dispersion * self$v(mu) * self$grad_link(mu)^2) * x * x
      f1
    },
    #' @description Initialise Normal Data Model
    #'
    #' @param y Numeric vector length n. Dependent variables y
    #' @param x Matrix of dimension n x d.  Response variables x
    initialize = function(y, x) {
      self$data_y <- y
      self$data_x <- as.matrix(x)



      self$n <- dim(self$data_x)[1]
      self$d <- dim(self$data_x)[2]


      beta_hat <- coef(lm(y ~ x - 1))[1:self$d]
      self$dispersion <- summary(lm(y ~ x - 1))$sigma^2
      names(beta_hat) <- NULL
      self$x_hat <- beta_hat


      super$initialize(self$n, self$d)
      invisible(self)
    }
  ),
  private = list(
    dispersion = 1,
    pi_actual = function(x) {

    },
    pi_observed = function(x) {
    }
  )
)

#' Binomial R6 Class for Scale
#' @description A logistic regression model using a GLM framework for the scale model.
#' @field dispersion The dispersion perameter of the Binomial Model
#' @examples
#' list2env(small_logistic_example, rlang::current_env())
#' model_glm <- glm(y ~ x - 1, family = binomial)
#' binomial_data <- Binomial$new(y = y, x = x)
#' @export
Binomial <- R6::R6Class("Binomial", inherit = GLM,
    public = list(
    dispersion = 1,

    #' @description The unit deviance of a binomial model
    #'
    #' @param mu Mean value mu
    #' @param y Dependent variable y
    unit_deviance = function(mu, y) {
      2*(log(y/mu) + (1-y)*log((1-y)/(1-mu)))
    },
    #' @description Variance Function
    #'
    #' @param mu Mean value mu
    v = function(mu) { # Variance Function
      mu * (1 - mu)
    },
    #' @description Link Function g(mu) = eta
    #'
    #' @param mu Mean value mu
    link = function(mu) { # g(mu) = eta
      log(mu / (1 - mu))
    },
    #' @description Derivative of Link Function g(mu) = eta
    #'
    #' @param mu Mean value mu
    grad_link = function(mu) {
      1 / (mu - mu^2)
    },
    #' @description Second derivative of Link Function
    #'
    #' @param mu Mean value mu
    grad2_link = function(mu) {
      (2*mu - 1) / ((mu - 1)^2 * mu^2)
    },
    #' @description Logistic Function - Inverse of Link Function
    #'
    #' @param eta Eta. Inverse of value mu
    logistic = function(eta) { # mu = g^-1(eta)
      exp(eta) / (1 + exp(eta))
    },
    #' @description Laplacian of log likelihood at datum i
    #'
    #' @param beta - Unkown parameter $\beta$
    #' @param i - Index of data matrix
    lap_ll = function(beta, i) {
      if (i == 0) return(rep(0, self$d))
      y <- self$data_y[i]
      x <- self$data_x[i, ]
      eta_i <- self$eta(beta, x)
      mu <- self$logistic(eta_i)
      if (mu == 1 || mu == 0) return(rep(0, self$d))
      else f1 <- -self$prior_weights[i]  / (self$dispersion * self$v(mu) * self$grad_link(mu)^2) * x * x
      # f2 <- (y - mu) * self$prior_weights[i] / (self$dispersion * self$v(mu) * self$grad2_link(mu)) * x * x
      f2 <- 0
      f1 + f2
    },
    #' @description Gradient of log likelihood at datum i
    #'
    #' @param beta - Unkown parameter $\beta$
    #' @param i - Index of data matrix
    grad_ll = function(beta, i) {
      if (i == 0) return(rep(0, self$d))
      y <- self$data_y[i]
      x <- self$data_x[i, ]
      eta_i <- self$eta(beta, x)

      mu <- self$logistic(eta_i)
      if (mu == 1 || mu == 0) return(rep(0, self$d))
      else (y - mu) * (self$prior_weights[i] * x) / (self$dispersion * self$v(mu) * self$grad_link(mu))
    },
    #' @description Estimated Fisher Information for datum i
    #'
    #' @param i - Index of data matrix
    fisher_information = function(i) { # Fisher information
      if (i == 0) return(rep(0, self$d))
      y <- self$data_y[i]
      x <- self$data_x[i, ]
      eta_i <- self$eta(self$x_hat, x)
      mu <- self$logistic(eta_i)
      if (mu == 1 || mu == 0) return(rep(0, self$d))
      else self$prior_weights[i] * x^2 / (self$dispersion * self$v(mu) * self$grad_link(mu)^2)
    },
    #' @description Initialise Binomial Data Model
    #'
    #' @param y Numeric vector length n. Dependent variables y
    #' @param x Matrix of dimension n x d.  Response variables x
    initialize = function(y, x) {
      self$data_y <- y
      self$data_x <- as.matrix(x)

      self$n <- dim(self$data_x)[1]
      self$d <- dim(self$data_x)[2]

      beta_hat <- coef(glm(self$data_y ~self$data_x - 1, family = binomial))
      names(beta_hat) <- NULL
      self$x_hat <- beta_hat

      super$initialize(self$n, self$d)
      invisible(self)
    }
  ),
  private = list(
    pi_actual = function(x) {

    },
    pi_observed = function(x) {
    }
  )
)

