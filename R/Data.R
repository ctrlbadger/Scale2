# Base Class for defining posterior distribution info and providing methods for subsampling the killing rate
Data <- R6::R6Class("Data",
                public = list(
                  n = 1L, # Data Size
                  d = 1L, # Dimension
                  lambda = NULL, # Preconditioning Matrix
                  inv_lambda = NULL, # Inverse of Preconditioning Matrix
                  hessian_bound = NULL, # Hessian Bound used in obtaining phi lower and upper bounds
                  x_hat = NULL, # Centering Control Variate
                  x_hat_total_grad_ll = NULL,
                  x_hat_total_lap_ll = NULL,
                  x_hat_C = NULL,
                  x_hat_grad_ll_i = NULL,
                  x_hat_lap_ll_i = NULL,
                  x_hat_norm_2grad_log = NULL,
                  initialize = function(n, d, x_hat = NULL, hessian_bound = NULL) {
                    self$n <- n
                    self$d <- d


                    if (is.numeric(x_hat) && is.numeric(hessian_bound)) {
                      self$x_hat <- x_hat
                      self$hessian_bound <- hessian_bound
                    }
                    self$precalc_control_var(self$x_hat)
                  },
                  grad_ll = function(x, idx) { }, # Grad log likelihood of posterior
                  lap_ll = function(x, idx) { }, # Laplacian log likelihood of posterior
                  total_grad_ll = function(x) { # Grad log likelihood over the full posterior returns an R^d vector
                    Reduce("+", map(0:self$n, ~ self$grad_ll(x, .x)))
                  },
                  total_lap_ll = function(x) { # Grad^2 log likelihood over the full posterior returns an R^d vector
                    Reduce("+", map(0:self$n, ~ self$lap_ll(x, .x)))
                  },
                  phi = function(x) { # Phi killing rate function over full posterior distribution
                    ## ||grad(log(pi(x))) || ^ 2 equiv to sum of squared first derivatives of log(pi(x))
                    grad_ll <- self$total_grad_ll(x)
                    sqrd_norm_grad_ll <- sum(grad_ll^2)

                    ## lap(log(pi(x))) equiv to sum of second derivatives of log(pi(x))
                    sum_lap_ll <- self$total_lap_ll(x) %>% sum()

                    0.5 * (sqrd_norm_grad_ll + sum_lap_ll)
                  },
                  phi_control_var = function(x) { # Alternative phi function using x_hat
                    t(self$alpha(x)) %*% as.matrix(2*self$total_grad_ll(self$x_hat) + self$alpha(x)) / 2 +
                      self$div_alpha(x) / 2 + self$x_hat_C |> as.double()
                  },
                  phi_estimator = function(x, i, j) { # Phi estimator for a subsampling mechanism of size two
                    a_i <- self$alpha_subsample(x, i)
                    a_j <- self$alpha_subsample(x, j)

                    # vector_step <- sum(a_i * (2 * self$total_grad_ll + a_j))

                    inner_prod_term <- 0.5 * t(a_i) %*% as.matrix(2 * self$x_hat_total_grad_ll + a_j) %>% as.double()

                    0.5 * (inner_prod_term + self$div_alpha_subsample(x, i)) + self$x_hat_C
                  },
                  alpha = function(x) { # Helper function difference between grad log likelihoods of centering value and x
                    self$total_grad_ll(x) - self$total_grad_ll(self$x_hat)
                  },
                  div_alpha = function(x) { # Helper function difference between laplacian log likelihoods of centering value and x
                    sum(self$total_lap_ll(x) - self$total_lap_ll(self$x_hat))
                  },
                  alpha_subsample = function(x, i) { # subsampling alpha at a single datum i
                    # Whether prior is included in calculating estimator
                    (self$n + 1) * (self$grad_ll(x, i) - self$x_hat_grad_ll_i[[i + 1]])
                  },
                  div_alpha_subsample = function(x, i) { # Subsampling div{alpha} at a single datum i
                    sum((self$n + 1) * (self$lap_ll(x, i) - self$x_hat_lap_ll_i[[i + 1]]))
                  },
                  precalc_control_var = function(x_hat = NULL) { # Precalculates control variate constants used in subsampling phi
                    if (is.numeric(x_hat)) self$x_hat <- x_hat

                    # Precalculates x_hat at each respective datum
                    # NOTE: For large n this will be memory intensive and may need to be specified whether to precalculate
                    self$x_hat_grad_ll_i <- map(0:self$n, ~ self$grad_ll(self$x_hat, .x))
                    self$x_hat_lap_ll_i <- map(0:self$n, ~ self$lap_ll(self$x_hat, .x))

                    ## ||grad(log(pi(x))) || ^ 2 equiv to sum of squared first derivatives of log(pi(x))
                    self$x_hat_total_grad_ll <- self$total_grad_ll(self$x_hat)
                    self$x_hat_total_lap_ll <- self$total_lap_ll(self$x_hat)

                    self$x_hat_norm_2grad_log <- sum((2*self$x_hat_total_grad_ll)^2)^(1/2)

                    sqrd_norm_grad_ll <- sum(self$x_hat_total_grad_ll^2)

                    ## lap(log(pi(x))) equiv to sum of second derivatives of log(pi(x))

                    sum_total_lap_ll <- sum(self$x_hat_total_lap_ll)

                    self$x_hat_C <- 0.5 * (sqrd_norm_grad_ll + sum_total_lap_ll)
                    invisible(self)
                  },
                  phi_estimator_bounds = function(maximal_distance) { # Calculate phi bounds from equation 22
                    absolute_bound <- (self$n + 1) * self$hessian_bound *
                      maximal_distance * (self$x_hat_norm_2grad_log + self$hessian_bound * (self$n + 1) * maximal_distance) +
                      self$hessian_bound * self$d * (self$n + 1)


                    phi_u <- 0.5 * absolute_bound + self$x_hat_C
                    phi_l <- -0.5 * absolute_bound + self$x_hat_C
                    intensity <- phi_u - phi_l

                    list(phi_u = phi_u, phi_l = phi_l, intensity = intensity)
                  },
                  x_unscale = function(x) { # Reverse of x_xuncale, transforms by inverse of preconditioning matrix and recenters around x_hat
                    c((self$lambda %*% as.matrix(x)) + self$x_hat)
                  },
                  x_scale = function(x) { # Centers around 0 and preconditions by inv_lambda
                    c(self$inv_lambda %*% as.matrix(x - self$x_hat))
                  }
                )
)


#' R6 Cauchy Distribution
#'
#' @description Distribution Info for sampling from cauchy distribution, inherits from Data
#' Defines grad and laplacian with a flat prior
#'
#' @field mu location parameter
#' @field gamma scale parameter
#'
#' @export
CauchyData <- R6::R6Class("CauchyData", inherit = Data,
    public = list(
    mu = 0,
    gamma = 1,
    #' @description Cauchy PDF
    #'
    #' @param x Parameters of pdf
    pi = function(x) {
        exp(pi_ll(x))
    },
    #' @description Log Likelihood of Cauchy Distribution
    #'
    #' @param x Parameters of ll
    #' @param i not needed
    pi_ll = function(x) {
        -log(pi * self$gamma) - log(1 + self$cauchy_scale(x, self$mu, self$gamma)^2)
    },
    #' @description Gradient of Log Likelihood of Cauchy Distribution
    #'
    #' @param x Parameters of grad ll
    #' @param i not needed
    grad_ll = function(x, i=NA) { # Grad ll for cauchy
        if (i == 0) return(0)

        - 2 * (x - self$mu) / ((x - self$mu)^2 + self$gamma^2)
    },
    #' @description Laplacian of Log Likelihood of Cauchy Distribution
    #'
    #' @param x Parameters of lap ll
    #' @param i not needed
    lap_ll = function(x, i=NA) { # Laplacian ll for cauchy
        if (i == 0) return(0)
        t <- self$cauchy_scale(x, self$mu, self$gamma)
        2 * (t^2 - 1) / ((t^2 + 1) * ((x - self$mu)^2 + self$gamma^2))
    },
    #' @description Cauchy Constructor
    #'
    #' @param mu Location Parameter
    #' @param gamma Positive Scale Parameter
    initialize = function(mu, gamma) { # Initialising cauchy dist class
        self$mu <- mu
        self$gamma <- gamma


        self$lambda <- 1 / (4 * self$gamma^2)
        self$inv_lambda <- solve(self$lambda)

        super$initialize(n = 1, d = 1, x_hat = self$mu, hessian_bound = 1 / (4 * gamma^2))
        invisible(self)
    },
    #' @description Cauchy Scale Helper Function
    #'
    #' @param x Parameter
    #' @param mu Location Parameter
    #' @param gamma Positive Scale Parameter
    cauchy_scale = function(x, mu, gamma)  { # Scaling helper function
        (x - mu) / gamma
    },
    #' @description Exact bounds on phi
    #'
    #' @param x.l Lower Bound
    #' @param x.u Upper Bounds
    phi_bounds_exact = function(x.l, x.u) { # Exact Phi Bounds
        b1 <- self$phi(x.l)
        b2 <- self$phi(x.u)

        phiU <- max(c(b1, b2))
        phiL <- min(c(b1, b2))

        if ((x.l < self$mu) && (self$mu < x.u)) {
            phiL <- self$phi(self$mu)
        }

        phi.max <- -sqrt( 5/3)*self$gamma + self$mu

        if ((x.l < phi.max) && (phi.max < x.u)) {
            phiU <- self$phi(phi.max)
        }

        phi.max <- phi.max <- sqrt(5/3)*self$gamma + self$mu
        if ((x.l < phi.max) && (phi.max < x.u)) {
            phiU <- self$phi(phi.max)
        }

        return(list(phi_u = phiU, phi_l = phiL, intensity = phiU - phiL))
    }
    )
)

# Only mu normal data example
MeanNormalData <- R6::R6Class("MeanNormalData", inherit = Data,
    public = list(
    data_x = 0, # Vector of means for each normal datum i
    sigma_est = 0, # Estimated standard deviation
    constant_pi_observed = 1, # Normalising Constants for log posterior
    constant_pi_actual = 1,
    mu_true = NULL, # Optional true mean
    sigma_true = NULL, # Optional true sd
    log_pi_actual = function(x) {
        log(self$pi_actual(x))
    },
    log_pi_observed = function(x) {
        log(self$pi_observed(x))
    },
    pi_actual = function(x) {
        exp(self$n * dnorm(x, mean = self$mu_true, sd=self$sigma_true, log = TRUE)) / self$constant_pi_actual
    },
    pi_observed = function(x) {
        exp(sum(map_dbl(self$data_x, ~ dnorm(.x, mean = x, sd = self$sigma_est, log = TRUE))) %>% sum()) / self$constant_pi_observed
    },
    grad_ll = function(x, i) {
        if (i == 0) return(0)

        mu_comp <- (self$data_x[i] - x[1]) / self$sigma_est^2
        sigma_comp <- 1 / x[2] + (self$data_x[i] - x[1])^2 * self$sigma_est^(-3)

        return(c(mu_comp))
    },
    lap_ll = function(x, i) {
        if (i == 0) return(0)

        mu_comp <- - self$sigma_est^-2
        sigma_comp <- -self$sigma_est^-2 - 3 * (self$data_x[i] - self$sigma_est)^2 * x[2]^-4

        return(c(mu_comp))
    },
    fisher_information = function() { # Fisher information
        diag(c(self$sigma_est^-2, 4 * self$sigma_est^-2))
    },
    hessian_ll_i = function(i) { # Hessian log likelihood for datum i
        diag(c(-self$sigma_est^-2, -self$sigma_est^-2 - 3 * (self$data_x[i] - self$x_hat)^2*self$sigma_est^(-4)))
    },
    get_hessian_bound = function() { # Find hessian bound
        get_radial_bound <- function(i) { # Function for finding maximum eigenvalue for the ith hessian log likelihood
            (eigen(self$hessian_ll_i(i))$values %>% abs() %>% max())
        }

        # Upper bound based on n * fisher information
        hessian_bound <- length(self$data_x) * (eigen(self$fisher_information())$values %>% abs() %>% max())
        # Bound based on direct computation of n hessians at each datum
        hessian_bound <- map_dbl(seq_along(self$data_x), ~ get_radial_bound(.x)) %>% max()

        return(hessian_bound)
    },
    print_debug = function() {
      paste("Hessian Bound is: ", self$hessian_bound) %>% print()
      paste("Lambda & Inv_lambda", self$lambda, self$inv_lambda) %>% print()
      paste("Mean & SD:", mu_est, sigma_est) %>% print()
    },
    initialize = function(data_x, mu_true = NULL, sigma_true = NULL) {
        self$data_x <- data_x

        # MAP estimates and MLE Estimates
        mu_est <- mean(data_x)
        sigma_est <- sd(data_x)

        self$x_hat <- mu_est

        self$mu_true <- mu_true
        self$sigma_true <- sigma_true


        self$n <- length(data_x)
        self$sigma_est <- sigma_est

        self$inv_lambda <- as.matrix(c(sqrt(self$n * self$sigma_est^-2)))
        self$lambda <- solve(self$inv_lambda)



        # Normalising Constants
        self$constant_pi_observed <- integrate(function(x) {map_dbl(x, ~ self$pi_observed(.x))}, lower = -Inf, upper = Inf)$value
        self$constant_pi_actual <- integrate(function(x) {map_dbl(x, ~ self$pi_actual(.x))}, lower = -Inf, upper = Inf)$value

        self$hessian_bound <- self$get_hessian_bound()

        super$initialize(n = length(data_x), d = 1, x_hat = mu_est, hessian_bound = self$hessian_bound)

        invisible(self)
    }
    )
)
