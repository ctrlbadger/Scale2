# Gets subsampling functions, phi and lower and upper bounds for given hypercube
get_phi_info <- function(data, path_l, path_u, rescale = TRUE, subsample = FALSE) {
  if (c("CauchyData" %in% class(data))) {
    if (rescale) {
      path_l <- data$x_unscale(path_l)
      path_u <- data$x_unscale(path_u)
    }
    phi_bounds_info <- data$phi_bounds_exact(path_l, path_u)

  } else {
    maximal_distance <- sum(pmax(abs(path_l), abs(path_u))^2)^(1/2)

    if (rescale) maximal_distance <- c(data$lambda %*% as.matrix(maximal_distance))

    phi_bounds_info <- data$phi_estimator_bounds(maximal_distance)
  }

  # Return maximal distance away from origin of hypercube, used in calculating phi bounds
  phi_u <- phi_bounds_info$phi_u
  phi_l <- phi_bounds_info$phi_l
  intensity <- phi_bounds_info$intensity

  if (rescale) x_transf <- function(x) data$x_unscale(x)
  else x_transf <- function(x) x

  A = NULL
  if (subsample) {
    # Get Sample Points
    A <- sample(0:data$n, size = 2, replace = TRUE)

    phi_func <- function(x) {
      force(A)

      data$phi_estimator(x_transf(x), A[1], A[2])
    }
  } else {
    phi_func <- function(x) data$phi(x_transf(x))
  }

  return(list(phi_estimator_hyp = phi_func, phi_u = phi_u, phi_l = phi_l, intensity = intensity))
}

