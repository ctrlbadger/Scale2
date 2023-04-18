library(Scale2)


create_cauchy_data <- function() {

    # Centered Cauchy Example
    mu <- 0
    gamma <- 1
    dist_data <- CauchyData$new(mu, gamma)

    cauchy_centered <- list(
        gamma = gamma, mu = mu, dist_data = dist_data
    )

    use_data(cauchy_centered)
    

    # Deviated Cauchy Example
    mu <- -2
    gamma <- 5
    dist_data <- CauchyData$new(mu, gamma)

    cauchy_mu_neg2_gamma_5 <- list(
        gamma = gamma, mu = mu, dist_data = dist_data
    )

    use_data(cauchy_mu_neg2_gamma_5, overwrite = TRUE)
}

