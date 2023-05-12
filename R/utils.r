trap_mean_sd <- function(x, y) {
  # Compute the mean
  mean <- sum((x[-1] - x[-length(x)]) * (y[-1] + y[-length(y)]) / 2) / sum(y)
  
  # Compute the variance using the trapezoidal rule
  variance <- sum((x[-1] - x[-length(x)]) * ((y[-1] * (x[-1] - mean)^2) + (y[-length(y)] * (x[-length(x)] - mean)^2)) / 2) / sum(y)
  
  # Return the standard deviation
  sd <- sqrt(variance)
  return(list(mean=mean, sd=sd))
}

vv_transformer <- function(text, envir) {
  regex <- "=$"
  if (!grepl(regex, text)) {
    return(glue::identity_transformer(text, envir))
  }
  
  text <- sub(regex, "", text)
  res <- glue::identity_transformer(text, envir)
  n <- length(res)
  res <- glue::glue_collapse(res, sep = ", ")
  if (n > 1) {
    res <- c("[", res, "]")
  }
  glue::glue_collapse(c(text, " = ", res))
}


get_aliases <- function(treatment_tibble, generator) {
  library(purrr)
  treatment_tibble %>% 
    as.list() %>% 
    list_transpose() %>% # Operate Row_Wise
    map_dfr( ~ dot_operation(.x, generator))
}

dot_operation <- function(treat, gen) {
  (treat + gen) %% 2
}
