library(tidyverse)

clamp <- function(x, min = 0, max = 1) {
    x[x < min] <- min
    x[x > max] <- max
    x
}

create_data <- function(seed, n_set) {
    set.seed(seed)
    n_set <- 5
    seeds <- round(runif(n_set, 0, 1000))
    n_obs <- sample(5:10, n_set)
    offset <- runif(n_set, -7, 7)
    sdev   <- runif(n_set, 0.5, 6)
    res <- bind_rows(lapply(seq_len(n_set), function(i) {
        tibble(label = as.character(seeds[i]),
               time = sort(runif(n_obs[i], -4, 4)),
               obs = clamp(pnorm(time, sd = sdev[i]) + rnorm(length(time), sd = 0.05))) %>%
        mutate(time = time + offset[i],
               time = time - min(time)) # To start with t = 0
    }))
    as.data.frame(res)
}

dummy_data <- create_data(12345076, 5)
devtools::use_data(dummy_data)
