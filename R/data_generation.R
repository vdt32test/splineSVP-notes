# ============================================================
# SplineSVP: change-in-slope simulation scenarios
# ============================================================


# ============================================================
# 1. Generate a continuous piecewise-linear signal
# ============================================================

generate_signal <- function(n,
                            changepoints = integer(0),
                            slopes,
                            intercept = 0) {
  
  if (n < 1) {
    stop("n must be positive.")
  }
  
  if (length(slopes) != length(changepoints) + 1) {
    stop("There must be one more slope than changepoints.")
  }
  
  if (length(changepoints) > 0) {
    
    if (any(changepoints <= 1) || any(changepoints >= n)) {
      stop("Changepoints must lie strictly inside the observation range.")
    }
    
    if (is.unsorted(changepoints, strictly = TRUE)) {
      stop("Changepoints must be strictly increasing.")
    }
  }
  
  t <- seq_len(n)
  
  # First affine piece
  signal <- intercept + slopes[1] * (t - 1)
  
  # Add changes in slope using (t-r)_+.
  # This preserves continuity at every changepoint.
  if (length(changepoints) > 0) {
    
    for (j in seq_along(changepoints)) {
      
      r <- changepoints[j]
      
      slope_change <- slopes[j + 1] - slopes[j]
      
      signal <- signal +
        slope_change * pmax(t - r, 0)
    }
  }
  
  signal
}


# ============================================================
# 2. Add Gaussian noise
# ============================================================

generate_data <- function(signal,
                          sigma = 1) {
  
  if (sigma < 0) {
    stop("sigma must be non-negative.")
  }
  
  n <- length(signal)
  
  noise <- rnorm(
    n = n,
    mean = 0,
    sd = sigma
  )
  
  y <- signal + noise
  
  data.frame(
    t = seq_len(n),
    signal = signal,
    noise = noise,
    y = y
  )
}


# ============================================================
# 3. Equally spaced segment lengths
# ============================================================

equal_segment_lengths <- function(n, K) {
  
  if (K < 1 || K > n) {
    stop("K must satisfy 1 <= K <= n.")
  }
  
  base_length <- n %/% K
  remainder <- n %% K
  
  lengths <- rep(base_length, K)
  
  if (remainder > 0) {
    lengths[seq_len(remainder)] <-
      lengths[seq_len(remainder)] + 1
  }
  
  lengths
}


# ============================================================
# 4. Random segment lengths using a Dirichlet distribution
# ============================================================

# A symmetric Dirichlet distribution is generated using independent
# Gamma random variables and normalization.
#
# dirichlet_alpha controls the variability:
#   alpha = 1  : uniform distribution over segment proportions
#   alpha > 1  : segment lengths tend to be more similar
#   alpha < 1  : more unequal segment lengths
#
# min_segment_length guarantees that every segment contains at least
# the specified number of observations.

random_segment_lengths_dirichlet <- function(
    n,
    K,
    min_segment_length = 2,
    dirichlet_alpha = 1) {
  
  if (K < 1) {
    stop("K must be positive.")
  }
  
  if (min_segment_length < 1) {
    stop("min_segment_length must be at least 1.")
  }
  
  if (dirichlet_alpha <= 0) {
    stop("dirichlet_alpha must be positive.")
  }
  
  minimum_total <- K * min_segment_length
  
  if (minimum_total > n) {
    stop(
      "n is too small for K segments with the requested minimum length."
    )
  }
  
  remaining <- n - minimum_total
  
  # If there is no remaining length, all segments have the
  # minimum length.
  if (remaining == 0) {
    return(rep(min_segment_length, K))
  }
  
  # Draw Dirichlet proportions using Gamma variables.
  g <- rgamma(
    K,
    shape = dirichlet_alpha,
    rate = 1
  )
  
  proportions <- g / sum(g)
  
  # Convert the continuous proportions into integer lengths.
  raw_extra <- remaining * proportions
  
  extra <- floor(raw_extra)
  
  # floor() may leave a few observations unassigned.
  left <- remaining - sum(extra)
  
  if (left > 0) {
    
    fractional_parts <- raw_extra - extra
    
    order_fractional <- order(
      fractional_parts,
      decreasing = TRUE
    )
    
    extra[
      order_fractional[seq_len(left)]
    ] <- extra[
      order_fractional[seq_len(left)]
    ] + 1
  }
  
  segment_lengths <-
    min_segment_length + extra
  
  stopifnot(sum(segment_lengths) == n)
  
  segment_lengths
}


# ============================================================
# 5. Convert segment lengths to changepoints
# ============================================================

lengths_to_changepoints <- function(segment_lengths) {
  
  K <- length(segment_lengths)
  
  if (K == 1) {
    return(integer(0))
  }
  
  cumsum(segment_lengths)[seq_len(K - 1)]
}


# ============================================================
# 6. Generate one simulation scenario
# ============================================================

generate_scenario <- function(
    scenario,
    n = 200,
    K = 4,
    sigma = 1,
    intercept = 0,
    seed = NULL,
    min_segment_length = 2,
    dirichlet_alpha = 1) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # ----------------------------------------------------------
  # Scenario 1: no change in slope
  # ----------------------------------------------------------
  
  if (scenario == "no_change") {
    
    K_actual <- 1
    
    segment_lengths <- n
    
    changepoints <- integer(0)
    
    slopes <- 0.10
    
    
    # ----------------------------------------------------------
    # Scenario 2: one change in slope
    # ----------------------------------------------------------
    
  } else if (scenario == "single_change") {
    
    K_actual <- 2
    
    segment_lengths <- equal_segment_lengths(
      n = n,
      K = K_actual
    )
    
    changepoints <- lengths_to_changepoints(
      segment_lengths
    )
    
    slopes <- c(
      0.05,
      0.25
    )
    
    
    # ----------------------------------------------------------
    # Scenario 3: several changes in the same direction
    # ----------------------------------------------------------
    
  } else if (scenario == "same_direction") {
    
    if (K < 2) {
      stop("For same_direction, K must be at least 2.")
    }
    
    K_actual <- K
    
    # Equal segment lengths for this controlled scenario.
    segment_lengths <- equal_segment_lengths(
      n = n,
      K = K_actual
    )
    
    changepoints <- lengths_to_changepoints(
      segment_lengths
    )
    
    # All slopes are positive.
    # Their magnitudes vary across the K segments.
    slopes <- seq(
      from = 0.25,
      to = 1,
      length.out = K_actual
    )
    
    
    # ----------------------------------------------------------
    # Scenario 4: alternating +1 / -1 slopes
    # ----------------------------------------------------------
    
  } else if (scenario == "alternating") {
    
    if (K < 2) {
      stop("For alternating, K must be at least 2.")
    }
    
    K_actual <- K
    
    # Equal segment lengths make the alternating pattern
    # particularly easy to interpret.
    segment_lengths <- equal_segment_lengths(
      n = n,
      K = K_actual
    )
    
    changepoints <- lengths_to_changepoints(
      segment_lengths
    )
    
    # +1, -1, +1, -1, ...
    slopes <- rep(
      c(1, -1),
      length.out = K_actual
    )
    
    
    # ----------------------------------------------------------
    # Scenario 5: random slopes and random segment lengths
    # ----------------------------------------------------------
    
  } else if (scenario == "random") {
    
    if (K < 2) {
      stop("For random, K must be at least 2.")
    }
    
    K_actual <- K
    
    # Random segment lengths generated from Dirichlet
    # proportions.
    segment_lengths <- random_segment_lengths_dirichlet(
      n = n,
      K = K_actual,
      min_segment_length = min_segment_length,
      dirichlet_alpha = dirichlet_alpha
    )
    
    changepoints <- lengths_to_changepoints(
      segment_lengths
    )
    
    # Random slopes for the K segments.
    slopes <- runif(
      n = K_actual,
      min = -1,
      max = 1
    )
    
    
  } else {
    
    stop(
      paste(
        "Unknown scenario:",
        scenario
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Construct the continuous signal
  # ----------------------------------------------------------
  
  signal <- generate_signal(
    n = n,
    changepoints = changepoints,
    slopes = slopes,
    intercept = intercept
  )
  
  
  # ----------------------------------------------------------
  # Add Gaussian noise
  # ----------------------------------------------------------
  
  dat <- generate_data(
    signal = signal,
    sigma = sigma
  )
  
  
  # ----------------------------------------------------------
  # Return simulation information
  # ----------------------------------------------------------
  
  list(
    scenario = scenario,
    n = n,
    K = K_actual,
    sigma = sigma,
    intercept = intercept,
    segment_lengths = segment_lengths,
    changepoints = changepoints,
    slopes = slopes,
    data = dat
  )
}


# ============================================================
# 7. Plot one simulation
# ============================================================

plot_scenario <- function(sim,
                          main = NULL) {
  
  if (is.null(main)) {
    
    main <- paste(
      "Scenario:",
      sim$scenario
    )
  }
  
  plot(
    sim$data$t,
    sim$data$y,
    pch = 16,
    cex = 0.55,
    xlab = "t",
    ylab = "y",
    main = main
  )
  
  # True underlying signal
  lines(
    sim$data$t,
    sim$data$signal,
    lwd = 2
  )
  
  # True changepoints
  if (length(sim$changepoints) > 0) {
    
    abline(
      v = sim$changepoints,
      lty = 2
    )
  }
  
  invisible(sim)
}

# ============================================================
# 8. Generate the five scenarios
# ============================================================

n <- 200
sigma <- 1

simulations <- list(
  
  # Scenario 1: no change
  no_change = generate_scenario(
    scenario = "no_change",
    n = n,
    sigma = sigma,
    intercept = 0,
    seed = 101
  ),
  
  # Scenario 2: single change
  single_change = generate_scenario(
    scenario = "single_change",
    n = n,
    sigma = sigma,
    intercept = 0,
    seed = 102
  ),
  
  # Scenario 3: several changes in the same direction
  same_direction = generate_scenario(
    scenario = "same_direction",
    n = n,
    K = 4,
    sigma = sigma,
    intercept = 0,
    seed = 103
  ),
  
  # Scenario 4: alternating +1 / -1 slopes
  alternating = generate_scenario(
    scenario = "alternating",
    n = n,
    K = 6,
    sigma = sigma,
    intercept = 0,
    seed = 1
  ),
  
  # Scenario 5: random slopes and random segment lengths
  random = generate_scenario(
    scenario = "random",
    n = n,
    K = 6,
    sigma = sigma,
    intercept = 0,
    seed = 2,
    min_segment_length = 10,
    dirichlet_alpha = 1
  )
)


# ============================================================
# 9. Plot the five scenarios
# ============================================================

old_par <- par(
  mfrow = c(3, 2),
  mar = c(4, 4, 3, 1)
)

plot_scenario(
  simulations$no_change,
  main = "No change"
)

plot_scenario(
  simulations$single_change,
  main = "Single change"
)

plot_scenario(
  simulations$same_direction,
  main = "Same-direction changes, K = 4"
)

plot_scenario(
  simulations$alternating,
  main = "Alternating slopes, K = 6"
)

plot_scenario(
  simulations$random,
  main = "Random changes, K = 6"
)

# Empty sixth panel
plot.new()

par(old_par)


# ============================================================
# 10. Inspect the true simulation parameters
# ============================================================

for (scenario in names(simulations)) {
  
  sim <- simulations[[scenario]]
  
  cat("\n")
  cat("========================================\n")
  cat("Scenario:", sim$scenario, "\n")
  cat("========================================\n")
  
  cat(
    "Number of segments K:",
    sim$K,
    "\n"
  )
  
  cat(
    "Segment lengths:",
    paste(
      sim$segment_lengths,
      collapse = ", "
    ),
    "\n"
  )
  
  cat(
    "Changepoints:",
    if (length(sim$changepoints) == 0) {
      "none"
    } else {
      paste(
        sim$changepoints,
        collapse = ", "
      )
    },
    "\n"
  )
  
  cat(
    "Slopes:",
    paste(
      round(sim$slopes, 3),
      collapse = ", "
    ),
    "\n"
  )
  
  cat(
    "Sigma:",
    sim$sigma,
    "\n"
  )
}


# ============================================================
# 11. Generate repeated simulations
# ============================================================

B <- 100

repeated_simulations <- lapply(
  seq_len(B),
  function(b) {
    
    generate_scenario(
      scenario = "random",
      n = 200,
      K = 5,
      sigma = 1,
      seed = 1000 + b,
      min_segment_length = 10,
      dirichlet_alpha = 1
    )
  }
)
