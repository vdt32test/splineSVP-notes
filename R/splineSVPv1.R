# ============================================================
# Spline-SVP v1: exhaustive oracle for tiny examples
# ============================================================


# ------------------------------------------------------------
# 1. Validity statistic T_{s,t}
#
# R indices s,t correspond directly to mathematical indices.
# The candidate segment is (s,t], so observations s+1,...,t
# are used in the local validity statistic.
# ------------------------------------------------------------

segment_validity_stat <- function(x, y, s, t, sigma) {
  
  stopifnot(
    length(x) == length(y),
    1 <= s,
    s < t,
    t <= length(x),
    sigma > 0
  )
  
  m <- t - s
  
  # The interval (s,t] contains m observations.
  #
  # If m <= 2, an affine function can interpolate all observations
  # exactly, hence RSS0 = 0 and therefore T_{s,t} = 0.
  if (m <= 2) {
    return(0)
  }
  
  idx <- (s + 1):t
  
  xs <- x[idx]
  ys <- y[idx]
  
  # ----------------------------------------------------------
  # Null model: one affine line
  # ----------------------------------------------------------
  
  X0 <- cbind(
    intercept = 1,
    x = xs
  )
  
  fit0 <- lm.fit(X0, ys)
  
  rss0 <- sum(fit0$residuals^2)
  
  # ----------------------------------------------------------
  # Alternative: one continuous internal slope change
  #
  # a + b*x + c*(x - x_r)_+
  # ----------------------------------------------------------
  
  candidate_r <- (s + 1):(t - 1)
  
  rss1_values <- numeric(length(candidate_r))
  
  for (j in seq_along(candidate_r)) {
    
    r <- candidate_r[j]
    
    hinge <- pmax(xs - x[r], 0)
    
    X1 <- cbind(
      intercept = 1,
      x = xs,
      hinge = hinge
    )
    
    fit1 <- lm.fit(X1, ys)
    
    rss1_values[j] <- sum(fit1$residuals^2)
  }
  
  best_rss1 <- min(rss1_values)
  
  # Theoretically best_rss1 <= rss0.
  # max(..., 0) protects against tiny numerical roundoff.
  improvement <- max(rss0 - best_rss1, 0)
  
  improvement / sigma^2
}


# ------------------------------------------------------------
# 2. Precompute all T_{s,t}
# ------------------------------------------------------------

compute_T_matrix <- function(x, y, sigma) {
  
  n <- length(y)
  
  Tmat <- matrix(
    NA_real_,
    nrow = n,
    ncol = n,
    dimnames = list(1:n, 1:n)
  )
  
  for (s in 1:(n - 1)) {
    for (t in (s + 1):n) {
      Tmat[s, t] <- segment_validity_stat(
        x, y, s, t, sigma
      )
    }
  }
  
  Tmat
}


# ------------------------------------------------------------
# 3. Construct the continuous piecewise-linear design matrix
#    for a fixed partition tau.
#
# tau = c(tau_0, ..., tau_K)
#
# Columns correspond to boundary states
# z_0, ..., z_K.
#
# Each observation is represented by linear interpolation
# between the two adjacent boundary states.
# ------------------------------------------------------------

spline_design_matrix <- function(x, tau) {
  
  n <- length(x)
  K <- length(tau) - 1
  
  X <- matrix(0, nrow = n, ncol = K + 1)
  
  # First observation corresponds to z_0
  X[1, 1] <- 1
  
  for (k in 1:K) {
    
    s <- tau[k]
    t <- tau[k + 1]
    
    idx <- (s + 1):t
    
    w <- (x[idx] - x[s]) / (x[t] - x[s])
    
    # l_{s,t}^{u,v}(x_i)
    # = (1-w_i) u + w_i v
    
    X[idx, k]     <- 1 - w
    X[idx, k + 1] <- w
  }
  
  X
}


# ------------------------------------------------------------
# 4. Solve the continuous least-squares problem for one
#    fixed partition.
#
# Returns z*(tau) and C*(tau).
# ------------------------------------------------------------

fit_fixed_partition <- function(x, y, tau) {
  
  X <- spline_design_matrix(x, tau)
  
  # QR-based least squares is preferable to explicitly computing
  # solve(t(X) %*% X) %*% t(X) %*% y.
  fit <- lm.fit(x = X, y = y)
  
  z_star <- as.numeric(fit$coefficients)
  fitted <- as.numeric(X %*% z_star)
  residuals <- y - fitted
  
  C_star <- sum(residuals^2)
  
  list(
    tau = tau,
    K = length(tau) - 1,
    z = z_star,
    cost = C_star,
    fitted = fitted,
    residuals = residuals
  )
}


# ------------------------------------------------------------
# 5. Enumerate all partitions
#
# Endpoints 1 and n are always included.
# Each index 2,...,n-1 is either an interior boundary or not.
#
# Number of partitions = 2^(n-2).
# ------------------------------------------------------------

all_partitions <- function(n) {
  
  stopifnot(n >= 2)
  
  interior <- if (n > 2) 2:(n - 1) else integer(0)
  m <- length(interior)
  
  partitions <- vector("list", 2^m)
  
  for (mask in 0:(2^m - 1)) {
    
    if (m == 0) {
      chosen <- integer(0)
    } else {
      bits <- as.logical(
        intToBits(mask)[seq_len(m)]
      )
      
      chosen <- interior[bits]
    }
    
    partitions[[mask + 1]] <- c(1, chosen, n)
  }
  
  partitions
}


# ------------------------------------------------------------
# 6. Check whether a partition belongs to V_gamma
# ------------------------------------------------------------

partition_is_valid <- function(tau, Tmat, gamma) {
  
  K <- length(tau) - 1
  
  for (k in 1:K) {
    
    s <- tau[k]
    t <- tau[k + 1]
    
    if (Tmat[s, t] > gamma) {
      return(FALSE)
    }
  }
  
  TRUE
}


# ------------------------------------------------------------
# 7. Exhaustive spline-SVP oracle
#
# Primary objective:
#       minimise K
#
# Secondary objective:
#       minimise C*(tau)
# ------------------------------------------------------------

spline_svp_oracle <- function(x, y, sigma, gamma,
                              verbose = FALSE) {
  
  n <- length(y)
  
  stopifnot(
    length(x) == n,
    n >= 2,
    all(diff(x) > 0),
    sigma > 0,
    gamma >= 0
  )
  
  # Precompute state-independent validity
  Tmat <- compute_T_matrix(x, y, sigma)
  
  # Enumerate all candidate partitions
  partitions <- all_partitions(n)
  
  best <- NULL
  number_valid <- 0L
  
  for (tau in partitions) {
    
    # Reject invalid partitions immediately
    if (!partition_is_valid(tau, Tmat, gamma)) {
      next
    }
    
    number_valid <- number_valid + 1L
    
    candidate <- fit_fixed_partition(x, y, tau)
    
    if (verbose) {
      cat(
        "tau =", paste(tau, collapse = ","),
        " | K =", candidate$K,
        " | cost =", candidate$cost,
        "\n"
      )
    }
    
    # Lexicographic comparison:
    # first K, then cost
    if (
      is.null(best) ||
      candidate$K < best$K ||
      (
        candidate$K == best$K &&
        candidate$cost < best$cost
      )
    ) {
      best <- candidate
    }
  }
  
  if (is.null(best)) {
    stop("No valid partition found.")
  }
  
  best$T <- Tmat
  best$gamma <- gamma
  best$sigma <- sigma
  best$number_partitions <- length(partitions)
  best$number_valid_partitions <- number_valid
  
  best
}


# ------------------------------------------------------------
# 8. Evaluate fitted spline at arbitrary x values
# ------------------------------------------------------------

predict_spline <- function(object, x, x_new) {
  
  tau <- object$tau
  z <- object$z
  K <- object$K
  
  out <- numeric(length(x_new))
  
  for (j in seq_along(x_new)) {
    
    xx <- x_new[j]
    
    if (xx < x[1] || xx > x[length(x)]) {
      out[j] <- NA_real_
      next
    }
    
    # Find segment containing xx
    k <- which(
      xx >= x[tau[1:K]] &
        xx <= x[tau[2:(K + 1)]]
    )[1]
    
    s <- tau[k]
    t <- tau[k + 1]
    
    w <- (xx - x[s]) / (x[t] - x[s])
    
    out[j] <- (1 - w) * z[k] + w * z[k + 1]
  }
  
  out
}


# ============================================================
# Tiny example
# ============================================================

set.seed(1)

n <- 10
x <- 1:n
sigma <- 0.25

# Continuous PWL true signal:
# slope changes around x = 5
f0 <- ifelse(
  x <= 5,
  0.5 * x,
  2.5 - 0.7 * (x - 5)
)

y <- f0 + rnorm(n, sd = sigma)

# Try a threshold
gamma <- 4

result <- spline_svp_oracle(
  x = x,
  y = y,
  sigma = sigma,
  gamma = gamma,
  verbose = TRUE
)

cat("\nSelected partition:\n")
print(result$tau)

cat("\nNumber of segments K:\n")
print(result$K)

cat("\nOptimal boundary states:\n")
print(result$z)

cat("\nOptimal cost C*(tau):\n")
print(result$cost)

cat("\nPartitions checked:\n")
print(result$number_partitions)

cat("\nValid partitions:\n")
print(result$number_valid_partitions)


# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

grid <- seq(min(x), max(x), length.out = 500)

fhat_grid <- predict_spline(
  result,
  x = x,
  x_new = grid
)

plot(
  x, y,
  pch = 19,
  xlab = "x",
  ylab = "y",
  main = "Exhaustive spline-SVP v1 oracle"
)

lines(
  x, f0,
  lty = 2,
  lwd = 2
)

lines(
  grid, fhat_grid,
  lwd = 2
)

points(
  x[result$tau],
  result$z,
  pch = 4,
  cex = 1.4,
  lwd = 2
)

legend(
  "topright",
  legend = c(
    "observations",
    "true signal",
    "spline-SVP fit",
    "fitted boundary states"
  ),
  pch = c(19, NA, NA, 4),
  lty = c(NA, 2, 1, NA),
  lwd = c(NA, 2, 2, 2)
)


result$tau
result$K
result$z
result$cost
result$number_partitions
result$number_valid_partitions

result$T


x_test <- 1:6
y_test <- c(1.1, 2.0, 3.2, 4.1, 5.0, 6.2)

stopifnot(
  segment_validity_stat(x_test, y_test, 1, 2, sigma = 1) == 0
)

stopifnot(
  segment_validity_stat(x_test, y_test, 1, 3, sigma = 1) == 0
)

cat("Short-interval tests passed.\n")


x_test <- 1:10
y_line <- 2 + 3 * x_test

T_line <- segment_validity_stat(
  x = x_test,
  y = y_line,
  s = 1,
  t = 10,
  sigma = 1
)

print(T_line)

stopifnot(abs(T_line) < 1e-10)

cat("Perfect-line test passed.\n")

x_test <- 1:10

y_knot <- ifelse(
  x_test <= 5,
  x_test,
  5 - 2 * (x_test - 5)
)

T_knot <- segment_validity_stat(
  x = x_test,
  y = y_knot,
  s = 1,
  t = 10,
  sigma = 1
)

print(T_knot)

stopifnot(T_knot > 0)

cat("Slope-change test passed.\n")

set.seed(123)

x_test <- 1:12
sigma_test <- 0.5

y_test <- ifelse(
  x_test <= 6,
  0.4 * x_test,
  2.4 - 0.8 * (x_test - 6)
) + rnorm(12, sd = sigma_test)

T1 <- segment_validity_stat(
  x_test, y_test,
  s = 1,
  t = 12,
  sigma = sigma_test
)

lambda <- 7

T2 <- segment_validity_stat(
  x_test,
  lambda * y_test,
  s = 1,
  t = 12,
  sigma = lambda * sigma_test
)

print(c(T1 = T1, T2 = T2))

stopifnot(
  abs(T1 - T2) < 1e-8
)

cat("Scale-invariance test passed.\n")


validity_details <- function(x, y, s, t, sigma) {
  
  m <- t - s
  
  if (m <= 2) {
    return(list(
      rss0 = 0,
      best_rss1 = 0,
      best_r = NA_integer_,
      T = 0
    ))
  }
  
  idx <- (s + 1):t
  xs <- x[idx]
  ys <- y[idx]
  
  X0 <- cbind(1, xs)
  fit0 <- lm.fit(X0, ys)
  rss0 <- sum(fit0$residuals^2)
  
  candidate_r <- (s + 1):(t - 1)
  
  rss1 <- numeric(length(candidate_r))
  
  for (j in seq_along(candidate_r)) {
    
    r <- candidate_r[j]
    
    hinge <- pmax(xs - x[r], 0)
    
    X1 <- cbind(1, xs, hinge)
    
    fit1 <- lm.fit(X1, ys)
    
    rss1[j] <- sum(fit1$residuals^2)
  }
  
  j_best <- which.min(rss1)
  
  best_rss1 <- rss1[j_best]
  best_r <- candidate_r[j_best]
  
  T <- max(rss0 - best_rss1, 0) / sigma^2
  
  list(
    rss0 = rss0,
    best_rss1 = best_rss1,
    best_r = best_r,
    T = T
  )
}
details <- validity_details(
  x = x_test,
  y = y_test,
  s = 1,
  t = 12,
  sigma = sigma_test
)

print(details)

stopifnot(
  details$best_rss1 <= details$rss0 + 1e-10
)

cat("Nested-model RSS test passed.\n")


# ============================================================
# Additional unit tests for spline-SVP v1 reference oracle
# ============================================================

cat("\n========================================\n")
cat("Additional spline-SVP v1 unit tests\n")
cat("========================================\n")


# ============================================================
# TEST 1. Fixed-partition fitting
# Compare one-segment spline fit against ordinary linear regression
# ============================================================

cat("\nTEST 1: Fixed-partition fitting\n")

set.seed(101)

x1 <- 1:12
y1 <- 1.5 + 0.7 * x1 + rnorm(length(x1), sd = 0.4)

tau1 <- c(1, length(x1))

fit_spline_1 <- fit_fixed_partition(
  x = x1,
  y = y1,
  tau = tau1
)

X_lm <- cbind(
  intercept = 1,
  x = x1
)

fit_lm_1 <- lm.fit(
  x = X_lm,
  y = y1
)

rss_lm_1 <- sum(fit_lm_1$residuals^2)

difference_cost <- abs(
  fit_spline_1$cost - rss_lm_1
)

cat("Spline cost     =", fit_spline_1$cost, "\n")
cat("Linear-reg RSS  =", rss_lm_1, "\n")
cat("Absolute diff.  =", difference_cost, "\n")

stopifnot(
  difference_cost < 1e-10
)

cat("RESULT: PASSED\n")


# ============================================================
# TEST 2. Continuity of fitted spline
# Explicitly compare left and right values at all interior boundaries
# ============================================================

cat("\nTEST 2: Continuity of fitted spline\n")

x2 <- 1:10

y2 <- c(
  0.4, 0.9, 1.5, 2.2, 2.6,
  2.0, 1.4, 0.9, 0.2, -0.5
)

tau2 <- c(1, 5, 8, 10)

fit2 <- fit_fixed_partition(
  x = x2,
  y = y2,
  tau = tau2
)

z2 <- fit2$z
K2 <- fit2$K

continuity_differences <- numeric(K2 - 1)

for (k in 1:(K2 - 1)) {
  
  s_left <- tau2[k]
  t_left <- tau2[k + 1]
  
  s_right <- tau2[k + 1]
  t_right <- tau2[k + 2]
  
  x_boundary <- x2[tau2[k + 1]]
  
  # Left fitted value at boundary
  left_value <-
    z2[k] +
    (z2[k + 1] - z2[k]) /
    (x2[t_left] - x2[s_left]) *
    (x_boundary - x2[s_left])
  
  # Right fitted value at boundary
  right_value <-
    z2[k + 1] +
    (z2[k + 2] - z2[k + 1]) /
    (x2[t_right] - x2[s_right]) *
    (x_boundary - x2[s_right])
  
  continuity_differences[k] <-
    abs(left_value - right_value)
  
  cat(
    "Boundary index", tau2[k + 1],
    ": left =", left_value,
    ", right =", right_value,
    ", diff =", continuity_differences[k],
    "\n"
  )
}

max_continuity_error <- max(continuity_differences)

cat("Maximum continuity error =", max_continuity_error, "\n")

stopifnot(
  max_continuity_error < 1e-12
)

cat("RESULT: PASSED\n")


# ============================================================
# TEST 3. Partition enumeration
# Verify exactly 2^(n-2) unique partitions
# ============================================================

cat("\nTEST 3: Partition enumeration\n")

n_values <- 2:10

for (n in n_values) {
  
  parts <- all_partitions(n)
  
  expected <- 2^(n - 2)
  observed <- length(parts)
  
  # Convert partitions to strings to check uniqueness
  keys <- vapply(
    parts,
    function(tau) paste(tau, collapse = "-"),
    character(1)
  )
  
  unique_count <- length(unique(keys))
  
  cat(
    "n =", n,
    "| expected =", expected,
    "| observed =", observed,
    "| unique =", unique_count,
    "\n"
  )
  
  stopifnot(
    observed == expected,
    unique_count == expected
  )
}

cat("RESULT: PASSED\n")


# ============================================================
# TEST 4. Partition validity
# Controlled valid and invalid partitions
# ============================================================

cat("\nTEST 4: Partition validity\n")

# Exact continuous one-knot signal
x4 <- 1:10

y4 <- ifelse(
  x4 <= 5,
  x4,
  5 - 2 * (x4 - 5)
)

sigma4 <- 1
gamma4 <- 4

Tmat4 <- compute_T_matrix(
  x = x4,
  y = y4,
  sigma = sigma4
)

tau_invalid <- c(1, 10)
tau_valid <- c(1, 5, 10)

is_invalid_partition_valid <- partition_is_valid(
  tau_invalid,
  Tmat4,
  gamma4
)

is_valid_partition_valid <- partition_is_valid(
  tau_valid,
  Tmat4,
  gamma4
)

cat(
  "Partition (1,10):",
  is_invalid_partition_valid,
  "| T(1,10) =", Tmat4[1, 10],
  "\n"
)

cat(
  "Partition (1,5,10):",
  is_valid_partition_valid,
  "| T(1,5) =", Tmat4[1, 5],
  "| T(5,10) =", Tmat4[5, 10],
  "\n"
)

stopifnot(
  is_invalid_partition_valid == FALSE,
  is_valid_partition_valid == TRUE
)

cat("RESULT: PASSED\n")


# ============================================================
# TEST 5. Lexicographic oracle selection
# Verify that K is primary objective, cost secondary
# ============================================================

cat("\nTEST 5: Lexicographic oracle selection\n")

set.seed(1)

n5 <- 10
x5 <- 1:n5
sigma5 <- 0.25

f5 <- ifelse(
  x5 <= 5,
  0.5 * x5,
  2.5 - 0.7 * (x5 - 5)
)

y5 <- f5 + rnorm(n5, sd = sigma5)

gamma5 <- 4

oracle5 <- spline_svp_oracle(
  x = x5,
  y = y5,
  sigma = sigma5,
  gamma = gamma5
)

Tmat5 <- compute_T_matrix(
  x5, y5, sigma5
)

parts5 <- all_partitions(n5)

valid_results5 <- list()

for (tau in parts5) {
  
  if (partition_is_valid(tau, Tmat5, gamma5)) {
    
    fit_tau <- fit_fixed_partition(
      x5, y5, tau
    )
    
    valid_results5[[length(valid_results5) + 1]] <-
      data.frame(
        tau = paste(tau, collapse = ","),
        K = fit_tau$K,
        cost = fit_tau$cost,
        stringsAsFactors = FALSE
      )
  }
}

valid_table5 <- do.call(
  rbind,
  valid_results5
)

minimum_K5 <- min(valid_table5$K)

same_K5 <- valid_table5[
  valid_table5$K == minimum_K5,
]

best_cost5 <- min(same_K5$cost)

cat("Oracle K        =", oracle5$K, "\n")
cat("Minimum valid K =", minimum_K5, "\n")
cat("Oracle cost     =", oracle5$cost, "\n")
cat("Best cost at K* =", best_cost5, "\n")
cat(
  "Oracle tau      =",
  paste(oracle5$tau, collapse = ","),
  "\n"
)

stopifnot(
  oracle5$K == minimum_K5,
  abs(oracle5$cost - best_cost5) < 1e-10
)

cat("RESULT: PASSED\n")


# ============================================================
# TEST 6. Global oracle correctness
# Independently enumerate all valid objective values and compare
# ============================================================

cat("\nTEST 6: Global oracle correctness\n")

# Independent exhaustive table already constructed above.
# Sort lexicographically:
#   first by K,
#   then by cost.

ordered5 <- valid_table5[
  order(
    valid_table5$K,
    valid_table5$cost
  ),
]

independent_best <- ordered5[1, ]

cat("Independent best tau  =", independent_best$tau, "\n")
cat("Independent best K    =", independent_best$K, "\n")
cat("Independent best cost =", independent_best$cost, "\n")

cat(
  "Oracle best tau       =",
  paste(oracle5$tau, collapse = ","),
  "\n"
)
cat("Oracle best K         =", oracle5$K, "\n")
cat("Oracle best cost      =", oracle5$cost, "\n")

stopifnot(
  independent_best$K == oracle5$K,
  abs(
    independent_best$cost -
      oracle5$cost
  ) < 1e-10,
  independent_best$tau ==
    paste(oracle5$tau, collapse = ",")
)

cat("RESULT: PASSED\n")


# ============================================================
# SUMMARY
# ============================================================

cat("\n========================================\n")
cat("ALL SIX ADDITIONAL TESTS PASSED\n")
cat("========================================\n")




# ============================================================
# Spline-SVP v1
# Gaussian-null calibration of the validity statistic
#
# Current statistic:
#
# T_{s,t}
# =
# [ RSS0(s,t) - min_r RSS1(s,r,t) ] / sigma^2
#
# Null hypothesis on (s,t]:
#
# H0:
#   E[Y_i] = a + b x_i,
#   epsilon_i iid ~ N(0, sigma^2)
#
# Goal:
#   estimate null quantiles of T_{s,t}
#   and investigate dependence on
#       1. interval length m = t-s;
#       2. design geometry x.
# ============================================================


# ============================================================
# 1. Validity statistic
# ============================================================

segment_validity_stat <- function(x, y, s, t, sigma) {
  
  stopifnot(
    length(x) == length(y),
    1 <= s,
    s < t,
    t <= length(x),
    sigma > 0
  )
  
  m <- t - s
  
  # Under the present v1 convention, if the interval contains
  # at most two observations, an affine fit interpolates exactly.
  if (m <= 2) {
    return(0)
  }
  
  idx <- (s + 1):t
  
  xs <- x[idx]
  ys <- y[idx]
  
  # ----------------------------------------------------------
  # Null model:
  # y_i = a + b x_i + epsilon_i
  # ----------------------------------------------------------
  
  X0 <- cbind(
    intercept = 1,
    x = xs
  )
  
  fit0 <- lm.fit(
    x = X0,
    y = ys
  )
  
  rss0 <- sum(fit0$residuals^2)
  
  # ----------------------------------------------------------
  # Alternative:
  # y_i = a + b x_i + c (x_i - x_r)_+ + epsilon_i
  # ----------------------------------------------------------
  
  candidate_r <- (s + 1):(t - 1)
  
  rss1_values <- numeric(length(candidate_r))
  
  for (j in seq_along(candidate_r)) {
    
    r <- candidate_r[j]
    
    hinge <- pmax(xs - x[r], 0)
    
    X1 <- cbind(
      intercept = 1,
      x = xs,
      hinge = hinge
    )
    
    fit1 <- lm.fit(
      x = X1,
      y = ys
    )
    
    rss1_values[j] <- sum(fit1$residuals^2)
  }
  
  best_rss1 <- min(rss1_values)
  
  # Protection against tiny floating-point negatives.
  improvement <- max(
    rss0 - best_rss1,
    0
  )
  
  improvement / sigma^2
}


# ============================================================
# 2. Simulate T under H0 for one fixed design
# ============================================================

simulate_null_T_design <- function(
    x,
    B = 10000,
    sigma = 1,
    a = 0,
    b = 0,
    seed = NULL) {
  
  stopifnot(
    length(x) >= 4,
    all(diff(x) > 0),
    B >= 1,
    sigma > 0
  )
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # IMPORTANT:
  # segment_validity_stat() uses observations (s,t].
  #
  # With
  #   s = 1,
  #   t = length(x),
  #
  # the statistic uses observations 2,...,length(x).
  #
  # Therefore:
  #
  #   m = length(x) - 1
  #
  # observations enter the local validity calculation.
  
  n <- length(x)
  
  T_values <- numeric(B)
  
  for (sim in seq_len(B)) {
    
    eps <- rnorm(
      n,
      mean = 0,
      sd = sigma
    )
    
    y <- a + b * x + eps
    
    T_values[sim] <- segment_validity_stat(
      x = x,
      y = y,
      s = 1,
      t = n,
      sigma = sigma
    )
  }
  
  T_values
}


# ============================================================
# 3. Equally spaced design for a given m
# ============================================================

make_equal_design <- function(m) {
  
  stopifnot(m >= 3)
  
  # There are m+1 design points,
  # but observations 2,...,m+1 enter T.
  0:m
}


# ============================================================
# 4. Calibration for one segment length m
# ============================================================

calibrate_one_length <- function(
    m,
    B = 10000,
    seed = 123) {
  
  x <- make_equal_design(m)
  
  T_values <- simulate_null_T_design(
    x = x,
    B = B,
    sigma = 1,
    a = 0,
    b = 0,
    seed = seed
  )
  
  qs <- quantile(
    T_values,
    probs = c(
      0.90,
      0.95,
      0.99
    ),
    names = FALSE,
    type = 7
  )
  
  data.frame(
    m = m,
    B = B,
    mean_T = mean(T_values),
    sd_T = sd(T_values),
    q90 = qs[1],
    q95 = qs[2],
    q99 = qs[3]
  )
}


# ============================================================
# 5. Calibration over many segment lengths
# ============================================================

calibrate_by_length <- function(
    m_values,
    B = 10000,
    seed = 123) {
  
  stopifnot(
    all(m_values >= 3)
  )
  
  out <- vector(
    "list",
    length(m_values)
  )
  
  for (j in seq_along(m_values)) {
    
    m <- m_values[j]
    
    cat(
      "Calibrating m =",
      m,
      "(",
      j,
      "of",
      length(m_values),
      ")\n"
    )
    
    # Different deterministic seed for each m.
    current_seed <- seed + j - 1
    
    out[[j]] <- calibrate_one_length(
      m = m,
      B = B,
      seed = current_seed
    )
  }
  
  do.call(
    rbind,
    out
  )
}


# ============================================================
# 6. RUN 1:
# Calibration versus segment length
# ============================================================

# Quick exploratory choice.
#
# For a final run, you may use:
#
#   m_values <- 3:100
#   B <- 20000
#
# depending on runtime.

m_values <- c(
  3, 4, 5, 6, 8, 10,
  15, 20, 30, 40, 50
)

B_calibration <- 10000

calibration_by_length <- calibrate_by_length(
  m_values = m_values,
  B = B_calibration,
  seed = 1000
)

cat("\n========================================\n")
cat("Calibration by segment length\n")
cat("========================================\n")

print(
  calibration_by_length,
  row.names = FALSE
)


# ============================================================
# 7. Save calibration table
# ============================================================

write.csv(
  calibration_by_length,
  file = "gaussian_null_calibration_by_length.csv",
  row.names = FALSE
)


# ============================================================
# 8. Plot null quantiles versus m
# ============================================================

png(
  filename = "gaussian_null_quantiles_by_length.png",
  width = 1200,
  height = 800,
  res = 150
)

plot(
  calibration_by_length$m,
  calibration_by_length$q95,
  type = "b",
  pch = 19,
  lwd = 2,
  xlab = "Segment length m",
  ylab = "Null quantile of T",
  main = "Gaussian-null calibration by segment length"
)

lines(
  calibration_by_length$m,
  calibration_by_length$q90,
  type = "b",
  pch = 1,
  lty = 2
)

lines(
  calibration_by_length$m,
  calibration_by_length$q99,
  type = "b",
  pch = 2,
  lty = 3
)

legend(
  "topleft",
  legend = c(
    "90% quantile",
    "95% quantile",
    "99% quantile"
  ),
  lty = c(2, 1, 3),
  pch = c(1, 19, 2),
  lwd = 2
)

dev.off()


# ============================================================
# 9. Plot mean of T versus m
# ============================================================

png(
  filename = "gaussian_null_mean_T_by_length.png",
  width = 1200,
  height = 800,
  res = 150
)

plot(
  calibration_by_length$m,
  calibration_by_length$mean_T,
  type = "b",
  pch = 19,
  lwd = 2,
  xlab = "Segment length m",
  ylab = "Mean null value of T",
  main = "Mean Gaussian-null statistic by segment length"
)

dev.off()


# ============================================================
# 10. Independent validation of a calibrated threshold
#
# Example:
#   m = 20
#   alpha = 0.05
#
# Calibration sample estimates gamma.
# Independent sample checks P(T > gamma).
# ============================================================

m_validation <- 20

B_threshold <- 20000
B_validation <- 20000

x_validation <- make_equal_design(
  m_validation
)

# Calibration sample
T_cal <- simulate_null_T_design(
  x = x_validation,
  B = B_threshold,
  sigma = 1,
  a = 0,
  b = 0,
  seed = 2001
)

gamma_95 <- unname(
  quantile(
    T_cal,
    probs = 0.95
  )
)

# Independent validation sample
T_val <- simulate_null_T_design(
  x = x_validation,
  B = B_validation,
  sigma = 1,
  a = 0,
  b = 0,
  seed = 2002
)

rejection_rate <- mean(
  T_val > gamma_95
)

cat("\n========================================\n")
cat("Independent threshold validation\n")
cat("========================================\n")

cat(
  "m                  =",
  m_validation,
  "\n"
)

cat(
  "Estimated gamma95  =",
  gamma_95,
  "\n"
)

cat(
  "Target rejection   = 0.05\n"
)

cat(
  "Observed rejection =",
  rejection_rate,
  "\n"
)


# ============================================================
# 11. Histogram of one null distribution
# ============================================================

png(
  filename = "gaussian_null_T_m20_histogram.png",
  width = 1200,
  height = 800,
  res = 150
)

hist(
  T_val,
  breaks = 50,
  probability = TRUE,
  main = paste(
    "Gaussian-null distribution of T, m =",
    m_validation
  ),
  xlab = "T"
)

abline(
  v = gamma_95,
  lwd = 2,
  lty = 2
)

legend(
  "topright",
  legend = c(
    "Estimated 95% threshold"
  ),
  lty = 2,
  lwd = 2
)

dev.off()


# ============================================================
# 12. Check affine-shift invariance
#
# Under H0:
#
# y = epsilon
#
# and
#
# y = a + b x + epsilon
#
# should give the same T for the same epsilon realization.
# ============================================================

set.seed(3001)

m_inv <- 20
x_inv <- make_equal_design(m_inv)

eps_inv <- rnorm(
  length(x_inv),
  mean = 0,
  sd = 1
)

T_base <- segment_validity_stat(
  x = x_inv,
  y = eps_inv,
  s = 1,
  t = length(x_inv),
  sigma = 1
)

a_inv <- 7
b_inv <- -2.3

y_affine <- a_inv + b_inv * x_inv + eps_inv

T_affine <- segment_validity_stat(
  x = x_inv,
  y = y_affine,
  s = 1,
  t = length(x_inv),
  sigma = 1
)

cat("\n========================================\n")
cat("Affine-shift invariance check\n")
cat("========================================\n")

cat(
  "T without affine trend =",
  T_base,
  "\n"
)

cat(
  "T with affine trend    =",
  T_affine,
  "\n"
)

cat(
  "Absolute difference    =",
  abs(T_base - T_affine),
  "\n"
)

stopifnot(
  abs(T_base - T_affine) < 1e-8
)

cat("RESULT: PASSED\n")


# ============================================================
# 13. Check scale invariance
# ============================================================

lambda_inv <- 4.7

T_scaled <- segment_validity_stat(
  x = x_inv,
  y = lambda_inv * y_affine,
  s = 1,
  t = length(x_inv),
  sigma = lambda_inv
)

cat("\n========================================\n")
cat("Scale-invariance check\n")
cat("========================================\n")

cat(
  "Original T =",
  T_affine,
  "\n"
)

cat(
  "Scaled T   =",
  T_scaled,
  "\n"
)

cat(
  "Difference =",
  abs(T_affine - T_scaled),
  "\n"
)

stopifnot(
  abs(T_affine - T_scaled) < 1e-8
)

cat("RESULT: PASSED\n")


# ============================================================
# 14. Design-geometry sensitivity
#
# Hold m fixed and change the spacing of x.
# ============================================================

m_design <- 20
B_design <- 10000

# ------------------------------------------------------------
# Design A: equally spaced
# ------------------------------------------------------------

x_equal <- 0:m_design


# ------------------------------------------------------------
# Design B: alternating short/long gaps
# ------------------------------------------------------------

gaps_alt <- rep(
  c(0.5, 1.5),
  length.out = m_design
)

x_alternating <- c(
  0,
  cumsum(gaps_alt)
)


# ------------------------------------------------------------
# Design C: increasingly large gaps
# ------------------------------------------------------------

gaps_increasing <- seq(
  0.2,
  2,
  length.out = m_design
)

x_increasing <- c(
  0,
  cumsum(gaps_increasing)
)


# ------------------------------------------------------------
# Design D: clustered near the left side
# ------------------------------------------------------------

u_clustered <- seq(
  0,
  1,
  length.out = m_design + 1
)

x_clustered <- u_clustered^3


# ------------------------------------------------------------
# Simulate null distributions
# ------------------------------------------------------------

T_equal <- simulate_null_T_design(
  x = x_equal,
  B = B_design,
  seed = 4001
)

T_alternating <- simulate_null_T_design(
  x = x_alternating,
  B = B_design,
  seed = 4002
)

T_increasing <- simulate_null_T_design(
  x = x_increasing,
  B = B_design,
  seed = 4003
)

T_clustered <- simulate_null_T_design(
  x = x_clustered,
  B = B_design,
  seed = 4004
)


# ============================================================
# 15. Summarize design sensitivity
# ============================================================

summarize_T <- function(T_values) {
  
  qs <- quantile(
    T_values,
    probs = c(
      0.90,
      0.95,
      0.99
    ),
    names = FALSE
  )
  
  c(
    mean = mean(T_values),
    q90 = qs[1],
    q95 = qs[2],
    q99 = qs[3]
  )
}


design_calibration <- rbind(
  equally_spaced = summarize_T(T_equal),
  alternating_gaps = summarize_T(T_alternating),
  increasing_gaps = summarize_T(T_increasing),
  clustered = summarize_T(T_clustered)
)

cat("\n========================================\n")
cat("Design-geometry sensitivity\n")
cat("========================================\n")

print(
  design_calibration
)

write.csv(
  design_calibration,
  file = "gaussian_null_design_sensitivity.csv",
  row.names = TRUE
)


# ============================================================
# 16. Plot 95% quantile by design
# ============================================================

png(
  filename = "gaussian_null_design_q95.png",
  width = 1200,
  height = 800,
  res = 150
)

barplot(
  design_calibration[, "q95"],
  names.arg = c(
    "Equal",
    "Alternating",
    "Increasing",
    "Clustered"
  ),
  ylab = "95% null quantile of T",
  main = paste(
    "Design sensitivity at m =",
    m_design
  )
)

dev.off()


# ============================================================
# 17. Compare an equal-design threshold on other designs
#
# This asks:
#
# If gamma is calibrated under equally spaced x,
# what false-rejection probability results under
# the other designs?
# ============================================================

gamma_equal_95 <- unname(
  quantile(
    T_equal,
    0.95
  )
)

rejection_by_design <- data.frame(
  design = c(
    "equally_spaced",
    "alternating_gaps",
    "increasing_gaps",
    "clustered"
  ),
  rejection_rate = c(
    mean(T_equal > gamma_equal_95),
    mean(T_alternating > gamma_equal_95),
    mean(T_increasing > gamma_equal_95),
    mean(T_clustered > gamma_equal_95)
  )
)

cat("\n========================================\n")
cat("Using equal-design gamma95 on all designs\n")
cat("========================================\n")

cat(
  "Equal-design gamma95 =",
  gamma_equal_95,
  "\n\n"
)

print(
  rejection_by_design,
  row.names = FALSE
)

write.csv(
  rejection_by_design,
  file = "gaussian_null_cross_design_rejection_rates.csv",
  row.names = FALSE
)


# ============================================================
# 18. Optional: denser equally spaced calibration
#
# Uncomment for a larger final run.
# ============================================================

# calibration_dense <- calibrate_by_length(
#   m_values = 3:100,
#   B = 20000,
#   seed = 5000
# )
#
# write.csv(
#   calibration_dense,
#   "gaussian_null_calibration_dense.csv",
#   row.names = FALSE
# )


# ============================================================
# 19. Final console summary
# ============================================================

cat("\n========================================\n")
cat("GAUSSIAN-NULL CALIBRATION COMPLETED\n")
cat("========================================\n")

cat(
  "\nFiles written:\n",
  "- gaussian_null_calibration_by_length.csv\n",
  "- gaussian_null_quantiles_by_length.png\n",
  "- gaussian_null_mean_T_by_length.png\n",
  "- gaussian_null_T_m20_histogram.png\n",
  "- gaussian_null_design_sensitivity.csv\n",
  "- gaussian_null_design_q95.png\n",
  "- gaussian_null_cross_design_rejection_rates.csv\n"
)

plot(
  calibration_by_length$m,
  calibration_by_length$q95,
  type = "b",
  pch = 19,
  lwd = 2,
  ylim = range(
    calibration_by_length$q90,
    calibration_by_length$q95,
    calibration_by_length$q99
  ),
  xlab = "Segment length m",
  ylab = "Null quantile of T",
  main = "Gaussian-null calibration by segment length"
)

lines(
  calibration_by_length$m,
  calibration_by_length$q90,
  type = "b",
  pch = 1,
  lty = 2
)

lines(
  calibration_by_length$m,
  calibration_by_length$q99,
  type = "b",
  pch = 2,
  lty = 3
)

legend(
  "bottomright",
  legend = c("90% quantile", "95% quantile", "99% quantile"),
  lty = c(2, 1, 3),
  pch = c(1, 19, 2),
  lwd = 2
)