#!/usr/bin/env Rscript
# SPDX-License-Identifier: MIT
# SPDX-FileComment: Independent R reference driver for automatic time-series analysis.

# Independent R reference driver for the ts_auto Fortran program.

stop_cli <- function(message, status = 1L) {
  writeLines(message, con = stderr())
  quit(save = "no", status = status, runLast = FALSE)
}

print_usage <- function() {
  cat(paste0(
    "Usage: Rscript tools/ts_auto_reference.R FILE [options]\n\n",
    "FILE contains numeric columns with an optional header and leading index column.\n",
    "Options:\n",
    "  --frequency N       observations per seasonal cycle (default 1)\n",
    "  --horizon N         forecast horizon (default 10)\n",
    "  --max-lag N         maximum ACF and PACF lag (default 24)\n",
    "  --display-lags N    positive correlation lags to print (default 5)\n",
    "  --max-ar N          maximum autoregressive order (default 12)\n",
    "  --validation N      held-out tail length (default automatic)\n",
    "  --selection NAME    validation, aicc, or bic (default validation)\n",
    "  --target NAME       mean or volatility (default mean)\n",
    "  --print-parameters  print fitted parameters for every candidate\n",
    "  --print-param       alias for --print-parameters\n",
    "  --param             alias for --print-parameters\n",
    "  --time-fits         report elapsed time for each candidate fit\n",
    "  --print-all-ar      report every tested autoregressive order\n",
    "  --print-all-arma    report every tested ARMA order\n",
    "  --corr              print transformed-data correlation matrices\n",
    "  --resample [N]      analyze N IID row-bootstrap samples (default 1)\n",
    "  --seed N            bootstrap random-number seed (default 12345)\n",
    "  --transform NAME    none, log, diff, or log-diff (default none)\n",
    "  --log               alias for --transform log\n",
    "  --diff              alias for --transform diff\n",
    "  --log-diff          alias for --transform log-diff\n",
    "  --stride N...       analyze the series separately at each stride\n",
    "  --max-models N      maximum model summaries; zero displays all\n",
    "  --print-forecasts N maximum forecasts printed; zero prints none\n",
    "  --observations N    use at most the first N valid observations\n",
    "  --obs N             alias for --observations\n"
  ))
}

integer_value <- function(value, option, minimum = NULL) {
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || as.character(parsed) != value) {
    stop_cli(sprintf("Invalid integer for %s: %s", option, value))
  }
  if (!is.null(minimum) && parsed < minimum) {
    stop_cli(sprintf("Value for %s must be at least %d.", option, minimum))
  }
  parsed
}

set_transform <- function(options, requested) {
  requested <- tolower(requested)
  allowed <- c("none", "log", "diff", "log-diff")
  if (!(requested %in% allowed)) {
    stop_cli(sprintf(
      "Invalid transformation: %s. Expected none, log, diff, or log-diff.",
      requested
    ))
  }
  if (options$transform_set && options$transform != requested) {
    stop_cli(sprintf(
      "Conflicting transformations: %s and %s.",
      options$transform, requested
    ))
  }
  options$transform <- requested
  options$transform_set <- TRUE
  options
}

parse_options <- function(arguments) {
  if (length(arguments) == 0L) {
    print_usage()
    quit(save = "no", status = 1L, runLast = FALSE)
  }
  if (arguments[[1L]] %in% c("--help", "-h")) {
    print_usage()
    quit(save = "no", status = 0L, runLast = FALSE)
  }
  options <- list(
    file = arguments[[1L]],
    frequency = 1L,
    horizon = 10L,
    max_lag = 24L,
    display_lags = 5L,
    max_ar = 12L,
    validation = 0L,
    selection = "validation",
    target = "mean",
    print_parameters = FALSE,
    time_fits = FALSE,
    print_all_ar = FALSE,
    print_all_arma = FALSE,
    corr = FALSE,
    resamples = 0L,
    seed = 12345L,
    transform = "none",
    transform_set = FALSE,
    strides = integer(),
    max_models = 0L,
    print_forecasts = .Machine$integer.max,
    observations = 0L
  )
  value_options <- c(
    "--frequency", "--horizon", "--max-lag", "--display-lags",
    "--max-ar", "--validation", "--selection", "--target",
    "--seed", "--transform", "--max-models", "--print-forecasts",
    "--observations", "--obs"
  )
  i <- 2L
  while (i <= length(arguments)) {
    argument <- arguments[[i]]
    if (argument %in% c("--help", "-h")) {
      print_usage()
      quit(save = "no", status = 0L, runLast = FALSE)
    }
    if (argument %in% c("--print-parameters", "--print-param", "--param")) {
      options$print_parameters <- TRUE
      i <- i + 1L
      next
    }
    if (argument == "--time-fits") {
      options$time_fits <- TRUE
      i <- i + 1L
      next
    }
    if (argument == "--print-all-ar") {
      options$print_all_ar <- TRUE
      i <- i + 1L
      next
    }
    if (argument == "--print-all-arma") {
      options$print_all_arma <- TRUE
      i <- i + 1L
      next
    }
    if (argument == "--corr") {
      options$corr <- TRUE
      i <- i + 1L
      next
    }
    if (argument == "--resample") {
      options$resamples <- 1L
      if (i < length(arguments) && !startsWith(arguments[[i + 1L]], "--")) {
        options$resamples <- integer_value(
          arguments[[i + 1L]], "--resample", 1L
        )
        i <- i + 1L
      }
      i <- i + 1L
      next
    }
    if (argument %in% c("--log", "--diff", "--log-diff")) {
      options <- set_transform(options, substring(argument, 3L))
      i <- i + 1L
      next
    }
    if (argument == "--stride") {
      i <- i + 1L
      first <- i
      while (i <= length(arguments) && !startsWith(arguments[[i]], "--")) {
        options$strides <- c(
          options$strides,
          integer_value(arguments[[i]], "--stride", 1L)
        )
        i <- i + 1L
      }
      if (i == first) {
        stop_cli("Option --stride requires at least one value.")
      }
      next
    }
    if (!(argument %in% value_options)) {
      stop_cli(sprintf("Unrecognized option: %s", argument))
    }
    if (i == length(arguments)) {
      stop_cli(sprintf("Option %s requires a value.", argument))
    }
    value <- arguments[[i + 1L]]
    if (argument == "--frequency") {
      options$frequency <- integer_value(value, argument, 1L)
    } else if (argument == "--horizon") {
      options$horizon <- integer_value(value, argument, 1L)
    } else if (argument == "--max-lag") {
      options$max_lag <- integer_value(value, argument, 1L)
    } else if (argument == "--display-lags") {
      options$display_lags <- integer_value(value, argument, 1L)
    } else if (argument == "--max-ar") {
      options$max_ar <- integer_value(value, argument, 1L)
    } else if (argument == "--validation") {
      options$validation <- integer_value(value, argument, 0L)
    } else if (argument == "--selection") {
      options$selection <- tolower(value)
      if (!(options$selection %in% c("validation", "aicc", "bic"))) {
        stop_cli("Selection must be validation, aicc, or bic.")
      }
    } else if (argument == "--target") {
      options$target <- tolower(value)
      if (!(options$target %in% c("mean", "volatility"))) {
        stop_cli("Target must be mean or volatility.")
      }
    } else if (argument == "--seed") {
      options$seed <- integer_value(value, argument, 0L)
    } else if (argument == "--transform") {
      options <- set_transform(options, value)
    } else if (argument == "--max-models") {
      options$max_models <- integer_value(value, argument, 0L)
    } else if (argument == "--print-forecasts") {
      options$print_forecasts <- integer_value(value, argument, 0L)
    } else if (argument %in% c("--observations", "--obs")) {
      options$observations <- integer_value(value, argument, 0L)
    }
    i <- i + 2L
  }
  if (length(options$strides) == 0L) {
    options$strides <- 1L
  }
  options$max_lag <- max(options$max_lag, options$display_lags)
  options
}

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop_cli(sprintf(
      "R package '%s' is required. Install it with install.packages('%s').",
      package, package
    ))
  }
}

field_is_finite_number <- function(value) {
  parsed <- suppressWarnings(as.numeric(value))
  !is.na(parsed) && is.finite(parsed)
}

read_input_table <- function(path, observation_limit) {
  if (!file.exists(path)) {
    stop_cli(sprintf("Input file does not exist: %s", path))
  }
  lines <- readLines(path, warn = FALSE)
  nonempty <- which(nzchar(trimws(lines)))
  if (length(nonempty) == 0L) {
    stop_cli("The input file contains no data.")
  }
  first_fields <- strsplit(lines[[nonempty[[1L]]]], ",", fixed = TRUE)[[1L]]
  numeric_flags <- vapply(first_fields, field_is_finite_number, logical(1L))
  first_numeric <- match(TRUE, numeric_flags, nomatch = 0L)
  first_is_data <- first_numeric > 0L && all(numeric_flags[first_numeric:length(numeric_flags)])
  has_header <- !first_is_data
  table <- tryCatch(
    utils::read.csv(
      path,
      header = has_header,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      na.strings = c("NA", "NaN", "")
    ),
    error = function(error) stop_cli(conditionMessage(error))
  )
  if (nrow(table) == 0L || ncol(table) == 0L) {
    stop_cli("No data rows were found.")
  }
  numeric_columns <- vapply(table, is.numeric, logical(1L))
  first_data_column <- match(TRUE, numeric_columns, nomatch = 0L)
  if (first_data_column == 0L || !all(numeric_columns[first_data_column:ncol(table)])) {
    stop_cli("No contiguous numeric data columns were found after the index fields.")
  }
  numeric_table <- table[, first_data_column:ncol(table), drop = FALSE]
  values <- as.matrix(numeric_table)
  storage.mode(values) <- "double"
  invalid <- which(!is.finite(values), arr.ind = TRUE)
  if (nrow(invalid) > 0L) {
    header_rows <- if (has_header) 1L else 0L
    stop_cli(sprintf(
      "Numeric data are missing or nonfinite at input line %d, column %s.",
      invalid[1L, 1L] + header_rows,
      colnames(values)[invalid[1L, 2L]]
    ))
  }
  truncated <- FALSE
  if (observation_limit > 0L && nrow(values) > observation_limit) {
    values <- values[seq_len(observation_limit), , drop = FALSE]
    table <- table[seq_len(observation_limit), , drop = FALSE]
    truncated <- TRUE
  }
  row_labels <- if (first_data_column > 1L) {
    as.character(table[[1L]])
  } else {
    rep.int("", nrow(values))
  }
  list(values = values, names = colnames(values), labels = row_labels, truncated = truncated)
}

transform_vector <- function(values, transformation) {
  if (transformation %in% c("log", "log-diff") && any(values <= 0)) {
    stop("Logarithmic transformations require strictly positive observations.")
  }
  transformed <- switch(
    transformation,
    none = values,
    log = log(values),
    diff = diff(values),
    `log-diff` = diff(log(values))
  )
  if (length(transformed) < 8L) {
    stop("The transformation leaves fewer than eight observations for modeling.")
  }
  transformed
}

restore_forecasts <- function(forecast, original, transformation) {
  restored <- switch(
    transformation,
    none = forecast,
    log = exp(forecast),
    diff = tail(original, 1L) + cumsum(forecast),
    `log-diff` = exp(log(tail(original, 1L)) + cumsum(forecast))
  )
  if (any(!is.finite(restored))) {
    stop("Unable to restore forecasts to the original data scale.")
  }
  restored
}

profile_series <- function(values, frequency, max_lag) {
  lag_limit <- min(max_lag, length(values) - 1L)
  mean_value <- mean(values)
  level_acf <- as.numeric(stats::acf(values, lag.max = lag_limit, plot = FALSE)$acf)
  squared <- (values - mean_value)^2
  squared_acf <- as.numeric(stats::acf(squared, lag.max = lag_limit, plot = FALSE)$acf)
  partial <- as.numeric(stats::pacf(values, lag.max = lag_limit, plot = FALSE)$acf)
  centered_time <- seq_along(values) - 0.5 * (length(values) + 1)
  centered_value <- values - mean_value
  trend_strength <- abs(sum(centered_time * centered_value)) /
    sqrt(sum(centered_time^2) * max(sum(centered_value^2), .Machine$double.xmin))
  threshold <- 1.96 / sqrt(length(values))
  seasonal_strength <- 0
  seasonal_detected <- FALSE
  if (frequency > 1L && frequency <= lag_limit) {
    seasonal_strength <- abs(level_acf[frequency + 1L])
    seasonal_detected <- seasonal_strength > threshold
  }
  list(
    observations = length(values),
    frequency = frequency,
    max_lag = lag_limit,
    mean = mean_value,
    variance = stats::var(values),
    trend_strength = trend_strength,
    seasonal_strength = seasonal_strength,
    threshold = threshold,
    trend_detected = trend_strength > 0.35,
    seasonality_detected = seasonal_detected,
    autocorrelation_detected = max(abs(level_acf[-1L])) > threshold,
    conditional_variance_detected = max(abs(squared_acf[-1L])) > threshold,
    differencing_suggested = abs(level_acf[[2L]]) > 0.8,
    acf = level_acf,
    squared_acf = squared_acf,
    pacf = partial
  )
}

display_profile <- function(profile, display_lags) {
  cat("Univariate series profile\n")
  cat(sprintf("  observations: %d\n", profile$observations))
  cat(sprintf("  frequency: %d\n", profile$frequency))
  cat(sprintf("  mean: %14.6e\n", profile$mean))
  cat(sprintf("  variance: %14.6e\n", profile$variance))
  cat(sprintf("  trend strength: %8.4f\n", profile$trend_strength))
  cat(sprintf("  seasonal strength: %8.4f\n", profile$seasonal_strength))
  cat(sprintf("  trend detected: %s\n", profile$trend_detected))
  cat(sprintf("  seasonality detected: %s\n", profile$seasonality_detected))
  cat(sprintf("  level dependence detected: %s\n", profile$autocorrelation_detected))
  cat(sprintf(
    "  squared-value dependence detected: %s\n",
    profile$conditional_variance_detected
  ))
  cat(sprintf("  differencing suggested: %s\n\n", profile$differencing_suggested))
  cat("  lag          ACF     squared ACF          PACF\n")
  count <- min(display_lags, profile$max_lag)
  for (lag in seq_len(count)) {
    cat(sprintf(
      "  %3d%16.6f%16.6f%16.6f\n",
      lag, profile$acf[lag + 1L], profile$squared_acf[lag + 1L], profile$pacf[lag]
    ))
  }
  cat("\n")
}

gaussian_criteria <- function(residuals, parameter_count) {
  residuals <- residuals[is.finite(residuals)]
  n <- length(residuals)
  variance <- max(mean(residuals^2), .Machine$double.xmin)
  log_likelihood <- -0.5 * n * (log(2 * pi * variance) + 1)
  aic <- -2 * log_likelihood + 2 * parameter_count
  aicc <- if (n > parameter_count + 1L) {
    aic + 2 * parameter_count * (parameter_count + 1) /
      (n - parameter_count - 1)
  } else {
    Inf
  }
  bic <- -2 * log_likelihood + parameter_count * log(n)
  list(log_likelihood = log_likelihood, aic = aic, aicc = aicc, bic = bic,
       variance = variance, observations = n)
}

forecast_components <- function(object, horizon) {
  fitted_values <- tryCatch(as.numeric(stats::fitted(object)), error = function(e) numeric())
  residuals <- tryCatch(as.numeric(stats::residuals(object)), error = function(e) numeric())
  if (length(residuals) == 0L && length(fitted_values) > 0L) {
    residuals <- as.numeric(object$x) - fitted_values
  }
  list(
    forecast = as.numeric(object$mean)[seq_len(horizon)],
    fitted = fitted_values,
    residuals = residuals
  )
}

fit_ar_grid <- function(values, horizon, options) {
  fits <- list()
  selected <- NULL
  selected_score <- Inf
  maximum <- min(options$max_ar, max(1L, floor(length(values) / 3L)))
  for (order in seq_len(maximum)) {
    fit <- tryCatch(
      stats::ar.yw(values, order.max = order, aic = FALSE),
      error = function(error) NULL
    )
    if (is.null(fit)) {
      fits[[length(fits) + 1L]] <- list(order = order, converged = FALSE)
      next
    }
    prediction <- stats::predict(fit, n.ahead = horizon)$pred
    residuals <- as.numeric(fit$resid)
    criteria <- gaussian_criteria(residuals, order + 2L)
    score <- if (options$selection == "bic") criteria$bic else criteria$aicc
    named_coefficients <- setNames(fit$ar, paste0("phi(", seq_len(order), ")"))
    entry <- c(list(order = order, converged = TRUE,
                    coefficients = named_coefficients), criteria)
    fits[[length(fits) + 1L]] <- entry
    if (score < selected_score) {
      selected_score <- score
      selected <- list(
        forecast = as.numeric(prediction), residuals = residuals,
        parameters = c(mean = mean(values), named_coefficients),
        order = order, ma_order = 0L, order_fits = fits,
        parameter_count = order + 2L
      )
    }
  }
  if (is.null(selected)) stop("All autoregressive fits failed.")
  selected$order_fits <- fits
  selected
}

fit_arma_grid <- function(values, horizon, options) {
  fits <- list()
  selected <- NULL
  selected_score <- Inf
  for (p in 0:3) {
    for (q in 1:3) {
      if (p + q > 4L) next
      fit <- tryCatch(
        stats::arima(values, order = c(p, 0L, q), include.mean = TRUE, method = "ML"),
        error = function(error) NULL
      )
      if (is.null(fit)) {
        fits[[length(fits) + 1L]] <- list(p = p, q = q, converged = FALSE)
        next
      }
      prediction <- stats::predict(fit, n.ahead = horizon)$pred
      parameter_count <- length(fit$coef) + 1L
      criteria <- gaussian_criteria(as.numeric(stats::residuals(fit)), parameter_count)
      score <- if (options$selection == "bic") criteria$bic else criteria$aicc
      entry <- c(list(p = p, q = q, converged = fit$code == 0L,
                      coefficients = fit$coef), criteria)
      fits[[length(fits) + 1L]] <- entry
      if (score < selected_score) {
        selected_score <- score
        selected <- list(
          forecast = as.numeric(prediction), residuals = as.numeric(stats::residuals(fit)),
          parameters = fit$coef, order = p, ma_order = q,
          parameter_count = parameter_count, order_fits = fits
        )
      }
    }
  }
  if (is.null(selected)) stop("All ARMA fits failed.")
  selected$order_fits <- fits
  selected
}

fit_mean_candidate <- function(values, horizon, code, options) {
  require_package("forecast")
  if (code == "ar") return(fit_ar_grid(values, horizon, options))
  if (code == "arma") return(fit_arma_grid(values, horizon, options))
  series <- stats::ts(values, frequency = options$frequency)
  object <- switch(
    code,
    mean = forecast::meanf(series, h = horizon),
    naive = forecast::naive(series, h = horizon),
    drift = forecast::rwf(series, h = horizon, drift = TRUE),
    ses = forecast::ses(series, h = horizon),
    theta = forecast::thetaf(series, h = horizon),
    holt = forecast::holt(series, h = horizon),
    seasonal_naive = forecast::snaive(series, h = horizon),
    holt_winters = forecast::hw(series, h = horizon, seasonal = "additive")
  )
  components <- forecast_components(object, horizon)
  parameters <- numeric()
  if (!is.null(object$model$par)) parameters <- object$model$par
  if (code == "mean") parameters <- c(mean = mean(values))
  if (code == "naive") parameters <- c(last = tail(values, 1L))
  parameter_count <- max(2L, length(parameters) + 1L)
  c(components, list(
    parameters = parameters, order = 0L, ma_order = 0L,
    parameter_count = parameter_count, order_fits = list()
  ))
}

mean_model_name <- function(code) {
  switch(
    code,
    mean = "Mean",
    naive = "Naive",
    drift = "Random walk with drift",
    ses = "Simple exponential smoothing",
    theta = "Theta",
    holt = "Holt trend",
    seasonal_naive = "Seasonal naive",
    holt_winters = "Additive Holt-Winters",
    ar = "Yule-Walker autoregression",
    arma = "Gaussian ARMA"
  )
}

fit_mean_models <- function(values, options) {
  profile <- profile_series(values, options$frequency, options$max_lag)
  codes <- c("mean", "naive", "drift", "ses", "theta")
  if (profile$trend_detected) codes <- c(codes, "holt")
  if (profile$autocorrelation_detected) {
    codes <- c(codes, "ar")
    if (!profile$differencing_suggested) codes <- c(codes, "arma")
  }
  if (options$frequency > 1L && profile$seasonality_detected &&
      length(values) >= 2L * options$frequency) {
    codes <- c(codes, "seasonal_naive", "holt_winters")
  }
  validation_size <- 0L
  if (options$selection == "validation") {
    validation_size <- options$validation
    if (validation_size <= 0L) {
      validation_size <- max(4L, min(24L, floor(length(values) / 5L)))
    }
    validation_size <- min(validation_size, length(values) - 6L)
  }
  candidates <- vector("list", length(codes))
  full_fits <- vector("list", length(codes))
  for (index in seq_along(codes)) {
    code <- codes[[index]]
    start <- proc.time()[[3L]]
    fit <- tryCatch(
      fit_mean_candidate(values, options$horizon, code, options),
      error = function(error) structure(list(message = conditionMessage(error)), class = "fit_error")
    )
    elapsed <- proc.time()[[3L]] - start
    full_fits[[index]] <- fit
    if (inherits(fit, "fit_error")) {
      candidates[[index]] <- list(code = code, name = mean_model_name(code), info = 1L,
                                  message = fit$message, full_time = elapsed)
      next
    }
    criteria <- gaussian_criteria(fit$residuals, fit$parameter_count)
    candidate <- c(list(
      code = code, name = mean_model_name(code), info = 0L,
      parameters = fit$parameters, order = fit$order, ma_order = fit$ma_order,
      order_fits = fit$order_fits, full_time = elapsed,
      validation_time = 0, rmse = Inf, mae = Inf
    ), criteria)
    if (options$selection == "validation") {
      training <- values[seq_len(length(values) - validation_size)]
      actual <- tail(values, validation_size)
      start <- proc.time()[[3L]]
      validation_fit <- tryCatch(
        fit_mean_candidate(training, validation_size, code, options),
        error = function(error) NULL
      )
      candidate$validation_time <- proc.time()[[3L]] - start
      if (!is.null(validation_fit)) {
        errors <- actual - validation_fit$forecast
        candidate$rmse <- sqrt(mean(errors^2))
        candidate$mae <- mean(abs(errors))
      }
    }
    candidates[[index]] <- candidate
  }
  score <- vapply(candidates, function(candidate) {
    if (candidate$info != 0L) return(Inf)
    switch(options$selection, validation = candidate$rmse,
           aicc = candidate$aicc, bic = candidate$bic)
  }, numeric(1L))
  ordering <- order(score)
  list(
    profile = profile,
    candidates = candidates[ordering],
    selected_fit = full_fits[[ordering[[1L]]]],
    selected_code = candidates[[ordering[[1L]]]]$code,
    validation_size = validation_size
  )
}

fit_constant_variance <- function(values, horizon) {
  mean_value <- mean(values)
  variance <- max(mean((values - mean_value)^2), .Machine$double.xmin)
  log_likelihood <- -0.5 * length(values) * (log(2 * pi * variance) + 1)
  list(
    mean_forecast = rep(mean_value, horizon),
    variance_forecast = rep(variance, horizon),
    parameters = c(mean = mean_value, variance = variance),
    log_likelihood = log_likelihood,
    parameter_count = 2L,
    optimizer_converged = TRUE
  )
}

rugarch_specification <- function(code) {
  variance_order <- switch(code, arch = c(1L, 0L), c(1L, 1L))
  distribution <- if (code == "garch_student") "std" else "norm"
  rugarch::ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = variance_order),
    mean.model = list(armaOrder = c(0L, 0L), include.mean = TRUE),
    distribution.model = distribution
  )
}

fit_rugarch_candidate <- function(values, horizon, code) {
  require_package("rugarch")
  specification <- rugarch_specification(code)
  fit <- rugarch::ugarchfit(
    spec = specification,
    data = values,
    solver = "hybrid",
    solver.control = list(trace = 0)
  )
  forecast <- rugarch::ugarchforecast(fit, n.ahead = horizon)
  coefficients <- stats::coef(fit)
  list(
    mean_forecast = as.numeric(rugarch::fitted(forecast)),
    variance_forecast = as.numeric(rugarch::sigma(forecast))^2,
    parameters = coefficients,
    log_likelihood = as.numeric(rugarch::likelihood(fit)),
    parameter_count = length(coefficients),
    optimizer_converged = fit@fit$convergence == 0L
  )
}

fit_volatility_candidate <- function(values, horizon, code) {
  if (code == "constant") fit_constant_variance(values, horizon) else
    fit_rugarch_candidate(values, horizon, code)
}

volatility_model_name <- function(code) {
  switch(
    code,
    constant = "Gaussian constant variance",
    arch = "Gaussian ARCH(1)",
    garch_normal = "Gaussian GARCH(1,1)",
    garch_student = "Student-t GARCH(1,1)"
  )
}

fit_volatility_models <- function(values, options) {
  codes <- c("constant", "arch", "garch_normal", "garch_student")
  validation_size <- 0L
  if (options$selection == "validation") {
    validation_size <- options$validation
    if (validation_size <= 0L) {
      validation_size <- max(10L, min(50L, floor(length(values) / 5L)))
    }
    validation_size <- min(validation_size, length(values) - 20L)
  }
  candidates <- vector("list", length(codes))
  full_fits <- vector("list", length(codes))
  for (index in seq_along(codes)) {
    code <- codes[[index]]
    start <- proc.time()[[3L]]
    fit <- tryCatch(
      fit_volatility_candidate(values, options$horizon, code),
      error = function(error) structure(list(message = conditionMessage(error)), class = "fit_error")
    )
    elapsed <- proc.time()[[3L]] - start
    full_fits[[index]] <- fit
    if (inherits(fit, "fit_error")) {
      candidates[[index]] <- list(code = code, name = volatility_model_name(code),
                                  info = 1L, message = fit$message, full_time = elapsed)
      next
    }
    k <- fit$parameter_count
    n <- length(values)
    aic <- -2 * fit$log_likelihood + 2 * k
    aicc <- if (n > k + 1L) aic + 2 * k * (k + 1) / (n - k - 1L) else Inf
    bic <- -2 * fit$log_likelihood + k * log(n)
    candidate <- list(
      code = code, name = volatility_model_name(code), info = 0L,
      parameters = fit$parameters, optimizer_converged = fit$optimizer_converged,
      log_likelihood = fit$log_likelihood, aic = aic, aicc = aicc, bic = bic,
      full_time = elapsed, validation_time = 0, qlike = Inf
    )
    if (options$selection == "validation") {
      training <- values[seq_len(n - validation_size)]
      actual <- tail(values, validation_size)
      start <- proc.time()[[3L]]
      validation_fit <- tryCatch(
        fit_volatility_candidate(training, validation_size, code),
        error = function(error) NULL
      )
      candidate$validation_time <- proc.time()[[3L]] - start
      if (!is.null(validation_fit)) {
        variance <- pmax(validation_fit$variance_forecast, .Machine$double.xmin)
        errors <- actual - validation_fit$mean_forecast
        candidate$qlike <- mean(log(variance) + errors^2 / variance)
      }
    }
    candidates[[index]] <- candidate
  }
  score <- vapply(candidates, function(candidate) {
    if (candidate$info != 0L) return(Inf)
    switch(options$selection, validation = candidate$qlike,
           aicc = candidate$aicc, bic = candidate$bic)
  }, numeric(1L))
  ordering <- order(score)
  list(
    candidates = candidates[ordering],
    selected_fit = full_fits[[ordering[[1L]]]],
    selected_code = candidates[[ordering[[1L]]]]$code,
    validation_size = validation_size
  )
}

fixed_number <- function(value) {
  formatC(value, format = "f", digits = 8L, width = 16L)
}

display_parameters <- function(parameters, indent = "       ") {
  if (length(parameters) == 0L) return(invisible(NULL))
  for (name in names(parameters)) {
    cat(sprintf("%s%s %s\n", indent, name, fixed_number(parameters[[name]])))
  }
}

display_order_fits <- function(fits, type, print_parameters) {
  if (length(fits) == 0L) return(invisible(NULL))
  cat(sprintf("\n%s order search\n", toupper(type)))
  for (fit in fits) {
    if (!isTRUE(fit$converged)) {
      if (type == "ar") cat(sprintf("  %2d  failed\n", fit$order)) else
        cat(sprintf("  %2d %2d  failed\n", fit$p, fit$q))
      next
    }
    if (type == "ar") {
      cat(sprintf("  %2d  logLik %12.4e  AICc %12.4e  BIC %12.4e\n",
                  fit$order, fit$log_likelihood, fit$aicc, fit$bic))
    } else {
      cat(sprintf("  %2d %2d  logLik %12.4e  AICc %12.4e  BIC %12.4e\n",
                  fit$p, fit$q, fit$log_likelihood, fit$aicc, fit$bic))
    }
    if (print_parameters) display_parameters(fit$coefficients, "         ")
  }
}

display_mean_result <- function(result, options) {
  display_profile(result$profile, options$display_lags)
  cat("Automatic model comparison\n")
  cat(sprintf("  selection criterion: %s\n", options$selection))
  if (options$selection == "validation") {
    cat(sprintf("  validation observations: %d\n", result$validation_size))
  }
  display_count <- length(result$candidates)
  if (options$max_models > 0L) display_count <- min(display_count, options$max_models)
  for (rank in seq_len(display_count)) {
    candidate <- result$candidates[[rank]]
    cat(sprintf("  %3d  %s\n", rank, candidate$name))
    if (candidate$info != 0L) {
      cat(sprintf("       fit failed: %s\n\n", candidate$message))
      next
    }
    if (options$selection == "validation") {
      cat(sprintf("       RMSE %12.4e  MAE %12.4e\n", candidate$rmse, candidate$mae))
    }
    cat(sprintf(
      "       logLik %12.4e  AIC %12.4e  AICc %12.4e  BIC %12.4e\n",
      candidate$log_likelihood, candidate$aic, candidate$aicc, candidate$bic
    ))
    if (options$time_fits) {
      cat(sprintf("       full-data fit time %10.4f seconds\n", candidate$full_time))
      if (options$selection == "validation") {
        cat(sprintf("       validation fit time %10.4f seconds\n", candidate$validation_time))
      }
    }
    if (options$print_parameters) display_parameters(candidate$parameters)
    cat("\n")
  }
  if (options$print_all_ar) {
    candidate <- Filter(function(x) x$code == "ar", result$candidates)
    if (length(candidate) > 0L) {
      display_order_fits(candidate[[1L]]$order_fits, "ar", options$print_parameters)
    }
  }
  if (options$print_all_arma) {
    candidate <- Filter(function(x) x$code == "arma", result$candidates)
    if (length(candidate) > 0L) {
      display_order_fits(candidate[[1L]]$order_fits, "arma", options$print_parameters)
    }
  }
  cat(sprintf("Selected model: %s\n", mean_model_name(result$selected_code)))
  forecast <- result$selected_fit$forecast
  count <- min(length(forecast), options$print_forecasts)
  if (count > 0L) {
    cat("Forecasts\n")
    for (index in seq_len(count)) cat(sprintf("  %5d  %16.8e\n", index, forecast[[index]]))
  }
}

display_volatility_result <- function(result, options) {
  cat("Automatic volatility model comparison\n")
  cat(sprintf("  selection criterion: %s\n", options$selection))
  if (options$selection == "validation") {
    cat(sprintf("  validation observations: %d\n", result$validation_size))
  }
  display_count <- length(result$candidates)
  if (options$max_models > 0L) display_count <- min(display_count, options$max_models)
  for (rank in seq_len(display_count)) {
    candidate <- result$candidates[[rank]]
    cat(sprintf("  %3d  %s\n", rank, candidate$name))
    if (candidate$info != 0L) {
      cat(sprintf("       fit failed: %s\n\n", candidate$message))
      next
    }
    if (options$selection == "validation") {
      cat(sprintf("       validation QLIKE %12.4e\n", candidate$qlike))
    }
    cat(sprintf(
      "       logLik %12.4e  AIC %12.4e  AICc %12.4e  BIC %12.4e\n",
      candidate$log_likelihood, candidate$aic, candidate$aicc, candidate$bic
    ))
    cat(sprintf("       optimizer converged %s\n", candidate$optimizer_converged))
    if (options$time_fits) {
      cat(sprintf("       full-data fit time %10.4f seconds\n", candidate$full_time))
      if (options$selection == "validation") {
        cat(sprintf("       validation fit time %10.4f seconds\n", candidate$validation_time))
      }
    }
    if (options$print_parameters) display_parameters(candidate$parameters)
    cat("\n")
  }
  cat(sprintf(
    "Selected volatility model: %s\n",
    volatility_model_name(result$selected_code)
  ))
  forecast <- sqrt(pmax(result$selected_fit$variance_forecast, 0))
  count <- min(length(forecast), options$print_forecasts)
  if (count > 0L) {
    cat("Conditional standard-deviation forecasts\n")
    for (index in seq_len(count)) cat(sprintf("  %5d  %16.8e\n", index, forecast[[index]]))
  }
}

main <- function() {
  start <- proc.time()[[3L]]
  options <- parse_options(commandArgs(trailingOnly = TRUE))
  if (options$target == "mean") require_package("forecast") else require_package("rugarch")
  input <- read_input_table(options$file, options$observations)
  cat(sprintf(
    "Read %d valid rows and %d numeric columns.\n",
    nrow(input$values), ncol(input$values)
  ))
  if (input$truncated) {
    cat(sprintf(
      "Input truncated at the requested limit of %d observations.\n",
      options$observations
    ))
  }
  cat(sprintf("Modeling target: %s\n", options$target))
  if (options$resamples > 0L) {
    set.seed(options$seed)
    cat(sprintf("IID bootstrap resamples: %d\n", options$resamples))
    cat(sprintf("Random-number seed: %d\n", options$seed))
  }
  for (stride in options$strides) {
    cat("============================================================\n")
    cat(sprintf("Stride analysis: %d\n", stride))
    retained_rows <- seq.int(1L, nrow(input$values), by = stride)
    transformed_columns <- lapply(seq_len(ncol(input$values)), function(column) {
      transform_vector(input$values[retained_rows, column], options$transform)
    })
    transformed <- do.call(cbind, transformed_columns)
    colnames(transformed) <- input$names
    replicate_count <- max(1L, options$resamples)
    for (replicate in seq_len(replicate_count)) {
      analyzed <- transformed
      if (options$resamples > 0L) {
        cat(sprintf("\nIID resample: %d\n", replicate))
        indices <- sample.int(nrow(transformed), nrow(transformed), replace = TRUE)
        analyzed <- transformed[indices, , drop = FALSE]
      }
      if (options$corr) {
        cat("\nTransformed-data correlation matrix\n")
        print(round(stats::cor(analyzed), 6L))
      }
      for (column in seq_len(ncol(analyzed))) {
        cat("\n############################################################\n")
        cat(sprintf("Column analysis: %d (%s)\n", column, input$names[[column]]))
        cat(sprintf("Retained observations: %d\n", length(retained_rows)))
        cat(sprintf("Transformation: %s\n", options$transform))
        cat(sprintf("Modeling observations: %d\n", nrow(analyzed)))
        if (options$target == "mean") {
          cat("Forecast scale: original observations\n")
          if (options$transform %in% c("log", "log-diff")) {
            cat("Back-transformed log forecasts are median forecasts on the original scale.\n")
          }
          result <- fit_mean_models(analyzed[, column], options)
          result$selected_fit$forecast <- restore_forecasts(
            result$selected_fit$forecast,
            input$values[retained_rows, column],
            options$transform
          )
          display_mean_result(result, options)
        } else {
          cat("Forecast scale: transformed-series conditional standard deviation\n")
          result <- fit_volatility_models(analyzed[, column], options)
          display_volatility_result(result, options)
        }
      }
    }
  }
  cat(sprintf("\nElapsed time: %12.3f seconds\n", proc.time()[[3L]] - start))
}

tryCatch(
  main(),
  error = function(error) stop_cli(conditionMessage(error))
)
