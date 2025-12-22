#' Prepare acceleration–speed data for A–S profiling
#'
#' @param x A data frame/tibble containing at least columns `speed` and `acc`.
#'          Optional columns: `activity_id`, `athlete_id`.
#' @param speed_threshold Minimum speed to keep (m/s). Default 0.
#' @param bin_width Speed bin width (m/s). Default 0.1.
#' @param envelope_n Number of top acceleration points per bin. Default 2.
#' @param keep_frac Fraction of raw points to keep for plotting. Default 0.5.
#' @param seed Optional integer seed for reproducible downsampling.
#' @param print_plot Logical; return a raw plot. Default TRUE.
#'
#' @return An object of class `as_prep` (a list with as_data, prepared_data,
#'         as_initial_lm_data, plot_initial).
#' @export
as_prepare_data <- function(x,
                            speed_threshold = 0,
                            bin_width = 0.1,
                            envelope_n = 2,
                            keep_frac = 0.5,
                            seed = NULL,
                            print_plot = TRUE) {

  # ---- checks
  req_cols <- c("speed", "acc")
  missing <- setdiff(req_cols, names(x))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }

  # Ensure optional ID cols exist (keeps downstream code stable)
  if (!"activity_id" %in% names(x)) x$activity_id <- NA_character_
  if (!"athlete_id"  %in% names(x)) x$athlete_id  <- NA_character_

  # ---- core filtering
  as_data <- x %>%
    dplyr::select(speed, acc, activity_id, athlete_id) %>%
    dplyr::filter(!is.na(speed), !is.na(acc)) %>%
    dplyr::filter(acc >= 0)

  prepared_data <- as_data %>%
    dplyr::filter(speed > speed_threshold)

  if (nrow(prepared_data) == 0) {
    warning("No data left after filtering. Check thresholds.")
    return(NULL)
  }

  # ---- binning
  max_speed <- max(prepared_data$speed, na.rm = TRUE)
  if (!is.finite(max_speed) || max_speed <= speed_threshold) {
    warning("Speed range too small to create bins.")
    return(NULL)
  }

  breaks <- seq(from = speed_threshold, to = max_speed, by = bin_width)
  if (length(breaks) < 2) {
    warning("Not enough breaks to bin data. Increase date range or lower speed_threshold.")
    return(NULL)
  }

  prepared_data <- prepared_data %>%
    dplyr::mutate(
      cuts = base::cut(speed, breaks = breaks, include.lowest = TRUE, right = TRUE)
    )

  # ---- envelope: top N acc per bin
  as_initial_lm_data <- prepared_data %>%
    dplyr::group_by(.data$cuts) %>%
    dplyr::slice_max(acc, n = envelope_n, with_ties = FALSE) %>%
    dplyr::ungroup()

  # ---- optional raw plot
  plot_initial <- NULL
  if (isTRUE(print_plot)) {
    if (!is.null(seed)) set.seed(seed)

    keep_frac <- max(0, min(1, keep_frac))
    keep_n <- max(1L, ceiling(nrow(as_data) * keep_frac))
    idx <- sample.int(nrow(as_data), size = keep_n)

    as_reduced <- as_data[idx, , drop = FALSE]

    plot_initial <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = as_reduced,
        mapping = ggplot2::aes(x = speed, y = acc),
        alpha = 0.25,
        size = 0.4
      ) +
      ggplot2::geom_point(
        data = as_initial_lm_data,
        mapping = ggplot2::aes(x = speed, y = acc),
        color = "#F8766D",
        size = 3,
        alpha = 0.9
      ) +
      ggplot2::geom_smooth(
        data = as_initial_lm_data,
        mapping = ggplot2::aes(x = speed, y = acc),
        method = "lm",
        se = TRUE
      ) +
      ggpubr::theme_pubr() +
      ggplot2::coord_cartesian(expand = FALSE) +
      ggplot2::ggtitle("Raw Acceleration–Speed Profile") +
      ggplot2::ylim(0, NA)
  }

  out <- list(
    as_data            = tibble::as_tibble(as_data),
    prepared_data      = tibble::as_tibble(prepared_data),
    as_initial_lm_data = tibble::as_tibble(as_initial_lm_data),
    plot_initial       = plot_initial
  )
  class(out) <- c("as_prep", class(out))
  out
}


#' Build A–S profile from prepared data
#'
#' @param prep Output of `as_prepare_data()`.
#' @param method Cleaning method: "ci" (default) or "box".
#' @param ci_level Confidence level for CI filtering. Default 0.95.
#' @param ci_interval Interval type passed to `predict()`: "confidence" or "prediction".
#'        Default "confidence".
#' @param seed Optional integer seed for reproducible downsampling.
#' @param keep_frac Fraction of raw points to keep for plotting. Default 0.5.
#' @param print_plot Logical; return a profile plot. Default TRUE.
#'
#' @return A list with cleaned data, models, metrics, and plot (if requested).
#' @export
as_profile <- function(prep,
                       method = c("ci", "box"),
                       ci_level = 0.95,
                       ci_interval = c("confidence", "prediction"),
                       seed = NULL,
                       keep_frac = 0.5,
                       print_plot = TRUE) {

  if (is.null(prep) || !inherits(prep, "as_prep")) {
    stop("`prep` must be the output of `as_prepare_data()`.")
  }

  method <- match.arg(method)
  ci_interval <- match.arg(ci_interval)

  as_data <- prep$as_data
  prepared_data <- prep$prepared_data
  as_initial_lm_data <- prep$as_initial_lm_data

  if (nrow(as_initial_lm_data) < 2) {
    warning("Too few points in as_initial_lm_data to fit a model.")
    return(NULL)
  }

  # ---- initial LM (for CI method)
  lm_initial <- stats::lm(acc ~ speed, data = as_initial_lm_data)
  summary_lm_initial <- summary(lm_initial)

  # ---- clean envelope points
  if (method == "ci") {
    predicted <- tibble::as_tibble(
      stats::predict(
        lm_initial,
        newdata = as_initial_lm_data,
        interval = ci_interval,
        level = ci_level
      )
    )

    # predicted data has columns: fit, lwr, upr
    as_clean <- dplyr::bind_cols(as_initial_lm_data, predicted) %>%
      dplyr::filter(acc >= lwr, acc <= upr) %>%
      dplyr::select(names(as_initial_lm_data)) %>%
      tibble::as_tibble()

  } else {
    # box / quantile + IQR method, 1 row per bin
    as_clean <- prepared_data %>%
      dplyr::group_by(.data$cuts) %>%
      dplyr::filter(
        acc >= stats::quantile(acc, 0.98, na.rm = TRUE) |
          speed >= stats::quantile(speed, 0.98, na.rm = TRUE)
      ) %>%
      dplyr::mutate(
        q1 = stats::quantile(acc, 0.25, na.rm = TRUE),
        q3 = stats::quantile(acc, 0.75, na.rm = TRUE),
        iqr = .data$q3 - .data$q1,
        lower_bound = .data$q1 - 1.5 * .data$iqr,
        upper_bound = .data$q3 + 1.5 * .data$iqr
      ) %>%
      dplyr::filter(acc >= .data$lower_bound, acc <= .data$upper_bound) %>%
      dplyr::summarise(
        acc = max(acc, na.rm = TRUE),
        speed = max(speed, na.rm = TRUE),
        activity_id = dplyr::first(activity_id),
        athlete_id  = dplyr::first(athlete_id),
        .groups = "drop"
      ) %>%
      dplyr::distinct(acc, .keep_all = TRUE) %>%
      tibble::as_tibble()
  }

  if (nrow(as_clean) < 2) {
    warning("Too few points after cleaning to fit lm(acc ~ speed).")
    return(list(
      as_clean = as_clean,
      lm_initial = lm_initial,
      summary_lm_initial = summary_lm_initial,
      plot_profile = NULL
    ))
  }

  # ---- final LM on cleaned data
  lm_clean <- stats::lm(acc ~ speed, data = as_clean)
  summary_lm_clean <- summary(lm_clean)

  a0 <- unname(stats::coef(lm_clean)[["(Intercept)"]])
  slope <- unname(stats::coef(lm_clean)[["speed"]])
  v0 <- a0 / abs(slope)
  r2 <- summary_lm_clean$r.squared

  # ---- optional plot
  plot_profile <- NULL
  if (isTRUE(print_plot)) {
    if (!is.null(seed)) set.seed(seed)

    keep_frac <- max(0, min(1, keep_frac))
    keep_n <- max(1L, ceiling(nrow(as_data) * keep_frac))
    idx <- sample.int(nrow(as_data), size = keep_n)

    as_reduced <- as_data[idx, , drop = FALSE]

    x_limit <- v0 + 0.5
    y_limit <- a0 + 1

    plot_profile <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = as_reduced,
        mapping = ggplot2::aes(x = speed, y = acc),
        alpha = 0.25,
        size = 0.4
      ) +
      ggplot2::geom_abline(intercept = a0, slope = slope, linewidth = 1) +
      ggplot2::geom_point(
        data = as_clean,
        mapping = ggplot2::aes(x = speed, y = acc),
        color = "#F8766D",
        size = 2
      ) +
      ggpubr::theme_pubr() +
      ggplot2::coord_cartesian(expand = FALSE, xlim = c(0, x_limit), ylim = c(0, y_limit)) +
      ggplot2::labs(
        title = "In-situ Acceleration–Speed Profile",
        x = "Sprinting velocity (m/s)",
        y = expression(Acceleration~(m/s^2))
      ) +
      ggplot2::geom_label(
        x = v0 - 4.6, y = a0 - 0.6, hjust = 0,
        aes(label = paste0(
          "Max acceleration (A0): ", round(a0, 2), " m/s²\n\n",
          "Max speed (S0): ", round(v0, 2), " m/s\n\n",
          "Slope: ", round(slope, 2), "\n\n",
          "Fit (R²): ", round(r2, 4)
        )
        )
      )
  }

  list(
    as_clean = as_clean,
    lm_initial = lm_initial,
    summary_lm_initial = summary_lm_initial,
    lm_clean = lm_clean,
    summary_lm_clean = summary_lm_clean,
    metrics = list(
      a0 = round(a0, 2),
      v0 = round(v0, 2),
      slope = round(slope, 2),
      r2 = round(r2, 4)
    ),
    plot_profile = plot_profile
  )
}
