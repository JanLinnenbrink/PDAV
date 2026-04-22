#' Compute distance-weighted CV weights
#'
#' Constructs calibration weights that match the distribution of prediction
#' distance `d` between validation tasks and deployment tasks.
#'
#' @param sample_tasks Validation-task data frame containing a continuous `d`
#'   column.
#' @param grid_tasks Deployment-task data frame containing a continuous `d`
#'   column.
#' @param balance_by Quantile spacing used for discretizing `d`.
#' @param shrink_lambda Shrinkage parameter applied after calibration.
#' @param verbose Verbosity level.
#'
#' @return A list containing calibration weights and metadata.
#' @noRd
compute_dwcv_weights <- function(sample_tasks, grid_tasks, balance_by = 0.2, shrink_lambda = 0, verbose = 0) {
	bal <- prepare_balanced_tasks_cv(
		loss_df = sample_tasks,
		grid_tasks = grid_tasks,
		balancing_vars = "d",
		by = balance_by
	)

	balance_df <- data.frame(d = bal$sample_tasks_bal$d_cat)

	target_margins <- compute_target_margins_generic(
		grid_tasks_bal = bal$grid_tasks_bal,
		balancing_vars = "d"
	)

	out <- rake_weights(
		balance_df = balance_df,
		target_margins = target_margins,
		verbose = verbose
	)

	out$weights_raw <- normalize_weights(out$weights)
	out$weights <- shrink_weights(out$weights_raw, lambda = shrink_lambda)
	out$shrink_lambda <- shrink_lambda
	out$balancing_vars <- "d"

	out
}

#' Compute target-weighted CV weights
#'
#' Constructs calibration weights that match the empirical margins of selected
#' task descriptors between validation tasks and deployment tasks.
#'
#' @param sample_tasks_bal Validation-task data frame with discretized balancing
#'   variables.
#' @param grid_tasks_bal Deployment-task data frame with discretized balancing
#'   variables.
#' @param balancing_vars Character vector of balancing-variable names.
#' @param shrink_lambda Shrinkage parameter applied after calibration.
#' @param verbose Verbosity level.
#'
#' @return A list containing calibration weights and metadata.
#' @noRd
compute_twcv_weights <- function(sample_tasks_bal, grid_tasks_bal, balancing_vars, shrink_lambda = 0, verbose = 0) {
	balance_df <- as.data.frame(
		lapply(balancing_vars, function(v) sample_tasks_bal[[paste0(v, "_cat")]])
	)
	names(balance_df) <- balancing_vars

	target_margins <- compute_target_margins_generic(
		grid_tasks_bal = grid_tasks_bal,
		balancing_vars = balancing_vars
	)

	out <- rake_weights(
		balance_df = balance_df,
		target_margins = target_margins,
		verbose = verbose
	)

	out$weights_raw <- normalize_weights(out$weights)
	out$weights <- shrink_weights(out$weights_raw, lambda = shrink_lambda)
	out$shrink_lambda <- shrink_lambda
	out$balancing_vars <- balancing_vars

	out
}


#' Compute target margins for calibration weighting
#'
#' Extracts empirical target-domain margins for discretized balancing variables.
#'
#' @param grid_tasks_bal Data frame of discretized deployment-task descriptors.
#' @param balancing_vars Character vector of balancing-variable names. For each
#'   variable `v`, the function expects a column named `paste0(v, "_cat")`.
#'
#' @return Named list of empirical target margins.
#' @noRd
compute_target_margins_generic <- function(grid_tasks_bal, balancing_vars) {
	out <- vector("list", length(balancing_vars))
	names(out) <- balancing_vars

	for (v in balancing_vars) {
		vn <- paste0(v, "_cat")
		if (!(vn %in% names(grid_tasks_bal))) {
			stop("Missing column in grid_tasks_bal: ", vn, call. = FALSE)
		}
		out[[v]] <- tabulate_proportions(grid_tasks_bal[[vn]])
	}

	out
}


#' Tabulate empirical category proportions
#'
#' Computes empirical proportions for an integer-coded categorical variable.
#'
#' @param x Vector of category labels coercible to integer.
#' @param levels Optional vector of category levels. If `NULL`, the sorted
#'   unique values of `x` are used.
#'
#' @return Numeric vector of category proportions summing to one.
#' @noRd
tabulate_proportions <- function(x, levels = NULL) {
	x <- as.integer(x)
	if (is.null(levels)) {
		levels <- sort(unique(x))
	}
	out <- table(factor(x, levels = levels))
	as.numeric(out) / sum(out)
}


#' Compute calibration weights by iterative proportional fitting
#'
#' Reweights validation tasks so that weighted empirical margins match target
#' margins for a set of discretized balancing variables.
#'
#' @param balance_df Data frame of discretized balancing variables.
#' @param target_margins Named list of target proportions for each balancing
#'   variable.
#' @param base_weights Optional numeric vector of starting weights.
#' @param max_iter Maximum number of raking iterations.
#' @param tol Convergence tolerance based on relative weight change.
#' @param verbose Verbosity level.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{weights}{Final calibration weights.}
#'     \item{converged}{Logical indicating whether convergence was reached.}
#'     \item{iterations}{Number of iterations performed.}
#'   }
#' @noRd
rake_weights <- function(balance_df, target_margins, base_weights = NULL, max_iter = 500, tol = 1e-6, verbose = 0) {
	n <- nrow(balance_df)
	if (is.null(base_weights)) {
		w <- rep(1, n)
	} else {
		w <- as.numeric(base_weights)
	}

	margin_names <- names(target_margins)

	for (iter in seq_len(max_iter)) {
		w_old <- w

		for (m in margin_names) {
			x <- as.integer(balance_df[[m]])
			levs <- seq_along(target_margins[[m]])
			target_prop <- target_margins[[m]]

			current_totals <- tapply(w, factor(x, levels = levs), sum)
			current_totals[is.na(current_totals)] <- 0

			target_totals <- sum(w) * target_prop

			adj <- rep(1, length(levs))
			ok <- current_totals > 0

			adj[ok] <- target_totals[ok] / current_totals[ok]
			adj[!ok] <- NA_real_

			valid <- !is.na(adj[x])
			w[valid] <- w[valid] * adj[x[valid]]
		}

		rel_change <- max(abs(w - w_old) / pmax(abs(w_old), 1e-12))

		if (verbose >= 2) {
			log_message(verbose, 2, "  raking iter ", iter, ", max relative change = ", format(rel_change, digits = 3))
		}

		if (rel_change < tol) break
	}

	converged <- iter < max_iter || rel_change < tol

	list(
		weights = w,
		converged = converged,
		iterations = iter
	)
}


#' Normalize weights to unit mean
#'
#' @param w Numeric weight vector.
#'
#' @return Numeric vector with mean equal to one.
#' @noRd
normalize_weights <- function(w) {
	w / mean(w)
}


#' Shrink weights towards uniform weights
#'
#' Applies convex shrinkage to reduce weight dispersion.
#'
#' @param w Numeric weight vector.
#' @param lambda Shrinkage parameter in `[0, 1]`. `0` leaves weights unchanged
#'   after normalization; `1` yields uniform weights.
#'
#' @return Numeric vector of shrunk weights with mean equal to one.
#' @noRd
shrink_weights <- function(w, lambda = 0) {
	stopifnot(lambda >= 0, lambda <= 1)
	w <- normalize_weights(w)
	(1 - lambda) * w + lambda
}
