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
	# Extracts the balancing variables from the training points
	balance_df <- as.data.frame(
		lapply(balancing_vars, function(v) sample_tasks_bal[[paste0(v, "_cat")]])
	)
	names(balance_df) <- balancing_vars

	# Calculates proportion of predpoints in each quantile of each predictor used for weighting
	target_margins <- compute_target_margins_generic(
		grid_tasks_bal = grid_tasks_bal,
		balancing_vars = balancing_vars
	)

	# applies iterative proportional fitting ("raking")
	out <- rake_weights(
		balance_df = balance_df,
		target_margins = target_margins,
		verbose = verbose
	)

	# Weights are normalized by their mean and shrinked towards 1 to mitigate extreme values
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

			# calculate the number of training points currently in each quantile
			# uses the weights from the previous balancing variable
			# -> already weighted, but likely not ideally for this variable
			current_totals <- tapply(w, factor(x, levels = levs), sum)
			current_totals[is.na(current_totals)] <- 0

			# calculate the desired number of training points in each quantile that matches the target margins:
			# Number of training points * target proportion vector
			target_totals <- sum(w) * target_prop

			adj <- rep(1, length(levs))
			ok <- current_totals > 0

			# calculates the weight needed to adjust the training point distribution to the target margins
			# (weights > 1 are used to up-weigh underrepresented classes, weights < 1 to down-weight over-represented ones)
			adj[ok] <- target_totals[ok] / current_totals[ok]
			adj[!ok] <- NA_real_

			valid <- !is.na(adj[x])
			w[valid] <- w[valid] * adj[x[valid]]
		}

		# Compute the relative strength of the absolute change of weights from base (or previous) to new weights
		rel_change <- max(abs(w - w_old) / pmax(abs(w_old), 1e-12))

		if (verbose >= 2) {
			log_message(verbose, 2, "  raking iter ", iter, ", max relative change = ", format(rel_change, digits = 3))
		}

		# When the changes converge (i.e., when the weight difference from previous iteration to new one is small), stop and return weights
		# (when the weights for predictor B are changed, weights for predictor A might be off again, and another iteration starts,
		# until the diff between them is small)
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


#' Prepare balanced task representations for CV
#'
#' Constructs balanced representations of sample and grid task descriptors
#' for weighting-based estimators such as TWCV.
#'
#' @param loss_df Data frame of sample-level losses.
#' @param grid_tasks Data frame of grid-level task descriptors.
#' @param balancing_vars Character vector of variables used for balancing.
#' @param by Numeric step size for quantile binning.
#'
#' @return List with \code{sample_tasks_bal} and \code{grid_tasks_bal}.
#' @noRd
prepare_balanced_tasks_cv <- function(loss_df, grid_tasks, balancing_vars, by = 0.2) {
	stopifnot(length(balancing_vars) > 0)

	sample_tasks_bal <- prepare_for_balancing(
		df = loss_df,
		vars = balancing_vars,
		ref_df = grid_tasks,
		by = by
	)

	grid_tasks_bal <- prepare_for_balancing(
		df = grid_tasks,
		vars = balancing_vars,
		ref_df = grid_tasks,
		by = by
	)

	list(
		sample_tasks_bal = sample_tasks_bal,
		grid_tasks_bal = grid_tasks_bal
	)
}


#' Prepare variables for balancing via discretization
#'
#' Transforms variables into categorical representations suitable for
#' balancing. Numeric variables are discretized using quantiles of a
#' reference distribution; categorical variables are aligned to the
#' reference support.
#'
#' @param df Data frame to transform.
#' @param vars Variables to transform.
#' @param ref_df Reference data frame.
#' @param by Quantile step size.
#'
#' @return Modified data frame with additional \code{*_cat} variables.
#' @noRd
prepare_for_balancing <- function(df, vars, ref_df, by = 0.2) {
	df_out <- df

	for (v in vars) {
		if (!(v %in% names(df))) {
			stop("Variable '", v, "' not found in df.", call. = FALSE)
		}
		if (!(v %in% names(ref_df))) {
			stop("Variable '", v, "' not found in ref_df.", call. = FALSE)
		}

		x_ref <- ref_df[[v]]
		x <- df[[v]]
		out_name <- paste0(v, "_cat")

		# numeric variables
		if (is.numeric(x_ref) && length(unique(stats::na.omit(x_ref))) > 2) {
			probs <- seq(0, 1, by = by)
			qtiles <- stats::quantile(x_ref, probs = probs, na.rm = TRUE, names = FALSE)

			qtiles <- unique(qtiles)

			# degenerate case
			if (length(qtiles) < 2) {
				levs <- paste0(v, "_Q1")
				df_out[[out_name]] <- factor(
					ifelse(is.na(x), NA, levs),
					levels = levs,
					ordered = TRUE
				)
				next
			}

			qtiles[1] <- -Inf
			qtiles[length(qtiles)] <- Inf

			n_bins <- length(qtiles) - 1L
			levs <- paste0(v, "_Q", seq_len(n_bins))

			df_out[[out_name]] <- cut(
				x,
				breaks = qtiles,
				include.lowest = TRUE,
				labels = levs,
				ordered_result = TRUE
			)
		} else {
			ref_levels <- sort(unique(as.character(stats::na.omit(x_ref))))
			x_chr <- as.character(x)

			x_chr[!(x_chr %in% ref_levels) & !is.na(x_chr)] <- NA_character_

			df_out[[out_name]] <- factor(
				x_chr,
				levels = ref_levels,
				ordered = FALSE
			)
		}
	}

	df_out
}

#' Checks if all quintiles of the prediction points are supported by the training points.
#' Otherwise, raking is prone to errors.
#'
#' @param balance_df data frame containing the training margins.
#' @param target_margins data frame containing the prediction point margins.
#' @param eps Tolerance
#'
#' @return data frame containing the predictors, their quintiles and information if they are supported.
#' @noRd
check_balance_support <- function(balance_df, target_margins, eps = 1e-12) {
	out <- lapply(names(target_margins), function(m) {
		levs <- seq_along(target_margins[[m]])

		sample_counts <- table(
			factor(as.integer(balance_df[[m]]), levels = levs)
		)

		data.frame(
			var = m,
			level = levs,
			sample_n = as.numeric(sample_counts),
			target_prop = as.numeric(target_margins[[m]]),
			unsupported = as.numeric(sample_counts) == 0 &
				as.numeric(target_margins[[m]]) > eps
		)
	})

	dplyr::bind_rows(out)
}
