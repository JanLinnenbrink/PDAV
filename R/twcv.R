#' Tabulate empirical category proportions
#'
#' Computes empirical proportions for an integer-coded categorical variable.
#'
#' @param x Vector of category labels coercible to integer.
#' @param levels Optional vector of category levels. If `NULL`, the sorted
#'   unique values of `x` are used.
#'
#' @return Numeric vector of category proportions summing to one.
tabulate_proportions <- function(x, levels = NULL) {
	x <- as.integer(x)
	if (is.null(levels)) {
		levels <- sort(unique(x))
	}
	out <- table(factor(x, levels = levels))
	as.numeric(out) / sum(out)
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
compute_target_margins <- function(grid_tasks_bal, balancing_vars) {
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
#'
#' @return A list with elements:
#'   \describe{
#'     \item{weights}{Final calibration weights.}
#'     \item{converged}{Logical indicating whether convergence was reached.}
#'     \item{iterations}{Number of iterations performed.}
#'   }
rake_weights <- function(
	balance_df,
	target_margins,
	base_weights = NULL,
	max_iter = 500,
	tol = 1e-6
) {
	n <- nrow(balance_df)
	w <- if (is.null(base_weights)) rep(1, n) else as.numeric(base_weights)

	rel_change <- Inf

	for (iter in seq_len(max_iter)) {
		w_old <- w

		for (m in names(target_margins)) {
			x <- as.integer(balance_df[[m]])
			levs <- seq_along(target_margins[[m]])
			target_prop <- target_margins[[m]]

			current_totals <- tapply(w, factor(x, levels = levs), sum)
			current_totals[is.na(current_totals)] <- 0

			target_totals <- sum(w) * target_prop

			adj <- rep(NA_real_, length(levs))
			ok <- current_totals > 0
			adj[ok] <- target_totals[ok] / current_totals[ok]

			valid <- !is.na(adj[x])
			w[valid] <- w[valid] * adj[x[valid]]
		}

		rel_change <- max(abs(w - w_old) / pmax(abs(w_old), 1e-12))
		if (rel_change < tol) {
			break
		}
	}

	list(
		weights = w,
		converged = (iter < max_iter) || (rel_change < tol),
		iterations = iter
	)
}


#' Normalize weights to unit mean
#'
#' @param w Numeric weight vector.
#'
#' @return Numeric vector with mean equal to one.
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
shrink_weights <- function(w, lambda = 0) {
	stopifnot(lambda >= 0, lambda <= 1)
	w <- normalize_weights(w)
	(1 - lambda) * w + lambda
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
#'
#' @return A list containing calibration weights and metadata.
compute_twcv_weights <- function(
	sample_tasks_bal,
	grid_tasks_bal,
	balancing_vars,
	shrink_lambda = 0
) {
	sample_cols <- paste0(balancing_vars, "_cat")
	missing_sample <- setdiff(sample_cols, names(sample_tasks_bal))
	missing_grid <- setdiff(sample_cols, names(grid_tasks_bal))

	if (length(missing_sample) > 0) {
		stop(
			"Missing columns in sample_tasks_bal: ",
			paste(missing_sample, collapse = ", "),
			call. = FALSE
		)
	}
	if (length(missing_grid) > 0) {
		stop(
			"Missing columns in grid_tasks_bal: ",
			paste(missing_grid, collapse = ", "),
			call. = FALSE
		)
	}

	balance_df <- as.data.frame(
		lapply(sample_cols, function(v) sample_tasks_bal[[v]])
	)
	names(balance_df) <- balancing_vars

	target_margins <- compute_target_margins(
		grid_tasks_bal = grid_tasks_bal,
		balancing_vars = balancing_vars
	)

	out <- rake_weights(
		balance_df = balance_df,
		target_margins = target_margins
	)

	out$weights_raw <- normalize_weights(out$weights)
	out$weights <- shrink_weights(out$weights_raw, lambda = shrink_lambda)
	out$shrink_lambda <- shrink_lambda
	out$balancing_vars <- balancing_vars

	out
}


#' Create a TWCV specification
#'
#' Defines one TWCV weighting configuration, including balancing variables,
#' discretization level, and shrinkage strength.
#'
#' @param predictor_vars Character vector of predictor variables.
#' @param include_distance Logical; include prediction distance `d`.
#' @param balance_by Quantile spacing used for discretization.
#' @param shrink_lambda Shrinkage parameter for calibration weights.
#' @param name Name of the specification.
#'
#' @return Named list with one TWCV specification.
#' @export
make_twcv_specs <- function(
	predictor_vars,
	include_distance = TRUE,
	balance_by = 0.2,
	shrink_lambda = 0.2,
	name = "twcv_extended"
) {
	balancing_vars <- predictor_vars
	if (include_distance) {
		balancing_vars <- c(balancing_vars, "d")
	}

	specs <- list()
	specs[[name]] <- list(
		balancing_vars = balancing_vars,
		balance_by = balance_by,
		shrink_lambda = shrink_lambda
	)
	specs
}


#' Create a standard TWCV specification set
#'
#' @param predictor_vars Character vector of focal predictor variables.
#' @param env_vars Character vector of environmental variables for an extended
#'   specification.
#' @param include_distance Logical; include prediction distance `d`.
#' @param balance_by Quantile spacing used for discretization.
#' @param shrink_lambda Shrinkage parameter for calibration weights.
#'
#' @return Named list of TWCV specifications.
make_default_twcv_spec_set <- function(
	predictor_vars,
	env_vars = predictor_vars,
	include_distance = TRUE,
	balance_by = 0.2,
	shrink_lambda = 0.2
) {
	c(
		make_twcv_specs(
			predictor_vars = predictor_vars,
			include_distance = include_distance,
			balance_by = balance_by,
			shrink_lambda = shrink_lambda,
			name = "twcv"
		),
		make_twcv_specs(
			predictor_vars = env_vars,
			include_distance = include_distance,
			balance_by = balance_by,
			shrink_lambda = shrink_lambda,
			name = "twcv_extended"
		)
	)
}


#' Apply TWCV weighting to validation losses
#'
#' Applies one or more TWCV weighting specifications to a validation-loss data
#' frame and summarizes predictive performance under each weighted estimator.
#'
#' This function assumes that validation losses and task descriptors have
#' already been computed externally.
#'
#' @param loss_df Validation-loss data frame augmented with continuous task
#'   descriptors.
#' @param grid_tasks Deployment-task data frame with continuous task descriptors.
#' @param twcv_specs Named list of TWCV specifications.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{estimators}{Named list of performance summaries.}
#'     \item{weights}{Named list of TWCV weight objects.}
#'     \item{twcv_specs}{The TWCV specification set actually used.}
#'   }
apply_twcv_estimators <- function(loss_df, grid_tasks, twcv_specs) {
	est_list <- list(
		unweighted = summarize_losses(loss_df)
	)

	weight_objects <- list()

	for (nm in names(twcv_specs)) {
		spec <- twcv_specs[[nm]]

		if (
			is.null(spec$balancing_vars) ||
				is.null(spec$balance_by) ||
				is.null(spec$shrink_lambda)
		) {
			stop("Malformed twcv spec for '", nm, "'.", call. = FALSE)
		}

		bal <- prepare_balanced_tasks_cv(
			loss_df = loss_df,
			grid_tasks = grid_tasks,
			balancing_vars = spec$balancing_vars,
			by = spec$balance_by
		)

		tw <- compute_twcv_weights(
			sample_tasks_bal = bal$sample_tasks_bal,
			grid_tasks_bal = bal$grid_tasks_bal,
			balancing_vars = spec$balancing_vars,
			shrink_lambda = spec$shrink_lambda
		)

		est_list[[nm]] <- summarize_losses(loss_df, tw$weights)
		weight_objects[[nm]] <- tw
	}

	list(
		estimators = est_list,
		weights = weight_objects,
		twcv_specs = twcv_specs
	)
}


#' Gini coefficient of a weight vector
#'
#' Computes the Gini coefficient as a measure of inequality for a vector
#' of non-negative weights.
#'
#' @param w Numeric vector of weights.
#'
#' @return Numeric scalar giving the Gini coefficient, or `NA` if undefined.
gini_coefficient <- function(w) {
	w <- as.numeric(w)
	w <- w[is.finite(w) & w >= 0]
	if (length(w) == 0) {
		return(NA_real_)
	}
	if (sum(w) == 0) {
		return(NA_real_)
	}

	w <- sort(w)
	n <- length(w)
	idx <- seq_len(n)

	(2 * sum(idx * w) / (n * sum(w))) - (n + 1) / n
}


#' Compute nearest-neighbor distances between point sets
#'
#' Computes Euclidean nearest-neighbor distances from query points to reference
#' points based on coordinate columns `x` and `y`. When `exclude_self = TRUE`
#' and both point sets are identical, the distance to the second-nearest
#' neighbor is returned.
#'
#' @param query_coords Data frame, matrix, or similar object containing point
#'   coordinates in columns `x` and `y`.
#' @param ref_coords Data frame, matrix, or similar object containing reference
#'   coordinates in columns `x` and `y`.
#' @param exclude_self Logical; if `TRUE`, self-matches are excluded when query
#'   and reference coordinates are identical.
#'
#' @return A numeric vector of nearest-neighbor distances, one per query point.
nearest_neighbor_distance <- function(query_coords, ref_coords, exclude_self = FALSE) {
	q <- as.matrix(query_coords[, c("x", "y")])
	r <- as.matrix(ref_coords[, c("x", "y")])

	if (!exclude_self) {
		nn <- FNN::get.knnx(data = r, query = q, k = 1)
		return(as.numeric(nn$nn.dist[, 1]))
	}

	same_object <- nrow(q) == nrow(r) &&
		isTRUE(all.equal(q, r, check.attributes = FALSE))

	if (!same_object) {
		nn <- FNN::get.knnx(data = r, query = q, k = 1)
		return(as.numeric(nn$nn.dist[, 1]))
	}

	nn <- FNN::get.knnx(data = r, query = q, k = 2)
	as.numeric(nn$nn.dist[, 2])
}


#' Compute prediction task descriptors for sample and deployment locations
#'
#' Constructs task descriptors of the form \eqn{T = (x, d)}, where `d` is the
#' nearest-neighbor prediction distance and `x` denotes selected environmental
#' covariates.
#'
#' @param sample_dat Data frame containing sampled observations. Must include
#'   columns `id`, `x`, and `y`.
#' @param grid_dat Data frame containing deployment or prediction locations.
#'   Must include columns `id`, `x`, and `y`.
#' @param env_vars Optional character vector of environmental covariate names to
#'   include in the task descriptors. If `NULL`, all common non-coordinate,
#'   non-response variables are used.
#'
#' @return A named list with two data frames:
#'   \describe{
#'     \item{sample_tasks}{Task descriptors for sample locations.}
#'     \item{grid_tasks}{Task descriptors for deployment-grid locations.}
#'   }
compute_task_descriptors <- function(sample_dat, grid_dat, env_vars = NULL) {
	if (!all(c("x", "y") %in% names(sample_dat))) {
		stop("sample_dat must contain columns 'x' and 'y'.", call. = FALSE)
	}
	if (!all(c("x", "y") %in% names(grid_dat))) {
		stop("grid_dat must contain columns 'x' and 'y'.", call. = FALSE)
	}
	if (!("id" %in% names(sample_dat))) {
		stop("sample_dat must contain an 'id' column.", call. = FALSE)
	}
	if (!("id" %in% names(grid_dat))) {
		stop("grid_dat must contain an 'id' column.", call. = FALSE)
	}

	if (is.null(env_vars)) {
		env_vars <- intersect(
			setdiff(names(sample_dat), c("id", "x", "y", "z", "outcome")),
			setdiff(names(grid_dat), c("id", "x", "y", "z", "outcome"))
		)
	} else {
		env_vars <- intersect(env_vars, intersect(names(sample_dat), names(grid_dat)))
	}

	d_sample <- nearest_neighbor_distance(
		query_coords = sample_dat[, c("x", "y")],
		ref_coords = sample_dat[, c("x", "y")],
		exclude_self = TRUE
	)

	d_grid <- nearest_neighbor_distance(
		query_coords = grid_dat[, c("x", "y")],
		ref_coords = sample_dat[, c("x", "y")],
		exclude_self = FALSE
	)

	sample_tasks <- data.frame(
		id = sample_dat$id,
		d = d_sample,
		stringsAsFactors = FALSE
	)

	grid_tasks <- data.frame(
		id = grid_dat$id,
		d = d_grid,
		stringsAsFactors = FALSE
	)

	for (v in env_vars) {
		sample_tasks[[v]] <- sample_dat[[v]]
		grid_tasks[[v]] <- grid_dat[[v]]
	}

	list(
		sample_tasks = sample_tasks,
		grid_tasks = grid_tasks
	)
}


#' Compute pointwise prediction errors
#'
#' Computes residuals and derived error metrics for observed and predicted values.
#'
#' @param obs Numeric vector of observations.
#' @param pred Numeric vector of predictions.
#'
#' @return Data frame with columns:
#'   \describe{
#'     \item{obs}{Observed values}
#'     \item{pred}{Predicted values}
#'     \item{error}{Residuals (obs - pred)}
#'     \item{se}{Squared errors}
#'     \item{ae}{Absolute errors}
#'   }
compute_pointwise_errors <- function(obs, pred) {
	stopifnot(length(obs) == length(pred))

	err <- obs - pred

	data.frame(
		obs = obs,
		pred = pred,
		error = err,
		se = err^2,
		ae = abs(err)
	)
}


# ------------------------------------------------------------
# Weighted performance summaries
# ------------------------------------------------------------

#' Compute weighted mean with safeguards
#'
#' Computes a weighted mean while handling missing values and
#' degenerate weight vectors.
#'
#' @param x Numeric vector.
#' @param w Optional numeric weights.
#'
#' @return Numeric scalar.
weighted_mean_safe <- function(x, w = NULL) {
	if (is.null(w)) {
		return(mean(x, na.rm = TRUE))
	}

	ok <- is.finite(x) & is.finite(w)
	if (!any(ok)) {
		return(NA_real_)
	}

	x <- x[ok]
	w <- w[ok]

	w_sum <- sum(w)
	if (!is.finite(w_sum) || w_sum <= 0) {
		return(NA_real_)
	}

	sum(w * x) / w_sum
}


#' Summarize weighted validation losses
#'
#' Computes performance metrics from a loss data frame, optionally using weights.
#'
#' @param loss_df Data frame containing columns `error`, `se`, and `ae`.
#' @param weights Optional numeric vector of weights.
#'
#' @return Named numeric vector with elements:
#'   `bias`, `mse`, `rmse`, `mae`.
summarize_losses <- function(loss_df, weights = NULL) {
	if (is.null(weights)) {
		weights <- rep(1, nrow(loss_df))
	}

	stopifnot(length(weights) == nrow(loss_df))

	w_sum <- sum(weights)
	if (is.finite(w_sum) && w_sum > 0) {
		weights <- weights / w_sum
	}

	c(
		bias = weighted_mean_safe(loss_df$error, weights),
		mse = weighted_mean_safe(loss_df$se, weights),
		rmse = sqrt(weighted_mean_safe(loss_df$se, weights)),
		mae = weighted_mean_safe(loss_df$ae, weights)
	)
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

		if (is.numeric(x_ref) && length(unique(stats::na.omit(x_ref))) > 2) {
			probs <- seq(0, 1, by = by)
			qtiles <- stats::quantile(x_ref, probs = probs, na.rm = TRUE, names = FALSE)

			qtiles <- unique(qtiles)

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


#' Infer buffered loss engine from task list
#'
#' Chooses a buffered-loss computation strategy based on whether multiple tasks
#' share the same training configuration.
#'
#' @param task_list List of buffered CV tasks.
#'
#' @return Character scalar, either `"taskwise"` or `"reuse"`.
infer_buffered_loss_engine <- function(task_list) {
	if (length(task_list) <= 1) {
		return("taskwise")
	}

	keys <- vapply(
		task_list,
		function(task) {
			train_key <- paste(sort(as.integer(task$train_rows)), collapse = ",")
			excl_key <- if (!is.null(task$excluded_rows)) {
				paste(sort(as.integer(task$excluded_rows)), collapse = ",")
			} else {
				""
			}
			paste(train_key, excl_key, sep = "||")
		},
		character(1)
	)

	if (length(unique(keys)) < length(keys)) "reuse" else "taskwise"
}


#' Compute buffered losses taskwise
#'
#' Fits one model per buffered validation task and returns pointwise losses with
#' realized nearest-neighbor prediction distances.
#'
#' @param sample_dat Data frame with sampled observations.
#' @param task_list List of buffered CV tasks.
#' @param response Response column name.
#' @param predictors Predictor variable names.
#' @param fit_fun Model fitting function.
#' @param predict_fun Prediction function.
#' @param model Optional model label passed to `fit_fun`.
#' @param ... Additional arguments passed to `fit_fun`.
#'
#' @return Data frame of pointwise losses.
compute_buffered_task_losses_taskwise <- function(
	sample_dat,
	task_list,
	response,
	predictors,
	fit_fun,
	predict_fun,
	model = NULL,
	...
) {
	coords <- as.matrix(sample_dat[, c("x", "y")])
	out_list <- vector("list", length(task_list))

	for (i in seq_along(task_list)) {
		task <- task_list[[i]]

		if (is.null(task$test_rows) || is.null(task$train_rows)) {
			stop("Each task must contain 'test_rows' and 'train_rows'.", call. = FALSE)
		}

		test_idx <- as.integer(task$test_rows)
		train_idx <- as.integer(task$train_rows)

		train_dat <- sample_dat[train_idx, , drop = FALSE]
		test_dat <- sample_dat[test_idx, , drop = FALSE]

		fit <- fit_fun(
			train_dat = train_dat,
			model = model,
			response = response,
			predictors = predictors,
			...
		)

		pred <- predict_fun(fit, newdata = test_dat)

		pe <- compute_pointwise_errors(
			obs = test_dat[[response]],
			pred = pred
		)

		nn <- FNN::get.knnx(
			data = coords[train_idx, , drop = FALSE],
			query = coords[test_idx, , drop = FALSE],
			k = 1
		)
		d_realized <- nn$nn.dist[, 1]

		out_list[[i]] <- data.frame(
			id = sample_dat$id[test_idx],
			task_id = if (!is.null(task$task_id)) task$task_id else i,
			obs = pe$obs,
			pred = pe$pred,
			error = pe$error,
			se = pe$se,
			ae = pe$ae,
			d = d_realized
		)
	}

	loss_df <- do.call(rbind, out_list)
	rownames(loss_df) <- NULL
	loss_df
}


#' Compute buffered losses with training-set reuse
#'
#' Fits one model per unique buffered training configuration and reuses it
#' across all tasks sharing that configuration.
#'
#' @param sample_dat Data frame with sampled observations.
#' @param task_list List of buffered CV tasks.
#' @param response Response column name.
#' @param predictors Predictor variable names.
#' @param fit_fun Model fitting function.
#' @param predict_fun Prediction function.
#' @param model Optional model label passed to `fit_fun`.
#' @param ... Additional arguments passed to `fit_fun`.
#'
#' @return Data frame of pointwise losses.
compute_buffered_task_losses_reuse <- function(
	sample_dat,
	task_list,
	response,
	predictors,
	fit_fun,
	predict_fun,
	model = NULL,
	...
) {
	coords <- as.matrix(sample_dat[, c("x", "y")])

	keys <- vapply(
		task_list,
		function(task) {
			train_key <- paste(sort(as.integer(task$train_rows)), collapse = ",")
			excl_key <- if (!is.null(task$excluded_rows)) {
				paste(sort(as.integer(task$excluded_rows)), collapse = ",")
			} else {
				""
			}
			paste(train_key, excl_key, sep = "||")
		},
		character(1)
	)

	groups <- split(seq_along(task_list), keys)
	out_list <- vector("list", length(task_list))

	for (g in seq_along(groups)) {
		idx_group <- groups[[g]]
		first_task <- task_list[[idx_group[1]]]
		train_idx <- as.integer(first_task$train_rows)

		train_dat <- sample_dat[train_idx, , drop = FALSE]

		fit <- fit_fun(
			train_dat = train_dat,
			model = model,
			response = response,
			predictors = predictors,
			...
		)

		for (i in idx_group) {
			task <- task_list[[i]]
			test_idx <- as.integer(task$test_rows)
			test_dat <- sample_dat[test_idx, , drop = FALSE]

			pred <- predict_fun(fit, newdata = test_dat)

			pe <- compute_pointwise_errors(
				obs = test_dat[[response]],
				pred = pred
			)

			nn <- FNN::get.knnx(
				data = coords[train_idx, , drop = FALSE],
				query = coords[test_idx, , drop = FALSE],
				k = 1
			)
			d_realized <- nn$nn.dist[, 1]

			out_list[[i]] <- data.frame(
				id = sample_dat$id[test_idx],
				task_id = if (!is.null(task$task_id)) task$task_id else i,
				obs = pe$obs,
				pred = pe$pred,
				error = pe$error,
				se = pe$se,
				ae = pe$ae,
				d = d_realized
			)
		}
	}

	loss_df <- do.call(rbind, out_list)
	rownames(loss_df) <- NULL
	loss_df
}


#' Compute TWCV-weighted CV estimators from external folds
#'
#' Fits a model under externally supplied CV folds, computes pointwise losses,
#' augments them with task descriptors, and returns unweighted and TWCV-weighted
#' performance summaries.
#'
#' @param sample_dat Data frame with sampled observations.
#' @param grid_dat Data frame with deployment locations.
#' @param folds Integer vector of fold assignments for `sample_dat`.
#' @param response Response column name.
#' @param predictors Predictor variable names.
#' @param fit_fun Model fitting function.
#' @param predict_fun Prediction function.
#' @param model Optional model label passed to `fit_fun`.
#' @param twcv_specs Named list of TWCV specifications.
#' @param env_vars Optional task-descriptor variables. Defaults to `predictors`.
#' @param ... Additional arguments passed to `fit_fun`.
#'
#' @return List with losses, estimators, weights, and specs.
#' @export
compute_twcv_cv <- function(
	sample_dat,
	grid_dat,
	folds,
	response = "z",
	predictors,
	fit_fun,
	predict_fun,
	model = NULL,
	twcv_specs = NULL,
	env_vars = predictors,
	...
) {
	if (length(folds) != nrow(sample_dat)) {
		stop("'folds' must have length nrow(sample_dat).", call. = FALSE)
	}

	if (is.null(twcv_specs)) {
		twcv_specs <- make_default_twcv_spec_set(
			predictor_vars = predictors,
			env_vars = env_vars,
			include_distance = TRUE,
			balance_by = 0.2,
			shrink_lambda = 0.2
		)
	}

	pred <- rep(NA_real_, nrow(sample_dat))
	coords <- as.matrix(sample_dat[, c("x", "y")])
	d_cv <- rep(NA_real_, nrow(sample_dat))

	for (f in sort(unique(folds))) {
		test_idx <- which(folds == f)
		train_idx <- which(folds != f)

		fit <- fit_fun(
			train_dat = sample_dat[train_idx, , drop = FALSE],
			model = model,
			response = response,
			predictors = predictors,
			...
		)

		pred[test_idx] <- predict_fun(
			fit,
			newdata = sample_dat[test_idx, , drop = FALSE]
		)

		nn <- FNN::get.knnx(
			data = coords[train_idx, , drop = FALSE],
			query = coords[test_idx, , drop = FALSE],
			k = 1
		)
		d_cv[test_idx] <- nn$nn.dist[, 1]
	}

	pe <- compute_pointwise_errors(
		obs = sample_dat[[response]],
		pred = pred
	)

	loss_df <- data.frame(
		id = sample_dat$id,
		fold = folds,
		obs = pe$obs,
		pred = pe$pred,
		error = pe$error,
		se = pe$se,
		ae = pe$ae,
		d = d_cv
	)

	for (v in env_vars) {
		loss_df[[v]] <- sample_dat[[v]]
	}

	grid_tasks <- compute_task_descriptors(
		sample_dat = sample_dat,
		grid_dat = grid_dat,
		env_vars = env_vars
	)$grid_tasks

	twcv_out <- apply_twcv_estimators(
		loss_df = loss_df,
		grid_tasks = grid_tasks,
		twcv_specs = twcv_specs
	)

	list(
		losses = loss_df,
		estimators = twcv_out$estimators,
		weights = twcv_out$weights,
		twcv_specs = twcv_out$twcv_specs
	)
}


#' Compute TWCV-weighted estimators from buffered CV tasks
#'
#' Uses externally supplied buffered validation tasks, computes pointwise losses,
#' realized prediction distances, and returns unweighted and TWCV-weighted
#' performance summaries.
#'
#' Each task must provide at least `test_rows` and `train_rows`. Optionally it
#' may also contain `task_id`.
#'
#' @param sample_dat Data frame with sampled observations.
#' @param grid_dat Data frame with deployment locations.
#' @param task_list List of buffered CV tasks. Each element must contain
#'   `test_rows` and `train_rows`.
#' @param response Response column name.
#' @param predictors Predictor variable names used for model fitting.
#' @param fit_fun Model fitting function.
#' @param predict_fun Prediction function.
#' @param model Optional model label passed to `fit_fun`.
#' @param twcv_specs Named list of TWCV specifications. If `NULL`, a default
#'   set is constructed from `predictors`.
#' @param env_vars Optional task-descriptor variables. Defaults to `predictors`.
#' @param loss_engine Buffered loss engine. `"taskwise"` fits one model per
#'   task; `"reuse"` fits one model per unique training configuration; `"auto"`
#'   chooses between them. `"fold"` is accepted as an alias for `"reuse"`.
#' @param ... Additional arguments passed to `fit_fun`.
#'
#' @return List with losses, estimators, weights, specs, and chosen engine.
#' @export
compute_twcv_buffered <- function(
	sample_dat,
	grid_dat,
	task_list,
	response = "z",
	predictors,
	fit_fun,
	predict_fun,
	model = NULL,
	twcv_specs = NULL,
	env_vars = predictors,
	loss_engine = c("auto", "taskwise", "reuse", "fold"),
	...
) {
	loss_engine <- match.arg(loss_engine)

	if (is.null(sample_dat$id)) {
		stop("sample_dat must contain an 'id' column.", call. = FALSE)
	}
	if (!all(c("x", "y") %in% names(sample_dat))) {
		stop("sample_dat must contain columns 'x' and 'y'.", call. = FALSE)
	}
	if (length(task_list) == 0) {
		stop("'task_list' must not be empty.", call. = FALSE)
	}

	if (loss_engine == "fold") {
		loss_engine <- "reuse"
	}
	if (loss_engine == "auto") {
		loss_engine <- infer_buffered_loss_engine(task_list)
	}

	if (is.null(twcv_specs)) {
		twcv_specs <- make_default_twcv_spec_set(
			predictor_vars = predictors,
			env_vars = env_vars,
			include_distance = TRUE,
			balance_by = 0.2,
			shrink_lambda = 0.2
		)
	}

	loss_df <- switch(
		loss_engine,
		taskwise = compute_buffered_task_losses_taskwise(
			sample_dat = sample_dat,
			task_list = task_list,
			response = response,
			predictors = predictors,
			fit_fun = fit_fun,
			predict_fun = predict_fun,
			model = model,
			...
		),
		reuse = compute_buffered_task_losses_reuse(
			sample_dat = sample_dat,
			task_list = task_list,
			response = response,
			predictors = predictors,
			fit_fun = fit_fun,
			predict_fun = predict_fun,
			model = model,
			...
		)
	)

	idx <- match(loss_df$id, sample_dat$id)
	for (v in env_vars) {
		loss_df[[v]] <- sample_dat[[v]][idx]
	}

	grid_tasks <- compute_task_descriptors(
		sample_dat = sample_dat,
		grid_dat = grid_dat,
		env_vars = env_vars
	)$grid_tasks

	twcv_out <- apply_twcv_estimators(
		loss_df = loss_df,
		grid_tasks = grid_tasks,
		twcv_specs = twcv_specs
	)

	list(
		losses = loss_df,
		estimators = twcv_out$estimators,
		weights = twcv_out$weights,
		twcv_specs = twcv_out$twcv_specs,
		loss_engine = loss_engine
	)
}


#' Convert a buffered CAST::knndm object to TWCV task list
#'
#' Creates a task list suitable for `compute_twcv_buffered()` from a buffered
#' `CAST::knndm` object containing `indx_train`, `indx_test`, and optionally
#' `indx_exclude`.
#'
#' @param knn_obj A buffered `CAST::knndm` object.
#'
#' @return List of tasks with elements `task_id`, `test_rows`, `train_rows`,
#'   and `excluded_rows`.
#' @export
knndm_to_task <- function(knn_obj) {
	if (is.null(knn_obj$indx_train) || is.null(knn_obj$indx_test)) {
		stop(
			"knn_obj must contain 'indx_train' and 'indx_test'.",
			call. = FALSE
		)
	}

	n_tasks <- length(knn_obj$indx_test)

	if (length(knn_obj$indx_train) != n_tasks) {
		stop(
			"'indx_train' and 'indx_test' must have the same length.",
			call. = FALSE
		)
	}

	has_exclude <- !is.null(knn_obj$indx_exclude)
	if (has_exclude && length(knn_obj$indx_exclude) != n_tasks) {
		stop(
			"'indx_exclude' must have the same length as 'indx_test' if present.",
			call. = FALSE
		)
	}

	task_list <- vector("list", n_tasks)

	for (i in seq_len(n_tasks)) {
		task_list[[i]] <- list(
			task_id = i,
			test_rows = as.integer(knn_obj$indx_test[[i]]),
			train_rows = as.integer(knn_obj$indx_train[[i]]),
			excluded_rows = if (has_exclude) as.integer(knn_obj$indx_exclude[[i]]) else integer(0)
		)
	}

	task_list
}
