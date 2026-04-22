#' Compute cross-validation losses
#'
#' Computes prediction errors for cross-validation predictions.
#'
#' Model fitting and prediction are supplied explicitly through `fit_fun` and
#' `predict_fun`, allowing study-specific model adapters to be passed into the
#' shared CV engine.
#'
#' @param sample_dat Data frame with response and an `id` column.
#' @param folds Integer vector of fold assignments.
#' @param model Model identifier passed to `fit_fun`.
#' @param response Optional response variable name. If `NULL`, the function
#'   tries `z` and then `outcome`.
#' @param fit_fun Model-fitting function passed to
#'   [compute_cv_predictions()].
#' @param predict_fun Prediction function passed to
#'   [compute_cv_predictions()].
#' @param verbose Verbosity level.
#' @param ... Additional arguments passed to `fit_fun`.
#'
#' @return Data frame with IDs, fold assignments, predictions, and pointwise
#'   error metrics.
#' @noRd
compute_cv_losses <- function(
	sample_dat,
	folds,
	model,
	response = NULL,
	fit_fun = fit_model,
	predict_fun = predict_model,
	verbose = 0,
	...
) {
	if (is.null(sample_dat$id)) {
		stop("sample_dat must contain an 'id' column.", call. = FALSE)
	}

	response <- response %||%
		if ("z" %in% names(sample_dat)) {
			"z"
		} else if ("outcome" %in% names(sample_dat)) {
			"outcome"
		} else {
			stop("Could not infer response column.", call. = FALSE)
		}

	pred <- compute_cv_predictions(
		sample_dat = sample_dat,
		folds = folds,
		model = model,
		response = response,
		fit_fun = fit_fun,
		predict_fun = predict_fun,
		verbose = verbose,
		...
	)

	pe <- compute_pointwise_errors(
		obs = sample_dat[[response]],
		pred = pred
	)

	data.frame(
		id = sample_dat$id,
		fold = folds,
		obs = pe$obs,
		pred = pe$pred,
		error = pe$error,
		se = pe$se,
		ae = pe$ae
	)
}


#' Compute buffered task losses
#'
#' Evaluates prediction errors for buffered validation tasks.
#'
#' Model fitting and prediction are supplied explicitly through `fit_fun` and
#' `predict_fun`, allowing the shared buffered-task engine to be reused across
#' different studies without redefining global adapter functions.
#'
#' @param sample_dat Data frame containing sampled observations.
#' @param task_obj Buffered task object.
#' @param model Model identifier passed to `fit_fun`.
#' @param response Optional response variable name. If `NULL`, the function tries `z` and then `outcome`.
#' @param fit_fun Model-fitting function. It must accept the training data as its first argument, and must also accept `model`, and `response`.
#' @param predict_fun Prediction function. It must accept a fitted model object and `newdata`.
#' @param verbose Verbosity level.
#' @param ... Additional arguments passed to `fit_fun`.
#'
#' @return Data frame of task-level prediction errors.
#' @noRd
compute_buffered_task_losses <- function(
	sample_dat,
	task_obj,
	model,
	response = NULL,
	fit_fun = fit_model,
	predict_fun = predict_model,
	verbose = 0,
	...
) {
	check_xy_columns(sample_dat)

	response <- response %||%
		if ("z" %in% names(sample_dat)) {
			"z"
		} else if ("outcome" %in% names(sample_dat)) {
			"outcome"
		} else {
			stop("Could not infer response column.", call. = FALSE)
		}

	ti <- task_obj$task_info
	tl <- task_obj$task_list

	if (is.null(ti) || nrow(ti) == 0 || !("feasible" %in% names(ti))) {
		warning("Buffered task object contains no feasibility information.")
		return(data.frame())
	}

	if (!is.logical(ti$feasible)) {
		ti$feasible <- as.logical(ti$feasible)
	}

	keep <- which(!is.na(ti$feasible) & ti$feasible)

	if (length(keep) == 0) {
		warning("No feasible buffered tasks available for loss computation.")
		return(data.frame())
	}

	log_message(verbose, 1, "Computing losses for buffered tasks...")
	log_message(verbose, 2, "  feasible tasks = ", length(keep))

	out_list <- vector("list", length(keep))

	for (jj in seq_along(keep)) {
		j <- keep[jj]
		task <- tl[[j]]

		if (verbose >= 2 && (jj %% 25 == 0 || jj == length(keep))) {
			log_message(verbose, 2, "  task ", jj, "/", length(keep))
		}

		train_dat <- sample_dat[task$train_rows, , drop = FALSE]
		test_dat <- sample_dat[task$test_row, , drop = FALSE]

		fit <- fit_fun(
			train_dat,
			model = model,
			response = response,
			...
		)
		pred <- predict_fun(fit, test_dat)

		pe <- compute_pointwise_errors(obs = test_dat[[response]], pred = pred)

		out_list[[jj]] <- data.frame(
			task_id = task$task_id,
			scheme = task$scheme,
			test_id = task$test_id,
			test_row = task$test_row,
			buffer_radius = task$buffer_radius,
			d_realized = task$d_realized,
			n_train = task$n_train,
			n_removed = task$n_removed,
			removed_frac = task$removed_frac,
			obs = pe$obs,
			pred = pe$pred,
			error = pe$error,
			se = pe$se,
			ae = pe$ae
		)
	}

	do.call(rbind, out_list)
}


#' Compute cross-validation predictions
#'
#' Fits models on training folds and predicts on held-out data.
#'
#' Model fitting and prediction are supplied explicitly through `fit_fun` and
#' `predict_fun`. This avoids reliance on globally defined adapter functions and
#' allows the same CV engine to be reused across different case studies.
#'
#' @param sample_dat Data frame containing the sampled observations.
#' @param folds Integer vector of fold assignments.
#' @param model Model identifier passed to `fit_fun`.
#' @param response Optional response variable name. If `NULL`, the function tries `z` and then `outcome`.
#' @param fit_fun Model-fitting function. It must accept the training data as its first argument, and must also accept `model`, and `response`.
#' @param predict_fun Prediction function. It must accept a fitted model object and `newdata`.
#' @param verbose Verbosity level.
#' @param ... Additional arguments passed to `fit_fun`.
#'
#' @return Numeric vector of cross-validation predictions.
#' @noRd
compute_cv_predictions <- function(
	sample_dat,
	folds,
	model,
	response = NULL,
	fit_fun = fit_model,
	predict_fun = predict_model,
	verbose = 0,
	...
) {
	response <- response %||%
		if ("z" %in% names(sample_dat)) {
			"z"
		} else if ("outcome" %in% names(sample_dat)) {
			"outcome"
		} else {
			stop("Could not infer response column.", call. = FALSE)
		}

	n <- nrow(sample_dat)
	pred <- rep(NA_real_, n)

	unique_folds <- sort(unique(folds))
	n_folds <- length(unique_folds)

	for (ii in seq_along(unique_folds)) {
		f <- unique_folds[ii]

		if (verbose >= 2) {
			if ((ii %% 25 == 0) || (ii == n_folds)) {
				log_message(verbose, 2, "  CV fold ", ii, "/", n_folds)
			}
		}

		test_idx <- which(folds == f)
		train_idx <- which(folds != f)

		train_dat <- sample_dat[train_idx, , drop = FALSE]
		test_dat <- sample_dat[test_idx, , drop = FALSE]

		fit <- fit_fun(
			train_dat,
			model = model,
			response = response,
			...
		)
		pred[test_idx] <- predict_fun(fit, test_dat)
	}

	pred
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
#' @noRd
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


#' Summarize weighted validation losses
#'
#' Computes performance metrics from a loss data frame, optionally using weights.
#'
#' @param loss_df Data frame containing columns `error`, `se`, and `ae`.
#' @param weights Optional numeric vector of weights.
#'
#' @return Named numeric vector with elements:
#'   `bias`, `mse`, `rmse`, `mae`.
#' @noRd
summarize_losses <- function(loss_df, weights = NULL) {
	if (is.null(weights)) {
		weights <- rep(1, nrow(loss_df))
	}

	stopifnot(length(weights) == nrow(loss_df))

	# normalize weights for numerical stability (optional but useful)
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


#' Compute weighted mean with safeguards
#'
#' Computes a weighted mean while handling missing values and
#' degenerate weight vectors.
#'
#' @param x Numeric vector.
#' @param w Optional numeric weights.
#'
#' @return Numeric scalar.
#' @noRd
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
