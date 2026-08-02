#' Compute environmental similarity weighted CV metric from caret results
#'
#' Compute a deployment-oriented CV metric by reweighting caret out-of-fold
#' predictions to match the environmental distribution of a raster prediction grid.
#'
#' This is a simplified variant of the package's TWCV workflow that ignores
#' explicit spatial distance and uses only environmental predictor space.
#'
#' @param train_obj A `caret::train` object created with
#'   `trainControl(savePredictions = "all")`.
#' @param predictors A `terra::SpatRaster` containing environmental predictor
#'   layers for the prediction domain.
#' @param samples An `sf` object containing the training locations and the
#'   original response values used to fit `train_obj`.
#' @param response Optional response column name. If `NULL`, the response is
#'   inferred from `train_obj` or from `samples`.
#' @param env_vars Optional character vector of raster layer names to use for
#'   environmental similarity. If `NULL`, all layers in `predictors` are used.
#' @param metric_fun Function with signature `function(obs, pred, weights = NULL)`
#'   that returns a scalar performance metric. The default is RMSE.
#' @param balance_by Numeric quantile bin width used for discretizing environmental
#'   variables during weighting.
#' @param shrink_lambda Numeric shrinkage parameter for calibration weights.
#' @param verbose Integer verbosity level.
#' @param ... Additional arguments passed to `terra::extract()` when extracting
#'   raster values for the sample locations.
#'
#' @return A list with components:
#'   * `unweighted_metric`: the CV metric computed without weighting.
#'   * `weighted_metric`: the CV metric after environmental reweighting.
#'   * `weights`: the normalized calibration weights used for reweighting.
#'   * `losses`: data frame of out-of-fold observations, predictions, and metrics.
#'   * `grid_tasks`: environmental descriptors for the prediction grid.
#' @export
compute_cv_estimators_caret <- function(
	train_obj,
	predictors,
	samples,
	response = NULL,
	env_vars = NULL,
	metric_fun = rmse_metric,
	balance_by = 0.2,
	shrink_lambda = 0.2,
	verbose = 0,
	...
) {
	if (!inherits(train_obj, "train")) {
		stop("train_obj must be a caret::train object.", call. = FALSE)
	}
	if (!inherits(predictors, "SpatRaster")) {
		stop("predictors must be a terra::SpatRaster.", call. = FALSE)
	}
	if (!inherits(samples, "sf")) {
		stop("samples must be an sf object.", call. = FALSE)
	}

	if (!is.null(response)) {
		if (!response %in% names(samples)) {
			stop("response not found in samples.", call. = FALSE)
		}
	} else {
		response <- infer_caret_response(train_obj, samples)
	}

	if (!is.null(env_vars)) {
		bad <- setdiff(env_vars, names(predictors))
		if (length(bad) > 0) {
			stop("env_vars not found in predictors: ", paste(bad, collapse = ", "), call. = FALSE)
		}
	} else {
		env_vars <- names(predictors)
	}

	pred_df <- extract_caret_oof_predictions(train_obj)
	if (!"rowIndex" %in% names(pred_df)) {
		stop("caret predictions must include rowIndex. Use savePredictions = 'all'.", call. = FALSE)
	}
	if (max(pred_df$rowIndex, na.rm = TRUE) > nrow(samples)) {
		stop("rowIndex values exceed number of rows in samples.", call. = FALSE)
	}

	sample_env <- terra::extract(predictors[[env_vars]], terra::vect(samples), ID = FALSE, ...)
	if (nrow(sample_env) != nrow(samples)) {
		stop("Unable to extract environmental values for all sample locations.", call. = FALSE)
	}

	if (any(!complete.cases(sample_env))) {
		stop("Some sample locations contain NA predictor values. Remove or mask those locations.", call. = FALSE)
	}

	loss_df <- data.frame(
		rowIndex = pred_df$rowIndex,
		obs = as.numeric(pred_df$obs),
		pred = as.numeric(pred_df$pred),
		stringsAsFactors = FALSE
	)
	loss_df$sq_error <- (loss_df$obs - loss_df$pred)^2
	loss_df$ae <- abs(loss_df$obs - loss_df$pred)

	env_data <- sample_env[loss_df$rowIndex, , drop = FALSE]
	loss_df <- cbind(loss_df, env_data)

	grid_values <- terra::as.data.frame(predictors[[env_vars]], na.rm = FALSE)
	valid <- stats::complete.cases(grid_values)
	if (sum(valid) == 0) {
		stop("No valid prediction cells found in predictors.", call. = FALSE)
	}
	grid_tasks <- data.frame(grid_values[valid, , drop = FALSE], stringsAsFactors = FALSE)

	balance_vars <- env_vars
	bal <- prepare_balanced_tasks_cv(
		loss_df = loss_df,
		grid_tasks = grid_tasks,
		balancing_vars = balance_vars,
		by = balance_by
	)

	tw <- compute_twcv_weights(
		sample_tasks_bal = bal$sample_tasks_bal,
		grid_tasks_bal = bal$grid_tasks_bal,
		balancing_vars = balance_vars,
		shrink_lambda = shrink_lambda,
		verbose = max(0, verbose - 1)
	)

	unweighted_metric <- metric_fun(loss_df$obs, loss_df$pred)
	weighted_metric <- metric_fun(loss_df$obs, loss_df$pred, tw$weights)

	list(
		unweighted_metric = unweighted_metric,
		weighted_metric = weighted_metric,
		weights = tw$weights,
		losses = loss_df,
		grid_tasks = grid_tasks,
		twcv_spec = list(balancing_vars = balance_vars, balance_by = balance_by, shrink_lambda = shrink_lambda)
	)
}

infer_caret_response <- function(train_obj, samples) {
	if (!is.null(train_obj$terms)) {
		terms_vars <- all.vars(train_obj$terms)
		if (length(terms_vars) >= 1) {
			return(terms_vars[1])
		}
	}
	if (!is.null(train_obj$call$formula)) {
		formula_vars <- all.vars(as.formula(train_obj$call$formula))
		if (length(formula_vars) >= 1) {
			return(formula_vars[1])
		}
	}
	if ("z" %in% names(samples)) {
		return("z")
	}
	if ("outcome" %in% names(samples)) {
		return("outcome")
	}
	stop("Could not infer response column. Specify response explicitly.", call. = FALSE)
}

extract_caret_oof_predictions <- function(train_obj) {
	if (is.null(train_obj$pred) || nrow(train_obj$pred) == 0) {
		stop(
			"caret object does not contain saved predictions. Set savePredictions = 'all' in trainControl().",
			call. = FALSE
		)
	}

	pred_df <- train_obj$pred
	if (!is.null(train_obj$bestTune) && nrow(train_obj$bestTune) > 0) {
		tune_vars <- names(train_obj$bestTune)
		if (!all(tune_vars %in% names(pred_df))) {
			stop("caret object pred frame does not contain tuning columns.", call. = FALSE)
		}
		for (v in tune_vars) {
			pred_df <- pred_df[pred_df[[v]] == train_obj$bestTune[[v]], , drop = FALSE]
		}
	}

	if (!"rowIndex" %in% names(pred_df)) {
		stop("caret predictions do not include rowIndex. savePredictions must be set to 'all'.", call. = FALSE)
	}
	pred_df
}

#' Root mean squared error metric
#'
#' @param obs Numeric vector of observed values.
#' @param pred Numeric vector of predicted values.
#' @param weights Optional numeric weights.
#' @return Numeric RMSE.
#' @export
rmse_metric <- function(obs, pred, weights = NULL) {
	obs <- as.numeric(obs)
	pred <- as.numeric(pred)
	if (length(obs) != length(pred)) {
		stop("obs and pred must have equal length.", call. = FALSE)
	}
	if (is.null(weights)) {
		sqrt(mean((obs - pred)^2, na.rm = TRUE))
	} else {
		weights <- as.numeric(weights)
		if (length(weights) != length(obs)) {
			stop("weights must have the same length as obs and pred.", call. = FALSE)
		}
		sqrt(sum(weights * (obs - pred)^2, na.rm = TRUE) / sum(weights, na.rm = TRUE))
	}
}
