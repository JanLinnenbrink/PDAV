#' Estimate deployment-oriented predictive performance from CV tasks
#'
#' Computes validation losses for a given cross-validation design, augments them
#' with realized task descriptors, and summarizes predictive performance under
#' several estimators targeting different task distributions. In addition to the
#' unweighted estimator, the function computes distance-weighted CV (DWCV),
#' target-weighted CV (TWCV), and optionally importance-weighted CV (IWCV).
#'
#' Model fitting and prediction are supplied explicitly through `fit_fun` and
#' `predict_fun`, allowing the shared estimator engine to be reused across case
#' studies without redefining global adapter functions.
#'
#' @param sample_dat Data frame of sampled observations used for validation.
#'   Must contain an ID column and the response variable, as well as "x" and "y" columns for the coordinates.
#' @param grid_dat Data frame representing deployment or prediction locations.
#' @param folds Integer vector of fold assignments for `sample_dat`.
#' @param model Character string identifying the prediction model.
#' @param response Optional response variable name. If `NULL`, the function
#'   tries `z` and then `outcome`.
#' @param fit_fun Model-fitting function passed to compute_cv_losses(). It
#'   must accept at least `train_dat`, `model`, and `response`.
#' @param predict_fun Prediction function passed to compute_cv_losses(). It
#'   must accept a fitted model object and `newdata`.
#' @param verbose Verbosity level.
#' @param twcv_specs Optional named list of TWCV specifications. Each
#'   specification must contain `balancing_vars`, `balance_by`, and
#'   `shrink_lambda`. If `NULL`, a default extended TWCV specification is used.
#' @param predictor_vars Optional character vector of predictor variables used
#'   by the predictive model. If `NULL`, common non-coordinate, non-response
#'   variables are inferred from `sample_dat` and `grid_dat`.
#' @param env_vars Optional character vector of environmental variables used as
#'   task descriptors. Defaults to `predictor_vars`.
#' @param use_dist Should geographical distances be accounted for during raking? Defaults to TRUE.
#' @param ... Additional arguments passed to `fit_fun` during validation.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{losses}{Validation-loss data frame augmented with task descriptors.}
#'     \item{estimators}{Named list of performance summaries for the unweighted TWCV.}
#'     \item{weights}{Named list of weight object for TWCV.}
#'     \item{twcv_specs}{The TWCV specification set actually used.}
#'   }
#'
#' @seealso twcv_buffered(), compute_twcv_weights()]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example data:
#' # sample_dat contains observed responses at sampled locations,
#' # grid_dat contains the deployment locations where predictions are intended.
#'
#' set.seed(1)
#'
#' n_sample <- 120
#' n_grid <- 500
#'
#' sample_dat <- data.frame(
#'   id = seq_len(n_sample),
#'   x = runif(n_sample),
#'   y = runif(n_sample),
#'   x1 = rnorm(n_sample),
#'   x2 = runif(n_sample),
#'   z = NA_real_
#' )
#'
#' sample_dat$z <- 1 + 2 * sample_dat$x1 - sample_dat$x2 + rnorm(n_sample, sd = 0.3)
#'
#' grid_dat <- data.frame(
#'   id = seq_len(n_grid),
#'   x = runif(n_grid),
#'   y = runif(n_grid),
#'   x1 = rnorm(n_grid),
#'   x2 = runif(n_grid)
#' )
#'
#' # Five-fold CV assignment
#' folds <- sample(rep(1:5, length.out = n_sample))
#'
#' # Simple model adapters used by twcv()
#' fit_lm <- function(train_dat, model, response, ...) {
#'   stats::lm(stats::as.formula(
#'     paste(response, "~", paste(c("x1", "x2"), collapse = " + "))
#'   ), data = train_dat)
#' }
#'
#' predict_lm <- function(object, newdata) {
#'   stats::predict(object, newdata = newdata)
#' }
#'
#' res <- twcv(
#'   sample_dat = sample_dat,
#'   grid_dat = grid_dat,
#'   folds = folds,
#'   model = "lm",
#'   response = "z",
#'   fit_fun = fit_lm,
#'   predict_fun = predict_lm,
#'   predictor_vars = c("x1", "x2"),
#'   env_vars = c("x1", "x2"),
#'   verbose = 1
#' )
#'
#' plot(res)
#'
#' # Augmented validation losses
#' head(res$losses)
#'
#' # Performance summaries
#' res$estimators
#'
#' # TWCV weight object(s)
#' names(res$weights)
#' }
twcv <- function(
	sample_dat,
	grid_dat,
	folds,
	model = c("rf", "lm", "ked_het_x1", "ked", "ked_het_pop"),
	response = NULL,
	fit_fun = fit_model,
	predict_fun = predict_model,
	verbose = 0,
	twcv_specs = NULL,
	predictor_vars = NULL,
	env_vars = NULL,
	use_dist = TRUE,
	...
) {
	model <- match.arg(model)

	response <- response %||%
		if ("z" %in% names(sample_dat)) {
			"z"
		} else if ("outcome" %in% names(sample_dat)) {
			"outcome"
		} else {
			stop("Could not infer response column.", call. = FALSE)
		}

	predictor_vars <- predictor_vars %||%
		intersect(
			setdiff(names(sample_dat), c("id", "set", "x", "y", "z", "outcome")),
			setdiff(names(grid_dat), c("id", "set", "x", "y", "z", "outcome"))
		)

	env_vars <- env_vars %||% predictor_vars

	if (is.null(twcv_specs)) {
		twcv_specs <- make_twcv_specs(
			predictor_vars = predictor_vars,
			include_distance = TRUE,
			balance_by = 0.2,
			shrink_lambda = 0.2,
			name = "twcv_extended"
		)
	}

	# Computes CV losses by calling compute_cv_predictions and compute_pointwise_errors
	log_message(verbose, 1, "Computing CV predictions...")
	cv_losses <- compute_cv_losses(
		sample_dat = sample_dat,
		folds = folds,
		model = model,
		response = response,
		fit_fun = fit_fun,
		predict_fun = predict_fun,
		predictors = predictor_vars,
		verbose = verbose,
		...
	)

	# Augments the CV validation task by calling compute_task_descriptors to prepare predictor values
	# and NNDs between samples (not returned), and between predpoints and samples
	# (compute_task_descriptors calls nearest_neighbor_distance which uses FNN).
	# Furthermore, it calls compute_cv_prediction_distance to calculate NNDs between folds (pairwise_distance_matrix).
	# Returns `losses` with sample task (NNDs between folds, CV error and predictor values) and
	# `grid_task` with grid task (NNDs between predpoints and samples, predictor values)
	log_message(verbose, 1, "Computing realized task descriptors...")
	aug <- augment_cv_task_descriptors(
		cv_losses = cv_losses,
		sample_dat = sample_dat,
		grid_dat = grid_dat,
		task_vars = env_vars,
		use_dist = use_dist
	)

	cv_losses <- aug$losses
	grid_tasks <- aug$grid_tasks

	est_list <- list(
		unweighted = summarize_losses(cv_losses)
	)

	weight_objects <- list()

	for (nm in names(twcv_specs)) {
		spec <- twcv_specs[[nm]]

		if (use_dist) {
			balancing_vars <- spec$balancing_vars
		} else {
			balancing_vars <- spec$balancing_vars[spec$balancing_vars != "d"]
		}

		if (
			is.null(balancing_vars) ||
				is.null(spec$balance_by) ||
				is.null(spec$shrink_lambda)
		) {
			stop("Malformed twcv spec for '", nm, "'.", call. = FALSE)
		}

		# Constructs balanced representations of sample and grid task descriptors for weighting.
		# Calls `prepare_for_balancing` which discretizes the variables in quantiles.
		# The returned data.frames have additional columns "x1_cat" etc. for every predictor that is used for weighting and
		# that specifies in which quantile the values of the predictor fall into (also for the NNDs, which is just another weighting variable here).
		bal <- prepare_balanced_tasks_cv(
			loss_df = cv_losses,
			grid_tasks = grid_tasks,
			balancing_vars = balancing_vars,
			by = spec$balance_by
		)

		# Calculates weights so that CV matches the margins of different task descriptors (predictors, NND) of the prediction task.
		# Calls compute_target_margins_generic which iterates over every discretized predictor ("marginal weighting"), calculates
		# the proportion of points falling into each quantile for the predpoints (tabulate_proportions).
		# Then calls `rake_weights` to apply proportional fitting ("raking").
		# Returns `weights_raw` (which are normalized by their mean) and `weights` which are normalized and shrunk towards 1.
		tw <- compute_twcv_weights(
			sample_tasks_bal = bal$sample_tasks_bal,
			grid_tasks_bal = bal$grid_tasks_bal,
			balancing_vars = balancing_vars,
			shrink_lambda = spec$shrink_lambda,
			verbose = max(0, verbose - 1)
		)

		# Check for any unsupported quintiles and issue a warning
		balance_df <- bal$sample_tasks_bal[grepl("_cat", names(bal$sample_tasks_bal))]
		names(balance_df) <- sub("_cat", "", names(balance_df))

		target_margins <- compute_target_margins_generic(
			grid_tasks_bal = bal$grid_tasks_bal,
			balancing_vars = balancing_vars
		)

		support_check <- check_balance_support(balance_df, target_margins)
		if (any(support_check$unsupported)) {
			unsupported_vars <- support_check[support_check$unsupported == TRUE, "var"]
			unsupported_flag <- 1
			warning(paste0(
				"The predictor(s) ",
				paste0(unsupported_vars, collapse = ","),
				" have quintiles that are not supported by the training data.
				Raking is likely to fail in this context, and limiting the prediction area to avoid extrapolation is recommended."
			))
		} else {
			unsupported_vars <- NA
			unsupported_flag <- 0
		}

		est_list[[nm]] <- summarize_losses(cv_losses, tw$weights)
		weight_objects[[nm]] <- tw
	}

	res <- list(
		losses = cv_losses,
		estimators = est_list,
		weights = weight_objects,
		twcv_specs = twcv_specs,
		unsupported_flag = unsupported_flag,
		unsupported_vars = unsupported_vars,
		use_dist = use_dist,
		balance_df = balance_df, # added for plotting purpose
		target_margins = target_margins # added for plotting purpose
	)
	class(res) <- "twcv"
	return(res)
}


#' Estimate deployment-oriented predictive performance from buffered tasks (Not used/exported at the momement)
#'
#' Analogous to [twcv()], but based on externally generated
#' buffered validation tasks, typically from buffered leave-one-out resampling.
#' Validation losses are computed for the selected buffered tasks, augmented
#' with realized task descriptors, and summarized using unweighted, DWCV, TWCV,
#' and optionally IWCV estimators.
#'
#' Model fitting and prediction are supplied explicitly through `fit_fun` and
#' `predict_fun`, allowing the shared buffered-task estimator engine to be
#' reused across case studies without redefining global adapter functions.
#'
#' @param sample_dat Data frame of sampled observations used for validation.
#' @param grid_dat Data frame representing deployment or prediction locations.
#' @param task_obj Buffered task object, typically created by
#'   generate_buffered_loo_tasks().
#' @param model Character string identifying the prediction model.
#' @param response Optional response variable name. If `NULL`, the function
#'   tries `z` and then `outcome`.
#' @param fit_fun Model-fitting function passed to
#'   compute_buffered_task_losses(). It must accept at least `train_dat`,
#'   `model`, and `response`.
#' @param predict_fun Prediction function passed to
#'   compute_buffered_task_losses(). It must accept a fitted model object and
#'   `newdata`.
#' @param verbose Verbosity level.
#' @param twcv_specs Optional named list of TWCV specifications.
#' @param predictor_vars Optional character vector of predictor variables used
#'   by the predictive model.
#' @param env_vars Optional character vector of environmental variables used as
#'   task descriptors. Defaults to `predictor_vars`.
#' @param iwcv_vars Optional character vector of variables for IWCV density-ratio
#'   estimation. Prediction distance `d` is appended internally when IWCV is
#'   enabled.
#' @param run_iwcv Logical; if `TRUE`, also compute IWCV.
#' @param iwcv_shrink_lambda Shrinkage parameter for IWCV weights.
#' @param ... Additional arguments passed to `fit_fun`.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{losses}{Buffered validation-loss data frame augmented with task
#'       descriptors.}
#'     \item{sample_tasks_bal}{Currently `NULL`; retained for backward
#'       compatibility.}
#'     \item{grid_tasks_bal}{Currently `NULL`; retained for backward
#'       compatibility.}
#'     \item{estimators}{Named list of performance summaries.}
#'     \item{weights}{Named list of weight objects for DWCV, TWCV, and
#'       optionally IWCV.}
#'     \item{twcv_specs}{The TWCV specification set actually used.}
#'   }
#'
#' @seealso [twcv()], [generate_buffered_loo_tasks()]
twcv_buffered <- function(
	sample_dat,
	grid_dat,
	task_obj,
	model = c("rf", "lm", "ked_het_x1", "ked", "ked_het_pop"),
	response = NULL,
	fit_fun = fit_model,
	predict_fun = predict_model,
	verbose = 0,
	twcv_specs = NULL,
	predictor_vars = NULL,
	env_vars = NULL,
	iwcv_vars = NULL,
	run_iwcv = FALSE,
	iwcv_shrink_lambda = 0,
	...
) {
	model <- match.arg(model)

	response <- response %||%
		if ("z" %in% names(sample_dat)) {
			"z"
		} else if ("outcome" %in% names(sample_dat)) {
			"outcome"
		} else {
			stop("Could not infer response column.", call. = FALSE)
		}

	predictor_vars <- predictor_vars %||%
		intersect(
			setdiff(names(sample_dat), c("id", "set", "x", "y", "z", "outcome")),
			setdiff(names(grid_dat), c("id", "set", "x", "y", "z", "outcome"))
		)

	env_vars <- env_vars %||% predictor_vars

	if (is.null(twcv_specs)) {
		twcv_specs <- make_twcv_specs(
			predictor_vars = predictor_vars,
			include_distance = TRUE,
			balance_by = 0.2,
			shrink_lambda = 0.2,
			name = "twcv_extended"
		)
	}

	# Computes CV losses by iterating over all training points, fitting a model, making predictions for the left-out data,
	# and then calling compute_pointwise_errors to calculate the error metrics
	raw_losses <- compute_buffered_task_losses(
		sample_dat = sample_dat,
		task_obj = task_obj,
		model = model,
		response = response,
		fit_fun = fit_fun,
		predict_fun = predict_fun,
		predictors = predictor_vars,
		verbose = verbose,
		...
	)

	if (nrow(raw_losses) == 0) {
		return(list(
			losses = data.frame(),
			sample_tasks_bal = NULL,
			grid_tasks_bal = NULL,
			estimators = NULL,
			weights = NULL
		))
	}

	# Augments the CV validation task by calling compute_task_descriptors to prepare predictor values
	# and NNDs between samples, and between predpoints and samples (compute_task_descriptors calls nearest_neighbor_distance which uses FNN).
	# It does not call compute_cv_prediction_distance to calculate NNDs between folds.
	# twcv_buffered rather needs those NNDs as input attached to task_obj, and calculated by generate_buffered_loo_tasks
	aug <- augment_buffered_task_descriptors(
		task_losses = raw_losses,
		sample_dat = sample_dat,
		grid_dat = grid_dat,
		task_vars = env_vars
	)

	loss_df <- aug$losses
	grid_tasks <- aug$grid_tasks

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
			shrink_lambda = spec$shrink_lambda,
			verbose = max(0, verbose - 1)
		)

		est_list[[nm]] <- summarize_losses(loss_df, tw$weights)
		weight_objects[[nm]] <- tw
	}

	list(
		losses = loss_df,
		sample_tasks_bal = NULL,
		grid_tasks_bal = NULL,
		estimators = est_list,
		weights = weight_objects,
		twcv_specs = twcv_specs
	)
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

#' Conditional logging with verbosity levels
#'
#' Prints a timestamped message when the requested verbosity threshold is met.
#'
#' @param verbose Integer verbosity level supplied by the caller.
#' @param level Minimum verbosity level required for printing.
#' @param ... Objects to be concatenated into the message text.
#'
#' @return Invisibly `NULL`.
log_message <- function(verbose = 0, level = 1, ...) {
	if (verbose >= level) {
		timestamp <- format(Sys.time(), "%H:%M:%S")
		message(sprintf("[%s] %s", timestamp, paste(..., collapse = "")))
	}
}
