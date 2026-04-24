#' Augment CV losses with task descriptors
#'
#' Adds task-descriptor variables and realized prediction distances to a
#' cross-validation loss table by matching sample locations to precomputed
#' task descriptors.
#'
#' @param cv_losses Data frame of CV losses containing at least \code{id} and
#'   \code{fold}.
#' @param sample_dat Sample data.
#' @param grid_dat Grid or population data.
#' @param task_vars Character vector of task-descriptor variables.
#'
#' @return List with components \code{losses} (augmented loss table) and
#'   \code{grid_tasks}.
#' @noRd
augment_cv_task_descriptors <- function(cv_losses, sample_dat, grid_dat, task_vars) {
	stopifnot(!is.null(task_vars), length(task_vars) > 0)
	assert_required_columns(cv_losses, c("id", "fold"))
	assert_required_columns(sample_dat, c("id"))

	task_vars <- intersect(task_vars, names(sample_dat))
	task_vars <- intersect(task_vars, names(grid_dat))

	if (length(task_vars) == 0) {
		stop("No valid task_vars found in sample_dat and grid_dat.", call. = FALSE)
	}

	# Calculates NNDs between predpoints and samples, and between samples using FNN, as well as the predictor values
	tasks <- compute_task_descriptors(
		sample_dat = sample_dat,
		grid_dat = grid_dat,
		env_vars = task_vars
	)

	sample_desc <- tasks$sample_tasks

	# enforce uniqueness of id
	if (anyDuplicated(sample_desc$id)) {
		stop("sample_tasks$id must be unique.", call. = FALSE)
	}

	idx <- match(cv_losses$id, sample_desc$id)

	if (anyNA(idx)) {
		stop("Some ids in cv_losses not found in sample_tasks.", call. = FALSE)
	}

	out <- cv_losses
	for (v in task_vars) {
		if (!(v %in% names(sample_desc))) {
			stop("Variable '", v, "' missing in sample_tasks.", call. = FALSE)
		}
		out[[v]] <- sample_desc[[v]][idx]
	}

	# Calculates the NND between folds based on distance matrix / matrices
	d_realized <- compute_cv_prediction_distances(
		sample_dat = sample_dat,
		folds = cv_losses$fold
	)

	idx_d <- match(out$id, sample_dat$id)
	out$d <- d_realized[idx_d]

	list(
		losses = out,
		grid_tasks = tasks$grid_tasks
	)
}


#' Augment buffered CV task losses with descriptors
#'
#' Adds task-descriptor variables to buffered cross-validation task losses,
#' using \code{test_id} to match tasks to sample locations. The balancing
#' variable \code{d} is set to the realized prediction distance under the
#' buffered training configuration.
#'
#' @param task_losses Data frame of buffered CV losses containing
#'   \code{test_id}.
#' @param sample_dat Sample data.
#' @param grid_dat Grid or population data.
#' @param task_vars Character vector of task-descriptor variables.
#'
#' @return List with components \code{losses} and \code{grid_tasks}.
#' @noRd
augment_buffered_task_descriptors <- function(task_losses, sample_dat, grid_dat, task_vars = NULL) {
	if (is.null(task_vars)) {
		stop("task_vars must be provided explicitly.", call. = FALSE)
	}

	task_vars <- intersect(task_vars, names(sample_dat))
	task_vars <- intersect(task_vars, names(grid_dat))

	if (length(task_vars) == 0) {
		stop("No valid task_vars found in sample_dat and grid_dat.", call. = FALSE)
	}

	# Calculates NNDs between predpoints and samples, and between samples using FNN, as well as the predictor values
	tasks <- compute_task_descriptors(
		sample_dat = sample_dat,
		grid_dat = grid_dat,
		env_vars = task_vars
	)

	task_vars <- intersect(task_vars, names(tasks$sample_tasks))
	if (length(task_vars) == 0) {
		stop("No valid task_vars found in sample_tasks.", call. = FALSE)
	}

	sample_desc <- tasks$sample_tasks[, c("id", task_vars), drop = FALSE]
	idx_desc <- match(task_losses$test_id, sample_desc$id)

	if (anyNA(idx_desc)) {
		stop("Some test_id values in task_losses were not found in sample_dat.", call. = FALSE)
	}

	out <- task_losses
	for (v in task_vars) {
		out[[v]] <- sample_desc[[v]][idx_desc]
	}

	if (!("d_realized" %in% names(out))) {
		stop("task_losses must contain 'd_realized'.", call. = FALSE)
	}

	# Canonical balancing variable: realized distance under buffering
	out$d <- out$d_realized

	# Drop redundant copy to avoid ambiguity downstream
	out$d_realized <- NULL

	list(
		losses = out,
		grid_tasks = tasks$grid_tasks
	)
}


#' Compute prediction task descriptors for sample and deployment locations
#'
#' Constructs task descriptors of the form \eqn{T = (x, d)}, where `d` is the
#' nearest-neighbor prediction distance and `x` denotes selected environmental
#' covariates. For sample locations, `d` is computed relative to the remaining
#' sample locations. For deployment-grid locations, `d` is computed relative to
#' the sample locations.
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
#'
#' @examples
#' \dontrun{
#' sample_dat <- data.frame(
#'   id = 1:3,
#'   x = c(0, 1, 2),
#'   y = c(0, 0, 0),
#'   z = c(10, 11, 12),
#'   elev = c(100, 110, 120)
#' )
#' grid_dat <- data.frame(
#'   id = 1:2,
#'   x = c(0.5, 1.5),
#'   y = c(0, 0),
#'   elev = c(105, 115)
#' )
#' compute_task_descriptors(sample_dat, grid_dat, env_vars = "elev")
#' }
#' @noRd
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


#' Compute CV prediction distances
#'
#' Computes nearest-neighbor distances between test points and their respective
#' training sets for each fold.
#'
#' @param sample_dat Data with coordinates.
#' @param folds Fold assignments.
#'
#' @return Numeric vector of distances.
#' @noRd
compute_cv_prediction_distances <- function(sample_dat, folds) {
	n <- nrow(sample_dat)
	d <- rep(NA_real_, n)

	check_xy_columns(sample_dat)
	coords <- as.matrix(sample_dat[, c("x", "y")])

	for (f in unique(folds)) {
		test_idx <- which(folds == f)
		train_idx <- which(folds != f)

		if (length(train_idx) == 0) {
			d[test_idx] <- NA_real_
			next
		}

		dmat <- pairwise_distance_matrix(
			coords1 = coords[test_idx, , drop = FALSE],
			coords2 = coords[train_idx, , drop = FALSE]
		)

		d[test_idx] <- apply(dmat, 1, min)
	}

	d
}

#' Check for numeric coordinate columns
#'
#' Verifies that a data object contains numeric columns named `x` and `y`.
#'
#' @param dat A data frame or similar object with named columns.
#'
#' @return Invisibly `TRUE` if the check succeeds.
#' @examples
#' \dontrun{
#' dat <- data.frame(x = 1:3, y = 4:6)
#' check_xy_columns(dat)
#' }
#' @noRd
check_xy_columns <- function(dat) {
	if (!all(c("x", "y") %in% names(dat))) {
		stop("Data must contain columns named 'x' and 'y'.", call. = FALSE)
	}
	if (!is.numeric(dat$x) || !is.numeric(dat$y)) {
		stop("Columns 'x' and 'y' must be numeric.", call. = FALSE)
	}
	invisible(TRUE)
}

#' Compute a Euclidean distance matrix
#'
#' Computes pairwise Euclidean distances either within one coordinate set or
#' between two coordinate sets.
#'
#' @param coords1 Matrix or data frame of coordinates with x/y columns or two
#'   columns.
#' @param coords2 Optional second matrix or data frame of coordinates. If
#'   omitted, distances within `coords1` are returned.
#'
#' @return A numeric matrix of Euclidean distances.
#' @examples
#' \dontrun{
#' xy <- data.frame(x = c(0, 1), y = c(0, 1))
#' pairwise_distance_matrix(xy)
#'
#' xy2 <- data.frame(x = c(1, 2), y = c(1, 2))
#' pairwise_distance_matrix(xy, xy2)
#' }
#' @noRd
pairwise_distance_matrix <- function(coords1, coords2 = NULL) {
	coords1 <- as.matrix(coords1)

	if (is.null(coords2)) {
		return(as.matrix(stats::dist(coords1)))
	}

	coords2 <- as.matrix(coords2)

	x1 <- coords1[, 1]
	y1 <- coords1[, 2]
	x2 <- coords2[, 1]
	y2 <- coords2[, 2]

	sqrt(outer(x1, x2, "-")^2 + outer(y1, y2, "-")^2)
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
#'
#' @examples
#' \dontrun{
#' pts <- data.frame(x = c(0, 1, 3), y = c(0, 0, 0))
#' nearest_neighbor_distance(pts, pts, exclude_self = TRUE)
#' }
#' @noRd
nearest_neighbor_distance <- function(query_coords, ref_coords, exclude_self = FALSE) {
	q <- as.matrix(query_coords[, c("x", "y")])
	r <- as.matrix(ref_coords[, c("x", "y")])

	if (!exclude_self) {
		nn <- FNN::get.knnx(data = r, query = q, k = 1)
		return(as.numeric(nn$nn.dist[, 1]))
	}

	# Self-exclusion only matters when query and reference coordinates coincide.
	same_object <- nrow(q) == nrow(r) &&
		isTRUE(all.equal(q, r, check.attributes = FALSE))

	if (!same_object) {
		nn <- FNN::get.knnx(data = r, query = q, k = 1)
		return(as.numeric(nn$nn.dist[, 1]))
	}

	# For identical point sets, use the second-nearest neighbor.
	nn <- FNN::get.knnx(data = r, query = q, k = 2)
	as.numeric(nn$nn.dist[, 2])
}


#' Assert presence of required columns
#'
#' Checks whether a data object contains all required columns and throws an
#' error otherwise.
#'
#' @param x Data frame or similar object with named columns.
#' @param required Character vector of required column names.
#'
#' @return Invisibly returns \code{x} if all required columns are present.
#' @noRd
assert_required_columns <- function(x, required) {
	miss <- setdiff(required, names(x))
	if (length(miss) > 0) {
		stop("Missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)
	}
	invisible(x)
}
