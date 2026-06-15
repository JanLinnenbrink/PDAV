#' Generate buffered leave-one-out tasks
#'
#' Constructs buffered LOO validation tasks by excluding observations within
#' a radius around each test point.
#'
#' @param sample_dat Data with coordinates and IDs.
#' @param buffer_radii Optional vector of radii.
#' @param target_d Optional target distance distribution.
#' @param n_candidates Number of candidate radii per point.
#' @param n_tasks Optional number of tasks to retain.
#' @param radius_correction Scaling factor for derived radii.
#' @param max_holdout_frac Maximum fraction of excluded data.
#' @param min_train_n Minimum training size.
#' @param max_dist Maximum allowed prediction distance.
#' @param include_zero Logical; include zero radius.
#' @param verbose Verbosity level.
#' @param seed Optional random seed.
#'
#' @return Object of class `"twcv_buffered_tasks"`.
generate_buffered_loo_tasks <- function(
	sample_dat,
	buffer_radii = NULL,
	target_d = NULL,
	n_candidates = 30,
	n_tasks = NULL,
	radius_correction = 0.6,
	max_holdout_frac = 0.2,
	min_train_n = floor(nrow(sample_dat) / 2),
	max_dist = NULL,
	include_zero = TRUE,
	verbose = 0,
	seed = NULL
) {
	if (!is.null(seed)) {
		set.seed(seed)
	}

	log_message(verbose, 1, "Generating buffered LOO tasks...")

	if (is.null(sample_dat$id)) {
		stop("sample_dat must contain an 'id' column.")
	}

	n <- nrow(sample_dat)
	coords <- as.matrix(sample_dat[, c("x", "y")])
	dmat <- pairwise_distance_matrix(coords)

	# ----------------------------------------------------------
	# Derive candidate radii
	# ----------------------------------------------------------
	if (is.null(buffer_radii)) {
		if (is.null(target_d)) {
			stop("Either 'buffer_radii' or 'target_d' must be provided.")
		}

		target_d <- as.numeric(target_d)
		target_d <- target_d[is.finite(target_d) & target_d > 0]

		if (length(target_d) == 0) {
			stop("'target_d' must contain at least one finite positive value.")
		}

		buffer_radii <- sort(unique(radius_correction * target_d))

		log_message(
			verbose,
			1,
			"Derived ",
			length(buffer_radii),
			" buffer radii from target distance distribution"
		)
	} else {
		buffer_radii <- sort(unique(as.numeric(buffer_radii)))
		buffer_radii <- buffer_radii[is.finite(buffer_radii) & buffer_radii >= 0]
	}

	if (length(buffer_radii) == 0) {
		stop("No valid buffer radii available.")
	}

	n_candidates <- min(n_candidates, length(buffer_radii))

	# ----------------------------------------------------------
	# Generate candidate tasks
	# ----------------------------------------------------------
	task_list <- vector("list", length = n * n_candidates)
	k <- 1L

	for (i in seq_len(n)) {
		if (verbose >= 2 && (i %% 25 == 0 || i == n)) {
			log_message(verbose, 2, "  focal point ", i, "/", n)
		}

		sel_buffer_radii <- if (include_zero) {
			n_nonzero <- min(n_candidates - 1L, sum(buffer_radii > 0))
			c(0, if (n_nonzero > 0) sample(buffer_radii[buffer_radii > 0], size = n_nonzero, replace = FALSE) else numeric(0))
		} else {
			sample(buffer_radii, size = n_candidates, replace = FALSE)
		}

		sel_buffer_radii <- unique(sel_buffer_radii)

		for (r in sel_buffer_radii) {
			excluded_idx <- which(dmat[i, ] <= r)
			excluded_idx <- sort(unique(c(i, excluded_idx)))

			train_idx <- setdiff(seq_len(n), excluded_idx)
			test_idx <- i

			n_removed <- length(excluded_idx)
			removed_frac <- n_removed / n
			n_train <- length(train_idx)

			feasible <- TRUE
			reason <- NA_character_

			if (removed_frac > max_holdout_frac) {
				feasible <- FALSE
				reason <- "max_holdout_frac"
			}

			if (n_train < min_train_n) {
				feasible <- FALSE
				reason <- if (is.na(reason)) "min_train_n" else paste(reason, "min_train_n", sep = ";")
			}

			d_realized <- if (n_train > 0) min(dmat[i, train_idx]) else NA_real_

			if (!is.null(max_dist) && is.finite(d_realized) && d_realized > max_dist) {
				feasible <- FALSE
				reason <- if (is.na(reason)) "max_dist" else paste(reason, "max_dist", sep = ";")
			}

			task_list[[k]] <- list(
				task_id = k,
				scheme = "buffered_loo",
				test_id = sample_dat$id[test_idx],
				test_row = test_idx,
				train_rows = train_idx,
				excluded_rows = excluded_idx,
				buffer_radius = r,
				d_realized = d_realized,
				n_train = n_train,
				n_removed = n_removed,
				removed_frac = removed_frac,
				feasible = feasible,
				infeasible_reason = reason
			)

			k <- k + 1L
		}
	}

	task_list <- task_list[seq_len(k - 1L)]

	log_message(verbose, 1, "Collecting task information")

	task_info <- do.call(
		rbind,
		lapply(task_list, function(x) {
			data.frame(
				task_id = x$task_id,
				scheme = x$scheme,
				test_id = x$test_id,
				test_row = x$test_row,
				buffer_radius = x$buffer_radius,
				d_realized = x$d_realized,
				n_train = x$n_train,
				n_removed = x$n_removed,
				removed_frac = x$removed_frac,
				feasible = x$feasible,
				infeasible_reason = if (is.na(x$infeasible_reason)) "" else x$infeasible_reason,
				stringsAsFactors = FALSE
			)
		})
	)

	rownames(task_info) <- NULL

	# ----------------------------------------------------------
	# Restrict to feasible candidates
	# ----------------------------------------------------------
	feasible_idx <- which(task_info$feasible & is.finite(task_info$d_realized))

	if (length(feasible_idx) == 0) {
		warning("No feasible buffered tasks were generated.")
		return(structure(
			list(
				candidate_task_info = task_info,
				candidate_task_list = task_list,
				task_info = NULL,
				task_list = NULL,
				buffer_radii = buffer_radii,
				constraints = list(
					max_holdout_frac = max_holdout_frac,
					min_train_n = min_train_n,
					max_dist = max_dist
				)
			),
			class = "twcv_buffered_tasks"
		))
	}

	log_message(
		verbose,
		1,
		"Generated ",
		nrow(task_info),
		" candidate tasks; ",
		length(feasible_idx),
		" feasible"
	)

	# ----------------------------------------------------------
	# Optional subsampling to match target distance distribution
	# ----------------------------------------------------------
	selected_idx <- feasible_idx

	if (!is.null(n_tasks)) {
		n_tasks <- min(as.integer(n_tasks), length(feasible_idx))

		if (is.null(target_d)) {
			selected_idx <- sample(feasible_idx, size = n_tasks)
			log_message(
				verbose,
				1,
				"Selected ",
				n_tasks,
				" feasible tasks by simple random sampling"
			)
		} else {
			target_d <- as.numeric(target_d)
			target_d <- target_d[is.finite(target_d) & target_d > 0]

			if (length(target_d) == 0) {
				stop("'target_d' must contain at least one finite positive value.")
			}

			target_draw <- if (length(target_d) >= n_tasks) {
				sample(target_d, size = n_tasks)
			} else {
				sample(target_d, size = n_tasks, replace = TRUE)
			}

			cand_d <- task_info$d_realized[feasible_idx]
			available <- seq_along(feasible_idx)
			chosen_local <- integer(0)

			for (td in target_draw) {
				if (length(available) == 0) {
					break
				}

				j <- available[which.min(abs(cand_d[available] - td))]
				chosen_local <- c(chosen_local, j)
				available <- setdiff(available, j)
			}

			selected_idx <- feasible_idx[chosen_local]

			log_message(
				verbose,
				1,
				"Selected ",
				length(selected_idx),
				" feasible tasks by approximate distance-distribution matching"
			)
		}
	}

	selected_task_info <- task_info[selected_idx, , drop = FALSE]
	selected_task_list <- task_list[selected_idx]
	rownames(selected_task_info) <- NULL

	structure(
		list(
			candidate_task_info = task_info,
			candidate_task_list = task_list,
			task_info = selected_task_info,
			task_list = selected_task_list,
			buffer_radii = buffer_radii,
			constraints = list(
				max_holdout_frac = max_holdout_frac,
				min_train_n = min_train_n,
				max_dist = max_dist
			)
		),
		class = "twcv_buffered_tasks"
	)
}
