#' @keywords internal
#' @noRd
plot_task_density_bands <- function(
	sample_tasks,
	grid_tasks,
	var,
	cat_var = paste0(var, "_cat"),
	sample_title = "Training points",
	grid_title = "Prediction points",
	x_lab = var,
	fill_lab = paste0(var, " quintile"),
	palette = "Set2",
	panel_labels = c("A", "B"),
	return_components = FALSE
) {
	#-----------------------------
	# Checks
	#-----------------------------

	if (!(var %in% names(sample_tasks))) {
		stop("Variable '", var, "' not found in sample_tasks.", call. = FALSE)
	}

	if (!(var %in% names(grid_tasks))) {
		stop("Variable '", var, "' not found in grid_tasks.", call. = FALSE)
	}

	if (!(cat_var %in% names(sample_tasks))) {
		stop("Categorical variable '", cat_var, "' not found in sample_tasks.", call. = FALSE)
	}

	if (!(cat_var %in% names(grid_tasks))) {
		stop("Categorical variable '", cat_var, "' not found in grid_tasks.", call. = FALSE)
	}

	#-----------------------------
	# Helper: order Q-levels
	#-----------------------------

	get_q_levels <- function(x) {
		if (is.factor(x)) {
			levs <- levels(x)
		} else {
			levs <- unique(as.character(stats::na.omit(x)))
		}

		qnum <- suppressWarnings(as.integer(sub(".*_Q", "", levs)))

		if (all(!is.na(qnum))) {
			levs[order(qnum)]
		} else {
			sort(levs)
		}
	}

	#-----------------------------
	# Helper: make bands from existing bins
	#-----------------------------

	make_bands_from_existing_bins <- function(ref_dat, x_var, cat_var, xlim = NULL) {
		levs <- get_q_levels(ref_dat[[cat_var]])

		ranges <- ref_dat |>
			dplyr::transmute(
				x = .data[[x_var]],
				cat = factor(as.character(.data[[cat_var]]), levels = levs)
			) |>
			dplyr::filter(!is.na(x), !is.na(cat)) |>
			dplyr::group_by(cat) |>
			dplyr::summarise(
				lo = min(x),
				hi = max(x),
				.groups = "drop"
			) |>
			dplyr::arrange(cat)

		if (nrow(ranges) == 0) {
			stop("No non-missing values found for '", x_var, "' and '", cat_var, "'.", call. = FALSE)
		}

		if (is.null(xlim)) {
			xlim <- range(ref_dat[[x_var]], na.rm = TRUE)
		}

		if (nrow(ranges) == 1) {
			return(
				tibble::tibble(
					cat = ranges$cat,
					xmin = xlim[1],
					xmax = xlim[2]
				)
			)
		}

		# Boundaries between adjacent observed bins.
		# Exact quantile cut points are not stored, so these are reconstructed.
		cuts <- (head(ranges$hi, -1) + tail(ranges$lo, -1)) / 2

		tibble::tibble(
			cat = ranges$cat,
			xmin = c(xlim[1], cuts),
			xmax = c(cuts, xlim[2])
		)
	}

	#-----------------------------
	# Helper: individual plot
	#-----------------------------

	make_density_band_plot <- function(dat, bands, title, x_var, x_lab, fill_lab, palette) {
		p <- ggplot(dat, aes(x = .data[[x_var]])) +
			geom_rect(
				data = bands,
				aes(
					xmin = xmin,
					xmax = xmax,
					ymin = -Inf,
					ymax = Inf,
					fill = cat
				),
				inherit.aes = FALSE,
				alpha = 0.18
			) +
			geom_density(
				fill = "grey70",
				color = "grey20",
				alpha = 0.35,
				linewidth = 1,
				na.rm = TRUE
			) +
			coord_cartesian(
				xlim = range(c(bands$xmin, bands$xmax), na.rm = TRUE)
			) +
			scale_fill_brewer(
				palette = palette,
				drop = FALSE
			) +
			labs(
				title = title,
				x = x_lab,
				y = "Density",
				fill = fill_lab
			) +
			theme_minimal() +
			theme(
				legend.position = "bottom"
			)

		if (nrow(bands) > 1) {
			p <- p +
				geom_vline(
					data = bands[-nrow(bands), , drop = FALSE],
					aes(xintercept = xmax),
					inherit.aes = FALSE,
					linetype = "dashed",
					color = "grey30"
				)
		}

		p
	}

	#-----------------------------
	# Common x-range
	#-----------------------------

	x_rng <- c(
		sample_tasks[[var]],
		grid_tasks[[var]]
	) |>
		range(na.rm = TRUE)

	#-----------------------------
	# Build bands from grid categories
	#-----------------------------

	bands <- make_bands_from_existing_bins(
		ref_dat = grid_tasks,
		x_var = var,
		cat_var = cat_var,
		xlim = x_rng
	)

	#-----------------------------
	# Create plots
	#-----------------------------

	p_sample <- make_density_band_plot(
		dat = sample_tasks,
		bands = bands,
		title = sample_title,
		x_var = var,
		x_lab = x_lab,
		fill_lab = fill_lab,
		palette = palette
	)

	p_grid <- make_density_band_plot(
		dat = grid_tasks,
		bands = bands,
		title = grid_title,
		x_var = var,
		x_lab = x_lab,
		fill_lab = fill_lab,
		palette = palette
	)

	#-----------------------------
	# Shared legend
	#-----------------------------

	shared_legend <- cowplot::get_legend(
		p_sample +
			theme(
				legend.position = "bottom",
				legend.box = "horizontal"
			)
	)

	p_sample_noleg <- p_sample + theme(legend.position = "none")
	p_grid_noleg <- p_grid + theme(legend.position = "none")

	combined_plot <- cowplot::plot_grid(
		cowplot::plot_grid(
			p_sample_noleg,
			p_grid_noleg,
			nrow = 1,
			align = "hv",
			axis = "tblr",
			labels = panel_labels
		),
		shared_legend,
		ncol = 1,
		rel_heights = c(1, 0.12)
	)

	if (return_components) {
		return(
			list(
				plot = combined_plot,
				sample_plot = p_sample,
				grid_plot = p_grid,
				bands = bands,
				legend = shared_legend
			)
		)
	}

	combined_plot
}
