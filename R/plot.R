#' Plot PDAV classes
#' @description Generic plot function for prediction-domain adaptive validation results Classes
#'
#' @name plot
#' @param x An object of type \emph{nndm}.
#' @param type String, defaults to "strict" to show the original nearest neighbour distance definitions in the legend.
#' Alternatively, set to "simple" to have more intuitive labels.
#' @param ... other arguments.
#' @author Carles Milà
#'
#' @export
plot.nndm <- function(x, type = "strict", stat = "ecdf", ...) {
	# Prepare data for plotting: Gij function
	Gij_df <- data.frame(r = x$Gij[order(x$Gij)])
	Gij_df$val <- 1:nrow(Gij_df) / nrow(Gij_df)
	Gij_df <- Gij_df[Gij_df$r <= x$phi, ]
	Gij_df <- rbind(Gij_df, data.frame(r = 0, val = 0))
	Gij_df <- rbind(Gij_df, data.frame(r = x$phi, val = sum(x$Gij <= x$phi) / length(x$Gij)))
	Gij_df$Function <- "1_Gij(r)"

	# Prepare data for plotting: Gjstar function
	Gjstar_df <- data.frame(r = x$Gjstar[order(x$Gjstar)])
	Gjstar_df$val <- 1:nrow(Gjstar_df) / nrow(Gjstar_df)
	Gjstar_df <- Gjstar_df[Gjstar_df$r <= x$phi, ]
	Gjstar_df <- rbind(Gjstar_df, data.frame(r = 0, val = 0))
	Gjstar_df <- rbind(Gjstar_df, data.frame(r = x$phi, val = sum(x$Gjstar <= x$phi) / length(x$Gjstar)))
	Gjstar_df$Function <- "2_Gjstar(r)"

	# Prepare data for plotting: G function
	Gj_df <- data.frame(r = x$Gj[order(x$Gj)])
	Gj_df$val <- 1:nrow(Gj_df) / nrow(Gj_df)
	Gj_df <- Gj_df[Gj_df$r <= x$phi, ]
	Gj_df <- rbind(Gj_df, data.frame(r = 0, val = 0))
	Gj_df <- rbind(Gj_df, data.frame(r = x$phi, val = sum(x$Gj <= x$phi) / length(x$Gj)))
	Gj_df$Function <- "3_Gj(r)"

	# Merge data for plotting, get maxdist relevant for plotting
	if (any(Gj_df$val == 1) & any(Gjstar_df$val == 1) & any(Gij_df$val == 1)) {
		Gplot <- rbind(Gij_df, Gjstar_df, Gj_df)
		maxdist <- max(Gplot$r[Gplot$val != 1]) + 1e-9
		Gplot <- Gplot[Gplot$r <= maxdist, ]
		Gplot <- rbind(Gplot, data.frame(r = maxdist, val = 1, Function = c("1_Gij(r)", "2_Gjstar(r)", "3_Gj(r)")))
	} else {
		Gplot <- rbind(Gij_df, Gjstar_df, Gj_df)
	}

	# Define colours matching those of geodist
	myColors <- RColorBrewer::brewer.pal(3, "Dark2")

	# Plot
	if (stat == "ecdf") {
		p <- ggplot2::ggplot(data = Gplot, ggplot2::aes_string(x = "r", group = "Function", col = "Function")) +
			ggplot2::geom_vline(xintercept = 0, lwd = 0.1) +
			ggplot2::geom_hline(yintercept = 0, lwd = 0.1) +
			ggplot2::geom_hline(yintercept = 1, lwd = 0.1) +
			ggplot2::stat_ecdf(geom = "step", lwd = 0.8) +
			ggplot2::theme_bw() +
			ggplot2::ylab("ECDF") +
			ggplot2::labs(group = "Distance function", col = "Distance function") +
			ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 10))

		if (type == "strict") {
			p <- p +
				ggplot2::scale_colour_manual(
					values = c(myColors[2], myColors[3], myColors[1]),
					labels = c(expression(hat(G)[ij](r)), expression(hat(G)[j]^"*" * "(r,L)"), expression(hat(G)[j](r)))
				)
		} else if (type == "simple") {
			p <- p +
				ggplot2::scale_colour_manual(
					values = c(myColors[2], myColors[3], myColors[1]),
					labels = c("prediction-to-sample", "CV-distances", "sample-to-sample")
				)
		}
	} else if (stat == "density") {
		p <- ggplot2::ggplot(data = Gplot, ggplot2::aes_string(x = "r", group = "Function", fill = "Function")) +
			ggplot2::geom_density(adjust = 1.5, alpha = .5, stat = stat, lwd = 0.3) +
			ggplot2::theme_bw() +
			ggplot2::ylab("Density") +
			ggplot2::labs(group = "Distance function", col = "Distance function") +
			ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 10))

		if (type == "strict") {
			p <- p +
				ggplot2::scale_fill_manual(
					values = c(myColors[2], myColors[3], myColors[1]),
					labels = c(expression(hat(G)[ij](r)), expression(hat(G)[j]^"*" * "(r,L)"), expression(hat(G)[j](r)))
				)
		} else if (type == "simple") {
			p <- p +
				ggplot2::scale_fill_manual(
					values = c(myColors[2], myColors[3], myColors[1]),
					labels = c("prediction-to-sample", "CV-distances", "sample-to-sample")
				)
		}
	}

	p
}


#' @name plot
#' @param x An object of type \emph{knndm}.
#' @param type String, defaults to "strict" to show the original nearest neighbour distance definitions in the legend.
#' Alternatively, set to "simple" to have more intuitive labels.
#' @param stat String, defaults to "ecdf" but can be set to "density" to estimate density functions.
#' @param ... other arguments.
#' @author Carles Milà
#'
#' @export
plot.knndm <- function(x, type = "strict", stat = "ecdf", ...) {
	# Prepare data for plotting: Gij function
	Gij_df <- data.frame(r = x$Gij[order(x$Gij)])
	Gij_df$Function <- "1_Gij(r)"

	# Prepare data for plotting: Gjstar function
	Gjstar_df <- data.frame(r = x$Gjstar[order(x$Gjstar)])
	Gjstar_df$Function <- "2_Gjstar(r)"

	# Prepare data for plotting: G function
	Gj_df <- data.frame(r = x$Gj[order(x$Gj)])
	Gj_df$Function <- "3_Gj(r)"

	# Merge data for plotting
	Gplot <- rbind(Gij_df, Gjstar_df, Gj_df)

	# Define colours matching those of geodist
	myColors <- RColorBrewer::brewer.pal(3, "Dark2")

	# Plot
	if (stat == "ecdf") {
		p <- ggplot2::ggplot(data = Gplot, ggplot2::aes_string(x = "r", group = "Function", col = "Function")) +
			ggplot2::geom_vline(xintercept = 0, lwd = 0.1) +
			ggplot2::geom_hline(yintercept = 0, lwd = 0.1) +
			ggplot2::geom_hline(yintercept = 1, lwd = 0.1) +
			ggplot2::stat_ecdf(geom = "step", lwd = 0.8) +
			ggplot2::theme_bw() +
			ggplot2::ylab("ECDF") +
			ggplot2::labs(group = "Distance function", col = "Distance function") +
			ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 10))

		if (type == "strict") {
			p <- p +
				ggplot2::scale_colour_manual(
					values = c(myColors[2], myColors[3], myColors[1]),
					labels = c(expression(hat(G)[ij](r)), expression(hat(G)[j]^"*" * "(r,L)"), expression(hat(G)[j](r)))
				)
		} else if (type == "simple") {
			p <- p +
				ggplot2::scale_colour_manual(
					values = c(myColors[2], myColors[3], myColors[1]),
					labels = c("prediction-to-sample", "CV-distances", "sample-to-sample")
				)
		}
	} else if (stat == "density") {
		p <- ggplot2::ggplot(data = Gplot, ggplot2::aes_string(x = "r", group = "Function", fill = "Function")) +
			ggplot2::geom_density(adjust = 1.5, alpha = .5, stat = stat, lwd = 0.3) +
			ggplot2::theme_bw() +
			ggplot2::ylab("Density") +
			ggplot2::labs(group = "Distance function", col = "Distance function") +
			ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 10))

		if (type == "strict") {
			p <- p +
				ggplot2::scale_fill_manual(
					values = c(myColors[2], myColors[3], myColors[1]),
					labels = c(expression(hat(G)[ij](r)), expression(hat(G)[j]^"*" * "(r,L)"), expression(hat(G)[j](r)))
				)
		} else if (type == "simple") {
			p <- p +
				ggplot2::scale_fill_manual(
					values = c(myColors[2], myColors[3], myColors[1]),
					labels = c("prediction-to-sample", "CV-distances", "sample-to-sample")
				)
		}
	}

	p
}


#' @name plot
#' @param x An object of type \emph{da_cv}.
#' @param ... other arguments.
#' @author Jan Linnenbrink
#'
#' @export
plot.da_cv <- function(x, ...) {
	sim_levels <- terra::levels(x$category_raster)[[1]][, 2]

	if (length(sim_levels) == 2) {
		labs <- c("", "")
		vals <- c("grey", "transparent")
		legend_title <- "Dissimilar areas"
		cat_rast <- x$category_raster
		sim_rast <- x$similarity_raster
	} else {
		labs <- c("")
		vals <- NA
		legend_title <- "No dissimilar areas"
		cat_rast <- terra::as.factor(x$category_raster)
		sim_rast <- x$similarity_raster
	}

	ggplot2::ggplot() +
		# 1) Similarity raster first
		tidyterra::geom_spatraster(data = sim_rast) +
		ggplot2::scale_fill_viridis_c(
			"Similarity",
			option = "inferno",
			na.value = NA,
			guide = ggplot2::guide_colorbar(order = 1)
		) +

		# 2) New fill scale for categorical overlay
		ggnewscale::new_scale_fill() +
		tidyterra::geom_spatraster(data = cat_rast) +
		# Only show title, not factor levels
		ggplot2::scale_fill_manual(
			name = legend_title,
			values = vals,
			na.value = NA,
			labels = labs,
			na.translate = FALSE,
			guide = ggplot2::guide_legend(order = 2)
		) +

		ggplot2::theme_bw() +
		ggplot2::theme(
			legend.position = "bottom",
			legend.text = ggplot2::element_text(size = 10),
			panel.grid = ggplot2::element_blank()
		)
}
