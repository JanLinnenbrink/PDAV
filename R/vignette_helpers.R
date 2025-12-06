#' @keywords internal
#' @noRd
generate_rast <- function() {
	rast_grid <- rast(xmin = 0, xmax = 200, ymin = 0, ymax = 200, ncols = 200, nrows = 200)

	grad_predictors <- sim_covariates(rast_grid, vgm = gstat::vgm(psill = 1, model = "Exp", range = 50), n = 8)
	landcover <- sim_covariates(
		rast_grid,
		vgm = gstat::vgm(psill = 5, model = "Exp", range = 100),
		n = 2,
		beta = 10,
		indicators = TRUE
	)

	predictors <- c(grad_predictors, landcover)
	names(predictors) <- c("temp", "moisture", "ph", "elev", "slope", "solar", "dist_road", "prod", "forest", "grass")

	outcome <- blend_rasters(
		predictors,
		~
			# species/habitat suitability score (unscaled)
			1.2 *
				exp(-((ph - 6.5)^2) / (2 * 0.6^2)) *
				(1 - exp(-moisture / 0.25)) + # pH optimum × moisture saturation
				0.5 * exp(-((temp - 15)^2) / (2 * 5^2)) * exp(-elev / 150) + # mild temp optimum, penalised by elevation
				0.3 * (forest / (1 + exp(-2 * (forest - 0.5)))) -
				0.2 * (grass / (1 + exp(-2 * (grass - 0.5)))) +
				0.05 * slope
		# forest boosts, grass reduces, small slope term
	)

	r <- c(predictors, outcome)
	terra::crs(r) <- "EPSG:3857"
	return(r)
}

#' @keywords internal
#' @noRd
generate_samples <- function(r, n_samples) {
	sample_random <- sam_field(
		x = r,
		size = n_samples,
		type = "random",
		na.rm = TRUE
	) |>
		terra::extract(x = r, bind = TRUE) |>
		st_as_sf() |>
		mutate(sampling = "random")

	sample_biased <- sam_field(
		x = r,
		size = n_samples,
		type = "clustered",
		type_opts = list(nclusters = 10, radius = 30),
		na.rm = TRUE
	) |>
		terra::extract(x = r, bind = TRUE) |>
		st_as_sf() |>
		mutate(sampling = "biased")

	sample_clustered <- sam_field(
		x = r,
		size = n_samples,
		type = "clustered",
		type_opts = list(nclusters = 5, radius = 15),
		na.rm = TRUE
	) |>
		terra::extract(x = r, bind = TRUE) |>
		st_as_sf() |>
		mutate(sampling = "clustered")

	samples <- rbind(sample_random, sample_biased, sample_clustered) |>
		mutate(sampling = factor(sampling, levels = c("random", "biased", "clustered")))

	return(samples)
}
