#' Dissimilarity-Adaptive Cross-Validation (DA-CV)
#'
#' This function implements dissimilarity-adaptive cross-validation (DA-CV),
#' which combines random CV (RDM-CV) and spatial CV (SP-CV) based on
#' adversarial validation (AV). The method follows the framework described in
#' [Wang et al., 2025](https://doi.org/10.1016/j.ecoinf.2025.103287).
#'
#' @param samples An `sf` object with point samples. Must include the response
#'   variable and predictor variables as attributes.
#' @param predictors A `terra::rast` object containing predictor rasters
#'   aligned to the study area.
#' @param response Name of the response column in `samples`.
#' @param folds_k Integer. Number of folds for cross-validation.
#' @param autoc_threshold Numeric. Spatial autocorrelation threshold for spatial+ CV.
#' @param cate_num Integer. Number of spatial clusters (used in spatial+ CV).
#' @param seed Optional integer. Random seed for reproducibility (passed to \code{set.seed()}).
#' @return A list with components:
#' \describe{
#'   \item{dissimilarity}{Dissimilarity index D; ranging from 0 to 1}
#'   \item{threshold}{Threshold T(D) = D * 0.5}
#' 	 \item{similarity_raster}{Similarity raster}
#' 	 \item{category_raster}{Binarized similarity raster. 1 - Different, 2 - Similar}
#'   \item{weights}{Relative weights of "similar" and "different" areas}
#'   \item{folds_RDM}{folds from RDM-CV}
#'   \item{folds_SP}{folds from SP-CV}
#' }
#' @details
#' Steps:
#' \enumerate{
#'   \item Generate adversarial validation dataset by combining sample points
#'   with randomly sampled prediction locations from `predictors`.
#'   \item Train a classifier (random forest) to distinguish samples from
#'   prediction locations.
#'   \item Evaluate classifier via AUC → dissimilarity index D.
#'   \item Apply classifier to all raster cells → similarity map.
#'   \item Threshold map at T(D) = D * 0.5 to classify cells into
#'   \{similar, different\}.
#'   \item Calculate random and spatial+ cross-validation splits
#'   \item Combine RMSEs via weighted average by area proportion. (Needs to be done manually after running \code{DA_CV()})
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(terra)
#'
#' # Sample points
#' pts <- st_as_sf(data.frame(
#'   x = runif(50, 0, 100),
#'   y = runif(50, 0, 100),
#'   yvar = rnorm(50)
#' ), coords = c("x", "y"), crs = 32633)
#'
#' # Predictor rasters
#' r <- rast(ncols=100, nrows=100, xmin=0, xmax=100, ymin=0, ymax=100)
#' values(r) <- runif(ncell(r))
#'
#' # Run DA-CV (RDM_CV and spatial_plus_cv (must be available)
#' result <- DA_CV(
#'   samples = pts,
#'   predictors = r,
#'   response = "yvar",
#'   folds_k = 5,
#'   autoc_threshold = 0.2,
#'   cate_num = 5
#' )
#' print(result$weights)
#' }
#'
#' @export
DA_CV <- function(
	samples,
	predictors,
	response,
	folds_k = 5,
	autoc_threshold = 0.2,
	cate_num = 5,
	seed = NULL
) {
	# checks
	stopifnot(inherits(samples, "sf"))
	stopifnot(inherits(predictors, "SpatRaster"))
	if (!response %in% names(samples)) {
		stop("response column not found in samples")
	}
	if (!is.null(seed)) {
		set.seed(seed)
	}

	# number of samples
	n_samp <- nrow(samples)

	# 1) Build a data.frame of all raster cells with cell id + predictors
	all_cells_df <- terra::as.data.frame(predictors, cells = TRUE, na.rm = FALSE)
	# at this point 'all_cells_df' has a 'cell' column and one column per layer (layer names preserved)
	cell_col <- "cell"
	pred_cols <- setdiff(names(all_cells_df), cell_col)

	# 2) determine valid prediction cells: exclude any cell that has NA in any predictor
	na_row <- apply(all_cells_df[, pred_cols, drop = FALSE], 1, function(x) any(is.na(x)))
	valid_cells <- all_cells_df[[cell_col]][!na_row] # integer cell numbers that are valid
	if (length(valid_cells) < n_samp) {
		stop("Not enough valid prediction cells to sample same number as samples")
	}

	# 3) exclude cells that contain sample points
	samp_xy <- sf::st_coordinates(samples)
	samp_cells <- terra::cellFromXY(predictors, samp_xy) # may contain NA if some samples outside raster
	samp_cells <- unique(samp_cells[!is.na(samp_cells)])
	allowed_cells <- setdiff(valid_cells, samp_cells)
	if (length(allowed_cells) < n_samp) {
		stop("Not enough allowed cells after excluding sample cells")
	}

	# 4) sample exactly n_samp random prediction cells (no repeats)
	random_cells <- sample(allowed_cells, size = n_samp, replace = FALSE)

	# 5) extract predictor values for samples and for the selected random cells
	#    use terra::extract for samples (returns rows in order of samples)
	sample_vals <- terra::extract(predictors, samples, ID = FALSE)
	# for random cells, pull rows from all_cells_df
	rand_rows_idx <- match(random_cells, all_cells_df[[cell_col]])
	random_vals <- all_cells_df[rand_rows_idx, pred_cols, drop = FALSE]

	# 6) build AV dataset and labels (1 = sample, 0 = random pred)
	av_samples <- as.data.frame(sample_vals)
	av_samples$label <- factor(1, levels = c(0, 1))
	av_preds <- as.data.frame(random_vals)
	av_preds$label <- factor(0, levels = c(0, 1))
	av_df <- rbind(av_samples, av_preds)
	# ensure no rows with NA remain (shouldn't if we used valid cells); but be safe
	av_df <- av_df[stats::complete.cases(av_df), , drop = FALSE]
	if (nrow(av_df) < 2) {
		stop("AV dataset too small after removing NA rows")
	}

	# 7) split: randomly select exactly n_samp rows for training, remainder for testing
	if (nrow(av_df) < 2 * n_samp) {
		# if some rows removed due to NA, re-check sizes (we sampled allowed_cells so this rarely happens)
		warning("AV dataset has fewer rows than 2*n_samp; proceeding with available rows.")
	}
	train_idx <- sample(seq_len(nrow(av_df)), size = min(n_samp, nrow(av_df)), replace = FALSE)
	train_df <- av_df[train_idx, , drop = FALSE]
	test_df <- av_df[-train_idx, , drop = FALSE]

	# make sure label is a factor with levels "0","1"
	train_df$label <- factor(as.character(train_df$label), levels = c("0", "1"))
	test_df$label <- factor(as.character(test_df$label), levels = c("0", "1"))

	# 8) train classifier (ranger)
	rf <- ranger::ranger(
		formula = label ~ .,
		data = train_df,
		num.trees = 500,
		probability = TRUE
	)

	# 9) evaluate classifier: get probability of class "1" on test set
	pred_test <- stats::predict(rf, data = test_df)$predictions
	# predictions is a matrix with column names equal to factor levels ("0","1")
	if (!"1" %in% colnames(pred_test)) {
		# fallback: take second column if labels are ordered differently
		test_prob <- as.numeric(pred_test[, 2])
	} else {
		test_prob <- as.numeric(pred_test[, "1"])
	}
	# compute AUC
	auc_val <- as.numeric(pROC::auc(response = as.numeric(as.character(test_df$label)), predictor = test_prob))
	D_val <- if (auc_val <= 0.5) 0 else floor((auc_val - 0.5) * 2 * 100) / 100

	# 10) apply AV classifier to ALL valid prediction cells and build similarity map
	prob_raster <- terra::predict(
		predictors,
		rf,
		fun = function(model, data) {
			p <- stats::predict(model, data = data)$predictions
			if ("1" %in% colnames(p)) {
				return(p[, "1"])
			} else {
				return(p[, ncol(p)]) # fallback
			}
		},
		na.rm = TRUE,
		filename = "", # keep in memory; or set a filename for on-disk output
		cores = 1 # or >1 if parallel is safe in your setup
	)

	# 11) threshold and classify (T(D) = 0.5 * D)
	threshold <- D_val * 0.5
	category_raster <- prob_raster
	terra::values(category_raster) <- ifelse(
		terra::values(prob_raster) >= threshold,
		2L,
		1L
	)

	# count ratios directly on raster
	sim_count <- sum(terra::values(category_raster) == 2L, na.rm = TRUE)
	diff_count <- sum(terra::values(category_raster) == 1L, na.rm = TRUE)
	allpreds_count <- sim_count + diff_count

	if (allpreds_count == 0) {
		stop("No valid prediction cells after masking")
	}
	sim_ratio <- sim_count / allpreds_count
	diff_ratio <- diff_count / allpreds_count

	# Mask out sample cells explicitly (just like Python)
	# samp_cells <- terra::cellFromXY(predictors, sf::st_coordinates(samples))
	# category_raster[samp_cells] <- NA
	# prob_raster[samp_cells] <- NA

	# 12) Obtain random and spatial+ folds (to calculate weighted RMSE afterwards)
	folds_rdm <- RDM_CV(samples = samples, folds_k = folds_k)
	folds_sp <- spatial_plus_cv(
		samples = samples,
		response_name = response,
		k = folds_k,
		cate_col_start = 0,
		cate_col_end = cate_num,
		sp_threshold = autoc_threshold
	)

	# return
	return(list(
		dissimilarity = D_val,
		auc = auc_val,
		threshold = threshold,
		similarity_raster = prob_raster,
		category_raster = category_raster,
		weights = c(similar = sim_ratio, different = diff_ratio),
		folds_RDM = folds_rdm$fold,
		folds_SP = folds_sp$fold
	))
}
