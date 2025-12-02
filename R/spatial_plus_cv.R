#' Spatially-Aware Cross-Validation (SP-CV) Split
#'
#' Performs spatially-aware cross-validation (SP-CV) by splitting samples into folds.
#' Supports three methods:
#' \itemize{
#'   \item \strong{"SP1"}: Stage 1 spatial CV (hierarchical clustering, clusters assigned randomly to folds).
#'   \item \strong{"SP2"}: Each point treated as a cluster (finest granularity).
#'   \item \strong{"SP"}: General SP-CV with hierarchical clustering, k-means/k-modes,
#'         and ensemble majority voting over coordinates, environment, and target values.
#' }
#'
#' @param samples An \code{sf} object containing point geometries and a unique ID column.
#' @param response_name A data frame or matrix with target variable(s) (first column must match point IDs,
#'   last column treated as the target).
#' @param cate_col_start Integer, index of first categorical column in \code{env} (1-based).
#' @param cate_col_end Integer, index of last categorical column in \code{env} (1-based).
#' @param k Integer, number of folds/clusters.
#' @param sp_threshold Numeric, spatial distance threshold for hierarchical clustering (default 1).
#'   If set to 0, each point is its own cluster (SP2).
#' @param method Character, one of \code{"SP"}, \code{"SP1"}, or \code{"SP2"}.
#'
#' @return A data frame with columns:
#'   \itemize{
#'     \item \code{ID} - Original point ID
#'     \item \code{fold} - Assigned fold label
#'     \item \code{cluster} - Assigned cluster ID
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' library(sf)
#' pts <- st_as_sf(data.frame(ID = 1:10, x = runif(10), y = runif(10)),
#'                 coords = c("x", "y"), crs = 4326)
#' env <- data.frame(ID = 1:10, var1 = rnorm(10), var2 = runif(10), target = runif(10))
#' pts <- cbind(pts, env)
#' response_name <- "target"
#'
#' folds <- spatial_plus_cv(samples = pts, response_name = response_name,
#'                         cate_col_start = 0, cate_col_end = 0,
#'                         k = 3, sp_threshold = 1, method = "SP")
#' }
spatial_plus_cv <- function(
	samples,
	response_name,
	cate_col_start = 0,
	cate_col_end = 0,
	k = 5,
	sp_threshold = 1,
	method = c("SP", "SP1", "SP2")
) {
	method <- match.arg(method)

	# --- Helper: normalize continuous vars ---
	normalize_cols <- function(mat, start_cat, end_cat) {
		if (start_cat > 0) {
			mat[, seq_len(start_cat)] <- apply(mat[, seq_len(start_cat), drop = FALSE], 2, function(col) {
				(col - min(col)) / (max(col) - min(col))
			})
		}
		if (end_cat < ncol(mat)) {
			cols <- seq(end_cat + 1, ncol(mat))
			mat[, cols] <- apply(mat[, cols, drop = FALSE], 2, function(col) (col - min(col)) / (max(col) - min(col)))
		}
		mat
	}

	# --- Extract coordinates from sf ---
	coords <- as.data.frame(sf::st_coordinates(samples))
	coords$id <- 1:nrow(coords)
	pts_df <- sf::st_drop_geometry(samples)

	# --- Convert sf to input matrix for clustering ---
	if (method == "SP1") {
		sp_threshold <- sp_threshold
	}
	if (method == "SP2") {
		sp_threshold <- 0
	}

	# --- Spatial clustering ---
	if (sp_threshold == 0) {
		clusters <- lapply(seq_len(nrow(pts_df)), function(i) list(ids = i))
	} else {
		model <- cluster::agnes(coords, diss = FALSE, metric = "euclidean", method = "complete")
		cluster_labels <- stats::cutree(stats::as.hclust(model), h = sp_threshold)
		clusters <- lapply(unique(cluster_labels), function(cl) {
			list(ids = which(cluster_labels == cl))
		})
	}

	# --- Cluster averages ---
	norm_env <- normalize_cols(
		as.matrix(pts_df[, setdiff(names(pts_df), response_name), drop = FALSE]),
		cate_col_start,
		cate_col_end
	)

	coords_avg <- do.call(rbind, lapply(clusters, function(cl) colMeans(coords[cl$ids, , drop = FALSE])))
	env_avg <- do.call(rbind, lapply(clusters, function(cl) colMeans(norm_env[cl$ids, , drop = FALSE])))
	response_name_avg <- matrix(unlist(lapply(clusters, function(cl) mean(pts_df[cl$ids, response_name]))), ncol = 1)

	# --- K-means / K-modes on averages ---
	coords_k <- stats::kmeans(coords_avg, centers = k)$cluster
	env_k <- if (cate_col_start < cate_col_end) {
		klaR::kmodes(env_avg, modes = k)$cluster
	} else {
		stats::kmeans(env_avg, centers = k)$cluster
	}
	response_name_k <- stats::kmeans(response_name_avg, centers = k)$cluster

	# --- Ensemble majority vote ---
	ensemble_labels <- apply(cbind(coords_k, env_k, response_name_k), 1, function(x) {
		as.numeric(names(sort(table(x), decreasing = TRUE)[1]))
	})

	# --- Assign folds back to samples ---
	res <- do.call(
		rbind,
		lapply(seq_along(clusters), function(cid) {
			fold <- ensemble_labels[cid]
			cbind(ID = pts_df$ID[clusters[[cid]]$ids], fold = fold, cluster = cid)
		})
	)

	as.data.frame(res)
}
