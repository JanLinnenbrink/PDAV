#' Random k-fold CV (RDM-CV)
#'
#' @param samples sf object with POINT geometries. Must have a column "ID".
#' @param folds_k integer, number of folds
#' @return Data frame with fold assignments
#' @export
RDM_CV <- function(samples, folds_k) {
	stopifnot(inherits(samples, "sf"))

	coords <- sf::st_coordinates(samples)
	ids <- if ("ID" %in% names(samples)) samples$ID else seq_len(nrow(samples))

	df <- data.frame(ID = ids, x = coords[, 1], y = coords[, 2])

	# Shuffle
	set.seed(123)
	df <- df[sample(nrow(df)), ]

	# Assign folds
	df$fold <- rep(0:(folds_k - 1), length.out = nrow(df))

	# Reorder back by ID for consistency
	df <- df[order(df$ID), ]

	return(df)
}
