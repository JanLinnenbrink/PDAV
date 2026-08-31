# k-NNDM

## Introduction

k-fold nearest-neighbour distance matching (kNNDM) cross-validation was
developed by Linnenbrink et al. (2024). It is the k-fold version of NNDM
(Milà et al. (2022)), and as such much faster and better suited towards
medium to large sized datasets. As NNDM, kNNDM also matches the
distribution of nearest-neighbour distances (NND) between prediction
locations and training samples during CV, but opposed to NNDM uses a
clustering algorithm instead of excluding training samples.

The more detailed workflow is:

1.  Compute the distribution of NNDs between training samples, G_(j), as
    well as between prediction locations and training samples, G_(ij).

2.  Use the KS-test to evaluate whether G_(j) ≥ G_(ij), i.e. whether the
    samples are randomly distributed in the prediction area.

3.  Depending on the KS-test result:

    1.  If G_(j) ≥ G_(ij), return a random CV split and stop the
        algorithm here.
    2.  Otherwise, continue by clustering the training points into q ∈
        Q, where Q = \[k, …, N\].

4.  The q clusters are then merged along the first principal component
    of the training-point coordinates until k folds are reached.

5.  For each of the q fold configurations, the NND distribution between
    the unique folds is calculated, G_(j), and the distance between
    G_(j) and G_(ij) is calculated using the Wasserstein statistic, W.

6.  The fold configuration yielding the lowest W is returned.

## Setup

``` r

library(PDAV)
library(terra)
#> terra 1.9.46
library(sf)
#> Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.4.0; sf_use_s2() is TRUE
library(simsam)
library(ggplot2)
library(cowplot)
library(tidyterra)
#> 
#> Attaching package: 'tidyterra'
#> The following object is masked from 'package:stats':
#> 
#>     filter
library(ggnewscale)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:terra':
#> 
#>     intersect, union
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union

set.seed(100)

k <- 2
```

## Example data

### Simulate predictors and response

``` r

r <- PDAV:::generate_rast()
predictor_stack <- r[[setdiff(names(r), "outcome")]]
cate_rasters <- which(names(r) %in% c("forest", "grass"))

samples <- PDAV:::generate_samples(r, 100) |>
    filter(sampling == "biased")
samples$point_id <- 1:nrow(samples)
samples <- select(samples, point_id)
```

![Figure 1: The predictor
stack](k-NNDM_files/figure-html/unnamed-chunk-3-1.png)

Figure 1: The predictor stack

![Figure 2: The simulated outcome with training
locations.](k-NNDM_files/figure-html/unnamed-chunk-4-1.png)

Figure 2: The simulated outcome with training locations.

## kNNDM in geographical space

### 1. Compute G_(j) and G_(ij)

Firstly, the kNNDM algorithm calcuates the distribution of NND between
samples (G_(j)), as well as the distribution of NNDs between prediction
points and samples (G_(ij)). The calculation of NNDs is shown in the
following figure: the nearest neighbour distance between samples is
shown in the left panel, while the right panel shows nearest neighbour
distances between prediction points and samples.

![Figure 3: Nearest-Neighbour Distances (NNDs) between training points
(left) and between prediction and training points
(right).](k-NNDM_files/figure-html/unnamed-chunk-5-1.png)

Figure 3: Nearest-Neighbour Distances (NNDs) between training points
(left) and between prediction and training points (right).

To characterize the distribution of these NNDs, their empirical
cumulative density functions (ECDFs) are calculated:

``` r

Gj <- c(FNN::knn.dist(sample_coords, k = 1))
Gij <- c(FNN::knnx.dist(
    query = pred_coords,
    data = sample_coords,
    k = 1
))
```

![Figure 4: The Empirical Cumulative Density Functions (ECDFs) of the
Nearest-Neighbour Distances between prediction and training points
(orange), as well as between training points
(green).](k-NNDM_files/figure-html/unnamed-chunk-7-1.png)

Figure 4: The Empirical Cumulative Density Functions (ECDFs) of the
Nearest-Neighbour Distances between prediction and training points
(orange), as well as between training points (green).

### 2. Test wether G_(j) \<= G_(ij)

Next, we test if G_(j) \<= G_(ij), which would imply a random sampling
design. In that case, the algorithm would stop and return a random fold
assignment. However, as can bee seen below (and also from the figure
above), this is not the case:

``` r

testks <- suppressWarnings(stats::ks.test(Gj, Gij, alternative = "great"))
testks$p.value >= 0.05
#> [1] FALSE
```

### 3. Cluster samples into *q* clusters

Next, we cluster the training points into *q* clusters, where *q* ∈ *Q*
and *Q* = \[*k*, …, *N*\]. This means that *q* is an integer sequence
ranging from the number of folds, *k*, to the number of training points,
*N*. The length of this sequence is 100 by default, but for
visualization purposes it is reduced to 3 here.

``` r

nk_length <- 3

clustgrid <- data.frame(
    nk = as.integer(round(exp(seq(
        log(k),
        log(nrow(samples) - 2),
        length.out = nk_length
    ))))
)
clustgrid$W <- NA
clustgrid <- clustgrid[!duplicated(clustgrid$nk), ]
clustgroups <- list()

nk_df <- list()
tabclust <- list()
clust_nk <- list()
for (nk in clustgrid$nk) {
    # Create nk clusters
    clust_nk[[nk]] <- stats::kmeans(sample_coords, nk)$cluster
    # tabclust: table containing the unique clust_nks and their frequency
    tabclust[[nk]] <- as.data.frame(table(clust_nk[[nk]]))
    names(tabclust[[nk]]) <- c("clust_nk", "Freq")
    tabclust[[nk]]$clust_k <- NA
}

clustgrid
#>   nk  W
#> 1  2 NA
#> 2 14 NA
#> 3 98 NA
```

The sampled sequence of *q* clusters is in this case 2, 14, 98. The
following figure shows the training points coloured by their cluster for
each of the *q* clusters:

![Figure 5: The training points assigned to 2 (= \*k\*), 14 and 98 (=
number of training points) clusters. Point colours and filled areas show
the clusters.](k-NNDM_files/figure-html/unnamed-chunk-10-1.png)

Figure 5: The training points assigned to 2 (= *k*), 14 and 98 (= number
of training points) clusters. Point colours and filled areas show the
clusters.

### 4. Merge the *q* clusters along their first PC into *k* folds

The *q* clusters are then merged along the first principal component of
the samples coordinates if *q* \> *k*. This prevents contigous clusters
to end up in the same fold.

Therefore, we first calculate the PCA over the samples coordinates to
capture the axis of largest spatial variation (first PC, red arrow):

``` r

pcacoords <- stats::prcomp(sample_coords, center = TRUE, scale. = FALSE, rank = 1)
```

![Figure 6: The training points coloured by the first Principal
Component (PC 1) calculated from their coordinates, which is also shown
as a red arrow.](k-NNDM_files/figure-html/unnamed-chunk-12-1.png)

Figure 6: The training points coloured by the first Principal Component
(PC 1) calculated from their coordinates, which is also shown as a red
arrow.

Then, we order the *q* clusters according to the first principal
component:

``` r

for (nk in clustgrid$nk) {
    # compute cluster centroids and project that centroid to the first PC
    centr_tpoints <- sapply(tabclust[[nk]]$clust_nk, function(x) {
        centrpca <- matrix(apply(sample_coords[clust_nk[[nk]] %in% x, , drop = FALSE], 2, mean), nrow = 1)
        colnames(centrpca) <- colnames(sample_coords)
        return(predict(pcacoords, centrpca))
    })
    # Order the clusters along this first PC to ensure spatial separation when merging them later
    tabclust[[nk]]$centrpca <- centr_tpoints
    tabclust[[nk]] <- tabclust[[nk]][order(tabclust[[nk]]$centrpca), ]
}
```

![Figure 7: The \*q\* clusters coloured by the projection of their
centroid onto the first principal
component.](k-NNDM_files/figure-html/unnamed-chunk-14-1.png)

Figure 7: The *q* clusters coloured by the projection of their centroid
onto the first principal component.

And lastly, we merge the *q* clusters cyclically in the order of the
first principal component:

``` r

clust_k <- list()
tabclust2 <- list()
for (nk in clustgrid$nk) {
    # We don't merge big clusters
    clust_i <- 1
    for (i in 1:nrow(tabclust[[nk]])) {
        if (tabclust[[nk]]$Freq[i] >= nrow(samples) / k) {
            tabclust[[nk]]$clust_k[i] <- clust_i
            clust_i <- clust_i + 1
        }
    }
    rm("clust_i")

    ## And we merge the remaining nk clusters into k groups
    clust_i <- setdiff(1:k, unique(tabclust[[nk]]$clust_k))
    # tabclust is ordered along PC1 and thus repeating a sequence
    # from 1:k to assign folds results in spatially separated merging of clusters
    tabclust[[nk]]$clust_k[is.na(tabclust[[nk]]$clust_k)] <- rep(clust_i, ceiling(nk / length(clust_i)))[
        1:sum(is.na(tabclust[[nk]]$clust_k))
    ]

    # Then creating a new data frame containing the raw nk cluster assignment for each point
    tabclust2[[nk]] <- data.frame(ID = 1:length(clust_nk[[nk]]), clust_nk = clust_nk[[nk]])
    # And add the assignment to k folds to them (clust_k)
    tabclust2[[nk]] <- merge(tabclust2[[nk]], tabclust[[nk]], by = "clust_nk")
    tabclust2[[nk]] <- tabclust2[[nk]][order(tabclust2[[nk]]$ID), ]
    clust_k[[nk]] <- tabclust2[[nk]]$clust_k
    tabclust2[[nk]]$nk <- nk # For plotting only
}
```

    #> Warning: `scale_fill_colorblind()` was deprecated in ggthemes 5.2.0.
    #> This warning is displayed once per session.
    #> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    #> generated.

![Figure 8: Assignment of the clusters to \*k\* folds, for each
candidate number of clusters \*q\* (panels). Clusters larger than
\*n\*/\*k\* are not merged; the rest are assigned cyclically to the
remaining folds in PC1 order, so that spatially adjacent clusters land
in different folds. Colours show the fold id, outlines the cluster
convex hulls.](k-NNDM_files/figure-html/unnamed-chunk-16-1.png)

Figure 8: Assignment of the clusters to *k* folds, for each candidate
number of clusters *q* (panels). Clusters larger than *n*/*k* are not
merged; the rest are assigned cyclically to the remaining folds in PC1
order, so that spatially adjacent clusters land in different folds.
Colours show the fold id, outlines the cluster convex hulls.

### 5. Calculate G_(j)^(\*) for each *q* configuration and calculate *W*

Lastly, we calculate the ECDFs of NNDs between CV folds for each of the
three configurations. The absolute value difference between the ECDF of
NNDs between CV folds and the ECDF of NNDs between prediction and
training points is then quantified as the Wasserstein statistic (W). The
configuration yielding the lowest W is then chosen as the final
configuration (in this case configuration 2 with *nk* = 14, see table
below).

``` r

distclust_euclidean <- function(tr_coords, folds) {
    alldist <- rep(NA, length(folds))
    for (f in unique(folds)) {
        alldist[f == folds] <- c(FNN::knnx.dist(
            query = tr_coords[f == folds, , drop = FALSE],
            data = tr_coords[f != folds, , drop = FALSE],
            k = 1
        ))
    }
    alldist
}

Gjstar_i <- list()
for (nk in clustgrid$nk) {
    # Compute W statistic if not exceeding maxp
    if (!any(table(clust_k[[nk]]) / length(clust_k[[nk]]) > 0.8)) {
        Gjstar_i[[nk]] <- distclust_euclidean(sample_coords, clust_k[[nk]])
        clustgrid$W[clustgrid$nk == nk] <- twosamples::wass_stat(Gjstar_i[[nk]], Gij)
        clustgroups[[paste0("nk", nk)]] <- clust_k[[nk]]
    }
}
```

![Figure 9: The resulting ECDFs from the three different configurations
(panels).](k-NNDM_files/figure-html/unnamed-chunk-18-1.png)

Figure 9: The resulting ECDFs from the three different configurations
(panels).

    #>   nk         W
    #> 1  2 52.160269
    #> 2 14  9.287419
    #> 3 98 15.573773

### 6. Return the configuration yielding the lowest *W*

``` r

k_final <- clustgrid$nk[which.min(clustgrid$W)]
W_final <- min(clustgrid$W, na.rm = T)
clust <- clustgroups[[paste0("nk", k_final)]]
Gjstar <- distclust_euclidean(sample_coords, clust)

print(k_final)
#> [1] 14
```

The following plot shows the training samples with their final fold
assignment:

![Figure 10: The final fold assignment based on the configuration \*q\*
= 14.](k-NNDM_files/figure-html/unnamed-chunk-20-1.png)

Figure 10: The final fold assignment based on the configuration *q* =
14.

## Summary

![Figure 11: An overview plot over the steps involved in kNNDM: starting
from the initial \*q\* clusters (left panel), over projecting the first
PC of the training point coordinates to these clusters' centroids, to
merging them cyclically along the first PC into \*k\* folds and
calculating the ECDF of the resulting CV configurations (right
panel).](k-NNDM_files/figure-html/unnamed-chunk-21-1.png)

Figure 11: An overview plot over the steps involved in kNNDM: starting
from the initial *q* clusters (left panel), over projecting the first PC
of the training point coordinates to these clusters’ centroids, to
merging them cyclically along the first PC into *k* folds and
calculating the ECDF of the resulting CV configurations (right panel).

Linnenbrink, Jan, Carles Milà, Marvin Ludwig, and Hanna Meyer. 2024.
“kNNDM CV: *K* -Fold Nearest-Neighbour Distance Matching
Cross-Validation for Map Accuracy Estimation.” *Geoscientific Model
Development* 17 (15): 5897–912.
<https://doi.org/10.5194/gmd-17-5897-2024>.

Milà, Carles, Jorge Mateu, Edzer Pebesma, and Hanna Meyer. 2022.
“Nearest neighbour distance matching leave-one-out cross-validation for
map validation.” *Methods in Ecology and Evolution* 13 (6): 1304–16.
<https://doi.org/10.1111/2041-210X.13851>.
