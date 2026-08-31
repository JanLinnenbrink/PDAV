# Prediction-domain adaptive evaluation

## Introduction

This vignette provides an overview over the different methods of
prediction-domain adaptive evaluation. First, I will briefly summarize
the (dis-)advantages of the different methods (based on my subjective
impression), then I will briefly show how they work using a simulated
case study with three different sampling designs. For a more detailed
description of the methods, please refer to the dedicated vignettes.

## Overview of the currently developed methods:

[TABLE]

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
library(caret)
#> Loading required package: lattice

set.seed(100)
```

## Case Study

``` r

r <- PDAV:::generate_rast()
predictor_stack <- r[[setdiff(names(r), "outcome")]]
cate_rasters <- which(names(r) %in% c("forest", "grass"))

samples <- PDAV:::generate_samples(r, 100)
```

![Figure 1: The predictor
stack](Prediction-domain-adaptive-validation_files/figure-html/unnamed-chunk-3-1.png)

Figure 1: The predictor stack

![Figure 2: The simulated outcome with training
locations.](Prediction-domain-adaptive-validation_files/figure-html/unnamed-chunk-4-1.png)

Figure 2: The simulated outcome with training locations.

### NNDM CV

Nearest-Neighbour Distance Matching (NNDM, Milà et al. (2022)) tries to
match the prediction situation during cross-validation using the concept
of nearest-neighbour distances (NNDs). It formalizes the prediction
situation as the empirical cumulative density function (ECDF) of NNDs
between prediction and training locations. It then calculates the ECDF
of NNDs between training points, and tries to match this distribution to
the one found during prediction by consecutively excluding points the
training set in LOO-CV. If excluding a training point results in overly
difficult prediction situation during CV (this happens when excluding a
point results in a drop in the ECDF of NNDs during CV below the ECDF of
NNDs during prediction) the algorithm does not exclude this point.

``` r

results_nndm <- lapply(unique(samples$sampling), function(smpling) {
    samples_smpl <- samples |>
        filter(sampling == smpling) |>
        select(!sampling)

    nndm(
        tpoints = samples_smpl,
        modeldomain = predictor_stack,
        space = "geographical"
    )
})
#> Warning: [as.polygons] only the first layer is polygonized when 'dissolve=TRUE'
#> 1000 prediction points are sampled from the modeldomain
#> Warning: [as.polygons] only the first layer is polygonized when 'dissolve=TRUE'
#> 1000 prediction points are sampled from the modeldomain
#> Warning: [as.polygons] only the first layer is polygonized when 'dissolve=TRUE'
#> 1000 prediction points are sampled from the modeldomain

names(results_nndm) <- unique(samples$sampling)

nndm_plots <- lapply(names(results_nndm), function(smpl) {
    results_nndm_i <- results_nndm[[smpl]]
    plot(results_nndm_i) +
        ggtitle(smpl) +
        theme(
            legend.position = "bottom",
            plot.title = element_text(face = "bold", hjust = 0.5)
        )
})
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the PDAV package.
#>   Please report the issue to the authors.
#> This warning is displayed once per session.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.

nndm_legends <- get_legend(nndm_plots[[1]])[[1]][[1]]


nndm_plots <- plot_grid(
    plotlist = lapply(nndm_plots, function(x) x + theme(legend.position = "None")),
    nrow = 1
)
plot_grid(nndm_plots, nndm_legends, ncol = 1, rel_heights = c(1, 0.1))
```

![Figure 3: Results of
NNDM.](Prediction-domain-adaptive-validation_files/figure-html/unnamed-chunk-5-1.png)

Figure 3: Results of NNDM.

### kNNDM CV

k-fold Nearest-Neighbour Distance Matching (kNNDM, Linnenbrink et al.
(2024)) is the k-fold variant of NNDM. It is much faster, but also less
flexible in matching the prediction situation. kNNDM splits the training
points into a sequence of initial clusters, where the cluster size
ranges between the number of folds *k* and the number of training
points. Then, these initial clusters are merged until the number of
desired folds *k* is reached (not that for the first configuration no
merging is needed). The clusters are merged cyclically along the first
principal component of the training points coordinates, ensuring that no
neighbouring initial clusters are merged into the same folds (which
would result in contiguous blocks). Thus, kNNDM generates a sequence of
fold configurations that starts with a spatial block cross-validation
(when the number of initial clusters matches the number of folds) over
spatially structured but not completely blocked configurations towards
randomly assigned folds (when the number of initial clusters matches the
number of training points). From this sequence of fold configurations,
the one that yields the closest match to the prediction situation is
selected. This is measured as the absolute value integral (the
Wasserstein statistic) between the ECDF of NNDs between prediction and
training points (the prediction situation) and the ECDF of NNDs between
CV folds (the CV situation).

``` r

results_knndm <- lapply(unique(samples$sampling), function(smpling) {
    samples_smpl <- samples |>
        filter(sampling == smpling) |>
        select(!sampling)

    knndm(
        tpoints = samples_smpl,
        modeldomain = predictor_stack,
        space = "geographical"
    )
})
#> Warning: [as.polygons] only the first layer is polygonized when 'dissolve=TRUE'
#> 1000 prediction points are sampled from the modeldomain
#> Gij <= Gj; a random CV assignment is returned
#> Registered S3 methods overwritten by 'CAST':
#>   method     from
#>   plot.knndm PDAV
#>   plot.nndm  PDAV
#> Warning: [as.polygons] only the first layer is polygonized when 'dissolve=TRUE'
#> 1000 prediction points are sampled from the modeldomain
#> Warning: [as.polygons] only the first layer is polygonized when 'dissolve=TRUE'
#> 1000 prediction points are sampled from the modeldomain

names(results_knndm) <- unique(samples$sampling)

knndm_plots <- lapply(names(results_knndm), function(smpl) {
    results_knndm_i <- results_knndm[[smpl]]
    plot(results_knndm_i) +
        ggtitle(smpl) +
        theme(
            legend.position = "bottom",
            plot.title = element_text(face = "bold", hjust = 0.5)
        )
})

knndm_legends <- get_legend(knndm_plots[[1]])[[1]][[1]]

knndm_plots <- plot_grid(plotlist = lapply(knndm_plots, function(x) x + theme(legend.position = "None")), nrow = 1)
plot_grid(knndm_plots, knndm_legends, ncol = 1, rel_heights = c(1, 0.1))
```

![Figure 4: Results of
kNNDM.](Prediction-domain-adaptive-validation_files/figure-html/unnamed-chunk-6-1.png)

Figure 4: Results of kNNDM.

### DA-CV

Dissimilarity-Adaptive Cross-Validation (DA-CV, Wang et al. (2025a))
uses adversial validation (AV) to predict the probability that a
prediction location is similar to the training samples. Then it
calculates a RMSE based on random CV, as well as spatial+ CV (Wang et
al. (2023)) and weights both of them according to the relative area of
similar cells (random CV) and dissimilar cells (spatial+ CV).

``` r

results <- lapply(unique(samples$sampling), function(smpling) {
    samples_smpl <- samples |>
        filter(sampling == smpling) |>
        select(!sampling)

    da_cv(
        samples = samples_smpl,
        predictors = r,
        response = "outcome",
        folds_k = 5,
        autoc_threshold = 1
        #cate_col_start = min(cate_rasters),
        #cate_col_end = max(cate_rasters)
    )
})
#> Setting levels: control = 0, case = 1
#> Setting direction: controls > cases
#> Setting levels: control = 0, case = 1
#> Setting direction: controls < cases
#> Setting levels: control = 0, case = 1
#> Setting direction: controls < cases

names(results) <- unique(samples$sampling)
```

For randomly distributed samples, the AV classifier has a performance of
AUC = 0.5. This is then normalized to D = 0, because the classifier is
not better than random guessing at distinguishing whether a prediction
location is similar or dissimilar to the training samples; see Wang et
al. (2025b), section 2.1, for more information. The threshold is then
calculated as T(D) = 0.5 × 0, which is 0. This means that all prediction
locations are considered similar to the training samples, and no
extrapolation correction is required, as shown in the following map.
Consequently, all weight is assigned to random CV.

For clustered samples, the AV classifier achieves a performance of AUC =
0.77. This leads to D = (0.77 - 0.5) / (1 - 0.5) = 0.54. The threshold
is then T(D) = 0.5 × 0.54 = 0.27. Hence, all prediction cells with a
similarity score lower than 0.27 are classified as dissimilar. The
relative fraction of prediction locations classified as similar to the
sampling locations is 0.79, while the fraction classified as dissimilar
is 0.21.

``` r

lapply(unique(samples$sampling), function(smpling) {
    plot(results[[smpling]]) +
        new_scale_fill() +
        geom_sf(data = samples[samples$sampling == smpling, ], shape = 21, fill = "white") +
        coord_sf(expand = FALSE, datum = st_crs(samples[samples$sampling == smpling, ], )) +
        ggtitle(smpling) +
        theme(
            legend.position = "right",
            plot.title = element_text(face = "bold", hjust = 0.5)
        )
}) |>
    plot_grid(plotlist = _, nrow = 1, align = "vh")
```

![Figure 5: Similarity maps resulting from the AV
classifier.](Prediction-domain-adaptive-validation_files/figure-html/unnamed-chunk-8-1.png)

Figure 5: Similarity maps resulting from the AV classifier.

Lastly, the resulting cross-validation fold assignments can be used to
calculate the weighted RMSE of a spatial predictive model.

``` r

form <- as.formula(paste0("outcome~", paste0(names(predictor_stack), collapse = "+")))
pgrid <- data.frame(mtry = 6, splitrule = "variance", min.node.size = 5)

weighted_RMSE <- lapply(unique(samples$sampling), function(smpling) {
    samples_smpl <- samples |>
        filter(sampling == smpling) |>
        mutate("randomCV" = results[[smpling]]$folds_RDM, "spatialCV" = results[[smpling]]$folds_SP) |>
        st_drop_geometry()

    folds_random <- CAST::CreateSpacetimeFolds(samples_smpl, spacevar = "randomCV", k = 5)
    folds_spatial <- CAST::CreateSpacetimeFolds(samples_smpl, spacevar = "spatialCV", k = 5)

    train_cntrl_random <- trainControl(
        method = "CV",
        index = folds_random$index,
        indexOut = folds_random$indexOut,
        savePredictions = TRUE
    )

    train_cntrl_SP <- trainControl(
        method = "CV",
        index = folds_spatial$index,
        indexOut = folds_spatial$indexOut,
        savePredictions = TRUE
    )

    rand_mod <- train(
        form,
        data = samples_smpl,
        method = "ranger",
        trControl = train_cntrl_random,
        tuneGrid = pgrid
    )

    spat_mod <- train(
        form,
        data = samples_smpl,
        method = "ranger",
        trControl = train_cntrl_SP,
        tuneGrid = pgrid
    )

    err_stats_rand <- CAST::global_validation(rand_mod)
    err_stats_SP <- CAST::global_validation(spat_mod)

    err_stats_weighted <- sqrt(
        results[[smpling]]$weights[["similar"]] *
            (err_stats_rand[["RMSE"]]^2) +
            results[[smpling]]$weights[["different"]] * (err_stats_SP[["RMSE"]]^2)
    )

    prediction <- predict(r, rand_mod)
    list(err_stats_weighted, prediction)
})
```

The RMSE values obtained by DA-CV are 6.014 for the random sampling
design, 5.041 for the biased sampling and 5.508 for the clustered
design.

### TWCV

Target-Weighted Cross Validation (TWCV, Brenning and Suesse (2026))
takes a different viewpoint than the cross-validation approaches above.
Instead of adjusting the resampling design (or weighting two contrasting
resampling designs), it assigns a weight to each training points based
on the distribution of predictor values (and spatial nearest-neighbour
distances) in the training and in the prediction data. It uses raking to
match the distribution of predictors observed at the training points to
those encountered in the prediction area by up-weighting training points
that show predictor values underrepresented in the training data (as
compared to the prediction area) and down-weighting those points that
show overrepresented conditions. These weights are then normalized and
multiplied with the (cross-validation) error calculated for each
training point. Therefore, it relies on some sort of resampling
beforehand to calculate the final resulting error, and is highly
sensitive to the choice of resampling (the weights stay the same
regardless of the resampling method , but it makes a huge difference if
one uses random CV or spatial CV, since the error that is then weighted
can vary widely based on the chosen resampling method).

We provide an experimental implementation of this algorithm, alongside
dedicated plotting functions, in the package
[PredictionMatching](https://github.com/JanLinnenbrink/PredictionMatching).
In our example, we apply TWCV to the results of kNNDM:

    #> Registered S3 method overwritten by 'PredictionMatching':
    #>   method    from
    #>   plot.twcv PDAV
    #> Warning: [as.polygons] only the first layer is polygonized when 'dissolve=TRUE'
    #> 1000 prediction points are sampled from the modeldomain
    #> predictor values are extracted for prediction points
    #> Warning: [as.polygons] only the first layer is polygonized when 'dissolve=TRUE'
    #> 1000 prediction points are sampled from the modeldomain
    #> predictor values are extracted for prediction points
    #> Warning: [as.polygons] only the first layer is polygonized when 'dissolve=TRUE'
    #> 1000 prediction points are sampled from the modeldomain
    #> predictor values are extracted for prediction points

Now that we obtained the weights for each training point, we can
calculate the cross-validation error using kNNDM and then multiply the
error with these weights:

``` r


pe_samples <- lapply(unique(samples$sampling), function(smpling) {
    samples_smpl <- samples |>
        filter(sampling == smpling) |>
        select(!sampling)

    knndm_folds <- knndm(
        tpoints = samples_smpl,
        modeldomain = predictor_stack,
        space = "geographical"
    )

    ctrl_knndm <- trainControl(
        method = "cv",
        index = knndm_folds$indx_train,
        indexOut = knndm_folds$indx_test,
        savePredictions = "final",
        verboseIter = FALSE
    )

    # train a random forest model based on kNNDM resampling
    rf_knndm <- train(
        x = st_drop_geometry(samples_smpl)[, names(predictor_stack)],
        y = st_drop_geometry(samples_smpl)[["outcome"]],
        method = "ranger",
        trControl = ctrl_knndm,
        metric = "RMSE",
        num.trees = 300,
        importance = "impurity"
    )

    # Standardize the format of the pointwise errors obtained by kNNDM CV and add an ID column
    tw_pointwise_error(
        obs = rf_knndm$pred$obs,
        pred = rf_knndm$pred$pred,
        id = rf_knndm$pred$rowIndex
    )
})
#> Warning: [as.polygons] only the first layer is polygonized when 'dissolve=TRUE'
#> 1000 prediction points are sampled from the modeldomain
#> Gij <= Gj; a random CV assignment is returned
#> Warning: [as.polygons] only the first layer is polygonized when 'dissolve=TRUE'
#> 1000 prediction points are sampled from the modeldomain
#> Warning: [as.polygons] only the first layer is polygonized when 'dissolve=TRUE'
#> 1000 prediction points are sampled from the modeldomain

names(pe_samples) <- unique(samples$sampling)

weighted_rmse <- lapply(unique(samples$sampling), function(smpling) {
    pe_smpl <- pe_samples[[smpling]]
    w_smpl <- weights_smpl[[smpling]]

    # Weigh the pointwise errors using the weights obtained from raking
    tw_weighted_error_stats(w_smpl, pe_smpl)[["rmse"]]
})
```

The RMSE values obtained by TWCV are 4.664 for the random sampling
design, 6.431 for the biased sampling and 7.24 for the clustered design.
The following Figure shows the weights assigned to each CV error, where
for biased and clustered sampling a higher weight is given to higher CV
errors (hence resulting in a higher weighted error than kNNDM CV along
would produce):

``` r

lapply(unique(samples$sampling), function(smpling) {
    pe_smpl <- pe_samples[[smpling]]
    w_smpl <- weights_smpl[[smpling]]

    PredictionMatching:::plot.twcv(w_smpl, pointwise_error = pe_smpl)[[2]] +
        ggtitle(smpling) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
}) |>
    plot_grid(plotlist = _, nrow = 1)
```

![Figure 6: The weights assigned to each CV error in the different
sampling
designs.](Prediction-domain-adaptive-validation_files/figure-html/unnamed-chunk-13-1.png)

Figure 6: The weights assigned to each CV error in the different
sampling designs.

Brenning, Alexander, and Thomas Suesse. 2026. *Aligning Validation with
Deployment: Target-Weighted Cross-Validation for Spatial Prediction*.
arXiv. <https://doi.org/10.48550/ARXIV.2603.29981>.

Linnenbrink, Jan, Carles Milà, Marvin Ludwig, and Hanna Meyer. 2024.
“kNNDM CV: *K* -Fold Nearest-Neighbour Distance Matching
Cross-Validation for Map Accuracy Estimation.” *Geoscientific Model
Development* 17 (15): 5897–912.
<https://doi.org/10.5194/gmd-17-5897-2024>.

Milà, Carles, Jorge Mateu, Edzer Pebesma, and Hanna Meyer. 2022.
“Nearest neighbour distance matching leave-one-out cross-validation for
map validation.” *Methods in Ecology and Evolution* 13 (6): 1304–16.
<https://doi.org/10.1111/2041-210X.13851>.

Wang, Yanwen, Mahdi Khodadadzadeh, and Raúl Zurita-Milla. 2023.
“Spatial+: A New Cross-Validation Method to Evaluate Geospatial Machine
Learning Models.” *International Journal of Applied Earth Observation
and Geoinformation* 121: 103364.
https://doi.org/<https://doi.org/10.1016/j.jag.2023.103364>.

Wang, Yanwen, Mahdi Khodadadzadeh, and Raúl Zurita-Milla. 2025a. “A
Dissimilarity-Adaptive Cross-Validation Method for Evaluating Geospatial
Machine Learning Predictions with Clustered Samples.” *Ecological
Informatics* 90: 103287.
https://doi.org/<https://doi.org/10.1016/j.ecoinf.2025.103287>.

Wang, Yanwen, Mahdi Khodadadzadeh, and Raúl Zurita-Milla. 2025b. “On the
Use of Adversarial Validation for Quantifying Dissimilarity in
Geospatial Machine Learning Prediction.” *GIScience & Remote Sensing* 62
(1): 2460513. <https://doi.org/10.1080/15481603.2025.2460513>.
