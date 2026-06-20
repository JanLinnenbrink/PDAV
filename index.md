# PDAV

Prediction-domain adaptive validation methods provide estimates of map accuracy that are adapted to the deployment or prediction domain ([Linnenbrink, Nowosad & Meyer, 2026](https://doi.org/10.48550/arXiv.2605.13689)).

The package aims to:

1. provide consistent implementations of prediction-domain adaptive validation methods for spatial predictive modelling;
2. provide vignettes that describe how the methods work and how they can be applied;
3. provide an [overview](articles/Prediction-domain-adaptive-validation.html) of currently available methods.

## Available methods

| Method | Main idea | Reference | Documentation |
|---|---|---|---|
| NNDM | Leave-one-out nearest-neighbour distance matching. Matches the nearest-neighbour distance distribution encountered during prediction by excluding training points from the LOO CV. | [Milà et al. (2022)](https://doi.org/10.1111/2041-210X.13851) | [Function](reference/nndm.html), [Article](articles/NNDM.html) |
| kNNDM | k-fold nearest-neighbour distance matching. Uses clustering to create fold configurations ranging from random to spatially clustered CV and selects the configuration that best approximates the prediction situation. | [Linnenbrink & Milà et al. (2024)](https://doi.org/10.5194/gmd-17-5897-2024) | [Function](reference/knndm.html), [Article](articles/k-NNDM.html) |
| DA-CV | Uses adversarial validation to classify the prediction area into (dis-)similar to the training data, and then weights random and spatial CV results based on the proportion of the two areas. | [Wang et al. (2025)](https://doi.org/10.1016/j.ecoinf.2025.103287) | [Function](reference/da_cv.html), [Article](articles/DA-CV.html) |
| TWCV | Uses raking to weight the error estimates at the training points to match the frequency distribution of the predictor variables at the prediction locations. | [Brenning & Suesse (2026)](https://doi.org/10.48550/arXiv.2603.29981) | [Function](reference/compute_cv_estimators.html), [Article](articles/twcv.html) |


## Detailed comparison

For a more detailed comparison of advantages, disadvantages, and assumptions, see the [overview article](articles/Prediction-domain-adaptive-validation.html).