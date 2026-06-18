# PDAV
This repository is a collection of prediction-domain adaptive validation methods. Due to the growing application of spatial predictive models in geoscientific fields, there is also a growing need for reliable validation of the resulting maps. Prediction-domain adaptive validation methods provide reliable proxies of map accuracies that can be used during model selection and also -- in the absence of an independent probability sample -- as a proxy of the final map accuracy ([Linnenbrink, Nowosad & Meyer, 2026](https://doi.org/10.48550/arXiv.2605.13689)). This repository aims at:

1) collecting and consistently [implementing](R/) prediction-domain adaptive validation methods applicable for spatial predictive modelling.
2) providing [vignettes](vignettes/) that describe their functioning and application.
3) providing an [overview](vignettes/Prediction-domain-adaptive-validation.Rmd) over the developed methods.

As such, the repository is expected to grow and include newly developed methods falling in the class of prediction-domain adaptive validation.
It is still work-in-progress, so please reach out if you find inaccuracies or want to add methods.


## Overview over the currently developed methods:

Below, you can find a technical comparison of the different approaches:

| Method | Authors | Short Description | Critical Parameters |
| ----------- | ----------- | ----------- | ----------- |
| [NNDM](R/nndm.R) | [Milà et al. (2022)](https://doi.org/10.1111/2041-210X.13851) | - Formalizes the prediction situation as the Nearest Neighbour Distance (**NND**) distribution between prediction locations and training samples.<br>- Then, for every k in N, it calculates the NND between the training samples and the hold-out sample.<br>- It **excludes the training point that is clostest to the held-out point** until the NND between the training samples and the held-out point matches the NND between prediction points and training samples. | - `phi`: Autocorrelation threshold up to which Nearest-Neighbour distances are matched<br>- `min_train`: fold balancing |
| [kNNDM](R/knndm.R) | [Linnenbrink & Milà et al. (2024)](https://doi.org/10.5194/gmd-17-5897-2024) | - Formalizes the prediction situation as the **NND** distribution between prediction locations and training samples.<br>- Uses clustering to create a **continuum of fold configurations** ranging from random resampling to a blocked split.<br>- Then selects the configuration that best approximates the prediction situation. | - `maxp`: maximum fold size allowed: higher numbers lead to more imbalanced splits, but also may resemble the prediction situation better |
| [DA-CV](R/da_cv.R) |  [Wang et al. (2025)](https://doi.org/10.1016/j.ecoinf.2025.103287) | - Uses **adversial validation** to classify the prediction area into locations that are similar or dissimilar to the training samples.<br>- Returns a random and a spatial resampling split.<br>- The validation statistics obtained by these two splits are then **weighted by the proportion of similar or dissimilar areas**, respectively. | - `autoc_threshold`: block size in spatial cross-validation |
| [TWCV](R/twcv.R) |  [Brenning & Suisse (2026)](https://doi.org/10.48550/arXiv.2603.29981) | - Uses **raking** to weight the error estimates at the training points to match the frequency distribution of the predictor variables at the prediction locations | -  |

The table below summarises some advantages and disadvantages of the different methods (based on my subjective impression):

| Method | Advantages | Disadvantages |
| ----------- | ----------- | ----------- |
| [NNDM](R/nndm.R) | - **direct quantification of the prediction situation** through nearest-neighbour distances<br>- **NND plots** reveal situations where prediction situation cannot be matched during CV<br>- works natively with **distances in the geographic or feature space**  | - **high computational costs** since it uses LOO-CV<br>- **May remove large fractions of the training samples** and thus lead to unstable models<br>- depends on the estimated autocorrelation range (though not so important as in spatial+ CV used in DA-CV) |
| [kNNDM](R/knndm.R) | - **direct quantification of the prediction situation** through nearest-neighbour distances<br>- **NND plots** reveal situations where prediction situation cannot be matched during CV<br>- more computational efficient than NNDM<br>- works natively with **distances in the geographic or feature space** | - **less flexible** in matching the NND distributions than NNDM |
| [DA-CV](R/da_cv.R) | - creates a **map that depicts areas of inter- vs extrapolation**, which can be useful for subsequent sampling or for uncertainty assessment | - **less direct** approximation of the prediction situation<br>- No inspection method available to **assess if the prediction situation can be matched** during CV<br>- Matches the prediction situation in the feature space, but **mixes spaces** in the spatial+ resampling method<br>- requires calculating the weighted **average of two validation statistics**, which might be problematic for metrics other than RMSE |
| [TWCV](R/twcv.R) | - applies **weighting** to match the prediction situation, which is more flexible than resampling-based matching of the prediction situation | - **integration of geographic distances** implicit and weighting of the predictor vs geographical space is automatically done by the raking algorithm |

## List of research papers

* Linnenbrink, J., Nowosad, J. and Meyer, H. (2026): Moving beyond spatial and random cross-validation in environmental modelling: a call for prediction-domain adaptive evaluation. Preprint. https://doi.org/10.48550/arXiv.2605.13689

* Milà, C., Mateu, J., Pebesma, E. and Meyer, H. (2022): Nearest Neighbour Distance Matching Leave-One-Out Cross-Validation for map validation. Methods in Ecology and Evolution 00, 1– 13.
https://doi.org/10.1111/2041-210X.13851

* Linnenbrink, J., Milà, C., Ludwig, M., and Meyer, H. (2024): kNNDM: k-fold Nearest Neighbour Distance Matching Cross-Validation for map accuracy estimation. GMD, 17, 5897–5912.
https://doi.org/10.5194/gmd-17-5897-2024

* Wang, Y., Khodadadzadeh, M. and Zurita-Milla, R. (2025): A dissimilarity-adaptive cross-validation method for evaluating geospatial machine learning predictions with clustered samples. Ecological Informatics 90, 1574-9541.
https://doi.org/10.1016/j.ecoinf.2025.103287

* Brenning, A., Suesse, T (2026): Aligning Validation with Deployment in Spatial Prediction: Target-Weighted Cross-Validation. Preprint.
https://arxiv.org/abs/2603.29981
