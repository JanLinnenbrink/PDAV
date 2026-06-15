# PDAV

This repository is a collection of prediction-domain adaptive validation
methods. Due to the growing application of spatial predictive models in
geoscientific fields, there is also a growing need for reliable
validation of the resulting maps. Prediction-domain adaptive validation
methods provide reliable proxies of map accuracies that can be used
during model selection and also – in the absence of an independent
probability sample – as a proxy of the final map accuracy. This
repository aims at:

1.  collecting and consistently
    [implementing](https://janlinnenbrink.github.io/PDAV/R/)
    prediction-domain adaptive validation methods applicable for spatial
    predictive modelling.
2.  providing
    [vignettes](https://janlinnenbrink.github.io/PDAV/vignettes/) that
    describe their functioning and application.
3.  providing an
    [overview](https://janlinnenbrink.github.io/PDAV/vignettes/Prediction-domain-adaptive-validation.Rmd)
    over the developed methods.

As such, the repository is expected to grow and include newly developed
methods falling in the class of prediction-domain adaptive validation.
It is still work-in-progress, so please reach out if you find
inaccuracies or want to add methods.

## Overview over the currently developed methods:

Below, you can find a technical comparison of the different approaches:

[TABLE]

The table below summarises some advantages and disadvantages of the
different methods (based on my subjective impression):

[TABLE]

## List of research papers

- Milà, C., Mateu, J., Pebesma, E., Meyer, H. (2022): Nearest Neighbour
  Distance Matching Leave-One-Out Cross-Validation for map validation.
  Methods in Ecology and Evolution 00, 1– 13.
  <https://doi.org/10.1111/2041-210X.13851>

- Linnenbrink, J., Milà, C., Ludwig, M., and Meyer, H. (2024): kNNDM:
  k-fold Nearest Neighbour Distance Matching Cross-Validation for map
  accuracy estimation. GMD, 17, 5897–5912.
  <https://doi.org/10.5194/gmd-17-5897-2024>

- Wang, Y., Khodadadzadeh, M. and Zurita-Milla, R. (2025): A
  dissimilarity-adaptive cross-validation method for evaluating
  geospatial machine learning predictions with clustered samples.
  Ecological Informatics 90, 1574-9541.
  <https://doi.org/10.1016/j.ecoinf.2025.103287>

- Brenning, A., Suesse, T (2026): Aligning Validation with Deployment in
  Spatial Prediction: Target-Weighted Cross-Validation. Preprint.
  <https://arxiv.org/abs/2603.29981>
