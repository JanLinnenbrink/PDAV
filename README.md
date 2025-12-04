# PDAV
This repository is a collection of prediction-domain adaptive validation methods. Due to the growing application of spatial predictive models in geoscientific fields, there is also a growing need for reliable validation of the resulting maps. Prediction-domain adaptive validation methods provide reliable proxies of map accuracies that can be used during model selection and also -- in the absence of an independent probability sample -- as a proxy of the final map accuracy. This repository aims at:

1) collecting and consistently implementing prediction-domain adaptive validation methods applicable for spatial predictive modelling.
2) providing vignettes that describe their functioning and application.
3) providing an overview over the developed methods.
4) comparing the performance of the validation methods as compared to random and spatial CV on a common benchmark dataset.

As such, the repository is expected to grow and include newly developed methods falling in the class of prediction-domain adaptive validation.


## Overview over the currently developed methods:

Below, you can find a technical comparison of the different approaches.

| Method      | Authors | Theoretical basis | Space | Critical Parameters |
| ----------- | ----------- | ----------- | ----------- | ----------- |
| NNDM      | [Milà et al. (2022)](https://doi.org/10.1111/2041-210X.13851) | Point Patterns | Geographical or Feature | - *phi*<br>- *min_train*: fold balancing |
| kNNDM   | [Linnenbrink & Milà et al. (2024)](https://doi.org/10.5194/gmd-17-5897-2024) | Point Patterns | Geographical or Feature | - *maxp*: fold balancing |
| DA-CV      |  [Wang et al. (2025)](https://doi.org/10.1016/j.ecoinf.2025.103287) | Adversial Validation | Feature for matching the prediction situation, but also mixes in geographical space during fold creation | - *autoc_threshold*: block size in spatial cross-validation |

As a guideline for choosing an appropriate method, the table below summarises some advantages and disadvantages of the different methods.

| Method      | Advantages | Disadvantages |
| ----------- | ----------- | ----------- |
| NNDM      | - direct quantification of the prediction situation through nearest-neighbour distances | - high computational costs since it uses LOO-CV |
| kNNDM   | - direct quantification of the prediction situation through nearest-neighbour distances<br>- more computational efficient than NNDM | ab |
| DA-CV      | - creates a map that depicts areas of inter- vs extrapolation, which can be useful for subsequent sampling or for uncertainty assessment<br>- | - less direct approximation of the prediction situation |


## Benchmarking

The following figure shows the results of the benchmarking study:

<img src="man/figures/rmse_comp.png" width="100%" />

This figure shows the computational costs of each method:
