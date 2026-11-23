# np-OE-inference

This repository contains R code for a Monte Carlo simulation assessing Overlap Excluded (OE) participation probabilities for Non-Probability Samples.  The simulation compares the performance of the OE approach with other methods in the literature.

Settings for the simulation can be changed as required:
* Lines 278-280: Bootstrap interations, population size, and number of simulation iterations
* Lines 300-302: Correlation between Y and X*beta, sample size of A, and sample size of B
* Line 306: Number of multiple imputes when designing optimal sample

The code is used in the simulation study described in:

Ang, L., R. Clark, B. Loong, and A. Holmberg. 2026. “Overlap-Excluded Estimation of Participation Probabilities for Non-Probability Samples”. In review.

Ang, L., R. Clark, B. Loong, and A. Holmberg. 2024. “Estimating Propensities of Selection for Big Datasets via Data Integration”. Presentation at the 7th International Conference on Establishment Statistics, Glasgow, United Kingdom (https://arxiv.org/abs/2501.04185)
