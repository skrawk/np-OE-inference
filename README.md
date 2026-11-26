# np-OE-inference

This repository contains R code for a Monte Carlo simulation assessing Overlap Excluded (OE) participation probabilities for Non-Probability Samples.  The simulation compares the performance of the OE approach with other methods in the literature.

## OE Simulation Parallel.R

Settings for the simulation can be changed as required:
* Lines 278-280: Bootstrap iterations, population size, and number of simulation iterations
* Lines 300-302: Correlation between $Y$ and $X\beta$, sample size of the probability and non-probability samples
* Line 306: Number of multiple imputes when producing results for an optimally designed probability sample

To obtain results for Tables 5.1, 5.2, A.1 and A.2 in the Ang et al (2026) paper, you will need to run the code multiple times under different settings:
* Set the correlation to 0.3 (for Tables 5.1 and 5.2) or 0.7 (for Tables A.1 and A.2) as needed
* Set the sample size of the non-probability sample to the values 750, 2,000, etc as needed
* Set the probability sample size to 500 or 5,000 as needed (Note: For the smaller probability sample, you will run into issues with insufficient sample for the KW approach, so you will need to comment out this method in the code)

For the paper, the population size was set to 200,000 and number of simulations was 2,000.  The bootstrap variance algorithm for the EM-based approach used 200 iterations.

The overlap results in Table 5.3 of the paper can be obtained for each simulation via the *overlap* variable.

The coverage probabilities and average lengths of the confidence intervals shown in Tables 5.4 and 5.5 of the paper are obtained from the *c_rate* and *al* variables in the code (lines 967 and 970).

## Reference

The above code is used in the simulation study described in:

Ang, L., R. Clark, B. Loong, and A. Holmberg. 2026. “Overlap-Excluded Estimation of Participation Probabilities for Non-Probability Samples”. In review. Pre-print available at http://arxiv.org/abs/2611.10962
