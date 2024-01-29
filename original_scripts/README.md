# R script files

This folder refers to the scripts used throughout the manuscript, for Greenlist creation 
and testing, made available for the purpose of reproducibility. The original commit presents
the scripts as they were used, while later commits are aimed at enhancing understandability
without changing the core content.

**norm+entropy.R** quantile normalizes pre-computed bin counts (outputted from [deeptools multiBamSummary](https://deeptools.readthedocs.io/en/develop/content/tools/multiBamSummary.html))
and calculates Shannon's entropy of each bin, using the Maximum Likelihood approach from the [entropy package](https://cran.r-project.org/web/packages/entropy/index.html)

**threshold_test.R** tests the performance of different normalizing sets (Figure 2c of the manuscript)

**montecarlo.R** runs a Monte-Carlo test using random bins as normalizers (Figure 3c of the manuscript)

**subsampling.R** runs a Monte-Carlo-like test to define *de novo* greenlists from subsamplings of our
total dataset, and compare the similarity of these *de novo* lists with their originals (Figure 3f of the manuscript)

**get_sizeFactors.R** calculates DESeq2 size factors based on a quantification of greenlist regions for a set of samples 
