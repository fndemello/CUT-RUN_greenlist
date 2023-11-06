# R script file

This folder refers to the scripts used throughout the manuscript, for Greenlist creation 
and testing, made available for the purpose of reproducibility. The original commit presents
the scripts as they were used, while later commits are aimed at enhancing understandability
without changing the core content (aka my code is a mess and I almost never use comments, so
I'm cleaning it all up now)

**norm+entropy.R** quantile normalizes pre-computed bin counts (outputted from [deeptools multiBamSummary](https://deeptools.readthedocs.io/en/develop/content/tools/multiBamSummary.html))
and calculates Shannon's entropy of each bin, using the Maximum Likelihood approach from the [entropy package](https://cran.r-project.org/web/packages/entropy/index.html)

**threshold_test.R** tests the performance of different normalizing sets (Figure 2c of the manuscript)

**montecarlo.R** runs a Monte-Carlo test using random bins as normalizers (Figure 3c of the manuscript)
