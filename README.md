This repository contains code and data relevant to the Dysfunctional Family Chinese Restaurant Process, as seen in "Clustering Craters on the
Moon with Dysfunctional Families".

The Introductory Plots folder contains code and PDFs for the data visualization plots presented in Section 1 of the paper (Figures 2 and 3).

The Perm Proposal Modification folder contains code and results for our analysis of the permutation proposal algorithm modification as seen in Section 2.3 of the supplemental material. It includes code and results measuring the speedup obtained, and also provides code to produce Table 1 in the supplement.

The Prior Specifications folder contains scripts used to specify the hyperparameters for the priors on cluster-specific means and covariances. Prior_Testing.R contains these calculations for the full dataset (as seen in Section 5.1 of the main paper), and Further_Priors.R contains similar calculations done for the subsets used in the paper (as detailed in Section 4 of the supplement).

The Radius Testing folder contains scripts used, draws obtained, and plots examined in determining the value of the radius parameter in the neighborhood-modification of the sampler (Section 6.1 of the main paper)

The Rand Index Code folder contains code used and draws obtained in the simulation   study presented in Section 7 of the paper.

The Samplers folder contains the DFCRP Gibbs sampler and supporting code. It includes code for the permutation and alpha Metropolis-Hastings samplers (Sections 1 and 2.2 in the supplement), R code to evaluate the DFCRP prior, and C++ code to evaluate the DFCRP prior (rewritten for efficiency).

The Timing folder contains the timing calculations used to determine the speedup of the neighborhood modification. It also contains Figure 4 in the paper.
