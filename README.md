# EPFL-CryoEM

This is an implementation of the [Relion]{https://www.sciencedirect.com/science/article/pii/S1047847712002481#s0085} algorithm in 2D. We frame the problem of estimating the structure in the bayesian framework and implement a expected maximization algorithm to correctly estimate the model. In particular, the reconstruction problem is formulated as finding the model that has the highest probability of being the correct one in the light of both the observed data and available prior information. Optimization of this posterior distribution is called maximum a posteriori (MAP), or regularized likelihood optimization.

## Mathematical details of the algorithm
TO have a look at the derivation of all the formulas, and how the expectation maximation algorithm has been formulated, you may have a look at [reports/math-relion.pdf]{}

## Implementation details
The main algorithm resides in [scripts/MAP2D.m]{https://github.com/Arunabh98/EPFL-CryoEM/blob/master/scripts/MAP2D.m} which, from the projections of the image, correctly estimates the orienation of each of those projections and then reconstructs the image using the formulae derived in the expectation maximization algorithm. To hva