# EPFL-CryoEM

This repository contains the implementation of two leading software in Cryo-EM - [Relion](https://www.sciencedirect.com/science/article/pii/S1047847712002481#s0085) and [CryoSPARC](https://www.nature.com/articles/nmeth.4169) in 2-D. We frame the entire problem in a Bayesian framework. In particular, the reconstruction problem is formulated as finding the model that has the highest probability of being the correct one in the light of both the observed data and available prior information. We optimize it using stochastic gradient descent to arrive at an initial estimate followed by an expected maximization algorithm to refine the model.  

## Mathematical details of the algorithm
TO have a look at the derivation of all the formulas, how the SGD is implemented and how the expectation maximation algorithm has been formulated, you may have a look at [reports/cryosparc-algorithm](https://github.com/Arunabh98/EPFL-CryoEM/blob/master/reports/cryosparc-algorithm.pdf) for CryoSPARC and [reports/math-relion.pdf](https://github.com/Arunabh98/EPFL-CryoEM/blob/master/reports/math-relion.pdf) for RELION. 

## Implementation details
The CryoSPARC algorithm resides in [scripts/cryoSPARC.m](https://github.com/Arunabh98/EPFL-CryoEM/blob/master/scripts/cryoSPARC.m) which produces an initial estimate for the image by iteratively gradient descending to the right structure.

The RELION algorithm resides in [scripts/MAP2D.m](https://github.com/Arunabh98/EPFL-CryoEM/blob/master/scripts/MAP2D.m) which, from the projections of the image, correctly estimates the orientation of each of those projections and then reconstructs the image using the formulae derived in the expectation maximization algorithm.

The back-projection and the projection algorithms algortihms have been implemented in [scripts/backproject_fourier_alternate.m](https://github.com/Arunabh98/EPFL-CryoEM/blob/master/scripts/backproject_fourier_alternate.m) and [scripts/project_fourier_alternate.m](https://github.com/Arunabh98/EPFL-CryoEM/blob/master/scripts/project_fourier_alternate.m). The probability of each projection having a particular orientation is calculated by [scripts/calc_prob_for_each_orientation.m](https://github.com/Arunabh98/EPFL-CryoEM/blob/master/scripts/calc_prob_for_each_orientation.m).

## Credits
This work would not have been possible without the guidance of professor [Victor Panaretos](http://smat.epfl.ch/victor/) and the resources provided to me by [EPFL](https://www.epfl.ch/). Also, some of the references used in this work are -
* [A Bayesian View on Cryo-EM Structure Determination](https://www.sciencedirect.com/science/article/pii/S0022283611012290)
* [RELION: Implementation of a Bayesian approach to cryo-EM structure determination](https://www.sciencedirect.com/science/article/pii/S1047847712002481#b0155)
* For an idea about how to go about the projection and back-projection operation - [Direct Fourier Reconstruction of a Tomographic Slice](https://uk.mathworks.com/matlabcentral/fileexchange/60257-direct-fourier-reconstruction-of-a-tomographic-slice)
* For the interpolation-scheme using least-squares approach - [inpaint_nans](https://uk.mathworks.com/matlabcentral/fileexchange/4551-inpaint-nans)
