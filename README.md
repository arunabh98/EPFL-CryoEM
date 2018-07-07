# EPFL-CryoEM

This is an implementation of the [Relion](https://www.sciencedirect.com/science/article/pii/S1047847712002481#s0085) algorithm in 2D. We frame the problem of estimating the structure in the bayesian framework and implement a expected maximization algorithm to correctly estimate the model. In particular, the reconstruction problem is formulated as finding the model that has the highest probability of being the correct one in the light of both the observed data and available prior information. Optimization of this posterior distribution is called maximum a posteriori (MAP), or regularized likelihood optimization.

## Mathematical details of the algorithm
TO have a look at the derivation of all the formulas, and how the expectation maximation algorithm has been formulated, you may have a look at [reports/math-relion.pdf](https://github.com/Arunabh98/EPFL-CryoEM/blob/master/reports/math-relion.pdf). 

## Implementation details
The main algorithm resides in [scripts/MAP2D.m](https://github.com/Arunabh98/EPFL-CryoEM/blob/master/scripts/MAP2D.m) which, from the projections of the image, correctly estimates the orienation of each of those projections and then reconstructs the image using the formulae derived in the expectation maximization algorithm. The back-projection and the projection algorithms algortihms have been implemented in [scripts/backproject_fourier_alternate.m](https://github.com/Arunabh98/EPFL-CryoEM/blob/master/scripts/backproject_fourier_alternate.m) and [scripts/project_fourier_alternate.m](https://github.com/Arunabh98/EPFL-CryoEM/blob/master/scripts/project_fourier_alternate.m). The probability of each projecttion having a particular orientation is calculated by [scripts/calc_prob_for_each_orientation.m](https://github.com/Arunabh98/EPFL-CryoEM/blob/master/scripts/calc_prob_for_each_orientation.m).

## Credits
This work would not have been possible without the guidance of professor [Victor Panaretos](http://smat.epfl.ch/victor/) and the resources provided to me by [EPFL](https://www.epfl.ch/). Also, some of the referenes used in this work are -
* [A Bayesian View on Cryo-EM Structure Determination](https://www.sciencedirect.com/science/article/pii/S0022283611012290)
* [RELION: Implementation of a Bayesian approach to cryo-EM structure determination](https://www.sciencedirect.com/science/article/pii/S1047847712002481#b0155)
* For an idea about how to go about the projection and back-projection operation - [Direct Fourier Reconstruction of a Tomographic Slice](https://uk.mathworks.com/matlabcentral/fileexchange/60257-direct-fourier-reconstruction-of-a-tomographic-slice)
* For the interpolation-scheme using least-squares approach - [inpaint_nans](https://uk.mathworks.com/matlabcentral/fileexchange/4551-inpaint-nans)