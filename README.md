# Reciprocity Toolbox

This Matlab toolbox provides a set of functions to implement EEG-informed optimization of the tDCS montage, as outlined in:
Dmochowski, J. P., Koessler, L., Norcia, A. M., Bikson, M., & Parra, L. C. (2017). Optimal use of EEG recordings to target active brain areas with transcranial electrical stimulation. NeuroImage.

The best place to start is generateFigure.m, which demonstrates the usage of the core functions:
  * reciprocate.m (unconstrained reciprocity)
  * tibshirani.m (L^1 constrained reciprocity)
while also generating Figure 3 of the manuscript.

Dependencies:
  * you must have the EEGLAB function topoplot, and its dependencies on your system, in order to visualize the resulting montage.  
  * NB: EEGLAB is not required to implement reciprocity!
  
To solve the L^1 constrained least squares problem, we implemented the technique described in:

Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. Journal of the Royal Statistical Society. Series B (Methodological), 267-288.

Any questions, comments, or bug reports should be directed to Jacek Dmochowski (jdmochowski@ccny.cuny.edu).
