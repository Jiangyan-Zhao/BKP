## R CMD check results

0 errors | 0 warnings | 0 note

* This is a new release.

## revdepcheck results

We checked 0 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

* We saw 0 new problems
* We failed to check 0 packages

## Additional comments

This is a maintenance and feature update for the BKP package.  
All changes are backward compatible and have been tested on major platforms.

We have added a fix in the vignettes to limit BLAS/OpenMP threads during vignette builds 
to avoid "OMP: Warning #96" on CRAN Linux platforms.  
All examples and vignettes run without errors or warnings on local checks.

