# BKP 0.2.1 (2025-09-17)

* CRAN Linux fix: limit threads to avoid OMP warnings

# BKP 0.2.0 (2025-09-16)

* Added `fitted()`, `parameter()`, and `quantile()` methods.
* Updated `predict()` and `simulate()` methods: both now return results for the training data by default when `Xnew` is not provided.  
* Extended `plot()` method with new `dims` argument for higher-dimensional inputs.
* Added Section 5 to the vignette, presenting a real-data application on Loa loa parasite infection in North Cameroon.  
* Added argument checking with informative error messages.

# BKP 0.1.1 (2025-08-19)

* Added a vignette introducing the package.  
* Fixed minor bugs and improved stability.  

# BKP 0.1.0 (2025-07-23)

* Initial release on CRAN.  
