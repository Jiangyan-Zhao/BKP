## Resubmission

This is a resubmission. I addressed the CRAN comments about examples/tests
using more than two cores by default.

Changes made:

* All package-controlled OpenMP loops are explicitly controlled by the
  existing `n_threads` argument, whose default is `1`.
* Examples and tests were revised to avoid unnecessary hyperparameter
  optimization where fixed kernel lengthscales are sufficient to demonstrate
  the user interface.
* Large kernel matrices are evaluated using the package's serial C++ loop
  engine to limit peak memory, while moderate-size kernel matrices retain the
  GEMM path for computational efficiency.
* Tests involving kernel matrix construction were reduced in size to avoid
  unnecessary CPU usage during CRAN checks.
* Users can still choose the number of package-controlled OpenMP threads
  through the `n_threads` argument.
* Added `inst/WORDLIST` for package-specific terms flagged by spell checking.
