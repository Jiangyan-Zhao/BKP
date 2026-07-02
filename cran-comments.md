## Resubmission

This is a resubmission. I addressed the CRAN comments about examples/tests
using more than two cores by default.

Changes made:

* All package-controlled OpenMP loops are explicitly controlled by the
  existing `n_threads` argument, whose default is `1`.
* Examples and tests were revised to avoid unnecessary hyperparameter
  optimization where fixed kernel lengthscales are sufficient to demonstrate
  the user interface.
* Small kernel matrices are evaluated using the package's serial C++ loop
  engine to avoid BLAS/GEMM thread startup overhead during CRAN
  examples/tests.
* Users can still choose the number of package-controlled OpenMP threads
  through the `n_threads` argument.
* Added `inst/WORDLIST` for package-specific terms flagged by spell checking.
