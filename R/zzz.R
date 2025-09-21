.onLoad <- function(libname, pkgname) {
  if (identical(Sys.getenv("_R_CHECK_LIMIT_CORES_"), "TRUE")) {
    # Unset problematic KMP variables
    Sys.unsetenv(c("KMP_DEVICE_THREAD_LIMIT", "KMP_ALL_THREADS", "KMP_TEAMS_THREAD_LIMIT"))

    # Limit threads for CRAN checks
    Sys.setenv(
      OMP_THREAD_LIMIT = "2",
      OMP_NUM_THREADS = "1",
      OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1",
      VECLIB_MAXIMUM_THREADS = "1",
      RCPP_PARALLEL_NUM_THREADS = "1"
    )
  }
}
