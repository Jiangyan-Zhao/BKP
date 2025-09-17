.onLoad <- function(libname, pkgname) {
  if (identical(Sys.getenv("_R_CHECK_LIMIT_CORES_"), "TRUE")) {
    Sys.setenv(
      "OMP_THREAD_LIMIT" = 2,
      "OMP_NUM_THREADS" = 1,
      "OPENBLAS_NUM_THREADS" = 1,
      "MKL_NUM_THREADS" = 1,
      "VECLIB_MAXIMUM_THREADS" = 1,
      "RCPP_PARALLEL_NUM_THREADS" = 1
    )
  }
}
