# Reuse compiled code across R sessions; the cache stays outside the repository.
load_fast_backend <- function(directory) {
    if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Install the Rcpp package first.")
    backend <- new.env(parent = globalenv())
    previous_flags <- Sys.getenv("PKG_CXXFLAGS", unset = NA_character_)
    on.exit({
        if (is.na(previous_flags)) {
            Sys.unsetenv("PKG_CXXFLAGS")
        } else {
            Sys.setenv(PKG_CXXFLAGS = previous_flags)
        }
    })
    flags <- if (is.na(previous_flags)) "" else previous_flags
    Sys.setenv(PKG_CXXFLAGS = paste(flags, "-ffp-contract=off"))
    cache_root <- Sys.getenv("ITR_CACHE_DIR", unset = tools::R_user_dir("integrative_itr", "cache"))
    makevars <- unique(c(
        Sys.getenv("R_MAKEVARS_USER"), path.expand("~/.R/Makevars"),
        Sys.getenv("R_MAKEVARS_SITE"), file.path(R.home("etc"), "Makevars.site"),
        file.path(R.home("etc"), "Makeconf")
    ))
    makevars <- makevars[file.exists(makevars)]
    signature <- list(
        R = R.version.string, platform = R.version$platform, R_home = R.home(),
        Rcpp = as.character(packageVersion("Rcpp")),
        source = unname(tools::md5sum(file.path(directory, "fast_kernels.cpp"))),
        flags = Sys.getenv(c("PKG_CXXFLAGS", "PKG_CPPFLAGS", "PKG_LIBS", "CXX", "CXXFLAGS")),
        makevars = tools::md5sum(makevars)
    )
    key_file <- tempfile("itr_build_key_")
    on.exit(unlink(key_file), add = TRUE)
    saveRDS(signature, key_file, compress = FALSE)
    cache_directory <- file.path(path.expand(cache_root), unname(tools::md5sum(key_file)))
    if (!dir.create(cache_directory, recursive = TRUE, showWarnings = FALSE) &&
        !dir.exists(cache_directory)) {
        stop("Cannot create the compiled-code cache: ", cache_directory)
    }
    Rcpp::sourceCpp(
        file.path(directory, "fast_kernels.cpp"), env = backend,
        cacheDir = cache_directory, showOutput = FALSE
    )
    backend
}
