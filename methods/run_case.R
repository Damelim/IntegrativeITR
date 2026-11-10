# Load shared code for exactly one setting/condition directory.
run_case <- function(setting, condition, project_directory, case_directory) {
    run_environment <- new.env(parent = globalenv())
    sys.source(file.path(project_directory, "setup.R"), envir = run_environment)
    code_directory <- file.path(project_directory, "methods")
    run_environment$.opt_dir <- code_directory
    for (filename in c("design.R", "methods.R", "optimizer.R", "fast_backend.R",
                       "fitting.R", "plots.R", "simulation.R")) {
        sys.source(file.path(code_directory, filename), envir = run_environment)
    }
    configuration <- mget(
        c("replications", "start_seed", "cores", "methods", "output_folder"),
        envir = run_environment, inherits = FALSE
    )
    folder <- configuration$output_folder
    if (!is.character(folder) || length(folder) != 1L || is.na(folder) ||
        !grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", folder)) {
        stop("output_folder must be a simple folder name, for example 'output'.")
    }
    expected_directory <- normalizePath(file.path(project_directory, setting, condition),
                                        mustWork = TRUE)
    if (!identical(normalizePath(case_directory), expected_directory)) {
        stop("The run.R setting/condition must match its directory.")
    }
    configuration$setting <- setting
    configuration$condition <- condition
    configuration$output_directory <- file.path(case_directory, folder)
    run_environment$run_simulation(configuration, project_directory)
}
