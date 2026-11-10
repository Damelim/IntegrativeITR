# Run only this setting and condition, using the shared setup.R.
local({
    setting <- "setting1"
    condition <- "strong_shift"

    # Locate this file when called with source() or Rscript.
    source_files <- Filter(
        function(path) !is.null(path) && basename(path) == "run.R",
        lapply(sys.frames(), function(frame) frame$ofile)
    )
    if (length(source_files)) {
        script_path <- tail(source_files, 1L)[[1L]]
    } else {
        argument <- grep("^--file=", commandArgs(), value = TRUE)
        if (length(argument) != 1L) stop("Use source('run.R') or Rscript run.R.")
        script_path <- sub("^--file=", "", argument)
    }
    # source(..., chdir = TRUE) may already have changed to this directory.
    if (!file.exists(script_path) && file.exists(basename(script_path))) {
        script_path <- basename(script_path)
    }
    case_directory <- dirname(normalizePath(script_path, mustWork = TRUE))
    project_directory <- normalizePath(file.path(case_directory, "../.."), mustWork = TRUE)

    runner <- new.env(parent = globalenv())
    sys.source(file.path(project_directory, "methods", "run_case.R"), envir = runner)
    runner$run_case(setting, condition, project_directory, case_directory)
})
