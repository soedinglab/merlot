#' @importFrom glue glue
.onLoad <- function(libname, pkgname) {
  if (!dir.exists(venv_path())) {
    reinstall()
  }
}

#' @importFrom glue glue
reinstall <- function() {
  is_checking <- find.package("merlot") != system.file(package = "merlot")

  # if the package is being run during an R CMD CHECK,
  # the venv needs to be created in a temporary folder
  if (is_checking) {
    venv_path <- tempfile("venv")
    on.exit(unlink(venv_path, recursive = TRUE, force = TRUE))
  } else {
    venv_path <- venv_path()
  }

  make_path <- system.file(package = "merlot", "make")
  system(glue::glue("bash {make_path} {venv_path}"))
}

venv_path <- function() {
  paste0(find.package("merlot"), "/venv")
}
