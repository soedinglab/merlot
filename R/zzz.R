#' @importFrom glue glue
.onLoad <- function(libname, pkgname) {
  if (!dir.exists(system.file(package = "merlot", "venv"))) {
    reinstall()
  }
}

#' @importFrom glue glue
reinstall <- function() {
  is_checking <- find.package("merlot") != system.file(package = "merlot")
  if (is_checking) {
    venv_path <- tempfile("venv")
    on.exit(unlink(venv_path, recursive = TRUE, force = TRUE))
  } else {
    venv_path <- system.file(package = "merlot", "venv")
  }

  make_path <- system.file(package = "merlot", "make")
  system(glue::glue("bash {make_path} {venv_path}"))
}
