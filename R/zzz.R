#' @importFrom glue glue
.onLoad <- function(libname, pkgname) {
  path <- find.package('merlot')
  if(!dir.exists(glue::glue("{path}/venv"))) {
    reinstall()
  }
}

#' @importFrom glue glue
reinstall <- function() {
  path <- find.package('merlot')
  system(glue::glue("bash {path}/make {path}/venv"))
}
