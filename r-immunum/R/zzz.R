.onLoad <- function(libname, pkgname) {
  invisible(NULL)
}

.onAttach <- function(libname, pkgname) {
  ver <- tryCatch(immunum_version(), error = function(e) NA_character_)
  if (!is.na(ver)) {
    packageStartupMessage(sprintf("immunum: linked against Rust crate %s", ver))
  }
}
