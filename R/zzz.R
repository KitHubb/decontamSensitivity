.onAttach <- function(libname, pkgname) {
  version <- as.character(utils::packageVersion(pkgname))
  packageStartupMessage(
    pkgname, " ", version,
    " loaded. Start with run_decontam_qc()."
  )
}
