.onAttach <- function(libname, pkgname) {
    version <- read.dcf(file = system.file("DESCRIPTION", package = pkgname), fields = "Version")
    packageStartupMessage(" ")
    packageStartupMessage("###################################################################")
    packageStartupMessage("This is ", paste(pkgname, version))
    packageStartupMessage("See the NEWS file for comparison with deprecated semTools features.")
    packageStartupMessage("###################################################################")
}


