.onAttach <- function(libname, pkgname) {
    version <- read.dcf(file = system.file("DESCRIPTION", package = pkgname), fields = "Version")
    packageStartupMessage(" ")
    packageStartupMessage("###################################################################")
    packageStartupMessage("This is ", paste(pkgname, version))
    packageStartupMessage("See the README file on github.com/TDJorgensen/lavaan.mi")
    packageStartupMessage("for a table comparing it with deprecated semTools features.")
    packageStartupMessage("###################################################################")
}


