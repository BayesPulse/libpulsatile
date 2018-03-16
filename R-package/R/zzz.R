
# when loading packages, load runs first, then attach

# Startup messages
.onAttach <- function(libname, pkgname) {

    msg <- paste0("bayespulse\n")

    packageStartupMessage(msg)

}

# Set default url options on load
.onUnload <- function (libpath) {
    library.dynam.unload("bayespulse", libpath)
}
