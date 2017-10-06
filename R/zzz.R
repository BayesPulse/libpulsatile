
# when loading packages, load runs first, then attach

# Startup messages
.onAttach <- function(libname, pkgname) {

    msg <- paste0("poppulsatile\n")

    packageStartupMessage(msg)

}

# Set default url options on load
.onUnload <- function (libpath) {
    library.dynam.unload("poppulsatile", libpath)
}
