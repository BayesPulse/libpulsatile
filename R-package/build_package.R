# cribbed from Microsoft's LightGBM package
unlink("./inst/include", recursive = TRUE)
unlink("./src/tests", recursive = TRUE)
unlink("./src/singlesubject.cpp", recursive = TRUE)
if (!file.copy("./../include", "./inst/", overwrite = TRUE, recursive = TRUE)) {
  stop("Cannot find folder libpulsatile/include")
}
if (!file.copy("./../tests", "./src/", overwrite = TRUE, recursive = TRUE)) {
  stop("Cannot find folder libpulsatile/tests")
}
if (!file.copy("./../src/singlesubject.cpp", "./src/", overwrite = TRUE, recursive = TRUE)) {
  stop("Cannot find folder libpulsatile/src")
}
#if (!file.copy("./../CMakeLists.txt", "./src/", overwrite = TRUE, recursive = TRUE)) {
#  stop("Cannot find file LightGBM/CMakeLists.txt")
#}
if (!file.exists("./src/_IS_FULL_PACKAGE")) {
  file.create("./src/_IS_FULL_PACKAGE")
}
#system("R CMD build --no-build-vignettes .")
file.remove("./src/_IS_FULL_PACKAGE")

