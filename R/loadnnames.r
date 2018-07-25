##' Run SS function
##' @param file file name to load
##' @author Jason Cope
##' @export

loadnnames <- function (file, ...)
{
  ls.ext <- function(file) {
    local({
      base::load(file)
      base::ls()
    })
  }
  base::load(file, .GlobalEnv, ...)
  ls.ext(file)
}

