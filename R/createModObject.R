#' This function checks a sim list for a specific object at a specific time, or
#' searches for it in an input folder (i.e. saved outputs). It simulates the
#' existence of a simList with specific objects at time t.
#'
#' When multiple files match, the most recently modified is chosen.
#'
#' @param data       Name of the object (pattern) to fetch.
#' @param sim        simList; if non-NULL and contains `data`, that object is returned.
#' @param pathInput  Directory to search for saved files (required if sim[[data]] is NULL).
#' @param currentTime  Numeric time stamp used to match file names (zero-padded).
#' @param fun        Loader function (default: `readRDS`).
#' @param returnNULL Logical; if TRUE, missing or load failures return NULL instead of stopping.
#'
#' @return The requested object (from `sim` or disk), or NULL if missing and `returnNULL = TRUE`.
#'
#' @author Tati Micheletti
#' @export
#' @include grepMulti.R
#' @importFrom crayon green magenta red
#' @importFrom raster raster
#' @importFrom SpaDES.core paddedFloatToChar
#' @rdname createModObject

createModObject <- function(data,
                            sim = NULL,
                            pathInput = NULL,
                            currentTime,
                            fun = readRDS,
                            returnNULL = FALSE) {
  # require a source
  if (all(is.null(sim), is.null(pathInput))) {
    stop("Either a simList or a folder containing the data need to be supplied")
  }
  
  # check in-memory
  if (!is.null(sim) && !is.null(sim[[data]])) {
    return(sim[[data]])
  }
  message(crayon::yellow(
    sprintf("%s not supplied by another module. Trying files in %s", data, pathInput)
  ))
  
  # list files
  files <- list.files(pathInput, recursive = TRUE)
  if (length(files) == 0) {
    if (returnNULL) {
      message(crayon::red(sprintf("The file for %s was not found. Returning NULL", data)))
      return(NULL)
    }
    stop(sprintf("Please place the data in the input folder %s", pathInput))
  }
  
  # validate currentTime
  if (!is.numeric(currentTime)) {
    stop("Current time needs to be numeric!")
  }
  
  # match pattern
  paddedTime <- SpaDES.core::paddedFloatToChar(currentTime, padL = 3)
  matches <- grepMulti(x = files, patterns = c(data, paddedTime))
  if (length(matches) == 0) {
    if (returnNULL) {
      message(crayon::red(sprintf("No file found for %s. Returning NULL", currentTime)))
      return(NULL)
    }
    message(crayon::red(sprintf("No file found for %s. Returning NULL", currentTime)))
    return(NULL)
  }
  
  # choose most recently modified when many
  if (length(matches) > 1) {
    paths <- file.path(pathInput, matches)
    infos <- file.info(paths)
    newest <- which.max(infos$mtime)
    chosen <- matches[newest]
    message(crayon::magenta(
      sprintf("Multiple files found for %s at time %s. Choosing most recent: %s",
              data, paddedTime, chosen)
    ))
  } else {
    chosen <- matches
  }
  
  # load chosen file
  fullPath <- file.path(pathInput, chosen)
  dt <- tryCatch(
    do.call(fun, list(fullPath)),
    error = function(e) {
      message(crayon::red(sprintf("Failed to load %s. Returning NULL", chosen)))
      return(NULL)
    }
  )
  if (is.null(dt)) {
    return(NULL)
  }
  message(crayon::green(
    sprintf("%s loaded from %s for time %s", data,
            crayon::magenta(fullPath), paddedTime)
  ))
  return(dt)
}
