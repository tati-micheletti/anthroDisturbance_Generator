# tests/testthat/test-reqdPkgs.R
# --- Helpers ---------------------------------------------------------------

# Package names declared in a module's reqdPkgs, without GitHub owner, branch or version
moduleReqdPkgs <- function(moduleFile) {
  findReqd <- function(e) {
    if (is.call(e)) {
      if (!is.null(names(e)) && "reqdPkgs" %in% names(e)) return(e[["reqdPkgs"]])
      for (a in as.list(e)[-1L]) {
        r <- findReqd(a)
        if (!is.null(r)) return(r)
      }
    }
    NULL
  }
  for (e in parse(moduleFile, keep.source = FALSE)) {
    reqd <- findReqd(e)
    if (!is.null(reqd)) break
  }
  pkgs <- unlist(eval(reqd, baseenv()))
  sub("@.*$", "", sub("^.*/", "", sub("[ (].*$", "", pkgs)))
}

# Walk the parsed code in `files` (comments and strings excluded; no parse data needed, which
# getParseData() would not give when the keep.parse.data option is FALSE, as it can be under covr).
# Returns the packages called as pkg::fun() or pkg:::fun(), without base R's own; the packages the
# code checks for with requireNamespace("pkg"), which it treats as optional; and the calls that
# attach packages (Require(), library(), require(), attachNamespace()), which SpaDES should do.
modulePkgCalls <- function(files) {
  basePkgs <- c("base", "utils", "stats", "methods", "tools", "parallel", "grDevices", "graphics")
  attachers <- c("Require", "library", "require", "attachNamespace")
  qualified <- character(0)
  optional <- character(0)
  attaching <- character(0)
  walk <- function(e, where) {
    if (is.pairlist(e)) e <- as.list(e) # function arguments and their defaults
    if (is.call(e)) {
      fn <- e[[1L]]
      if (identical(fn, as.name("::")) || identical(fn, as.name(":::"))) {
        qualified <<- c(qualified, as.character(e[[2L]]))
        if (as.character(e[[3L]]) %in% attachers) attaching <<- c(attaching, where)
      }
      if (is.name(fn) && as.character(fn) %in% attachers) attaching <<- c(attaching, where)
      if (identical(fn, as.name("requireNamespace")) && length(e) > 1L && is.character(e[[2L]])) {
        optional <<- c(optional, e[[2L]])
      }
    }
    if (is.call(e) || is.list(e)) {
      for (i in seq_along(e)) {
        if (!identical(e[[i]], quote(expr = ))) walk(e[[i]], where) # skip empty arguments, as in x[, 1]
      }
    }
  }
  for (f in files) for (e in parse(f, keep.source = FALSE)) walk(e, basename(f))
  list(qualified = setdiff(unique(qualified), basePkgs), optional = unique(optional),
       attaching = attaching)
}

moduleDir <- testthat::test_path("..", "..")
moduleFile <- file.path(moduleDir, "anthroDisturbance_Generator.R")
moduleCode <- c(moduleFile, list.files(file.path(moduleDir, "R"), "[.]R$", full.names = TRUE))

# --- Tests -----------------------------------------------------------------

testthat::test_that("reqdPkgs declares every package the module code calls with pkg::, unless it is optional", {
  calls <- modulePkgCalls(moduleCode)
  testthat::expect_true("terra" %in% calls$qualified) # the scan found the module's terra:: calls
  reqd <- moduleReqdPkgs(moduleFile)
  testthat::expect_equal(setdiff(calls$qualified, c(reqd, calls$optional)), character(0))
})

testthat::test_that("the module code never attaches packages itself; SpaDES attaches reqdPkgs", {
  calls <- modulePkgCalls(moduleCode)
  # files that call Require(), library(), require() or attachNamespace()
  testthat::expect_equal(unique(calls$attaching), character(0))
})
