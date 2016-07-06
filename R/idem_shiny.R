
#' Run Web-Based \code{idem} application
#'
#' Call Shiny to run \code{idem} as a web-based application. A web browser will
#' be brought up.
#'
#' @examples
#' \dontrun{
#' run.idem()}
#'
#'
#' @export
#'
run.idem <- function() {

    if (!requireNamespace("shiny", quietly = TRUE)) {
        stop("Shiny needed for this function to work. Please install it.",
             call. = FALSE)
    }

    appDir <- system.file("shiny", package = "idem")
    if (appDir == "") {
        stop("Could not find Shiny directory. Try re-installing `idem`.",
             call. = FALSE)
    }

    shiny::runApp(appDir, display.mode = "normal");
}
