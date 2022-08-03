
# Functions for setting the pinboard and getting data from it
#' @include haemdata.R

# Setup the pinboard for use in the current R session
# Get data from either the OneDrive or development pin board
# Assign the pinboard to the haemdata_env environment

#' @title Setup the pinboard for use
#' @description Sets up the pinboard to use for the current R session. Rather
#' than issuing this command each session, you can add the line
#' `if (require(haemdata)) haemdata::use_pinboard("onedrive")` to your project folder
#' [`.Rprofile`](https://support.rstudio.com/hc/en-us/articles/360047157094-Managing-R-with-Rprofile-Renviron-Rprofile-site-Renviron-site-rsession-conf-and-repos-conf) file.
#' @param pin_board Indicate which pinboard to use. One of c('onedrive,
#' 'devel'), Default: NULL
#' @param haemdata_folder The full path of the folder for devel
#' pins, Default: '/net/isi-dcnl/ifs/user_data/rrockne/MHO/haemdata'
#' @return Assigns the `pin_board` variable in the haemdata_env environment
#' @details Package releases always use the OneDrive pinboard with versioned pins.
#' Pin versions associated with each release can be retrieved with `data(published_pins)`.
#' Development pins are published to the `haemdata_folder`, but are not versioned.
#' @author Denis O'Meally
#' @examples
#' \dontrun{
#' if (interactive()) {
#'     # Use the OneDrive pinboard
#'     use_pinboard("onedrive")
#'     # Use the Development pinboard
#'     use_pinboard("devel")
#' }
#' }
#' @seealso
#'  \code{\link[pins]{board_ms365}}, \code{\link[pins]{board_folder}}
#'  \code{\link[Microsoft365R]{get_personal_onedrive}}
#' @rdname use_pinboard
#' @export
#' @importFrom pins board_ms365 board_folder
#' @importFrom Microsoft365R get_team
use_pinboard <- function(pin_board = NULL,
                         haemdata_folder = "/net/isi-dcnl/ifs/user_data/rrockne/MHO/haemdata") {
    # setup pin_board
    if (is.null(pin_board)) {
        message("`pin_board` is NULL: set the pinboard to `onedrive` or `devel`")
        assign("pin_board", NULL, haemdata_env)
    } else if (pin_board == "onedrive") {
        # OneDrive pin_board
        assign("pin_board", pins::board_ms365(
            drive = Microsoft365R::get_team("PSON AML State-Transition")$get_drive(),
            path = "haemdata",
            versioned = TRUE
        ), envir = haemdata_env)
    } else if (pin_board == "devel") {
        # check that haemdata_folder is accessible
        if (!dir.exists(haemdata_folder)) {
            stop(glue::glue("'{haemdata_folder}' is not is not accessible.
            Check the path"))
        }
        # devel pin_board
        assign("pin_board",
            pins::board_folder(
                haemdata_folder,
                versioned = FALSE
            ),
            envir = haemdata_env
        )
    } else {
        stop("Please set the `pin_board` parameter to either `onedrive` or `devel`")
    }
}

# get the latest version of a pin
#' @title Get latest pin version
#' @description Finds a pin with a given name and returns the latest
#' version from the pin_board. The pinboard must be set up for the current
#' R session using the [`use_pinboard()`] function.
#' @param pin_name the name of the pin for which the version is required.
#' @return the pin version
#' @details Use this function to get the latest version of a pin.
#' Uses the [`pins::pin_versions()`] function to get the list of pin
#' versions and `dplyr` to get the top element of the returned list.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'     metadata_mmu_version <- get_latest_pin_version("metadata_mmu.csv")
#' }
#' }
#' @seealso
#'  \code{\link[pins]{pin_versions}}
#'  \code{\link[dplyr]{pull}}, \code{\link[dplyr]{nth}}
#' @rdname get_latest_pin_version
#' @export
#' @importFrom pins pin_versions
#' @importFrom dplyr pull first
get_latest_pin_version <- function(pin_name) {
    if (is.null(pin_name)) {
        stop("`pin_name` is NULL")
    } else if (is.null(haemdata::haemdata_env$pin_board)) {
        stop(haemdata::haemdata_env$pin_board_msg)
    } else {
        haemdata::haemdata_env$pin_board |>
            pins::pin_versions(pin_name) |>
            dplyr::pull(version) |>
            dplyr::first()
    }
}

#' @title Get pin
#' @description Get an object from the pin_board
#' @param pin_name The name of the pin to retrieve
#' @param version The version of the pin to retrieve; if NULL, the latest version is used. Default: NULL
#' @return An object from the pin_board
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'     metadata_mmu <- get_pin("metadata_mmu.csv")
#' }
#' }
#' @seealso
#'  \code{\link[pins]{pin_read}}
#' @rdname get_pin
#' @export
#' @importFrom pins pin_read
get_pin <- function(pin_name, version = NULL) {
    if (is.null(haemdata::haemdata_env$pin_board)) {
        stop(haemdata::haemdata_env$pin_board_msg)
    } else {
        available_pins <- get_pin_list()
        if (pin_name %in% available_pins) {
            if (is.null(version)) {
                version <- get_latest_pin_version(pin_name)
            }
            hash <- gsub(".*-", "", version)

            pins::pin_read(haemdata::haemdata_env$pin_board, pin_name, version = version, hash = hash)
        } else {
            stop("Pin not found in pin_board; use get_pin_list() to see available pins")
        }
    }
}

# get names of all heamdata pins
#' @title Get pin list
#' @description Return a list of pin names from the pin_board.
#' If `onedrive` is TRUE, the pin is retrieved from Sharepoint.
#' If `onedrive` is FALSE, the pin is retrieved from the MHO
#' haemdata folder.
#' @return list of pin names
#' @examples
#' \dontrun{
#' if (interactive()) {
#'     get_pin_list()
#' }
#' }
#' @seealso
#'  \code{\link[pins]{pin_list}}
#' @rdname get_pin_list
#' @export
#' @importFrom pins pin_list
get_pin_list <- function() {
    if (is.null(haemdata::haemdata_env$pin_board)) {
        stop(haemdata::haemdata_env$pin_board_msg)
    } else {
        haemdata::haemdata_env$pin_board |>
            pins::pin_list()
    }
}