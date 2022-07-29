
# Functions for getting data from the Heamdata pin_board
#' @include haemdata.R

# Setup the pinboard to use
# Post data to OneDrive pin board (pinboard = "onedrive") or Isilon pin board (pinboard = "isilon")
# haemdata_folder = Location for posting data files, if not OneDrive

#' @title Setup the pinboard for use
#' @description Sets up the pinboard to use for the current R session.
#' @param pin_board Indicate which pinboard to use. One of c('onedrive,
#' 'isilon'), Default: NULL
#' @param haemdata_folder The full path of the Isilon folder holding Haemdata
#' pins, Default: '/net/isi-dcnl/ifs/user_data/rrockne/MHO/haemdata'
#' @return Nothing. Assigns the `pin_board` global variable accordingly
#' @details Package releases always use the OneDrive pinboard with versioned pins.
#' Pin versions can be retrieved from the `published_pins` variable.
#' Development versions of the package use the Isilon pinboard with unversioned pins.
#' Only the most recent pin is kept.
#' @author Denis O'Meally
#' @examples
#' \dontrun{
#' if (interactive()) {
#'     # Use the OneDrive pinboard (default)
#'     use_pinboard("onedrive")
#'     # Use the Isilon pinboard
#'     use_pinboard("isilon")
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
        message("`pin_board` is NULL: set the pinboard to `onedrive` or `isilon`")
        assign("pin_board", NULL, haemdata_env)
    } else if (pin_board == "onedrive") {
        # OneDrive pin_board
        pin_board <- pins::board_ms365(
            drive = Microsoft365R::get_team("PSON AML State-Transition")$get_drive(),
            path = "haemdata",
            versioned = TRUE
        )
        assign("pin_board", pin_board, envir = haemdata_env)
    } else if (pin_board == "isilon") {
        # Isilon pin_board
        pin_board <- pins::board_folder(
            haemdata_folder,
            versioned = FALSE
        )
        assign("pin_board", pin_board, envir = haemdata_env)
    } else {
        stop("Please set the `pin_board` parameter to either `onedrive` or `isilon`")
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
    if (is.null(haemdata_env$pin_board)) {
        rlang::inform(haemdata_env$pin_board_msg, .frequency = "always")
    } else {
        haemdata_env$pin_board |>
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
    if (is.null(haemdata_env$pin_board)) {
        rlang::inform(haemdata_env$pin_board_msg, .frequency = "always")
    } else {
        available_pins <- get_pin_list()
        if (pin_name %in% available_pins) {
            if (is.null(version)) {
                version <- get_latest_pin_version(pin_name)
            }
            hash <- gsub(".*-", "", version)

            pins::pin_read(haemdata_env$pin_board, pin_name, version = version, hash = hash)
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
    if (is.null(haemdata_env$pin_board)) {
        rlang::inform(haemdata_env$pin_board_msg, .frequency = "always")
    } else {
        haemdata_env$pin_board |>
            pins::pin_list()
    }
}