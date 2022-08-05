# Utility functions for the target pipeline

# cross platform hostname
get_hostname <- function() {
    return(as.character(Sys.info()["nodename"]))
}

# package logo -----------------------------------------------------------
#' Make the package logo
#'
#' makes a hex sticker given a url
#' @param url a URL to a picture
#' @param text the name on the sticker
#' @return a hex sticker
#' @export

# Hex sticker
# https://nelson-gon.github.io/12/06/2020/hex-sticker-creation-r/
make_logo <- function(url = NULL, text = NULL) {
    if (is.null(url)) {
        img <- magick::image_read("https://www.maxpixel.net/static/photo/2x/Blood-Group-0-Blood-Rh-factor-Positive-2781421.jpg")
    } else {
        img <- magick::image_read(url)
    }

    if (is.null(text)) {
        text <- "haemdata"
    }

    # make the logo 600dpi
    img |>
        magick::image_convert("png") |>
        magick::image_resize("900 x 900") |>
        hexSticker::sticker(
            package = text,
            p_size = 25, p_y = 1,
            s_x = 1, s_y = 0.8,
            s_width = 10, s_height = 10,
            white_around_sticker = TRUE,
            h_color = "#FE2B3F",
            filename = "data-raw/haemdata_icon.png",
            dpi = 600
        )

    # use the logo
    usethis::use_logo("data-raw/haemdata_icon.png", retina = TRUE)

    return(img)
}


#' Write rda file to `data/`
#'
#' This function writes an `rda` file in `data`
#' and names it `schemeName`. Makes saving objects easier
#' from within a function. Inspired by [this](https://stackoverflow.com/questions/56293910/create-r-data-with-a-dynamic-variable-name-from-function-for-package)
#' StackOverflow post. Uses 'xz' to compress the file.
#'
#' @param schemeName the name of the rda file, without the `.rda` extension
#' @param data the object to save
#' @examples
#' \dontrun{
#' if (interactive()) {
#'     write_rda("an_rda_file_named_that", an_object_named_this)
#' }
#' }
#' @return the path to the rda file
#' @export
#' @rdname write_data
write_data <- function(schemeName, data) {
    rdaFile <- paste0(schemeName, ".rda")
    fileLocation <- file.path(".", "data", rdaFile)
    varName <- paste0(schemeName)

    assign(varName, data)
    eval(parse(text = sprintf("save(%s, file = '%s', compress = 'xz')", varName, fileLocation)))
    return(paste0("data/", rdaFile))
}