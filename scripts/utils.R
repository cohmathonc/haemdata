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
#' @source https://nelson-gon.github.io/12/06/2020/hex-sticker-creation-r/
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
# get a file from the PSON Teams folder and copy it to
# the data-raw folder
get_teams_file <- function(filename) {
    drive <- Microsoft365R::get_team("PSON AML State-Transition")$get_drive()
    drive$download_file(filename)
    file.copy(
        from = basename(filename),
        to = "data-raw",
        overwrite = TRUE
    )
    file.remove(basename(filename))
}

# build and install the package -------------------------------------------
build_package <- function(latest_published_data) {

    # Check the package
    devtools::check(error_on = "error", vignettes = FALSE)

    # Make a release on GitHub
    ## put fields in standard order and alphabetises dependencies
    usethis::use_tidy_description()
    # use_tidy_eval()
    # usethis::use_version("dev")
    # usethis::use_dev_version("minor")

    # build the package
    tgz <- devtools::build(vignettes = FALSE)
    # Publish to cgt.coh.org
    drat::insertPackage(tgz, "/net/nfs-irwrsrchnas01/labs/rrockne/MHO")

    # install locally
    # devtools::install()
    install.packages("haemdata", repo = "http://cgt.coh.org/MHO")

    # build the pkgdown site
    pkgdown::clean_site()
    devtools::document()

    pkgdown::build_site()
}
