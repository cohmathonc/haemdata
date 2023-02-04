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
            filename = here::here("data-raw/haemdata_icon.png"),
            dpi = 600
        )

    # use the logo
    usethis::use_logo(here::here("data-raw/haemdata_icon.png"), retina = TRUE)

    return(img)
}
# get a file from the PSON Teams folder and copy it to
# the data-raw folder
get_teams_file <- function(filename) {
    drive <- Microsoft365R::get_team("PSON AML State-Transition")$get_drive()
    drive$download_file(filename)
    file.copy(
        from = basename(filename),
        to = here::here("data-raw"),
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

    # cp over target markdown
    #fs::file_copy(here::here("preprocessing/_targets.html"), "~/MHO/haemdata-www", overwrite = TRUE)

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

#' Add age and weeks columns to a data.frame
#'
#' @param sample_sheet A data frame containing mouse ID, date of birth (dob), and sample date
#' @return A data frame with additional columns: sample_weeks, age_at_end, age_at_start, age_at_sample
#' @export
#'
#' @examples
#' add_age_and_weeks_columns(sample_sheet)
#'
add_age_and_weeks_columns <- function(sample_sheet) {
    sample_sheet |>
        dplyr::group_by(mouse_id) |>
        dplyr::mutate(
            # add columns: sample_weeks, age_at_end, age_at_start, age_at_sample
            sample_weeks = difftime(as.Date(sample_date), min(as.Date(sample_date)), units = "weeks") |> round(1),
            age_at_end = difftime(max(as.Date(sample_date)), as.Date(dob), units = "weeks"),
            age_at_start = difftime(min(as.Date(sample_date)), as.Date(dob), units = "weeks"),
            age_at_sample = difftime(as.Date(sample_date), as.Date(dob), units = "weeks"),
            dplyr::across(dplyr::starts_with("age"), round, 1),
        ) |>
        dplyr::ungroup() |>
        purrr::modify_if(lubridate::is.timepoint, as.numeric)
}

#' Assign New Sample IDs
#'
#' ChatGPT explains: This code is a function that assigns new sample_ids to a given sample sheet.
#' It takes in the sample sheet as an argument and filters out any rows that don't have
#' a sample_id. It then creates a numeric version of both the library_id and sample_id
#' columns, arranges them in order, and selects only the mouse_id, tissue, sample_date,
#' sample_id_num, and library_id_num columns. It then creates a new column for the sample_id
#' by adding "MHO" to the beginning of each number in the sample_id column and another column
#' for the library id by adding "COHP" to each number in the library id column. Finally, it
#' filters out any rows with library ids that are not in the new samples list and patches them
#' into the original sample sheet.
#'
#' @param sample_sheet A data frame containing the columns `mouse_id`, `tissue`, `sample_date`, `sample_id` and `library_id`.
#' @return A data frame with the same columns as \code{sample_sheet}, but with new sample IDs assigned to rows containing missing values in the `sample_id` column.
#' @export
#'
#' @examples
#' assign_new_sample_ids(sample_sheet)
#'
assign_new_sample_ids <- function(sample_sheet) {
        new_samples <- sample_sheet |>
            dplyr::filter(is.na(sample_id)) |>
            pull(library_id)

        new_sample_ids <- sample_sheet |>
            dplyr::mutate(
                library_id_num = as.numeric(gsub("[^[:digit:]\\.]", "", library_id)),
                sample_id_num = as.numeric(gsub("[^[:digit:]\\.]", "", sample_id))
            ) |>
            dplyr::arrange(sample_id_num, library_id_num) |>
            dplyr::select(mouse_id, tissue, sample_date, sample_id_num, library_id_num) |>
            dplyr::distinct() |>
            dplyr::ungroup() |>
            dplyr::mutate(
                sample_id_num = ifelse(is.na(sample_id_num),
                    max(sample_id_num, na.rm = TRUE) + seq(1, length(which(is.na(sample_id_num))), by = 1),
                    sample_id_num
                ),
                sample_id = sprintf("MHO_%04d", sample_id_num),
                library_id = paste0("COHP_", library_id_num)
            ) |>
            dplyr::select(sample_id, library_id) |>
            dplyr::filter(library_id %in% new_samples)

        dplyr::rows_patch(sample_sheet, new_sample_ids, by = "library_id")
}
