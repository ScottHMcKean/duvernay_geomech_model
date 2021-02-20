#' Read a file into an ASCII text file
#' @param file_name File path to LAS
#' @return char string
#' @export
read_file_to_char <- function(file_name){
  fi <- file.info(file_name)
  if (!is.na(fi$size) & (!fi$isdir)) {
    conn <- file(file_name, open = "r")
    asciitxt <- readLines(conn)
  }
  else {
    print("Bad file : either does not exist or is a folder !")
    return(NA)
  }

  close(conn)
  rm(conn)
  # trim whitespace
  asciitxt <- unlist(lapply(asciitxt, trimws))
  # eliminate commented lines
  asciitxt <- asciitxt[which(substring(asciitxt, 1, 1) != "#")]
  return(asciitxt)
}

#' Check and pull a substring from LAS
#' @param asciitxt las ascii object
#' @param sub.str pattern to search for
#' @return character string
pull_substring <- function(asciitxt, sub.str) {
  print(sprintf("Looking for a substring  %s....", sub.str))

  item <-
    asciitxt[which(substring(asciitxt, 1, nchar(sub.str)) == sub.str)]

  if (nchar(item) > 1) {
    item <- stringr::str_replace_all(item, pattern = "_", "")
    item <- stringr::str_replace_all(item, pattern = "-", "")
    item <- stringr::str_replace_all(item, pattern = " ", "")
    item <- strsplit(item, split = ".", 2)[[1]][[2]]
    item <- strsplit(item, split = ":", 2)[[1]][[1]]
  }
  else {
    print(sprintf("The substring  %s was not found....", sub.str))
    item <- NA
  }
}

#' Make a numeric table by cutting the headers and unnecessary las characters
#' @param asciitxt las ascii object
#' @param headers headers to be set
#' @export
make_curves_table <- function(asciitxt, headers) {
  topdata <-which(
    substring(asciitxt, 1, 6) == "~Ascii" |
    substring(asciitxt, 1, 2) == "~A") + 1
  dataarray <- asciitxt[topdata:length(asciitxt)] %>%
    str_replace_all(., " ", ";") %>%
    str_replace_all(., ";+", ";") %>%
    str_replace_all(., "^;", "")

  dataarray <- unlist(strsplit(dataarray, split = ";"))

  if (length(dataarray) %% length(headers) != 0) {
    cat(paste0("Invalid number of data-points. Well ",uwi,
        " may not be parsed correctly!\r\n")
    )
  }
  dataarray <- as.numeric(dataarray)
  dataarray[dataarray == -999.25] <- NA
  return(dataarray)
}

#' Get Las Headers and return cleaned column
#' @description Retrieve the columns provided in the LAS header
#' @param asciitxt las ascii object
#' @export
get_las_headers <- function(asciitxt) {
  inds <- which(substring(asciitxt, 1, 1) == "~", arr.ind = TRUE)
  topheader <- which(substring(asciitxt, 1, 6) == "~Curve") + 1
  btmheader <- inds[inds > topheader][1] - 1
  headers <- asciitxt[topheader:btmheader]
  headers <- strsplit(headers, split = ":")
  headers <- suppressWarnings(
    as.character(do.call(rbind,headers)[, 1]))
  headers <- stringr::str_replace_all(headers,pattern = " ","")
  headers <- tolower(headers)
  headers <- gsub("[0-9]+", "", headers)
  return(headers)
}

#' Load LAS file into a table with a name ID and formatted headers
#' @param lasfile path to an las file
#' @param short id of that file (user-defined)
#' @export
load_las_file <- function(lasfile, short_uwi) {
  print(paste0("Reading ", lasfile))
  asciitxt <- read_file_to_char(lasfile)

  if (length(asciitxt) == 1) {
    print(sprintf("Could not read %s", lasfile))
    return(NA)
  }
  # Get UWI
  uwi <- pull_substring(asciitxt = asciitxt, sub.str = "UWI")

  # This section looks for the curve section that has headers for the data
  # it will ignore the tops section, if it is present
  headers <- get_las_headers(asciitxt)
  dataarray <- make_curves_table(asciitxt, headers = headers)

  result <- data.table(matrix(dataarray,
    ncol = length(headers),
    byrow = TRUE
  ))
  setnames(result, headers)
  result$uwi <- uwi[[1]]
  setkeyv(result, c("uwi"))
  result <- subset(result, select = c("uwi", headers))
  result$short <- short_uwi
  return(result)
}

#' Mappable function to read LAS files and subset to duvernay base and top
#' @param filepath full filepath
#' @param study_locs an sf dataframe with study locations and kelly bushing
#' @param dv_base a raster of Duvernay base
#' @param dv_top a raster of Duvernay top
#' @export
get_subset_log_data <- function(filepath, study_locs, dv_base, dv_top){
  uwi <- str_extract(filepath, "(?<=//).+(?=.las)") %>%
    str_to_upper()

  las_df = load_las_file(filepath, short_uwi=uwi)

  well_loc = study_locs[study_locs$uwi == uwi,]

  base = raster::extract(dv_base, well_loc %>% as(.,'Spatial'))
  top = raster::extract(dv_top, well_loc %>% as(.,'Spatial'))

  las_df %>%
    mutate(elev_masl = well_loc$kbe - depth.m) %>%
    filter(elev_masl > base & elev_masl < top) %>%
    janitor::clean_names()
}

#' Mappable function to subset logs into shale and carbonate based on tops
#' works on elev_masl (what the tops are picked in)
#' @param this_uwi the uwi for subsetting and facies classification
#' @param comb_las_df the entire combined las df
#' @param dv_carb_tops the tops of the Duvernay (has start and end of carbonate)
#' @return a dataframe of a single well, to be combined using purrr::map_df
#' @export
classify_log_intervals = function(this_uwi, comb_las_df, dv_carb_tops){
  carb_top = dv_carb_tops %>% filter(uwi == this_uwi) %>% pull(dv_b_start) *-1
  if (is.na(carb_top)){
    carb_top = comb_las_df %>%
      filter(uwi == this_uwi) %>%
      pull(elev_masl) %>%
      min() - 0.1
  }
  carb_bot = dv_carb_tops %>% filter(uwi == this_uwi) %>% pull(dv_b_end) *-1
  if (is.na(carb_bot)){
    carb_bot = comb_las_df %>%
      filter(uwi == this_uwi) %>%
      pull(elev_masl) %>%
      min() - 0.1
  }

  out_las_df = comb_las_df %>%
    filter(uwi == this_uwi) %>%
    mutate(facies = case_when(
      elev_masl > carb_top ~ 'dv_shale',
      elev_masl <= carb_top & elev_masl >= carb_bot ~ 'dv_carb',
      TRUE ~ 'na'
    ))
}
