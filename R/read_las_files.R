#' readFileIntoChar
#' @description reads LAS into ascii string
#' @param file_name File path to LAS
#'
#' @return
#' @export
#'
#' @examples
readFileIntoChar <- function(file_name){
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
  asciitxt <- unlist(lapply(asciitxt, trimws))
  # Get rid of the commented lines
  asciitxt <- asciitxt[which(substring(asciitxt, 1, 1) != "#")]
  return(asciitxt)
}

#' check_pull_SubString
#' @description check for a pattern in the LAS-as-ascii
#' @param asciitxt las ascii object
#' @param sub.str pattern to search for
#'
#' @return
#' @export
#'
#' @examples
check_pull_SubString <- function(asciitxt, sub.str) {
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

#' getNumbersTable
#' @description sets / cuts the headers, and unnecessary characters from las ascii, returns numeric valued table
#' @param asciitxt las ascii object
#' @param headers headers to be set
#'
#' @return
#' @export
#'
#' @examples
getNumbersTable <- function(asciitxt, headers) {
  topdata <-
    which(substring(asciitxt, 1, 6) == "~Ascii" |
      substring(asciitxt, 1, 2) == "~A") + 1
  dataarray <- asciitxt[topdata:length(asciitxt)]
  dataarray <- stringr::str_replace_all(dataarray,
    pattern = " ",
    replacement = ";"
  )
  dataarray <- stringr::str_replace_all(dataarray,
    pattern = ";+",
    ";"
  )
  dataarray <- stringr::str_replace_all(dataarray,
    pattern = "^;",
    ""
  )
  dataarray <- unlist(strsplit(dataarray, split = ";"))
  if (length(dataarray) %% length(headers) != 0) {
    cat(
      paste0(
        "Invalid number of data-points. Well ",
        uwi,
        " may not be parsed correctly!\r\n"
      )
    )
    cat(paste0(
      "Length. ",
      length(headers), " "
    ))
    print(headers)
    #  return(NA)
  }
  dataarray <- as.numeric(dataarray)
  dataarray[dataarray == -999.25] <- NA
  return(dataarray)
}

#' getCurveHeaders
#' @description Retrieve the columns provided in the LAS header
#' @param asciitxt las ascii object
#'
#' @return
#' @export
#'
#' @examples
getCurveHeaders <- function(asciitxt) {
  inds <- which(substring(asciitxt, 1, 1) == "~", arr.ind = TRUE)
  topheader <- which(substring(asciitxt, 1, 6) == "~Curve") + 1
  btmheader <- inds[inds > topheader][1] - 1


  headers <- asciitxt[topheader:btmheader]
  headers <- strsplit(headers, split = ":")
  headers <- suppressWarnings(as.character(do.call(
    rbind,
    headers
  )[, 1]))
  headers <- stringr::str_replace_all(headers,
    pattern = " ",
    ""
  )
  headers <- tolower(headers)
  headers <- gsub("[0-9]+", "", headers)
  return(headers)
}


#' getLASList
#' @description get a list of las files from the folder path provided
#' @param folder_path folderto look for las files
#'
#' @return
#' @export
#'
#' @examples
getLASList <- function(folder_path) {
  if (file.info(folder_path)$isdir) {
    lasfiles <-
      list.files(
        path = folder_path,
        full.names = TRUE,
        recursive = TRUE,
        pattern = "las"
      )
    lasfilesshort <-
      list.files(
        path = folder_path,
        full.names = FALSE,
        recursive = TRUE,
        pattern = "las"
      )
    nlasfiles <- length(lasfiles)
  }
  else {
    print("The folder does not exist, or the path is wrong !")
  }

  return(list(path = lasfiles, short = lasfilesshort, n = nlasfiles))
}

#' readLASFileIntoTable
#' @description Main function, reads lasfile provided, and adds a name ID as a 'short' variable (you can )
#' @param lasfile path to an las file
#' @param short id of that file (user-defined)
#'
#' @return
#' @export
#'
#' @examples
readLASFileIntoTable <- function(lasfile, short) {
  asciitxt <- readFileIntoChar(lasfile)

  if (length(asciitxt) == 1) {
    print(sprintf("Could not read  %s....", lasfile))
    return(NA)
  }
  # Get UWI
  uwi <- check_pull_SubString(asciitxt = asciitxt, sub.str = "UWI")
  # Get Well ID
  well_id <- check_pull_SubString(asciitxt = asciitxt, sub.str = "WELL")

  # This section looks for the curve section that has headers for the data
  # it will ignore the tops section, if it is present
  headers <- getCurveHeaders(asciitxt)
  dataarray <- getNumbersTable(asciitxt, headers = headers)

  result <- data.table(matrix(dataarray,
    ncol = length(headers),
    byrow = TRUE
  ))
  setnames(result, headers)
  result$uwi <- uwi[[1]]
  result$well_id <- well_id[[1]]
  setkeyv(result, c("uwi"))
  result <- subset(result, select = c("uwi", headers))
  result$file <- lasfile
  result$short <- short
  return(result)
}

fillWellLogs <- function(las_df) {
  cols <- names(las_df)

  names_to_find <-
    c(
      "rhob.k/m",
      "rhob.g/cc",
      "dtc.us/m",
      "dts.us/m",
      "dt.m/s",
      "sdt.m/s",
      "phie.v/v",
      "efac_shale."
    )
  names_right <- c("rho", "rho", "vp", "vs", "vp", "vs", "phie", "efac")
  unit_c <- c(1, 1e3, 1e6, 1e6, 1, 1, 1, 1)
  invert <- c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)

  for (i in seq_along(names_to_find))
  {
    print(paste0("i=", i))
    ind <- grep(sprintf("^%s", names_to_find[i]), cols)
    print(paste0("ind=", ind))
    # Index found:
    if (length(ind) > 0) {
      print(sprintf("Renaming %s to %s", cols[ind], names_right[i]))
      # Adjust Unit, invert if necessary
      if (invert[i]) {
        las_df[[ind]] <- 1 / las_df[[ind]]
      }
      # print(ind)
      las_df[[ind]] <- unit_c[i] * las_df[[ind]]
      cols[ind] <- names_right[i]
    }
  }
  names(las_df) <- cols
  print(names(las_df))
  las_df <-
    mutate(
      las_df,
      VpVs = vp / vs,
      Ip = rho * vp,
      Is = rho * vs,
      lambda = rho * (vp**2 - 2 * vs**2),
      mu = rho * vs**2
    )
  las_df <-
    mutate(las_df,
      lambda_rho = 1e-12 * lambda * rho,
      mu_rho = 1e-12 * mu * rho
    ) %>%
    filter(rho > 0 & vp > 0 & vs > 0 & rho < 8e3 & vp < 13e3 & vs < 13e3 & lambda_rho > 0, mu_rho > 0) %>%
    filter((is.na(phie)) | (phie > 0))
  return(las_df)
}

findClosestAlias <- function(x, dict.vars, type) {
  ind <- which.min(adist(x, dict.vars$CurveName, ignore.case = T))

  distance <- min(adist(x, dict.vars$CurveName, ignore.case = T))
  sub <- dict.vars$Alias[ind]
  coeff.sub <- dict.vars$Units[ind]
  inv.var <- dict.vars$Invert[ind]
  y <- switch(type,
    sub = sub,
    coeff.sub = coeff.sub,
    inv.var = inv.var,
    dist = distance
  )
  return(y)
}

GetLASColumnNames <- function(col.names, dict.vars) {
  col_names <- toupper(col.names)
  matches.sub <- sapply(col.names, function(x) findClosestAlias(x, dict.vars, type = "sub"))
  matches.inv <- sapply(col.names, function(x) findClosestAlias(x, dict.vars, type = "inv.var"))
  matches.coeffsub <- sapply(col.names, function(x) findClosestAlias(x, dict.vars, type = "coeff.sub"))
  matches.dist <- sapply(col.names, function(x) findClosestAlias(x, dict.vars, type = "dist"))

  matches.df <- tbl_df(data.frame(
    sub = matches.sub, inv = matches.inv, coeff.sub = matches.coeffsub, orig = col_names,
    dist = matches.dist
  ))
  # IGNORE columns that match with a large distance
  matches.df[matches.df$dist >= 3, ]$sub <- "IGNORE"
  cat("\nI replaced the columns according to the following matches:\n\n")
  print(matches.df)
  cat("\n\nBe careful, the match could be wrong... In that case add the entry into the dictionary with a proper Curve title and unit !")

  return(matches.df)
}

calculateElasticAttrs <- function(las_df) {
  las_df <- mutate(las_df,
    VpVs = vp / vs, Ip = rho * vp, Is = rho * vs, lambda = rho * (vp**2 - 2 * vs**2),
    mu = rho * vs**2
  )
  las_df <- mutate(las_df, lambda_rho = 1e-12 * lambda * rho, mu_rho = 1e-12 * mu * rho) %>%
    filter(
      rho > 0 & vp > 0 & vs > 0 & rho < 8e3 & vp < 13e3 & vs < 13e3 & lambda_rho > 0,
      mu_rho > 0
    )
  return(las_df)
}

fillWellLogs.dict <- function(LAS, dict.filename = "./Dict-LAS-clean.csv") {
  dict.LAS <- tbl_df(read.csv(dict.filename, stringsAsFactors = FALSE)) %>% na.omit() %>% mutate(CurveName = toupper(gsub(" ", "", CurveName)))
  matches <- GetLASColumnNames(col.names = colnames(LAS), dict.vars = dict.LAS)
  new.col.names <- matches$sub
  print(new.col.names)
  colnames(LAS) <- as.character(new.col.names)
  LAS <- tbl_df(LAS %>%
    setNames(make.names(names(.), unique = TRUE)) %>%
    select(-starts_with("IGNORE")))

  names(LAS) <- gsub(".", "", names(LAS), fixed = TRUE)

  for (name in colnames(LAS))
  {
    print(name)
    submatch <- matches %>% filter(sub == name)
    if (nrow(submatch) == 0) next
    if (submatch$inv) {
      LAS[, name] <- 1 / LAS[, name]
    }
    if (submatch$coeff.sub != 1) {
      LAS[, name] <- LAS[, name] * submatch$coeff.sub
    }
  }
  LAS <- calculateElasticAttrs(LAS)
  return(LAS)
}
