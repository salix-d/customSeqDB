#'=================================================================================
#' Misc functions
#'=================================================================================
#'
#'--------------------------------------------------------------------------------
#' Write a fasta file from sequences and their names
#'--------------------------------------------------------------------------------
#'
#' @param sequences vector or named list of character string. The sequences to be writen in the file.
#' @param ids       vector of character string. The ids for the sequences; Optional if a named list is provided
#' @param outDir    character string. The directory in which to write the file
#' @param outFile   character string. The name of the file to write
#' @param nbChar    integer. The number of sequence characters per line. Default is NULL, sequence will be one line
#' @param open      character string. A description of how to open the connection. Default is "w" for write, could be set to a to append the sequences.
#' @return returns a string with the path to the fasta file.
#' @export
write.fasta <- function(sequences,
                        ids = NULL,
                        outDir = ".",
                        outFile = paste0("sequences", Sys.Date(), ".fasta"),
                        nbChar = NULL,
                        open = "w"){
  if(!dir.exists(outDir)) dir.create(outDir)
  outFile <- file.path(outDir, outFile)
  fasta <- file(description = outFile, open = open)
  # if the sequences should be split on multiple lines,
  # makes a list of split sequences of nbChar length (or less for the last line)
  if(!missing(nbChar)){
    sequences <- lapply(sequences, function(txt){
      ss <- strsplit(txt, "")[[1]]
      ns <- 0:floor(length(ss)/nbChar) * nbChar + 1
      return(sapply(ns, function(i) paste(na.omit(ss[i:(i+59)]), collapse = "")))
    })
  }
  if(missing(ids)) ids <- names(sequences)
  # write pairs of ids and sequences
  sapply(seq_along(sequences), function(i){
    writeLines(paste0(">", ids[i]), fasta)
    writeLines(as.character(sequences[[i]]), fasta)
  })
  close(fasta)
  return(outFile)
}


#' Writes parsed data to a csv file
#'
#' @description
#' `write_parsed2csv` writes the parsed data to a csv file and returns a character string of its path.
#'
#' @details
#' Writes the parsed data to csv file and returns a character string of its path.
#' If no path provided for outFile, file will be saved to "./db_downloads/pasredRecords/pasredRecords.[Sys.Date()].csv"
#' Is used by  \code{\link{parse_EMBLxml}}, \code{\link{parse_flatFile}}, \code{\link{parse_INSDxml}} when save2csv is set to TRUE.
#'
#' @param parsedInfo  data.frame. The data.frame returned from \code{\link{parse_EMBLxml}}, \code{\link{parse_flatFile}} or \code{\link{parse_INSDxml}}.
#' @param outFile     character string. The path to the csv file to be written. Optional. If save2csv is set to TRUE and outFile is missing, the file will be saved to './db_downloads/parsedRecords/parsedRecords_[Sys.Date()]_[###].csv'
#' @export
write_parsed2csv <- function(parsedInfo, outFile){
  if(!outFile){
    mainDir <- "./db_downloads"
    subDir <- "pasredRecords"
    if(!dir.exists(mainDir)) dir.create(mainDir)
    if(!dir.exists(file.path(mainDir, subDir))) dir.create(file.path(mainDir, subDir))

    pattern <- paste0("pasredRecords.", Sys.Date(), ".")
    outFile <- file.path(mainDir, subDir, paste0(pattern, "001.csv"))
    if(file.exists(outFile)){
      n <- length(dir(subDir, pattern = pattern))
      outFile <- file.path(mainDir, subDir, paste0(pattern, sprintf("%03d", n+1), ".csv"))
    }
  } else {
    if(!grepl(".csv$", outFile)) outFile <- paste0(outFile, ".csv")
  }
  write.csv(pasredInfo, outFile)
  return(normalizePath(outFile))
}

#' Download the fetched records instead of reading them directly
#'
#' @description
#' `saveRecords` downloads a file from a URL and returns a character string of its path.
#'
#' @details
#' Downloads a file from a URL and returns a character string of its path.
#' If no path provided for outFile, file will be saved to "./db_downloads/records/records.[Sys.Date()].csv"
#' Is used by  \code{\link{parse_EMBLxml}}, \code{\link{parse_flatFile}}, \code{\link{parse_INSDxml}} when save2csv is set to TRUE.
#'
#' @param parsedInfo  data.frame. The data.frame returned from \code{\link{parse_EMBLxml}}, \code{\link{parse_flatFile}} or \code{\link{parse_INSDxml}}.
#' @param outFile     character string. The path to the csv file to be written. Optional. If missing, the file will be saved to './db_downloads/records/records_[Sys.Date()]_[###].csv'
#' @export
saveRecords <- function(URL, outFile, ext = c("xml", "txt")){
  if(missing(outFile)){
    mainDir <- "./db_downloads"
    subDir <- "records"
    if(!dir.exists(mainDir)) dir.create(mainDir)
    if(!dir.exists(file.path(mainDir, subDir))) dir.create(file.path(mainDir, subDir))

    pattern = paste0("records.", Sys.Date(), ".")
    outFile <- file.path(mainDir, subDir, paste0(pattern, "001.", ext))
    if(file.exists(outFile)){
      n <- length(dir(subDir, pattern = pattern))
      outFile <- file.path(mainDir, subDir, paste0(pattern, sprintf("%03d", n+1), ".", ext))
    }
  } else {
    if(!grepl(paste0(".",ext,"$"), outFile)) outFile <- paste0(outFile, ".", ext)
  }
  download.file(URL, destfile = outFile)
  return(outFile)
}


#' Makes the matrix with long row names (like RNA sequences) printable.
#'
#' @description
#' `mk_printDF` returns a data.frame without named rows. Optionally, also replaces NAs for a value.
#'
#' @details
#' Convert matrix to data.frame and remove row names. Can replace the NA by "NA" or custom value make comparison between tables easier.
#'
#' @param mat    matrix to be converted. Usually the one returned from dada2::assignedTaxonomy()
#' @param na2char Weither to convert NAs to a character string or not; Default TRUE
#' @param naString  String to use to replace the NA if na2char is set to TRUE; Default is "unidentified"
#' @export
mk_printDF <- function(mat, na2char = TRUE, naString = "unidentified"){
  df <- as.data.frame(mat)
  rownames(df) <- NULL
  if(na2char) df[is.na(df)] <- naString
  return(df)
}
