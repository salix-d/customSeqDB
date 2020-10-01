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
