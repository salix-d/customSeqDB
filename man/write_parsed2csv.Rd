% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/savingFiles.R
\name{write_parsed2csv}
\alias{write_parsed2csv}
\title{Writes parsed data to a csv file}
\usage{
write_parsed2csv(parsedInfo, outFile)
}
\arguments{
\item{parsedInfo}{data.frame. The data.frame returned from \code{\link{parse_EMBLxml}}, \code{\link{parse_flatFile}} or \code{\link{parse_INSDxml}}.}

\item{outFile}{character string. The path to the csv file to be written. Optional. If save2csv is set to TRUE and outFile is missing, the file will be saved to './db_downloads/parsedRecords/parsedRecords_[Sys.Date()]_[###].csv'}
}
\description{
`write_parsed2csv` writes the parsed data to a csv file and returns a character string of its path.
}
\details{
Writes the parsed data to csv file and returns a character string of its path.
If no path provided for outFile, file will be saved to "./db_downloads/pasredRecords/pasredRecords.[Sys.Date()].csv"
Is used by  \code{\link{parse_EMBLxml}}, \code{\link{parse_flatFile}}, \code{\link{parse_INSDxml}} when save2csv is set to TRUE.
}
