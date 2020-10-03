#' Get the taxIds of taxon names from the ENA database
#'
#' @description
#' `get_taxIds` returns a data.frame with taxId and taxon scientific name as columns.
#'
#' @details
#' Uses the ENA API to search the taxonomy database with taxon names and parse the data to return only the taxid and the scientific name for each taxon.
#'
#' @param names    vector of character string. Names of the taxa for which to find their taxId.
#' @param nameType character string. The type of name used. Default is 'scientific' Can be set to 'commun' if the commun names should be used.
#' @family ena-db
#' @export
#' @seealso \code{\link{https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/taxon-based-search.html}}
#' @references
get_enaTaxIds <- function(names, nameType = "scientific"){
  types <- list(scientific = 'scientific-name', commun = 'any-name')
  baseURL <- "https://www.ebi.ac.uk/ena/taxonomy/rest/"
  cols <- c("taxId", "scientificName")
  out <- lapply(names, function(name) jsonlite::fromJSON(readLines(paste0(baseURL, types[nameType], "/", name)))[,cols])
  return(do.call(rbind,out))
}
#' Get accession numbers and descriptions of records from the ENA databases
#'
#' @description
#' `get_enaDetails` returns a data.frame with accession numbers and descriptions (as columns) of all records of a specified database matching specified taxa tree and description.
#'
#' @details
#' Uses the ENA Portal API to look up the records from each taxIds taxonomic tree that have the result type and the description keywords specified and return their accession number and description.
#'
#'
#' @param taxIds      character string, integer or vector of character string or integer. The taxonomic ids used to include their sub-tree in the search.
#' @param description character string or vector of character string or integer. The keywords used to search in the "description" field of the records
#' @param result      character string. The result type (data set) to search against. Default is "sequence" (Nucleotide sequences). See \url{https://www.ebi.ac.uk/ena/portal/api/results?dataPortal=ena} for all possible results.
#' @param outFile     character string. The path to where the file should be downloaded. (With this, the URL can't always be read directly, so it's downloaded.)
#' @param keepFile    logical. Whether to keep or remove the downloaded file once it's been read.
#' @param mode        character. Whether to write  ("w") over or append ("a") to an existing outFile. Default is "w".
#' @param overwrite   logical. If the file exist and the mode isn't "a" (append), should the file be overwritten. If set to FALSE, will read the existing file instead of downloading a new one. Default is TRUE.
#' @family ena-db
#' @export
#' @references \url{https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/taxon-based-search.html}, \url{https://www.ebi.ac.uk/ena/portal/api/}
get_enaDetails <- function(taxIds,
                           description = NULL,
                           result = "sequence",
                           outFile="./db_downloads/ena_description.json",
                           keepFile=T,
                           mode = "w",
                           overwrite=T){

  if(file.exists(outFile) & mode!="a" & !overwrite){
    message(outFile, "exist. Details will be read from this file. \nIf that wasn't the goal, set the mode argument to 'a', the overwrite argument to TRUE or provide a different path for the outFile argument.")
  } else {
    # making the URL
    baseURL <- "https://www.ebi.ac.uk/ena/portal/api/search?"
    taxIds <- if(length(taxIds)>1) paste0("(", paste(paste0("tax_tree(", taxIds, ")"), collapse = " OR "), ")") else paste0("tax_tree(", taxIds, ")")
    if(!missing(description)) description <- if(length(description)>1) paste0(" AND (", paste(paste0('description="*', description, '*"'), collapse = " OR "), ")") else paste0(' AND description="*', description, '*"')
    params <- list(
      result = result,
      query = gsub(" ", "%20", paste0(taxIds, description)),
      format = "json"
    )
    URL <- paste0(baseURL, paste(paste0(names(params), "=", params), collapse = "&"))
    #making sure the directory exist, if not creates it
    if(!dir.exists(dirname(outFile))) dir.create(dirname(outFile))
    # download the file
    utils::download.file(URL, destfile = outFile, mode = mode)
    message("    File saved : ", outFile)
  }
  # read the file and convert to data.frame
  out <- jsonlite::read_json(outFile, simplifyVector = T)
  # remove file if needed
  if(!keepFile){
    unlink(outFile)
    message("    ", outFile, " has been removed")
  }
  return(out)
}
#' Fetch records and parse them for the required information.
#'
#' @description
#' `fetch_enaSeq` returns a data.frame with the accession number, the taxid and the region of the sequence for the specified gene (as columns) for the provided accession numbers (as rows).
#'
#' @details
#' Uses DBFetch to fetch the records for all provided accession numbers in insdxml format and parse them to return the accession number, the taxid and the region of the sequence for the specified gene.
#'
#' @param accList     character string, integer or vector of character string or integer. The accession numbers used to fetch the records.
#' @param gene        character string or vector of character string. The gene for which to find the sequence region.
#' @param codon_start logical. Must be set to T if the gene has reading frames so it can adjust the region start accordingly. Default is false.
#' @param full_seq    logical. Whether to fetch the full sequence or only the gene region. Default is FALSE, gets the gene region.
#' @param format      character string. The format of the fetched file. Supported choice for parsing : 'flatfile', 'emblxml-1.1' or 'insdseq'.
#' @family ena-db
#' @export
#' @references \url{http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/dbfetch.databases#ena_sequence}, \url{http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/dbfetch}
get_enaFasta <- function(accList, db="ena_sequence", format = c('embl', 'emblxml-1.1', 'insdseq'), style = "raw", saveRec = F, outRec = NULL, saveParsedRec = F, outParsedRec = NULL){
  baseURL <- "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?"
  if(missing(format)) format <- "embl"
  params <- list(
    db     = "ena_sequence",
    id     = paste(accList, collapse = ","),
    format = "indsxml",
    style  = "raw"
  )
  URL <- paste0(baseURL, paste(names(params), "=", params), collapse = "&")
  if(saveRec) URL <- save_records(URL, outFile = outRec, ext = iselse(format=="embl", "txt", "xml"))

  if(format == "emblxml-1.1") format <- "emblxml"
  parseFuns <- list(
    embl = parse_flatFile,
    indsxml = parse_INSDxml,
    emblxml = parse_EMBLxml
  )
  out <- parseFuns[[format]](URL = URL, gene = gene, codon_start = codon_start, full_seq = full_seq, save2csv = saveParsedRec, outFile = outParsedRec)
  return(out)
}
