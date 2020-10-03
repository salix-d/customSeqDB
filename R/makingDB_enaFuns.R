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
#'
#'
#' @description
#'
#'
#' @details
#'
#'
#' @param
#' @family ena-db
#' @export
#' @seealso \code{\link{}}
#' @references
get_enaDetails <- function(taxIds, description = NULL, outFile="./ena_description.json", keepFile=T, mode = "w", overwrite=T){

  if(file.exists(outFile) & mode!="a" & !overwrite){
    message(outFile, "exist. Details will be read from this file. \nIf that wasn't the goal, set the mode argument to 'a', the overwrite argument to TRUE or provide a different path for the outFile argument.")
  } else {
    baseURL <- "https://www.ebi.ac.uk/ena/portal/api/search?"
    taxIds <- if(length(taxIds)>1) paste0("(", paste(paste0("tax_tree(", taxIds, ")"), collapse = " OR "), ")") else paste0("tax_tree(", taxIds, ")")
    if(!missing(description)) description <- if(length(description)>1) paste0(" AND (", paste(paste0('description="*', description, '*"'), collapse = " OR "), ")") else paste0(' AND description="*', description, '*"')
    params <- list(
      result = "sequence",
      query = gsub(" ", "%20", paste0(taxIds, description)),
      format = "json"
    )
    URL <- paste0(baseURL, paste(paste0(names(params), "=", params), collapse = "&"))
    utils::download.file(URL, outFile)
  }
  out <- jsonlite::read_json(outFile, simplifyVector = T)
  if(!keepFile){
    unlink(outFile)
    message("    ", outFile, " has been removed")
  }
  return(out)
}
#'
#'
#' @description
#'
#'
#' @details
#'
#'
#' @param
#' @family ena-db
#' @export
#' @seealso \code{\link{}}
#' @references
fetch_enaSeq <- function(accList, gene, codon_start = F){
  URL <- paste0("http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ena_sequence;id=", paste(accList, collapse = ","),";format=emblxml-1.2")

  ena_xml <- xml2::read_xml(URL)
  INSDSeqs <- xml2::xml_contents(xml2::xml_children(ena_xml))

  # function to keep only wanted information from the xml entries
  parse_dbFetch_entries <- function(INSDSeq, gene){
    AC <- xml2::xml_contents(xml2::xml_contents(INSDSeq)[grep("accession-version", xml2::xml_contents(INSDSeq))])
    SQ <- gsub("\\n|\\t", "", xml2::xml_contents(xml2::xml_find_all(INSDSeq, ".//sequence")))
    FT <- xml2::xml_find_all(INSDSeq,"./feature")

    TX <- xml2::xml_find_all(FT[xml2::xml_attr(FT, "name") == "source"], "./taxon")
    taxId <- as.numeric(xml2::xml_attr(TX, "taxId"))
    sciName <- xml2::xml_attr(TX, "scientificName")
    CDS <- FT[xml2::xml_attr(FT, "name") == "CDS"]
    which.COI <- grep(gene, CDS, ignore.case = T)
    if(length(which.COI)==0){
      which.COI <- grep(gene, ge, ignore.case = T)
        warning("    Warning : ", AC, " (", sciName, ") ", ": ", DE, " has no COI gene and will be ignored.\n")
        return()
    } else if(length(which.COI)>1){
      warning("    Warning : ", AC, " has more than one COI gene, first one will be use. Record should be checked manually.\n")
      which.COI <- which.COI[1]
    }
    COI <- CDS[which.COI]
    COI_loc <- xml2::xml_attr(COI, "location")
    if(is.na(COI_loc)) COI_loc <- grep("complement", COI, value = T)
    COI_location <- as.numeric(strsplit(gsub("complement| |\\(|\\)|<|>", "", COI_loc), "\\.\\.")[[1]])
    COI_location <- list(start = COI_location[1], end = COI_location[2])
    if(codon_start){
      # fix the start location depending on the reading frame
      COI_start <- COI_location$start + as.numeric(gsub("\\n|\\t", "", xml2::xml_contents(xml2::xml_contents(xml2::xml_contents(COI)[which(xml2::xml_attr(xml2::xml_contents(COI), "name")=="codon_start")])))) - 1
      COI_SQ <- substring(SQ, COI_start, COI_location$end)
    } else {
      COI_SQ <- substring(SQ, COI_location$start, COI_location$end)
    }
    return(list(accession = AC, description = DE, sequence = toupper(COI_SQ), taxId = taxId, sciName = sciName))
  }
  parsed_entries <- lapply(ena_entries, function(INSDSeq) parse_dbFetch_entries(INSDSeq, gene))
  return(do.call(rbind,parsed_entries))
}
