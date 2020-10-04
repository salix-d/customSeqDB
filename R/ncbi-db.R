#' Get the id list of records matching a query with ESearch
#'
#' @description
#' `get_ncbiAcc` returns a list of the accession numbers and the history server information (WebEnv and Query_key)
#'
#' @details
#' ESearch a NCBI databases for records and returns the accession number.
#' Can get the details information with \code{\link{get_ncbiDetails}} with the returned list or use \code{\link{get_ncbiDetails}} directly.
#'
#' @note The ESearch can return up to 100k results. (Will update function to deal with bigger result sets)
#' @param taxa       list of taxa to search for
#' @param gene       list of the genes/variations of a gene names
#' @param customTerm string of terms to add to the eSearch query
#' @param db         ncbi db to search
#' @param retmode    return mode of the results; Default is 'json' since this function parse it with jsonlite (could be adapted to parse xml)
#' @param retmax     maximum result returns; Default is 10000 which is the maximum
#' @family ncbi-db
#' @export
#' @seealso \code{\link{get_ncbiDetails}}
#' @references \url{https://www.ncbi.nlm.nih.gov/books/NBK25499/}, \url{https://www.ncbi.nlm.nih.gov/books/NBK49540/#_chapter4_ESearch_}
get_ncbiAcc <- function(taxa,
                       gene = "",
                       customTerm = "",
                       db="nucleotide",
                       retmode = "json",
                       retmax = "100000",
                       WebEnv = ""){

  #' arguments check ------------------------------------------------------------------------------------------------
  if(missing(taxa)) stop("Error : taxa argument missing")
  if(missing(gene) & missing(customTerm)) warning("Only taxa was used for the search")

  #' making the eSearch URL : ---------------------------------------------------------------------------------------
  baseURL <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"

  gene <- if(length(gene)>1) paste(paste0(gene, "[GENE]"), collapse = " OR ") else paste0(gene, "[GENE]")
  taxa <- if(length(taxa)>1) paste(paste0(taxa, "[ORGN]"), collapse = " OR ") else paste0(taxa, "[ORGN]")
  term <- paste(paste0("(", taxa, ") AND (", gene, ")"), customTerm)
  params <- list(
    db = db,
    term = gsub(" ", "%20", term),
    retmode = retmode,
    idtype = "acc",
    retmax = retmax,
    usehistory = "y"
  )
  if(!missing(WebEnv)) params$WebEnv <- WebEnv
  URL <- paste0(baseURL, paste(paste0(names(params), "=", params), collapse = "&"))
  message("\n    searching for records matching :", term)

  #' parsing the data for the accList (accessions) -------------------------------------------------------------------
  dat <- jsonlite::fromJSON(URL)
  out <- list(accession = dat$esearchresult$idlist,
              WebEnv    = dat$esearchresult$webenv,
              Query_key = dat$esearchresult$querykey)

  #' printing number of results -------------------------------------------------------------------------------------
  cat("    ", length(out$accession), " records found\n")

  return(out)
}


#' Get the details (description and taxId) from a list of ids with ESummary
#'
#' @description
#' `get_ncbiDetails` returns a list containing a list of the history server information ($query) and a list with the fetch details ($details).
#'
#' @details
#' Uses ESummary and the history server to get the summary information of NCBI records and parse the data for the desired information (description and taxId).
#' If \code{\link{get_ncbiAcc}} was used previously, the returned list should be provided as the query argument.
#' If not, taxa, gene and customTerm should be use. The function will do the ESearch with those first and then get the details.
#'
#' @param taxa       list of taxa to search for; Optional if using the arguments WebEnv, Query_key and nIDs.
#' @param gene       list of the genes/variations of a gene names; Optional if using the arguments WebEnv, Query_key and nIDs.
#' @param customTerm string of terms to add to the eSearch query; Optional if using the arguments WebEnv, Query_key and nIDs.
#' @param db         ncbi db to search. Default is 'nucleotide'.
#' @param retmode    return mode of the results; Default is 'json' since this function parse it with jsonlite (could be adapted to parse xml)
#' @param retmax     maximum result returns; Default is 10000 which is the maximum
#' @param query      A list of the history server info and the ids associated with it($accession, $Query_key, $WebEnv); Object returned from \code{\link{get_ncbiAcc}}
#' @family ncbi-db
#' @export
#' @seealso \code{\link{get_ncbiAcc}}
#' @references \url{https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_ESummary_}
get_ncbiDetails <-function(taxa,
                       gene = "",
                       customTerm = "",
                       db = "nucleotide",
                       retmode = "json",
                       retmax = "100000",
                       query = NULL){

  #' arguments check --------------------------------
  if(missing(query)){
    if(missing(taxa)) stop("Error : taxa argument missing")
    if(missing(gene) & missing(customTerm)) warning("Only taxa was used for the search")
    #' getting the records and history server info
    query <- get_ncbiAcc(taxa = taxa,
                        gene = gene,
                        customTerm = customTerm,
                        db = db,
                        retmode = retmode,
                        retmax = retmax)
  } else {
    if(!missing(taxa)){
      #' getting the records and adding them to the same history server -----------------------------------------------
      query <- get_ncbiAcc(taxa = taxa,
                           gene = gene,
                           customTerm = customTerm,
                           db = db,
                           retmode = retmode,
                           retmax = retmax,
                           WebEnv = query$WebEnv)
    }
    #' make sure query is valid
    if(length(query)<3) stop("Error : missing values in query")
    if(! all(names(query) %in% c("WebEnv", "Query_key", "accession"))) stop("Invalid query list, must contain $WebEnv, $Query_key, $accession")
  }

  #' making the eSummary URL ------------------------------------------------------------------------------------------
  baseURL <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
  params <- list(
    db = db,
    retmax = "10000",
    WebEnv = query$WebEnv,
    Query_key = query$Query_key
  )
  URL <- paste0(baseURL, paste(paste0(names(params), "=", params), collapse = "&"))

  #' Function to parse the data ---------------------------------------------------------------------------------------
  parse_url <- function(URL){
    lines <- readLines(URL)
    values <- grep("AccessionVersion|TaxId|Title", lines, value = T)
    out <- as.data.frame(matrix(gsub("^.*>(.*)<.*$", "\\1", values), ncol=3, byrow = T, dimnames = list(NULL, c("description", "taxid", "accession")))[,c(3,2,1)])
    return(out)
  }
  nIDs <- length(query$accession)
  #' split if more than 5000 -------------------------------------------------------------------------------------------
  #' server can return up to 10000, but doing so I was missing some records.
  #' Could be tested with more than 5000.
  if(nIDs>5000){
    n <- floor(nIDs/5000)
    ranges <- list(start = c(0:n * 5000), end = c(1:n * 5000 - 1, nIDs))
    smmry <- list()
    for(i in seq_along(ranges$start)){
      message("    parsing records ", ranges$start[i], ":", ranges$end[i]," of ", nIDs, "\n")
      smmry[[i]] <- parse_url(paste0(URL, "&retstart=", ranges$start[i]))
    }
    out <- do.call(rbind, smmry)
  } else {
    message("    parsing records ", 0, ":", nIDs," of ", nIDs, "\n")
    out <- parse_url(URL)
  }
  return(list(query = query, details = out))
}

#' Get the fasta files from a list of ids with EFetch -------------------------------------------------------------------
#'
#' @description
#' `get_ncbiFasta` returns a vector of the paths to the downloaded fasta files.
#'
#' @details
#' Download EFetched fasta files of the records matching the accList argument and returns a vector of their paths.
#' History server information could also be use, but has not been implemented yet.
#'
#' @note This can only get full sequences or CDS region.
#' @param accList     list of taxa to search for; Optional if using the arguments WebEnv, Query_key and nIDs.
#' @param outDir     list of the genes/variations of a gene names; Optional if using the arguments WebEnv, Query_key and nIDs.
#' @param db         ncbi db to search. Default is 'nucleotide'.
#' @param rettype    return type of the results; Default is 'fasta_cds_na' which gets only the CDS. Could be set to fasta to get whole sequences. Note: to get only the gene region, the whole gb file should be fetched and parsed similarly to how it's done for the ENA db;
#' @param Query_key  integer specifying which of the id lists attached to the given Web Environment will be used as input to EFetch (to be implemented).
#' @param WebEnv     string indicating the Web Environment that contains the id list to be provided as input to EFetch (to be implemented).
#' @family ncbi-db
#' @export
#' @seealso \code{\link{get_ncbiDetails}}
#' @references \url{https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_EFetch_}
get_ncbiFasta <- function(accList,
                         outDir=".",
                         db="nucleotide",
                         rettype = "fasta_cds_na",
                         WebEnv = NULL,
                         Query_key = NULL){

  #' makes the ranges of sequences to download (200 at a time) -------------------------------------------------------------
  len <- length(accList)
  n <- floor(len/200)
  r = list(start = 0:n * 200 + 1, end = c(1:n * 200, len))
  ranges <- lapply(seq_len(n+1), function(i) c(r$start[i]:r$end[i]))

  #' makes file names ------------------------------------------------------------------------------------------------------
  fileNames <- paste(paste0("eFetch_", Sys.Date()), seq_len(n+1), sep = "_")
  if(!dir.exists(outDir)) dir.create(outDir)
  outFiles<-file.path(outDir, paste0(fileNames, ".fasta"))

  #' make a list of URL to EFetch the fasta files --------------------------------------------------------------------------
  URLs <- sapply(ranges, function(r) paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=", db,"&id=", paste(accList[r], collapse = ","), "&rettype=", rettype, "&retmode=text"))

  #' downloads the EFetched fasta files ------------------------------------------------------------------------------------
  mapply(utils::download.file, URLs, outFiles)
  return(outFiles)
}

#' Fetch records and parse them for the required information.
#'
#' @description
#' `fetch_ncbiSeq` returns a data.frame with the accession number, the taxid and the region of the sequence for the specified gene (as columns) for the provided accession numbers (as rows).
#'
#' @details
#' Uses EFetch to fetch the records for all provided accession numbers in specified format and parse them to return the accession number, the taxid and the region of the sequence for the specified gene(s).
#'
#' @param accList     character string, integer or vector of character string or integer. The accession numbers used to fetch the records.
#' @param gene        character string or vector of character string. The gene for which to find the sequence region.
#' @param codon_start logical. Must be set to T if the gene has reading frames so it can adjust the region start accordingly. Default is false.
#' @param full_seq    logical. Whether to fetch the full sequence or only the gene region. Default is FALSE, gets the gene region.
#' @param rettype     character string. The format of the fetched file. Supported choice for parsing : 'gb'(flat file) or 'gbc' (INSDSeq XML).
#' @family ncbi-db
#' @export
#' @references \url{https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_EFetch_}, \url{https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly}
fetch_ncbiSeq <- function(accList = NULL, query = NULL, gene, codon_start = F, full_seq = F, db="nucleotide", rettype = c('gb', 'gbc'), saveRec = F, outRec = NULL, saveParsedRec = F, outParsedRec = NULL){
  if(missing(accList) & missing(query)) stop("A list of accession numbers (accList) or a named list containing the WebEnv and Query_key to use with the history server (query) must be provided.")
  if(missing(rettype)) rettype <- "gbwithparts"
  if(rettype == "gb") rettype <- "gbwithparts"
  retmode = ifelse(rettype=="gbwithparts", "txt", "xml")
  baseURL <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
  if(!missing(accList) & missing(query)){

    # Separating the accLi in groups of 200 so the URLs can be opened.
    len <- length(accList)
    if(len > 200){
      n <- floor(len/200)
      r = list(start = 0:n * 200 + 1, end = c(1:n * 200, len))
      ranges <- lapply(seq_along(r$start), function(i) c(r$start[i]:r$end[i]))
    } else {
      ranges <- list(1:len)
    }

    # make a list of URL to EFetch the files
    URLs <- sapply(ranges, function(r){
      params <- list(
        db      = db,
        id      = paste(accList[r], collapse = ","),
        rettype = rettype,
        retmode = retmode
      )
      return(paste0(baseURL, paste(paste0(names(params), "=", params), collapse = "&")))
    })
  } else {
    if(length(query)<3) stop("Error : missing values in query")
    if(! all(names(query) %in% c("WebEnv", "Query_key", "accession"))) stop("Invalid query list, must contain $WebEnv, $Query_key, $accession")

    baseURL <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
    len <- length(query$accession)
    if(len > 10000){
      n <- floor(len/10000)
      retstart <- 0:n * 10000 + 1
    } else {
      retstart <- 0
    }

    # make a list of URL to EFetch the files
    URLs <- sapply(retstart, function(i){
      params <- list(
        db = db,
        rettype = rettype,
        retmode = retmode,
        WebEnv    = query$WebEnv,
        Query_key = query$Query_key,
        retmax = "10000",
        retstart = i
      )
      return(paste0(baseURL, paste(paste0(names(params), "=", params), collapse = "&")))
    })
  }

  if(saveRec){
    URLs <- sapply(URLs, save_records, outFile = outRec, ext = retmode)}
  parseFuns <- list(
    gbwithparts  = parse_flatFile,
    gbc = parse_INSDxml
  )

  out <- do.call(rbind, lapply(URLs, function(URL) parseFuns[[rettype]](URL = URL, gene = gene, codon_start = codon_start, full_seq = full_seq, save2csv = saveParsedRec, outCsv = outParsedRec)))
  return(out)
}
