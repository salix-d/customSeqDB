#' Get the id list of records matching a query with ESearch
#'
#' @description
#' `get_idList` returns a list of the accession numbers and the history server information (WebEnv and Query_key)
#'
#' @details
#' ESearch a NCBI databases for records and returns the accession number.
#' Uses the ncbi history server to fetch details faster with \code{\link{get_idListDetails}}.
#' Can also use \code{\link{get_idListDetails}} directly.
#' NOTE : The ESearch can return up to 100k results. (Will update function to deal with bigger result sets)
#'
#' @param taxa       list of taxa to search for
#' @param gene       list of the genes/variations of a gene names
#' @param customTerm string of terms to add to the eSearch query
#' @param db         ncbi db to search
#' @param retmode    return mode of the results; Default is 'json' since this function parse it with jsonlite (could be adapted to parse xml)
#' @param idtype     type of id; default is accession (UID are depricated for most case)
#' @param retmax     maximum result returns; Default is 10000 which is the maximum
#' @export
#' @seealso \code{\link{get_idListDetails}}
#' @references \url{https://www.ncbi.nlm.nih.gov/books/NBK25499/}, \url{https://www.ncbi.nlm.nih.gov/books/NBK49540/#_chapter4_ESearch_}
get_idList <- function(taxa,
                       gene = "",
                       customTerm = "",
                       db="nucleotide",
                       retmode = "json",
                       idtype = "acc",
                       retmax = "100000"
){

  #' arguments check --------------------------------
  if(missing(taxa)) stop("Error : taxa argument missing")
  if(missing(gene) & missing(customTerm)) warning("Only taxa was used for the search")

  #' making the eSearch URL : -----------------------
  baseURL <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"

  gene <- if(length(gene)>1) paste(paste0(gene, "[GENE]"), collapse = " OR ") else paste0(gene, "[GENE]")
  taxa <- if(length(taxa)>1) paste(paste0(taxa, "[ORGN]"), collapse = " OR ") else paste0(taxa, "[ORGN]")
  term <- paste(paste0("(", taxa, ") AND (", gene, ")"), customTerm)
  params <- list(
    db = db,
    term = gsub(" ", "%20", term),
    retmode = retmode,
    idtype = idtype,
    retmax = retmax,
    usehistory = "y"
  )
  URL <- paste0(baseURL, paste(paste0(names(params), "=", params), collapse = "&"))
  cat(paste0("\tsearching for records matching :", term, "\n"))

  #' parsing the data for the idList (accessions) ---
  dat <- jsonlite::fromJSON(URL)
  out <- list(accession  = dat$esearchresult$idlist,
              WebEnv = dat$esearchresult$webenv,
              Query_key = dat$esearchresult$querykey)

  #' printing number of results ---------------------
  cat(paste0("\t", length(out$accession), " records found\n"))

  return(out)
}


#' Get the details (description and taxId) from a list of ids with ESummary
#'
#' @description
#' `get_details` returns a list containing a list of the history server information ($query) and a list with the fetch details ($details).
#'
#' @details
#' Uses ESummary to get the summary information of NCBI records and parse the data for the desired information (description and taxId).
#' Uses the ncbi history server info returned from \code{\link{get_idList}} to fetch details faster.
#' If used like \code{\link{get_idList}}, will use that function to retrieved the history server info first.
#'
#' @param taxa       list of taxa to search for; Optional if using the arguments WebEnv, Query_key and nIDs.
#' @param gene       list of the genes/variations of a gene names; Optional if using the arguments WebEnv, Query_key and nIDs.
#' @param customTerm string of terms to add to the eSearch query; Optional if using the arguments WebEnv, Query_key and nIDs.
#' @param db         ncbi db to search. Default is 'nucleotide'.
#' @param retmode    return mode of the results; Default is 'json' since this function parse it with jsonlite (could be adapted to parse xml)
#' @param idtype     type of id; default is accession (UID are depricated for most case)
#' @param retmax     maximum result returns; Default is 10000 which is the maximum
#' @param query      A list of the history server info and the ids associated with it($accession, $Query_key, $WebEnv); Object returned from \code{\link{get_idList}}
#' @export
#' @seealso \code{\link{get_idList}}
#' @references \url{https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_ESummary_}
get_details <-function(taxa,
                       gene = "",
                       customTerm = "",
                       db="nucleotide",
                       retmode = "json",
                       idtype = "acc",
                       retmax = "100000",
                       query = NULL){

  #' arguments check --------------------------------
  if(missing(query)){
    if(missing(taxa)) stop("Error : taxa argument missing")
    if(missing(gene) & missing(customTerm)) warning("Only taxa was used for the search")

    #' getting the records and history server info ----
    query <- get_idList(taxa = taxa,
                        gene = gene,
                        customTerm = customTerm,
                        db = db,
                        retmode = retmode,
                        idtype = idtype,
                        retmax = retmax)
  } else {
    #' make sure query is valid ---------------------
    if(length(query)<3) stop("Error : missing values in query")
    if(! names(query) %in% c("WebEnv", "Query_key", "accession")) stop("Invalid query list, must contain $WebEnv, $Query_key, $accession")
  }
  #' making the eSummary URL ------------------------
  baseURL <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
  params <- list(
    db = db,
    retmax = "10000",
    WebEnv = query$WebEnv,
    Query_key = query$Query_key
  )
  URL <- paste0(baseURL, paste(paste0(names(params), "=", params), collapse = "&"))
  print(URL)

  #' Function to parse the data ---------------------
  parse_url <- function(URL){
    lines <- readLines(URL)
    values <- grep("AccessionVersion|TaxId|Title", lines, value = T)
    out <- as.data.frame(matrix(gsub("^.*>(.*)<.*$", "\\1", values), ncol=3, byrow = T, dimnames = list(NULL, c("description", "taxid", "accession")))[,c(3,2,1)])
    return(out)
  }

  nIDs <- length(query$accession)
  #' split if more than 5000
  #' server can return up to 10000, but doing so I was missing some records.
  #' Could be testd with more than 5000.
  if(nIDs>5000){
    n <- floor(nIDs/5000)
    ranges <- list(start = c(0:n * 5000), end = c(1:n * 5000 - 1, nIDs))
    #out <- do.call(rbind, lapply(URLs, parse_url))
    smmry <- list()
    for(i in seq_along(ranges$start)){
      cat(paste0("\tparsing records ", ranges$start[i], " to ", ranges$end[i]," of ", nIDs, "\n"))
      smmry[[i]] <- parse_url(paste0(URL, "&retstart=", ranges$start[i]))
    }
    out <- do.call(rbind, smmry)
  } else {
    out <- parse_url(URL)
  }
  return(list(query = query, details = out))
}

#' Get the fasta files from a list of ids with EFetch
#'
#' @description
#' `get_seqFasta` returns a vector of the paths to the downloaded fasta files.
#'
#' @details
#' Download EFetched fasta files and returns a vector of their paths.
#' Could use the ncbi history server info if all the sequences were wanted. (To be implemented)
#' NOTE : This can only get full sequences or CDS region, function for other gene region will be added.
#'
#' @param taxa       list of taxa to search for; Optional if using the arguments WebEnv, Query_key and nIDs.
#' @param gene       list of the genes/variations of a gene names; Optional if using the arguments WebEnv, Query_key and nIDs.
#' @param customTerm string of terms to add to the eSearch query; Optional if using the arguments WebEnv, Query_key and nIDs.
#' @param db         ncbi db to search. Default is 'nucleotide'.
#' @param rettype    return type of the results; Default is 'fasta_cds_na' which gets only the CDS. Could be set to fasta to get whole sequences. Note: to get only the gene region, the whole gb file should be fetched and parsed similarly to how it's done for the ENA db;
#' @param Query_key  integer specifying which of the id lists attached to the given Web Environment will be used as input to EFetch (to be implemented).
#' @param WebEnv     string indicating the Web Environment that contains the id list to be provided as input to EFetch (to be implemented).
#' @export
#' @seealso \code{\link{get_idListDetails}}
#' @references \url{https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_EFetch_}
get_seqFasta <- function(idList,
                         outDir=".",
                         db="nucleotide",
                         rettype = "fasta_cds_na",
                         WebEnv = NULL,
                         Query_key = NULL){
  #' makes the ranges of sequences to download (200 at a time)
  n <- floor(length(idList)/200)
  r = list(start = 0:n * 200 + 1, end = c(1:n * 200, len))
  ranges <- lapply(seq_len(n+1), function(i) c(r$start[i]:r$end[i]))
  #' makes file names
  fileNames <- paste(paste0("eFetch_", Sys.Date()), seq_len(n+1), sep = "_")
  if(!dir.exists(outDir)) dir.create(outDir)
  outFiles<-file.path(outDir, paste0(fileNames, ".fasta"))
  #' make a list of URL to EFetch the fasta files
  URLs <- sapply(ranges, function(r) paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=", db,"&id=", paste(idList[r], collapse = ","), "&rettype=", rettype, "&retmode=text"))
  #' downloads the EFetched fasta files
  mapply(utils::download.file, URLs, outFiles)
  return(outFiles)
}
