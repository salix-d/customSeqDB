#' Parse records of INSDxml format for the required information.
#'
#' @description
#' `parse_INSDxml` parses a INSDxml file and returns a data.frame with the details of the records (as columns) for the provided accession numbers (as rows).
#'
#' @details
#' Reads a INSDxml file from a file path or a URL and parse it for the accession number, taxid, description, gene, product, protein id, translation exception, note and sequence fields.
#' Is used by \code{\link{get_enaFasta}}
#'
#' @param accList     character string, integer or vector of character string or integer. The accession numbers used to fetch the records.
#' @param gene        character string or vector of character string. The gene for which to find the sequence region.
#' @param codon_start logical. Must be set to T if the gene has reading frames so it can adjust the region start accordingly. Default is false.
#' @param full_seq    logical. Whether to fetch the full sequence or only the gene region. Default is FALSE, gets the gene region.
#' @param save2csv    logical. Whether to save returned data.frame of the parsed records to a csv file or not. Done using \code\link{write_parsed2csv}.
#' @param outFile     character string. The path to the csv file to be written. Optional. If save2csv is set to TRUE and outFile is missing, the file will be saved to './db_downloads/parsedRecords/parsedRecords_[Sys.Date()]_[###].csv'
#' @export
#' @seealso \code{\link{parse_EMBLxml}}, \code{\link{parse_flatFile}}, \code{\link{get_enaFasta}}
parse_INSDxml <- function(URL, gene, codon_start = F, full_seq = F, save2csv = F, outCsv = NULL){

  INSDSeqs <- xml2::xml_children(xml2::read_xml(URL))
  # function to keep only wanted information from each record
  parse_INSDSeq <- function(INSDSeq, gene, codon_start, full_seq){
    AC <- xml2::xml_contents(xml2::xml_find_all(INSDSeq, "./INSDSeq_accession-version"))
    print(AC)
    DE <- xml2::xml_contents(xml2::xml_find_all(INSDSeq, "./INSDSeq_definition"))
    SQ <- xml2::xml_contents(xml2::xml_find_all(INSDSeq, "./INSDSeq_sequence"))
    FT <- xml2::xml_find_all(INSDSeq,"./INSDSeq_feature-table/INSDFeature")

    # get the feature 'source' which contain the taxon id and find all the children values, including the taxon id; displays only the content
    SRC <- xml2::xml_contents(xml2::xml_find_all(FT[grep("source", FT)], ".//INSDQualifier_value" ))
    # grep the line with the taxon and keep only the id
    TX <- gsub("taxon:([0-9]*)$", "\\1", grep("taxon", SRC, value = T))

    GENE_info <- data.frame(
      accession     = paste(AC),
      taxId         = paste(TX),
      description   = paste(DE)
    )

    if(full_seq){
      GENE_info$sequence <- toupper(SQ)
      return(GENE_info)
    }

    whichIsGene <- grep(gene, FT, ignore.case = T, perl = T)
    if(length(whichIsGene)==0){
      warning("    Warning : ", AC, " has no COI gene and will be ignored.\n")
      return()
    }
    GENE <- xml2::xml_children(xml2::xml_children(FT[whichIsGene]))

    fields <- c("gene","product","proteinId","transl_except","note")
    for(f in fields){
      i <- grep(f, GENE)
      out[,f] <- ifelse(length(i)==0, NA, gsub("^.*<INSDQualifier_value>([^<]*)<.*$", "\\1", GENE[i[1]]))
    }

    GENE_location <- lapply(c("from>", "to>"), function(x) as.numeric(gsub("^.*from>([0-9]*)<.*$", "\\1", GENE[grep(x, GENE)][1]))) #sometimes is in field gene and field cds. both are the same, use the first.
    # if needed, fix the start location depending on the reading frame
    if(codon_start){
      GENE_frame <- as.numeric(gsub("^.*value>([0-9])<.*$", "\\1", GENE[grep("codon_start", GENE)]))
      if(length(GENE_frame)==0) GENE_frame = 1
      if(GENE_frame>1) GENE_start <- GENE_start - 1 + GENE_frame
    }
    # get only the region of the gene from the sequence
    GENE_info$sequence <- toupper(substring(SQ, GENE_start, GENE_end))
    return(GENE_info)
  }
  parsed_INSDSeqs <- do.call(rbindlapply(INSDSeqs, function(INSDSeq) parse_INSDSeq(INSDSeq = INSDSeq, gene = gene, codon_start = codon_start, full_seq = full_seq)))
  if(save2csv) write_parsed2csv(parsed_INSDSeqs, outCsv)
  return(parsed_INSDSeqs)
}

#' Parse records of EMBLxml format for the required information.
#'
#' @description
#' `parse_EMBLxml` parses a EMBLxml file and returns a data.frame with the details of the records (as columns) for the provided accession numbers (as rows).
#'
#' @details
#' Reads a EMBLxml file from a file path or a URL and parse it for the accession number, taxid, description, gene, product, protein id, translation exception, note and sequence fields.
#' Is used by \code{\link{get_enaFasta}}
#'
#' @param accList     character string, integer or vector of character string or integer. The accession numbers used to fetch the records.
#' @param gene        character string or vector of character string. The gene for which to find the sequence region.
#' @param codon_start logical. Must be set to T if the gene has reading frames so it can adjust the region start accordingly. Default is false.
#' @param full_seq    logical. Whether to fetch the full sequence or only the gene region. Default is FALSE, gets the gene region.
#' @param save2csv    logical. Whether to save returned data.frame of the parsed records to a csv file or not. Done using \code\link{write_parsed2csv}.
#' @param outFile     character string. The path to the csv file to be written. Optional. If save2csv is set to TRUE and outFile is missing, the file will be saved to './db_downloads/parsedRecords/parsedRecords_[Sys.Date()]_[###].csv'
#' @export
#' @seealso \code{\link{parse_INSDxml}}, \code{\link{parse_flatFile}}, \code{\link{get_enaFasta}}
parse_EMBLxml <- function(URL, gene, codon_start = F, full_seq = F, save2csv = F, outCsv = NULL){

  entries <- xml2::xml_children(xml2::read_xml(URL))

  # function to keep only wanted information from the embl entries
  parse_entry <- function(entry, gene, full_seq){
    AC <- paste(xml2::xml_attr(entry, "accession"), xml2::xml_attr(entry, "version"), sep = ".")
    SQ <- toupper(gsub("\\n|\\t", "", xml2::xml_contents(xml2::xml_find_all(entry, ".//sequence"))))
    FT <- xml2::xml_find_all(entry,"./feature")

    TX <- xml2::xml_find_all(FT[xml2::xml_attr(FT, "name") == "source"], "./taxon")
    TX <- as.numeric(xml2::xml_attr(TX, "taxId"))
    GENE_info <- data.frame(
      accession     = paste(AC),
      taxId         = paste(TX),
      description   = paste(DE)
    )

    if(full_seq){
      GENE_info$sequence <- toupper(SQ)
      return(GENE_info)
    }

    whichIsGene <- grep(gene, FT, ignore.case = T)
    if(length(whichIsGene)==0){
      warning("    Warning : ", AC, " has no COI gene and will be ignored.\n")
      return()
    }
    GENE <- xml2::xml_children(FT[whichIsGene])
    fields <- c("gene","product","proteinId","transl_except","note")
    for(f in fields){
      i <- grep(f, GENE)
      out[,f] <- ifelse(length(i)==0, NA, gsub("^.*<value>\\s*(\\S*)\\s.*$", "\\1", GENE[i[1]]))
    }

    GENE_loc <- xml2::xml_attr(FT[whichIsGene][1], "location")
    if(any(is.na(GENE_loc))) GENE_loc <- grep("complement", FT[whichIsGene], value = T)
    GENE_location <- setNames(as.list(as.numeric(strsplit(gsub("complement| |\\(|\\)|<|>", "", GENE_loc), "\\.\\.")[[1]])), c("start","end"))
    if(codon_start){
      # fix the start location depending on the reading frame
      GENE_frame <- as.numeric(gsub("^.*<value>\\s*(\\S*)\\s.*$", "\\1", GENE[grep("codon_start", GENE)]))
      if(length(GENE_frame)==0) GENE_frame <-  1
      if(GENE_frame>1) GENE_location$start <- GENE_location$start + GENE_frame  - 1
    }
    GENE_info$sequence <- substring(SQ, GENE_location$start, GENE_location$end)
    return(GENE_info)
  }
  parsed_entries <- do.call(rbind,lapply(entries, function(entry) parse_entry(entry, gene, full_seq = full_seq)))
  if(save2csv) write_parsed2csv(parsed_entries, outCsv)
  return(parsed_entries)
}

#' Parse records of EMBLxml format for the required information.
#'
#' @description
#' `parse_flatFile` parses a EMBLxml file and returns a data.frame with the details of the records (as columns) for the provided accession numbers (as rows).
#'
#' @details
#' Reads a flat file from a file path or a URL and parse it for the accession number, taxid, description, gene, product, protein id, translation exception, note and sequence fields.
#' Is used by \code{\link{get_enaFasta}}
#'
#' @param accList     character string, integer or vector of character string or integer. The accession numbers used to fetch the records.
#' @param gene        character string or vector of character string. The gene for which to find the sequence region.
#' @param codon_start logical. Must be set to T if the gene has reading frames so it can adjust the region start accordingly. Default is false.
#' @param full_seq    logical. Whether to fetch the full sequence or only the gene region. Default is FALSE, gets the gene region.
#' @param save2csv    logical. Whether to save returned data.frame of the parsed records to a csv file or not. Done using \code\link{write_parsed2csv}.
#' @param outFile     character string. The path to the csv file to be written. Optional. If save2csv is set to TRUE and outFile is missing, the file will be saved to './db_downloads/parsedRecords/parsedRecords_[Sys.Date()]_[###].csv'
#' @export
#' @seealso \code{\link{parse_INSDxml}}, \code{\link{parse_EMBLxml}}, \code{\link{get_enaFasta}}
parse_flatFile <- function(URL, gene, codon_start = F, full_seq = F, save2csv = F, outCsv = NULL){

  records <- readLines(URL)
  rec_end <- grep("//", records)
  rec_start <- c(0, rec_end[-length(rec_end)])+1
  records <- lapply(seq_along(rec_end), function(i) records[rec_start[i]:rec_end[i]])

  # function to keep only wanted information from the xml entries
  parse_record <- function(record, gene, full_seq){
    AC <- gsub("VERSION|ACCESSION|\\s", "", grep("^VERSION\\s|^ACCESSION\\s",record, value=T))
    DE <- gsub("^DEFINITION\\s*|^DE\\s*", "", grep("^DEFINITION\\s*|^DE\\s*", record, value = T))

    SQ_start <- grep("^SQ\\s|^ORIGIN\\s", record)
    SQ <- record[(SQ_start+1):(length(record)-1)]
    SQ <- toupper(gsub("\\s|\\d", "", paste(SQ, collapse = "")))

    FT_start <- grep("^\\s{5}[[:alpha:]]", record)
    FT_end <- c(FT_start[-1], SQ_start)-1
    FT <- sapply(seq_along(FT_start), function(i) record[FT_start[i]:FT_end[i]])

    TX <- FT[[grep("source", FT)]]
    TX <- gsub("[^0-9]*", "", grep("taxon", TX, value = T))

    GENE_info <- data.frame(
      accession   = AC,
      taxId       = TX,
      description = DE
    )

    if(full_seq){
      GENE_info$sequence <- toupper(SQ)
      return(GENE_info)
    }

    whichIsGene <- sapply(FT, function(ft) any(grepl(gene, ft, ignore.case = T)))
    if(length(whichIsGene)==0){
      warning("    Warning : ", AC, " has no COI gene and will be ignored.\n")
      return()
    }
    fields <- c("gene","product","proteinId","transl_except","note")
    GENE <- unlist(FT[whichIsGene])
    for(f in fields){
      i <- grep(f, GENE)
      out[,f] <- ifelse(length(i)==0, NA, gsub(paste0('\\s*',f,'="(.*)"$'), "\\1", GENE[i[1]]))
    }
    GENE_loc <- grep("[0-9]\\.\\.>?[0-9]", GENE, value = T)[1]
    GENE_location <- setNames(as.list(as.numeric(gsub("\\D","", strsplit(GENE_loc, "\\.\\.")[[1]]))), c("start","end"))
    if(codon_start){
      # fix the start location depending on the reading frame
      GENE_frame <- as.numeric(gsub("\\D", "", grep("codon_start", GENE, value = T)))
      if(length(GENE_frame)==0) GENE_frame <-  1
      if(GENE_frame>1) GENE_location$start <- GENE_location$start + GENE_frame  - 1
    }
    GENE_info$sequence <- substring(SQ, GENE_location$start, GENE_location$end)
    return(GENE_info)
  }
  parsedRecords <- do.call(rbind, lapply(records, function(record) parse_record(record = record, gene = gene, full_seq = full_seq)))
  if(save2csv) write_parsed2csv(parsedRecords, outCsv)
  return(parsedRecords)
}
