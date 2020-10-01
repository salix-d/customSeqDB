#' Make a taxonomy SQL tables from NCBI nodes.dmp and names.dmp file
#'
#' @description
#' `mk_taxonomyTable` returns a character string of the path to the SQLite file in which the table was created.
#'
#' @details
#' If the files needed aren't already downloaded, it will download a taxdump.tar.gz file from NCBI servers and extract the names.dmp and nodes.dmp files from it.
#' The files will then be merged into one temporary .dmp file keeping only the required data (id, scientific name, rank and parentId).
#' This file is then used to create a SQLite table for easy taxonomy assignment with \code{\link{get_taxonomy}}.
#'
#' @param dmpDir    character string. The path to the directory where to find or download the names.dmp and nodes.dmp files. If directory doesn't exist, it will be created. If the directory is empty, it will ask if you want to download the files.
#' @param dbFile    character string. The path where the output SQLite file should be saved. Default is './ncbi-db/'. If directory doesn't exist, it will be created.
#' @param overwrite logical. If TRUE, overwrites the taxonomy table. Default is FALSE.
#' @param rmTaxDmp  logical. If FALSE, the merged file created (tax.dmp) won't be removed after the table is created. Default is TRUE.
#' @export
#' @references \url{ftp://ftp.ncbi.nih.gov/pub/taxonomy/}, \url{https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/}
#' @seealso \code{\link{get_namesNodes}}, \code{\link{get_taxonomy}}
mk_taxonomy.table <- function(dmpDir='.',
                              dbFile='./ncbi-db',
                              overwrite=FALSE,
                              rmTaxDmp=TRUE){
  if(!dir.exists(outDir)) dir.create(outDir)
  if(!all(file.exists(file.path(outDir, c('names.dmp','nodes.dmp'))))){
    download <- readline(paste0("The nodes.dmp and names.dmp files weren't found in ", outDir,". Should they be downloaded from the NCBI server? (y/N)"))
    if(download == "y"){
      dmp <- get_namesNodes(outDir = outDir)
      namesDmp <- dmp$names
      nodesDmp <- dmp$nodes
    } else {
      stop("Please provide the path to directory containing the names.dmp and nodes.dmp files.")
    }
  }
  # if the extension is wrong, add the correct one
  if(length(grep(".sqlite$", dbFile))==0) dbFile = paste0(dbFile, ".sqlite")
  # connect to db
  db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname=dbFile)
  on.exit(RSQLite::dbDisconnect(db),add=TRUE)
  # check if table already exist
  if('taxonomy' %in% RSQLite::dbListTables(db)){
    # check if table is to be updated
    if(overwrite){
      RSQLite::dbExecute(db,'DROP TABLE taxonomy')
    }else{
      stop(paste0(dbFile,' already contains table taxonomy. Set overwrite=TRUE to reload'))
    }
  }
  # make the merged tax.dmp file
  taxFile = mk_tax.dmp(nodesDmp = nodesDmp, namesDmp = namesDmp)
  # create SQL table
  RSQLite::dbExecute(db, 'CREATE TABLE taxonomy(id INTEGER, name TEXT, rank TEXT, parent INTEGER)')
  message('\twriting the taxonomy table\n')
  # create SQL table
  RSQLite::dbWriteTable(db, name = 'taxonomy', value = taxFile, append = TRUE, header = FALSE, sep = '\t')
  message('\tcreating index for the taxonomy table\n')
  RSQLite::dbExecute(db,"CREATE INDEX index_taxonomy_id ON taxonomy(id)")
  unlink(dbFile)
  unlink(taxFile)
  if(rmTaxDmp){
    file.remove(taxFile)
    message('\tMerged dump file removed')
  }
  message(paste0('\tDONE : ', basename(dbFile), ' now has a taxonomy table\n'))
  return(dbFile)
}
#' Download names and nodes files from NCBI
#'
#' @description
#' `get_namesNodes` returns a vector of file path strings of the locations of the output files
#'
#' @details
#' Download a taxdump.tar.gz file from NCBI servers and extract the names.dmp and nodes.dmp files from it.
#' These can then be used to create a SQLite database with \code{\link{mk_taxonomy.table}}.
#'
#' @param outDir the directory to put names.dmp and nodes.dmp in. Required. If directory doesn't exist, it will be created.
#' @keywords internal
get_namesNodes<-function(outDir='.'){
  outFiles<-setNames(file.path(outDir, c('names.dmp','nodes.dmp')), c("names", "nodes"))
  if(all(file.exists(outFiles))){
    message(paste(outFiles,collapse=', '),' already exist. Delete to redownload')
    return(outFiles)
  }
  url='ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz'
  base<-basename(url)
  tmp<-tempfile()
  dir.create(tmp)
  tarFile<-file.path(tmp,base)
  utils::download.file(url,tarFile,mode='wb')
  utils::untar(tarFile,fileNames,exdir=tmp,tar='internal',compressed='gzip')
  tmpFiles<-file.path(tmp,fileNames)
  if(!all(file.exists(tmpFiles)))stop("Problem finding files ",paste(tmpFiles[!file.exists(tmpFiles)],collapse=', '))
  mapply(file.copy,tmpFiles,outFiles)
  if(!all(file.exists(outFiles)))stop("Problem copying files ",paste(outFiles[!file.exists(outFiles)],collapse=', '))
  file.remove(c(tarFile,tmpFiles))
  return(outFiles)
}


#' Parse and merge names.dmp and nodes.dmp into one file (tax.dmp)
#'
#' @description
#' `mk_tax.dmp` returns a character string of the path to the merged file (tax.dmp).
#'
#' @details
#' Helper function for \code{\link{mk_taxonomy.table}}.
#' Merge names.dmp and nodes.dmp files from NCBI into one tax.dmp file keeping only the columns required (id, name, rank and parentId)
#' For names.dmp, also only keeps the rows containing the scientific names so there's only one name per id.
#'
#' @param nodesDmp character string. The path to the nodes.dmp file.
#' @param namesDmp character string. The path to the names.dmp file.
#' @keywords internal
mk_tax.dmp <- function(namesDmp, nodesDmp){
  # get taxdump path
  path <- dirname(nodesDmp)
  # path of tax.dmp file
  taxDmp <- paste0(path, 'tax.dmp')
  # function to make a tmp.dmp file formatted to be merged
  mk_tmp.dmp <- function(dmpFile, type){
    # setting fields according to type
    # keeping only the taxid and scientific name
    if(type == 'names') fields = c(1,3)
    # keeping only the rank and parent taxid
    # tested that when only keeping sciname, taxid are all in same order so no need to keep it twice
    # but it's field could be added here for double check
    else fields <- c(5,3)
    # reformatting fields for awk argument string
    fields <- paste0('$', fields, collapse = '"\t"')
    # reformatting in tmp files as precaution
    tmpFile <- gsub('.dmp', '.tmp', dmpFile)
    # awk arguments (if names, only print rows that have "scientific name" )
    argStr <- paste0('-F',"'\t' '", if(type=='names') '/scientific name/ ' else '',"{print ", fields,"}' ", dmpFile)
    # reformats files with awk
    # reformatting files with awk is faster than SQL or converting to dataframe first. Also, avoid loading big df in R for no reason.
    system2('awk', args = argStr, stdout = tmpFile)
    return(tmpFile)
  }
  # formatting names to be merged
  message('\tReading names dump file\n')
  namesTmp = mk_tmp.dmp(namesDmp, type = 'names')
  # formatting nodes to be merged
  message('\tReading nodes dump file\n')
  nodesTmp = mk_tmp.dmp(nodesDmp, type = 'nodes')
  # merging with paste using system2 (faster than SQL merge)
  argStr = paste('-d "\t"', namesTmp, nodesTmp)
  message('\tMerging names and nodes into one tax dump file\n')
  system2('paste', args = argStr, stdout = taxDmp)
  # removing tmp files
  file.remove(c(namesTmp, nodesTmp))
  return(taxDmp)
}

#' Parse and merge names.dmp and nodes.dmp into one file (tax.dmp)
#'
#' @description
#' `get_taxonomy` returns a data frame with the names of the desired ranks (columns) corresponding to each taxids (rows).
#'
#' @details
#' Uses the taxonomy table in a sqlite file to gets the the parentId of all parents of each taxIds and returns the names corresponding to the ranks desired a data frame with the taxids as row names and ranks as columns.
#'
#' @param ids       vector of character string or integer. The taxIds for which to find the taxonomic information.
#' @param sqlFile   character string. The path to the sqlite database containing the taxonomy table.
#' @param ranks     vector of character string. The desired ranks to get the taxonomic information of.
#' @export
get_taxonomy <- function (ids, sqlFile, ranks=c('superkingdom','phylum','class','order','family','genus','species')){
  if(length(ids)==0) return(NULL)
  workingIDs = unique(as.numeric(ids))
  tax.res = data.frame(matrix(as.character(NA), ncol = length(ranks), nrow = length(workingIDs), dimnames=list(workingIDs, ranks)))
  rep = 0

  get_parents<-function(ids, sqlFile){
    db  <- RSQLite::dbConnect(RSQLite::SQLite(), dbname=sqlFile)
    on.exit(RSQLite::dbDisconnect(db),add=TRUE)
    RSQLite::dbExecute(db, "CREATE TEMPORARY TABLE query(id INTEGER)")
    RSQLite::dbExecute(db, paste0("INSERT INTO query(id) VALUES(", paste0(as.integer(ids), collapse = '),('), ");"))
    taxaDf <- RSQLite::dbGetQuery(db, paste0("SELECT * FROM query LEFT JOIN taxonomy USING(id)"))
    if(!all(taxaDf$id %in% ids)) stop(simpleError('Problem finding ids'))
    return(taxaDf)
  }

  while(any(stillWorking <- !is.na(workingIDs) & workingIDs != 1)){
    parents = get_parents(workingIDs[stillWorking], sqlFile = sqlFile)
    for(rank in ranks[ranks %in% parents$rank]){
      m <- c(which(parents$rank %in% rank), grep("unclassified", parents$name))
      tax.res[which(stillWorking)[m], rank] <- parents$name[m]
      is.incertae = grep('incertae sedis', parents$name)
      if(length(is.incertae) > 0){
        for(i in is.incertae){
          rank <- colnames(tax.res)[which(!is.na(tax.res[i,]))[1] - 1]
          if(length(rank) != 0) if(rank != ranks[1]) tax.res[i, rank] <- parents$name[i]
        }
      }
    }
    workingIDs[stillWorking] = parents$parent
    rep = rep + 1
    if(rep>200)stop('Found cycle in taxonomy')
  }
  out<-tax.res[as.character(ids),,drop=FALSE]
  return(out)
}
