
<!-- README.md is generated from README.Rmd. Please edit that file -->

# customSeqDB

<!-- badges: start -->

<!-- badges: end -->

The goal of customSeqDB is to help make custom sequences databases from
open source database like NCBI and ENA. Functions allow to search and
fetch records on their API and parse these records for relevant
information. It can also make a taxonomy SQLite file to allow to fetch
taxonomic information for the records used to make the custom database.

## Installation

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("salix-d/customSeqDB")
```

## Example

### Getting the sequences details from NCBI

This gets the ESummary of all records in the NCBI nucleotide database
matching these taxa and gene which have the word ‘complete’ in the
‘title’ fields but don’t have ‘UNVERIFIED’ in the ‘keyword’ fields.
The query isn’t case sensitive. For the usable fields and their name see
: <https://www.ncbi.nlm.nih.gov/books/NBK49540/>

It then parse the data and returns a list with history server
information in case we want to look at the records again and a list of
details including the accession numbers, the taxonomic ids and the
description. The ESearch can return up to 100k results. This part may
take a while.

``` r
library(customSeqDB)
taxa <- c("Trombidiformes", "Holothyrida")
gene <- c("COI", "CO1", "COXI", "COX1") # ncbi includes 6 orders and the COI gene can be written 4 ways
ncbi <- get_ncbiDetails(taxa = taxa, gene = gene, customTerm = "AND complete [TITL] NOT unverfified[KYWD]")
#> 
#>     searching for records matching :(Trombidiformes[ORGN] OR Holothyrida[ORGN]) AND (COI[GENE] OR CO1[GENE] OR COXI[GENE] OR COX1[GENE]) AND complete [TITL] NOT unverfified[KYWD]
#>      195  records found
#>     parsing records 0:195 of 195
head(ncbi$details)
#>     accession   taxid
#> 1  LC552027.1 2759127
#> 2  LC552026.1 2740590
#> 3 NC_049059.1 2740590
#> 4  EU345430.1   32264
#> 5 NC_039813.1 2047715
#> 6  MG701313.1 2047715
#>                                                      description
#> 1 Hygrobates taniguchii NMsp2 mitochondrial DNA, complete genome
#> 2 Hygrobates longiporus NMsp1 mitochondrial DNA, complete genome
#> 3 Hygrobates longiporus NMsp1 mitochondrial DNA, complete genome
#> 4             Tetranychus urticae mitochondrion, complete genome
#> 5               Sperchon plumifer mitochondrion, complete genome
#> 6               Sperchon plumifer mitochondrion, complete genome
```

More examples will be added.
