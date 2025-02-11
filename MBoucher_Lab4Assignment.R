# DIY.transcriptomics Lab 4 ----
# Masen Boucher 
# February 5, 2025
# Submitted due to absence on date of lab session

# load packages----
library(rhdf5) #provides functions for handling hdf5 file formats (kallisto outputs bootstraps in this format)
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tximport) # package for getting Kallisto results into R
library(ensembldb) #helps deal with ensembl
library(datapasta) #great for copy/paste data into the R environment

# TASK ONE: Generating dataframe containing transcript information for the ferret genome ----
# get annotations using BiomaRt----
library(biomaRt) # an alternative for annotation

listMarts() #default host is ensembl.org, and most current release of mammalian genomes

#choose the 'mart' you want to work with
myMart <- useMart(biomart="ENSEMBL_MART_ENSEMBL")
#take a look at all available datasets within the selected mart
available.datasets <- listDatasets(myMart)
#now grab the ensembl annotations for ferret
ferret.anno <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "mpfuro_gene_ensembl")
ferret.attributes <- listAttributes(ferret.anno)

#For each transcript in the ferret genome, we want: 
#transcript ID, start position, end position, gene name, gene description, entrez gene ID, and pfam domains
Tx.ferret <- getBM(attributes=c('ensembl_transcript_id',
                         'start_position', 'end_position', 'external_gene_name', 'description',
                         'entrezgene_id', 'pfam'),
            mart = ferret.anno)

Tx.ferret <- as_tibble(Tx.ferret)

#we need to rename the columns we just retreived from biomart
Tx.ferret <- dplyr::rename(Tx.ferret, target_id = ensembl_transcript_id, start_position = start_position, end_position = end_position,
                    gene_name = external_gene_name, gene_description = description, entrez_gene_ID = entrezgene_id, PFAM_domains = pfam)
#Now we have our annotated data frame

# TASK TWO: Retrieving the ferret promoter sequences ---- 
# We will do this for our PI's 5 favorite antiviral genes -- IFIT2, OAS2, IRF1, IFNAR1, MX1
# To do this, we will use the getSequence() function from the BiomaRt package
IFIT2_promoter_seq <- getSequence(id = "IFIT2",  #putting in our gene ID
                                  type = "external_gene_name", #we are using the external_gene_name to search
                                  seqType = "gene_flank", #the promoter region would be a flanking sequence 
                                  upstream = 1000, #the flanking sequence is 1kb upstream, which is 1000 b upstream
                                  mart = ferret.anno) #using our ferret annotation

OAS2_promoter_seq <- getSequence(id = "OAS2",   #putting in our gene ID
                                  type = "external_gene_name", #we are using the external_gene_name to search
                                  seqType = "gene_flank", #the promoter region would be a flanking sequence
                                  upstream = 1000, #the flanking sequence is 1kb upstream, which is 1000 b upstream
                                  mart = ferret.anno) #using our ferret annotation

IRF1_promoter_seq <- getSequence(id = "IRF1",    #putting in our gene ID
                                 type = "external_gene_name", #we are using the external_gene_name to search
                                 seqType = "gene_flank",#the promoter region would be a flanking sequence
                                 upstream = 1000, #the flanking sequence is 1kb upstream, which is 1000 b upstream
                                 mart = ferret.anno) #using our ferret annotation

IFNAR1_promoter_seq <- getSequence(id = "IFNAR1",  #putting in our gene ID
                                 type = "external_gene_name", #we are using the external_gene_name to search
                                 seqType = "gene_flank", #the promoter region would be a flanking sequence
                                 upstream = 1000, #the flanking sequence is 1kb upstream, which is 1000 b upstream
                                 mart = ferret.anno) #using our ferret annotation

MX1_promoter_seq <- getSequence(id = "MX1",   #putting in our gene ID
                                   type = "external_gene_name", #we are using the external_gene_name to search
                                   seqType = "gene_flank", #the promoter region would be a flanking sequence
                                   upstream = 1000, #the flanking sequence is 1kb upstream, which is 1000 b upstream
                                   mart = ferret.anno) #using our ferret annotation
