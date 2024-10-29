# NOTE: They changed some stuff in the NCBI database and there are no more 34 mtDNA for Neanderthals but 31,
# to get the names of the sequences from the msa I was using and to change the tips label of the tree I naively used the indeces,
# so if you run the script you need to take this into account
# I will fix this as soon as possible

# install genebankr (not available for BioConductor 3.19)
# BiocManager::install("gmbecker/genbankr")

library(rentrez)
library(seqinr)
library(Biostrings)
library(msa)
library(genbankr)
library(GenomicFeatures)
library(tidyverse)
library(parallel)
library(ape)
library(phangorn)
library(ggtree)

setwd("~/UCPH/Year2/Block1/BioProject1/Code/")
load("~/UCPH/Year2/Block1/BioProject1/Code/workspace_mut.RData")

# explore options to search for mtDNA of specific species -------------------------------------

# A few resoures to read on fetching data from online databases:

# https://docs.ropensci.org/rentrez/articles/rentrez_tutorial.html
# https://docs.ropensci.org/rentrez/reference/entrez_fetch.html
# https://gtpb.github.io/MEVR16/getseq/example1.html

# entrez_db_summary("nuccore")
# entrez_db_searchable("nuccore")

# Search queries to search at https://www.ncbi.nlm.nih.gov/nuccore/ before proceeding to
# automated fetching in R:

# We want to fetch mtDNA sequences from:
# - all Neanderthals
# - all Denisovans
# - all Sima de los Huesos early Neanderthals (see paper by Meyer et al. in Nature)
# - as large a representative sample of present-day human mtDNA sequences of *healthy individuals*
#   (do some text processing on abstracts of respective papers or metadata of the mtDNA
#    sequences to filter down to healthy individuals?)
# - include as many ancient modern human mtDNA sequences as well
#   (is this going to be helpful? https://amtdb.org)

# query <- "(015400[SLEN]:016600[SLEN]) AND Homo[Organism] AND mitochondrion[FILT]" # 67128 results
# query <- "(015400[SLEN]:016600[SLEN]) AND txid9604[orgn] AND mitochondrion[FILT]" # 67530 results
# query <- "(015400[SLEN]:016600[SLEN]) AND txid207598[orgn] AND mitochondrion[FILT]" # 67522 results

# Neanderthals only?
# query <- "(015400[SLEN]:016600[SLEN]) AND txid63221[orgn] AND mitochondrion[FILT]" # 34 results

# Denisovans only?
# query <- "(015400[SLEN]:016600[SLEN]) AND txid741158[orgn] AND mitochondrion[FILT]" # 8 results

# Sima de los Huesos?
# query <- "(015400[SLEN]:016600[SLEN]) AND txid1425170[orgn] AND mitochondrion[FILT]" # 2 results

# all homo sapiens? https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=9606
# query <- "(015400[SLEN]:016600[SLEN]) AND txid9606[orgn] AND mitochondrion[FILT]" # 2 results

# all chimpz
# query <- "(015400[SLEN]:016600[SLEN]) AND txid9598[orgn] AND mitochondrion[FILT])" # 185 results

# function to fetch sequences from the Nucletotide database and save them in gb and fasta file


fetch_seq <- function(query, dir_path, n=NULL){
  
  search <- entrez_search(db = "nuccore",term = query, use_history = TRUE)
  
  unlink(dir_path, recursive = TRUE, force = TRUE)
  dir.create(dir_path, recursive = TRUE)
  
  if (is.null(n)) {
    n=search$count
  }
  
  else if(search$count<n){
    n=search$count
  }
  
  # this condition it's needed in the case where we wanted to download a lot of sequences,
  # it's much more reasonable and faster to download them in chunks of fasta file due to limitations of the server,
  # we won't download the genebank file though
  
  if (n>=5000) {
    for( seq_start in seq(1,n,5000)){
      max=5000
      if(n-(seq_start-1)<5000){
        max=n-(seq_start-1)
      }
      recs <- entrez_fetch(db="nuccore", web_history=snail_coi$web_history,
                           rettype="fasta", retmax=max, retstart=seq_start-1)
      write(recs, file=paste0(path, seq_start, "-", seq_start+max-1, ".fa"))
      cat(seq_start+max-1, "sequences downloaded\r")
    }
  }
  
  else{
    for (i in seq(n)) {
      cat(sprintf("Processing sequence [%d/%d]\n", i, n))
      
      # fetch a GenBank file in a string format
      gb_txt <- entrez_fetch(db = "nuccore", web_history = search$web_history,
                             rettype = "gb", retmax = 1, retstart = i-1)
      
      cat("  - Fetched GenBank file from the internet\n")
      
      # extract accession # from the file for saving the data
      accession_id <-
        gb_txt %>%
        strsplit("\n") %>%
        { grep("ACCESSION", .[[1]], value = TRUE, ignore.case = TRUE) } %>%
        strsplit(" +") %>%
        { .[[1]][2] } # "ACCESSION <accession ID> <potentially other stuff>
      
      cat(sprintf("  - Extracted accession code %s\n", accession_id))
      
      gb_file <- file.path(dir_path, paste0(accession_id, ".gb"))
      fa_file <- file.path(dir_path, paste0(accession_id, ".fa"))
      
      cat(sprintf("  - Saving GenBank file to %s\n", gb_file))
      
      # write the GenBank to a file
      write(gb_txt, file = gb_file)
      
      cat(sprintf("  - Converting GenBank to FASTA file %s\n", fa_file))
      
      # convert to FASTA
      tryCatch(
        gb2fasta(gb_file, fa_file),
        error = function(e) {
          if (e$message == "argument of length 0") {
            cat("  - Empty FASTA string, skipping this record\n")
            unlink(gb_file)
            unlink(fa_file)
          } else {
            stop("Unknown error!", .call = FALSE)
          }
        }
      )
    }
  }
}

# fetch mtDNA Cambridge Reference Sequence: NC_012920
# the CRS has some rare variants and sequencing error (check link)
# https://www.mitomap.org/foswiki/bin/view/MITOMAP/CambridgeReanalysis
 
ref_query <- "NC_012920"
ref_path <- "data/mtdna/ref"
fetch_seq(ref_query, ref_path)

# fetch Neanderthals sequences
nean_query <- "(015400[SLEN]:016600[SLEN]) AND txid63221[orgn] AND mitochondrion[FILT]" # 34 results
nean_path <- "data/mtdna/neanderthals"
fetch_seq(nean_query, nean_path)


# fetch Denisovans sequences
deni_query <- "(015400[SLEN]:016600[SLEN]) AND txid741158[orgn] AND mitochondrion[FILT]" # 8 results
deni_path <- "data/mtdna/denisovans"
fetch_seq(deni_query, deni_path)

# fetch Sima de los Huesos sequences
sima_query <- "(015400[SLEN]:016600[SLEN]) AND txid1425170[orgn] AND mitochondrion[FILT]" # 2 results (repicate)
sima_path <- "data/mtdna/sima"
# n=1 because the two results are the same
fetch_seq(sima_query, sima_path, n=1)

# fetch all homo sapiens sequences
hs_query <- "(015400[SLEN]:016600[SLEN]) AND txid9606[orgn] AND mitochondrion[FILT] NOT txid1425170[orgn] NOT txid741158[orgn] NOT txid63221[orgn] NOT NC_012920[accn]" 
hs_path <- "data/mtdna/modern_humans"
# fetch only 100 sequences
fetch_seq(hs_query, hs_path, n=100)

# fetch all chimpz
chimp_query <- "(015400[SLEN]:016600[SLEN]) AND txid9598[orgn] AND mitochondrion[FILT])"
chimp_path <- "data/mtdna/chimpz"
# fetch only 10 sequences
fetch_seq(chimp_query, chimp_path, n=10)

# install msa packagage for multiple sequence alignment
# https://bioconductor.org/packages/release/bioc/vignettes/msa/inst/doc/msa.pdf

# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("msa")


# function to create set with all sequences
# pass the list of the files path

createStringSet <- function(files){
  combined_sequences <- DNAStringSet()
  for (file in files) {
    sequences <- readDNAStringSet(file)
    combined_sequences <- c(combined_sequences, sequences)
  }
  combined_sequences
}

# create a dnaStringSet with all the sequences we fetched

# get the files paths
all_paths <- list.files(path = "data/mtdna", pattern = "\\.fa$", full.names = TRUE, recursive = TRUE)
my_seqs <- createStringSet(all_paths)



# create a msa using all the sequences we fetched, method used is ClustalOmega
# keep the input order so to have different "taxonomic groups closer"

# run it in background with parallel
bg_msa <- mcparallel(msa(my_seqs, order = "input", method = "ClustalOmega"))
# collect the result, NULL if not done yet
my_msa <- mccollect(bg_msa, wait = FALSE)
my_msa <- my_msa[[1]]

# unmask the MsaDNAAlignment object to work in biostrings with DNAstringset object

unmask_msa <- unmasked(my_msa)

# without parallelization

# my_msa <- msa(my_seqs, order = "input", method = "ClustalOmega")
# my_msa

# alternative to run clustal omega in the terminal
# ADVANTAGES: more options and possibility of parallelization

# writeXStringSet(my_seqs, "data/")
# system("conda activate && clustalo -i data/my_seqs -o data/my_msa")
# my_msa2 <- readDNAStringSet("data/my_msa")


# 
# 
# 
# build the tree
# this is just for a visual idea
# to redo better


# this is temporary
# my_msa <- new("MsaDNAMultipleAlignment", unmasked = unmask_msa)
# this works if my_msa is a MsaDnaMultipleAlignment from the msa package
msaAln <- msaConvert(my_msa, type="seqinr::alignment")
d <- dist.alignment(msaAln, "identity")

mtdnaTree <- nj(d)

mtdnaTree <- midpoint(mtdnaTree)

mtdnaTree$tip.label[1:8]="d"
mtdnaTree$tip.label[9:108]="h"
mtdnaTree$tip.label[109:142]="n"
mtdnaTree$tip.label[143]="r"
mtdnaTree$tip.label[144]="s"

modernh_node <- getMRCA(mtdnaTree, which(mtdnaTree$tip.label== "r"|mtdnaTree$tip.label== "h"))
deni_node <- getMRCA(mtdnaTree, which(mtdnaTree$tip.label== "d"))
nean_node <- getMRCA(mtdnaTree, which(mtdnaTree$tip.label== "n"))
sima_node <- getMRCA(mtdnaTree, which(mtdnaTree$tip.label== "s"))

p2 <- ggtree(mtdnaTree, branch.length = "none") %>% scaleClade(modernh_node, 0.05) + 
  geom_hilight(node = deni_node, fill="red", to.bottom = TRUE)+
  geom_hilight(node = modernh_node, fill="darkgreen", to.bottom = TRUE, alpha=1, type = "gradient")+geom_hilight(node = nean_node, fill="blue", to.bottom = TRUE)+
  geom_cladelabel(node = deni_node, "Denisovans", angle = 90, hjust = "center", offset.text = 0.4, fontsize = 3.2)+
  geom_cladelabel(node=nean_node, "Neanderthals", angle = 90, hjust = "center", offset.text = 0.4, fontsize = 3.2)+
  geom_tiplab(aes(label="Sima",  subset=(label=="s")), size=3.2, colour="purple")+
  geom_cladelabel(node=modernh_node, "Present\nhumans", angle = 90, hjust = "center", offset.text = 0.5, fontsize = 3.2)
p2

ggsave("mtdna_tree.png", p2)



# extract annotated genes of the reference from the msa
# msa needs to be as msa object, not DNAStringSet
# !!! If the gap is within the gene the width is not changed
extract_genes <- function(msa){
  # read gb file of the reference mtDNA
  gb_ref <- readGenBank("data/mtdna/ref/NC_012920.gb")
  # extract gene details as dataframe
  genes <- data.frame(genes(gb_ref))
  # unmask the msa object to work on a DNAStringSet
  unmask_msa <- unmasked(msa)
  # identify locations of the gaps in the reference to adjust the indices of the genes
  gaps <- unlist(gregexpr("-",as.character(unmask_msa$`NC_012920 16569 bp`)))
  # adjust the indices according to the position of the gaps in the reference
  for (i in gaps) {
    genes$start <- genes$start+as.integer(genes$start>=i)
    genes$end <- genes$end+as.integer(genes$end>=i)
  }
  genes
}

genes_set <- extract_genes(my_msa)
genes_set

# write to fasta all the msa of each gene



# create different directories for cds, trna, rrna
cds_path <- "data/genes/cds"
trna_path <- "data/genes/trna"
rrna_path <- "data/genes/rrna"

unlink(cds_path, recursive = TRUE, force = TRUE)
dir.create(cds_path, recursive = TRUE)
unlink(trna_path, recursive = TRUE, force = TRUE)
dir.create(trna_path, recursive = TRUE)
unlink(rrna_path, recursive = TRUE, force = TRUE)
dir.create(rrna_path, recursive = TRUE)

for (i in rownames(genes_set)) {
  # save in fasta the subsets of the alignment containing only the genes
  if (startsWith(genes_set[i,]$gene_id, "R")) {
    writeXStringSet(subseq(unmask_msa, start = genes_set[i,]$start, end = genes_set[i,]$end), paste0(rrna_path, "/", genes_set[i,]$gene_id ,'.fa')) 
  }else if (startsWith(genes_set[i,]$gene_id, "T")) {
    writeXStringSet(subseq(unmask_msa, start = genes_set[i,]$start, end = genes_set[i,]$end), paste0(trna_path, "/", genes_set[i,]$gene_id ,'.fa')) 
  }else{
    writeXStringSet(subseq(unmask_msa, start = genes_set[i,]$start, end = genes_set[i,]$end), paste0(cds_path, "/", genes_set[i,]$gene_id ,'.fa'))
  }
}


# let's work only on the cds

cds_files <- list.files("data/genes/cds", full.names = FALSE)
cds_names <- substr(cds_files,1,nchar(cds_files)-3)
cds_list <- DNAStringSetList()

# create a list with the msa of the cds
for (file in cds_files) {
  # ND6 is the only gene on the - strand, since later I will translate them I reversed it now
  if(!startsWith(file, "ND6")){
    cds_list <- append(cds_list, DNAStringSetList(readDNAStringSet(paste0("data/genes/cds/", file))))
  }
  else{
    cds_list <- append(cds_list, DNAStringSetList(reverseComplement(readDNAStringSet(paste0("data/genes/cds/", file)))))
  }
}

names(cds_list) <- cds_names

# function to calculate the information content of a vector of probabilities (A,C,G,T frequencies per position in msa)
# based on Shannon entropy
# NOTE: I'm not gonna calcualte the information content but just the entropy

# function to calculate Shannon entropy in log base 2
# NOTE: the function receives as argument a column produced by the function Biostrings::consensusMatrix
# with the absolute frequencies (counts) of a column of a msa, by default the latter functions
# returns all possible IUPAC character for DNA, RNA seq (https://www.bioinformatics.org/sms/iupac.html).
# Since wherever I have ambiguities (N, R, etc.) I actually don't have information regarding that site I exlcuded those counts
# from the calculation of the entropy, on the other hand I included the gaps ("-")

# the argument prob is a column of probabilities
calc_entropy <- function(prob){
  # consider only non zero frequencies
  prob <- prob[prob>0]
  abs(sum(prob*log2(prob)))
}

# given a msa of dna sequences it calculates the entropy at each position
position_entropy_dna <- function(msa){
  # find consensus matrix as counts and calc the entropy for dna seq
  cons_matrix <- consensusMatrix(msa)
  freq <- cons_matrix[c("A", "C", "G", "T", "-"),]
  freq <- sweep(freq, 2, colSums(freq), "/")
  entropy <- apply(freq, 2, calc_entropy)
  entropy
}

# given a msa of protein sequences it calculates the entropy at each position
position_entropy_prot <- function(msa){
  # find consensus matrix as counts and calc the entropy for proteins seq
  cons_matrix <- consensusMatrix(msa)
  # consider only the 22 aminoacid, no ambiguities
  freq <- cons_matrix[1:22,]
  freq <- sweep(freq, 2, colSums(freq), "/")
  entropy <- apply(freq, 2, calc_entropy)
  entropy
}

# extract entropies for codon positions
# the function takes as argument the entropy for each position for a cds
extract_codons <- function(ent_pos){
  # number of codons excluding stop codon, check if integer because of polyA stop codon
  if(length(ent_pos)%%3==0){
    n_codons <- length(ent_pos)/3-1 
  }
  else{
    n_codons <- as.integer(length(ent_pos)/3)
  }
  df_codons <- data.frame(matrix(ent_pos[1:(n_codons*3)], ncol = 3, byrow = TRUE))
}

# entropies over the list of proteins' CDS
cds_entropies <- lapply(cds_list, position_entropy_dna)
codons_entropies <- lapply(cds_entropies, extract_codons)

# create a matrix with all the entropies of the proteins' CDS together 
# to study the distribution according to the position

# bind all the df for each protein
tot_entropies <- do.call(rbind, codons_entropies)
colnames(tot_entropies) <- c("1st position", "2nd position", "3rd position")

# function to translate a cds into the aa sequence
translation <- function(msa){
  # number of codons excluding stop codon, check if integer because of polyA stop codon
  if(width(msa)[1]%%3==0){
    n_codons <- width(msa)[1]/3-1 
  }
  else{
    n_codons <- as.integer(width(msa)[1]/3)
  }
  
  # extract cds without stop codon
  prot2trans <- subseq(msa, start = 1, end = n_codons*3)
  
  # translate
  Biostrings::translate(prot2trans, genetic.code = getGeneticCode("2"), if.fuzzy.codon = "X")
}

# translate the proteins
proteins_list <- lapply(cds_list, translation)
proteins_list

# calculate the entropies
prot_entropies <- lapply(proteins_list, position_entropy_prot)
prot_entropies_bind <- unlist(prot_entropies)

# let's work on rrna

rrna_files <- list.files("data/genes/rrna", full.names = FALSE)
rrna_names <- substr(rrna_files,1,nchar(rrna_files)-3)
rrna_list <- DNAStringSetList()

for (file in rrna_files) {
  rrna_list <- append(rrna_list, DNAStringSetList(readDNAStringSet(paste0("data/genes/rrna/", file))))
}

names(rrna_list) <- rrna_names

# calculate the entropies

rrna_entropies <- lapply(rrna_list, position_entropy_dna)
rrna_entropies_bind <- unlist(rrna_entropies)


# let's work on trna

trna_files <- list.files("data/genes/trna", full.names = FALSE)
trna_names <- substr(trna_files,1,nchar(trna_files)-3)
trna_list <- DNAStringSetList()

for (file in trna_files) {
  trna_list <- append(trna_list, DNAStringSetList(readDNAStringSet(paste0("data/genes/trna/", file))))
}

names(trna_list) <- trna_names

# calc the entropies

trna_entropies <- lapply(trna_list, position_entropy_dna)
trna_entropies_bind <- unlist(trna_entropies)

# plot distribution of entropy per position

# a bit messy but it is to set the dataframe for the plot
tot_entropies_long <- tot_entropies %>% pivot_longer(cols = everything(), names_to = "Position", values_to = "Entropy")
df_protein <- data.frame(Position = rep("Protein", length(prot_entropies_bind)), Entropy=prot_entropies_bind)
df_rrna <- data.frame(Position = rep("rRNA", length(rrna_entropies_bind)), Entropy=rrna_entropies_bind)
df_trna <- data.frame(Position = rep("tRNA", length(trna_entropies_bind)), Entropy=trna_entropies_bind)
entropies_def <- rbind(tot_entropies_long, df_protein, df_rrna, df_trna)

plot_entropies <- ggplot(entropies_def, aes(x = Entropy)) +
  geom_density() +
  facet_wrap(~Position, scales = "fixed") +
  labs(x = "Entropy", y = "Density") +
  theme_minimal()

ggsave("plot_entropies.png", plot_entropies)

# 
# 
# 
# 
# This section is everything about finding specific mutations in the MSA

d_names <- names(proteins_list$ATP6)[1:8]
mh_names <- names(proteins_list$ATP6)[9:108]
n_names <- names(proteins_list$ATP6)[109:142]
r_names <- names(proteins_list$ATP6)[143]
s_names <- names(proteins_list$ATP6)[144]


# 
# 
# 
# function to identify the non conserved columns in a msa

find_var_pos <- function(myMsa, type="AA"){
  myMsa_matrix <- as.matrix(myMsa)
  
  # find positions of non conserved columns
  if(type=="AA"){
    mismatch_positions <- which(apply(myMsa_matrix, 2, function(col) length(unique(col)) > 1 & !(any(unique(col)=="X") & length(unique(col))==2)))  
  }
  else if(type=="DNA"){
    mismatch_positions <- which(apply(myMsa_matrix, 2, function(col) length(unique(col)) > 1 & !(any(unique(col)=="N") & length(unique(col))==2)))
  }
  
  # return something only if the msa is not totally conserved, otherwise is NULL
  if(!identical(mismatch_positions, integer(0))){
    # extract only non conserved columns
    non_cons_cols <- myMsa_matrix[, mismatch_positions, drop=FALSE]
    colnames(non_cons_cols) <- mismatch_positions
    non_cons_cols
  }
  
}

# try to identify common known diseases in the msa

# function to create a fasta file to submit in http://mutpred.mutdb.org/index.html
# in order to see if there are any potential deeterious mutation

# the functions are a bit a mess

# create_names takes as input a matrix with as colnames the position of the substituion
# the first row corresponds to the aminoacid in the reference while the second the aminoacid in the target group
# densiovans+sima in this case. It converts 
# example:
#      7   45  51  59  112 185
# [1,] "A" "T" "K" "T" "T" "N"
# [2,] "T" "A" "Q" "A" "A" "S"
# The function creates a file containing all the reference sequence of the proteins containing amino acid substitutions
# with a header of this type:
# ATP6 A7T T45A etc.
# where ATP6 is the name of the protein and the rest are the aa substitution that we want to investigate


create_names <- function(x){
  subst <- unlist(lapply(colnames(x), function(y) paste0(x[1, y], y, x[2, y])))
  paste(subst, collapse = " ")
}

create_mutpred_file <- function(prot_list, target_names, ref_name, path){
  
  target_consseq <- lapply(prot_list, function(x) unlist(strsplit(consensusString(x[target_names]), "")))
  ref_prot_list <- lapply(prot_list, function(x) unlist(strsplit(as.character(x[[ref_name]]), "")))
  prot_names <- names(ref_prot_list)
  mutations <- lapply(prot_names, function(x) find_var_pos(matrix(c(ref_prot_list[[x]], target_consseq[[x]]), nrow = 2, byrow = TRUE)))
  names(mutations) <- prot_names
  mutations <- Filter(Negate(is.null), mutations)
  prot_names <- names(mutations)
  ref_prot_list <- lapply(ref_prot_list, function(x) paste(x, collapse = ""))
  ref_prot_list <- ref_prot_list[prot_names] 
  names(ref_prot_list) <- NULL
  ref_prot_list <- do.call(c, AAStringSetList(ref_prot_list))
  mutpred_names <- lapply(mutations, function(x) create_names(x))
  def_names <- sapply(1:length(mutations), function(i) paste(prot_names[[i]], mutpred_names[[i]], sep = " "))
  names(ref_prot_list) <- def_names
  writeXStringSet(ref_prot_list, path)
}


target_names <- c(d_names, s_names)
create_mutpred_file(proteins_list, target_names, r_names, "data/mutpred.fa")

# 
# 
# function to find variations between two sequences (target and refernce) in a msa
# it returns a vector with all the variations, I followed this nomenclature https://www.hgmd.cf.ac.uk/docs/mut_nom.html
# NOTE: it's pretty common to find different nomenclature in different papers
find_variants <- function(my_msa, target, ref="NC_012920 16569 bp"){
  # ref seq goes in front
  seqs <- as.matrix(my_msa[c(ref,target)])
  variants_list <- apply(seqs, 2, unique)
  gaps <- unlist(gregexpr("-",as.character(unmask_msa$`NC_012920 16569 bp`)))
  mutations <- c()
  i <- 1
  # variants_list
  while (i<=length(variants_list)) {
    # if in the position there is only 1 character (no difference) skip the iteration
    if(length(variants_list[[i]])==1 |  "N" %in% variants_list[[i]]){
      i <- i+1
    }
    # Insertion case: if there is a gap in the reference
    else if(variants_list[[i]][1]=="-"){
      # position after the beginning of the gap
      if(i+1<=length(variants_list)){
        c <- i+1  
      }
      # if the insertion is longer than one base we need to check how long it is
      while (variants_list[[c]][1]=="-" & c+1<=length(variants_list) & length(variants_list[[c]])==2) {
        c <- c+1
      }
      # adjust the insertion indices according to the position in the ref seq not in the msa, as always
      ins_start <- i-sum(gaps<i)-1
      ins_end <- ins_start+1
      # extract the bases inserted in the target seq
      ins_bases <- paste0(sapply(variants_list[i:(c-1)], function(x) x[2]), collapse = "")
      mutations <- c(mutations, paste0(ins_start, "_", ins_end, "ins", ins_bases))
      # this is needed in case we have a gap in the last position, otherwise if don't do this we get stuck in the loop
      if(!is.na(variants_list[[c]][2]) & variants_list[[c]][2]=="-" & c+1>length(variants_list)){
        i <- c+1  
      }
      else{
        i <- c
      }
    }
    # Deletion case: if there is a gap in the target
    else if(variants_list[[i]][2]=="-"){
      # position after the beginning of the gap
      if(i+1<=length(variants_list)){
        c <- i+1  
      }
      # if the deletion is longer than one base we need to check how long it is
      # the other conditions are to avoid to go out of bound when if at the end of the alignment
      while (!is.na(variants_list[[c]][2]) & variants_list[[c]][2]=="-" & c+1<=length(variants_list) & length(variants_list[[c]])==2) {
        c <- c+1
      }
      # adjust the insertion indices according to the position in the ref seq not in the msa, as always
      del_start <- i-sum(gaps<i)
      del_end <- del_start+(c-i)-1
      # one base del
      if(del_start==del_end){
        del_bases <- variants_list[[i]][1]
        mutations <- c(mutations, paste0(del_start, "del", del_bases))
      }
      # more than one base deletion
      else{
        # extract the bases inserted in the target seq
        del_bases <- paste0(sapply(variants_list[i:(c-1)], function(x) x[1]), collapse = "")
        mutations <- c(mutations, paste0(del_start, "_", del_end, "del", del_bases))
      }
      # this is needed in case we have a gap in the last position, otherwise if don't do this we get stuck in the loop
      if(!is.na(variants_list[[c]][2]) & variants_list[[c]][2]=="-" & c+1>length(variants_list)){
        i <- c+1  
      }
      else{
        i <- c
      }
    }
    # subsituition case
    else{
      # adjust index according to the msa
      sub_i <- i-sum(gaps<i)
      mutations <- c(mutations, paste0(sub_i, variants_list[[i]][1], ">", variants_list[[i]][2]))
      i <- i+1
    }
  }
  mutations
}


mutations_list <- sapply(c(d_names, mh_names, n_names, s_names), function(x) find_variants(unmask_msa, x, r_names))


# given the mutation recorded in https://www.mitomap.org/foswiki/bin/view/MITOMAP/MutationsCodingControl
# let's see if they are present in our msa
# this function crates a dataframe containing the counts of how many potentially deleterious mutations are contained
# in each sequence of the msa
# NOTE: it's not perfect when working with indel, I need more time to think about how to face that, cause the nomencalture
# changes from db to db and it's a bit a mess

scan_mutations <- function(my_msa, file_path){
  disease_df <- read_csv(file_path)
  
  # let's look only at the point mutations and point gap
  # the nomenclature can be a bit different sometimes, I'm then considering only the one base deletion
  point_disease_df <- disease_df %>% filter(grepl( "^[A-Za-z]\\d+[A-Za-z]$|^m\\.\\d+[A-Za-z]?>[A-Za-z]$|[A-Za-z]\\d+del|[A-Za-z]\\.\\d+del
",Allele))
  
  # adjust position according to ref seq in msa, as always
  gaps <- unlist(gregexpr("-",as.character(unmask_msa$`NC_012920 16569 bp`)))
  for (i in gaps) {
    point_disease_df$Position<- point_disease_df$Position+as.integer(point_disease_df$Position>i)
  }
  
  # extract potentially deleterious alleles
  
  del_allels <- str_sub(point_disease_df$Allele, -1,-1)
  # change the l (last letter of del) into the gap symbol
  del_allels[del_allels=="l"] <- "-"
  
  # scan our msa to see which of those is present
  
  species=c(rep("Denisovans", 8), rep("Present Humans", 100), rep("Neanderthals", 34), "Sima")
  mut_count <- data.frame(sample=names(my_msa), count=numeric(length(names(my_msa))), species=species)
  
  for (i in 1:length(my_msa)) {
    # vector that tells if the mutation is present or no
    bool_mut <- unlist(str_split(as.character(my_msa[[i]][point_disease_df$Position]), ""))==del_allels
    # count of how many mutations are present
    count <- sum(bool_mut)
    # identification of the mutation position
    # df_pos is the index in the df of the positions, i'll do this to make it easier to have a look at the position
    df_pos <- which(bool_mut)
    # msa_pos is the position of the mutation in the msa
    msa_pos <- point_disease_df$Position[df_pos]
    mut_count$count[i] <- count
    mut_count$df_positions[i] <- list(df_pos)
    mut_count$msa_positions[i] <- list(msa_pos)
  }
  list(df_count=mut_count,
       df_disease=point_disease_df)
}

# on reported cds and control region substitutions
cds_ctrl_subs <- scan_mutations(my_msa = unmask_msa[c(d_names, mh_names, n_names, s_names)], file_path = "data/patho_mut/MutationsCodingControl MITOMAP Foswiki.csv")

# cds_ctrl_subs$df_disease[cds_ctrl_subs$df_count$df_positions[[6]],]$Status

# on reported rnas substitutions
rna_trna_subs <- scan_mutations(my_msa = unmask_msa[c(d_names, mh_names, n_names, s_names)], file_path = "data/patho_mut/MutationsRNA MITOMAP Foswiki.csv")


cds_ctrl_plot <- ggplot(cds_ctrl_subs$df_count, aes(x=reorder(sample, species, FUN=first), y=count, fill=species))+
  geom_bar(stat = "identity",  width = 1)+ theme_minimal() + theme(axis.text.x=element_blank(),
                                      axis.ticks.x=element_blank())+
  xlab("Samples")+ylab("Count of potentially deleterious alleles")+ggtitle("Coding and control regions")

cds_ctrl_plot

t_rrna_plot <- ggplot(rna_trna_subs$df_count, aes(x=reorder(sample, species, FUN=first), y=count, fill=species))+
  geom_bar(stat = "identity",  position = "dodge", width = 1) + theme_minimal()  + theme(axis.text.x=element_blank(),
                                      axis.ticks.x=element_blank())+
  xlab("Samples")+ylab("Number of potentially deleterious alleles")+ggtitle("tRNA and rRNA")

t_rrna_plot

# normalize version
# create normalizing vector that contains the distance in number of mutations from the rCRS

normalizing_vec <- sapply(mutations_list, length)

cds_ctrl_plot_norm <- ggplot(cds_ctrl_subs$df_count, aes(x=reorder(sample, species, FUN=first), y=count/normalizing_vec, fill=species))+
  geom_bar(stat = "identity",  width = 1)+ theme_minimal() + theme(axis.text.x=element_blank(),
                                                                   axis.ticks.x=element_blank())+
  xlab("Samples")+ylab("Count of potentially deleterious alleles")+ggtitle("coding and control regions - normalized")

cds_ctrl_plot_norm

t_rrna_plot_norm <- ggplot(rna_trna_subs$df_count, aes(x=reorder(sample, species, FUN=first), y=count/normalizing_vec, fill=species))+
  geom_bar(stat = "identity",  position = "dodge", width = 1) + theme_minimal()  + theme(axis.text.x=element_blank(),
                                                                                         axis.ticks.x=element_blank())+
  xlab("Samples")+ylab("Number of potentially deleterious alleles")+ggtitle("tRNA and rRNA - normalized")

t_rrna_plot_norm

# normalizing the total number of variants, both in cds and rnas
tot_subs <- rna_trna_subs
tot_subs$df_count$count <- rna_trna_subs$df_count$count + cds_ctrl_subs$df_count$count 


tot_plot_norm <- ggplot(tot_subs$df_count, aes(x=reorder(sample, species, FUN=first), y=count/normalizing_vec, fill=species))+
  geom_bar(stat = "identity", width = 1)+ theme_minimal() + theme(axis.text.x=element_blank(),
                                      axis.ticks.x=element_blank())+
  xlab("Samples")+ylab("Number of potentially deleterious alleles")+ggtitle("tRNA, rRNA, coding, control regions - normalized ")

tot_plot_norm


ggsave("cds_ctrl_mitomap.png", cds_ctrl_plot)
ggsave("trna_rrna_mitomap.png", t_rrna_plot)
ggsave("cds_ctrl_norm_mitomap.png", cds_ctrl_plot_norm)
ggsave("trna_rrna_mitomap_norm.png", t_rrna_plot_norm)
ggsave("tot_norm_mitomap.png", tot_plot_norm)


# other normalization
# normalize according to the number of mutation in that category

gb_ref <- readGenBank("data/mtdna/ref/NC_012920.gb")
genes <- data.frame(genes(gb_ref))

cds_ranges <- genes %>% filter(!(startsWith(gene_id, "R") | startsWith(gene_id, "T"))) %>% select(start, end)
rnas_ranges <- genes %>% filter((startsWith(gene_id, "R") | startsWith(gene_id, "T"))) %>% select(start, end)

mutations_positions <- lapply(mutations_list, function(x) sub("^([0-9]+).*", "\\1", x))
norm_vec_cds <- lapply(mutations_positions, function(y) apply(cds_ranges, 1, function(x) sum(y>x[1] & y<x[2])))
norm_vec_cds <-  sapply(norm_vec_cds, sum)

norm_vec_rna <- lapply(mutations_positions, function(y) apply(rnas_ranges, 1, function(x) sum(y>=x[1] & y<=x[2])))
norm_vec_rna <-  sapply(norm_vec_rna, sum)

cds_ctrl_plot_norm <- ggplot(cds_ctrl_subs$df_count, aes(x=reorder(sample, species, FUN=first), y=count/norm_vec_cds, fill=species))+
  geom_bar(stat = "identity",  width = 1)+ theme_minimal() + theme(axis.text.x=element_blank(),
                                                                   axis.ticks.x=element_blank())+
  xlab("Samples")+ylab("Count of potentially deleterious alleles")+ggtitle("coding and control regions - normalized")

cds_ctrl_plot_norm

t_rrna_plot_norm <- ggplot(rna_trna_subs$df_count, aes(x=reorder(sample, species, FUN=first), y=count/norm_vec_rna, fill=species))+
  geom_bar(stat = "identity",  position = "dodge", width = 1) + theme_minimal()  + theme(axis.text.x=element_blank(),
                                                                                         axis.ticks.x=element_blank())+
  xlab("Samples")+ylab("Number of potentially deleterious alleles")+ggtitle("tRNA and rRNA - normalized")

t_rrna_plot_norm

ggsave("cds_plot_norm.png", cds_ctrl_plot_norm)
ggsave("rna_plot_norm.png", t_rrna_plot_norm)
# 
# 
# 
# 
# 
# Check if completely conserved column in protein for modern humans are not so conserved in Denisovans and Neanderthals


# function to count the number of non conserved sites that are completely conserved in modern humans
# target is what we want to check (ex. denisovans, neanderthals etc), reference is modern humans
count_non_conserved <- function(target, reference){
  vec_target <- unlist(target)
  vec_ref <- unlist(reference)
  if(length(vec_target)!=length(vec_ref)) {
    stop("Something is wrong in the length")
  }
  sum(vec_target[vec_ref!=0]!=0)
}

# function to build the plot, not well done
build_count_plot <- function(targ){
  # proteins
  prot_entropies_mh <- lapply(proteins_list, function(prot_msa) position_entropy_prot(prot_msa[mh_names]))
  prot_entropies_den <- lapply(proteins_list, function(prot_msa) position_entropy_prot(prot_msa[targ]))
  den_count_prot <- count_non_conserved(prot_entropies_den, prot_entropies_mh)
  # same but on trna
  trna_entropies_mh <- lapply(trna_list, function(trna_msa) position_entropy_dna(trna_msa[mh_names]))
  trna_entropies_den <- lapply(trna_list, function(trna_msa) position_entropy_dna(trna_msa[targ]))
  den_count_trna <- count_non_conserved(trna_entropies_den, trna_entropies_mh)
  # rrna
  rrna_entropies_mh <- lapply(rrna_list, function(rrna_msa) position_entropy_dna(rrna_msa[mh_names]))
  rrna_entropies_den <- lapply(rrna_list, function(rrna_msa) position_entropy_dna(rrna_msa[targ]))
  den_count_rrna <- count_non_conserved(rrna_entropies_den, rrna_entropies_mh)
  # neanderthals
  # set counter to 0
  nean_count_rrna <- 0
  nean_count_trna <- 0
  nean_count_prot <- 0
  for (i in 1:50) {
    n_names_samp <- sample(n_names, length(targ))
    trna_entropies_nean <- lapply(trna_list, function(trna_msa) position_entropy_dna(trna_msa[c(n_names_samp)]))
    nean_count_trna <- nean_count_trna + count_non_conserved(trna_entropies_nean, trna_entropies_mh)
    rrna_entropies_nean <- lapply(rrna_list, function(rrna_msa) position_entropy_dna(rrna_msa[c(n_names_samp)]))
    nean_count_rrna <- nean_count_rrna + count_non_conserved(rrna_entropies_nean, rrna_entropies_mh)
    prot_entropies_nean <- lapply(proteins_list, function(prot_msa) position_entropy_prot(prot_msa[n_names_samp]))  
    nean_count_prot <- nean_count_prot + count_non_conserved(prot_entropies_nean, prot_entropies_mh)
  }
  nean_count_prot <- nean_count_prot/50
  nean_count_trna <- nean_count_trna/50
  nean_count_rrna <- nean_count_rrna/50

  count_df <- tibble(Count=c(den_count_prot, nean_count_prot, den_count_rrna, nean_count_rrna, den_count_trna, nean_count_trna), 
                     Species=rep(c("Denisovans", "Neanderthals"), 3), Type=c("Protein", "Protein", "rRNA",  "rRNA",  "tRNA",  "tRNA"))
  
  count_plot <- ggplot(count_df, aes(Species, Count))+geom_bar(stat = "identity")+facet_grid(~Type) + theme_minimal()
  count_plot
}

# densiovans + sima
count_plot_den_sim <- build_count_plot(c(d_names, s_names))
count_plot_den_sim

# Only densiovans
count_plot_den <- build_count_plot(d_names)
count_plot_den

# Whatever it comes after this point is just trying stuff


# anticodon test


extract_anticodons <- function(msa){
  # read gb file of the reference mtDNA
  gb_ref <- readGenBank("data/mtdna/ref/NC_012920.gb")
  # extract anticodons detaails
  trna_df <- data.frame(gb_ref@other_features) %>% filter(type=="tRNA")
  trna_product <- trna_df$product
  # extract only the anticodon start-end to adjust according to msa
  numbers <- gregexpr("[0-9]+", trna_df$anticodon)
  result <- regmatches(trna_df$anticodon, numbers)
  ant_pos <- data.frame(matrix(as.numeric(unlist(result)), ncol = 2, byrow = TRUE))
  colnames(ant_pos) <- c("start", "end")
  # unmask the msa object to work on a DNAStringSet
  unmask_msa <- unmasked(msa)
  # identify locations of the gaps in the reference to adjust the indexes of the genes
  gaps <- unlist(gregexpr("-",as.character(unmask_msa$`NC_012920 16569 bp`)))
  # add the shift do to the gaps inserted in the reference by the msa
  for (i in gaps) {
    ant_pos$start <- ant_pos$start+as.integer(ant_pos$start>i)
    ant_pos$end <- ant_pos$end+as.integer(ant_pos$end>i)
  }
  cbind(ant_pos, trna_product)
}

ant_pos <- extract_anticodons(my_msa)

# create a list of msa of anticodons
anticodon_list <- apply(ant_pos[c("start", "end")], 1, function(pos) subseq(unmask_msa, start = pos[1], end = pos[2]))
names(anticodon_list) <- ant_pos$trna_product

#check entropy in anticodons 
anticodon_entropies <- lapply(anticodon_list, position_entropy_dna)
anticodon_entropies




# 
# 


# alternative way to count mutations, not really perfect 

rna_disease_df <- read_csv("data/patho_mut/MutationsRNA MITOMAP Foswiki.csv")
cds_disease_df <- read_csv("data/patho_mut/MutationsCodingControl MITOMAP Foswiki.csv")

reformat_mutations_rna <- function(mutations) {
  sapply(mutations, function(mutation) {
    gsub("([A-Z])(\\d+)([A-Z])", "\\2\\1>\\3", mutation)
  })
}

reformat_mutations_cds <- function(mutations) {
  sapply(mutations, function(mutation) {
    gsub("^m\\.", "", mutation)
  })
}

trial <- reformat_mutations_rna(rna_disease_df$Allele)
trial_2 <- reformat_mutations_cds(cds_disease_df$Allele)

length(intersect(mutations_list[[1]], trial_2))



# save aa sequence of protein to predict with AlphaFold: https://alphafoldserver.com/

# save all the consensus sequences of Denisovans + Sima lineage and rCRS sequences

for (prot in names(proteins_list)) {
  unlink(paste0("data/proteins2fold/", prot), recursive = TRUE, force = TRUE)
  dir.create(paste0("data/proteins2fold/", prot), recursive = TRUE)
  deni_seq <- consensusString(proteins_list[[prot]][c(d_names, s_names)]) 
  crs_seq <- as.character(proteins_list[[prot]][r_names])
  write(deni_seq, file = paste0("data/proteins2fold/", prot, "/", prot,  "_deni.txt"))
  write(crs_seq, file = paste0("data/proteins2fold/", prot, "/", prot,  "_crs.txt"))
}

# write a file to store the number of the residues with the mutation
# this is needed to highlight them in pymol

for (prot in names(proteins_list)) {
  deni_seq <- consensusString(proteins_list[[prot]][c(d_names, s_names)]) 
  crs_seq <- as.character(proteins_list[[prot]][r_names])
  # get the position of the mutation
  res_mutated <- colnames(find_var_pos(AAStringSet(c(crs_seq, deni_seq))))
  cat(res_mutated, file = paste0("data/proteins2fold/", prot, "/", prot,  "_mut_resi.txt"), sep = " ")
}





