#!/usr/bin/env Rscript

rm(list=ls())

suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(seqinr))

args = commandArgs(trailingOnly=TRUE)
threads = as.numeric(args[4])
pident = as.numeric(args[5])

query <- read.fasta(args[2], as.string = T, forceDNAtolower = F, set.attributes = F)

set.seed(100)
taxann <- assignTaxonomy(seqs =  unlist(query), refFasta = args[1], 
          tryRC = T, outputBootstraps=F, minBoot = pident, multithread = threads, verbose = F, 
          taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
          )

# remove lengthy seqs from row names
queryid <- names(strsplit(rownames(taxann), " "))
taxodf <- as.data.frame(taxann)
for (i in 1:length(taxodf)) taxodf[,i] <- gsub('_', ' ', taxodf[,i])
names(taxodf)[1] <- "Domain"
taxodf$"#readid" <- queryid
write.table(taxodf[,c(7,1:6)], file = args[3], row.names=F, col.names=T, quote=F, sep='\t')


