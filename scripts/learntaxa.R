#!/usr/bin/env Rscript

rm(list=ls())

suppressPackageStartupMessages(library("DECIPHER"))

args = commandArgs(trailingOnly=TRUE)

# train classifyer
refseqs <- readDNAStringSet(args[1])
taxo = read.table(args[2], sep = "\t", comment.char = "", quote="" )

taxonomy <- paste0(taxo$V1, " Root;", taxo$V2)
trained_classifiyer <- LearnTaxa(train = refseqs, taxonomy = taxonomy, verbose = T)

save(trained_classifiyer, file = args[3])
