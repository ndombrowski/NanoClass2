#!/usr/bin/env Rscript

rm(list=ls())

suppressPackageStartupMessages(library("DECIPHER"))

args = commandArgs(trailingOnly=TRUE)

threads = as.numeric(args[4])

pctth = as.numeric(args[5])

query <- readDNAStringSet(args[2])
db <- load(args[1])

taxann <- IdTaxa(query, trained_classifiyer, strand = "both", threshold = pctth, processors = threads, verbose = F)

assignment <- sapply(taxann, function(x) paste(x$taxon, collapse="\t"))
names(assignment) <- sapply(strsplit(names(assignment)," "), `[`, 1)

dl <- lapply(assignment, function (x) {gsub("D_[0-6]__", "", x)})
dl <- lapply(dl, function (x) {gsub("Root\t", "", x)})

df <- setNames(as.data.frame(unlist(dl)),"Domain\tPhylum\tClass\tOrder\tFamily\tGenus")
df$'#readid' <- rownames(df)

write.table(df[,c(2,1)], file=args[3], row.names=F, col.names=T, quote=F, sep='\t')
