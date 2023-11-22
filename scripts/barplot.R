#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(tidyr))

args = commandArgs(trailingOnly=TRUE)

taxmat <- as.data.frame(unique(vroom(delim = '\t', args)))
taxmat <- taxmat[!is.na(taxmat$Domain),]

taxmat <- unique(taxmat)
names(taxmat)[1] <- "taxid"
write.table(taxmat, file="tables/taxonomy-table.tsv", row.names=F, col.names=T, sep='\t', quote=F)

rownames(taxmat) <- taxmat$taxid
taxmat$taxid <- NULL
taxmat <- as.matrix(taxmat)

file <- gsub(".taxmat$", ".otumat", args)
sam <- data.frame(run = rep(NA, length(file)), sample = rep(NA, length(file)), method = rep(NA, length(file)))

for (i in 1:length(file)){
  otumatje <- read.table(file[i], header = T, sep = '\t', comment = "")
  names(otumatje)[1] <- "taxid"
  names(otumatje)[2] <- paste0(strsplit(file[i], "/")[[1]][2],"_",names(otumatje)[2])
  ifelse (i == 1, otumat <- otumatje, otumat <- merge(otumat, otumatje, all=TRUE))

  lab = strsplit(file[i], "/")[[1]][4]
  sam$run[i] = strsplit(file[i], "/")[[1]][2]
  sam$sample[i] = strsplit(lab, "[.]")[[1]][1]
  sam$method[i] = strsplit(lab, "[.]")[[1]][2]
} 

otumat[is.na(otumat)] <- 0

## for some custom DB like BOLD local copy, the same taxonomic lineage may occur 2x.
## This causes an error when they are used as rownames, which should be unique --> fix: aggregate
otumat <- aggregate(otumat[2:length(otumat)], by=list(otumat$taxid), sum)
names(otumat)[1] <- "taxid"
write.table(otumat, file="tables/otu-table.tsv", row.names=F, col.names=T, sep='\t', quote=F)

rownames(otumat) <- otumat$taxid
otumat$taxid <- NULL
otumat <- as.matrix(otumat)

rownames(sam) <- colnames(otumat)


OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
SAM = sample_data(sam)

physeq = phyloseq(OTU, TAX, SAM)
pphyseq  = transform_sample_counts(physeq, function(x) x / sum(x) )

theme_set(theme_bw())

for (level in colnames(taxmat)){
    top.taxa <- tax_glom(physeq, level)
    TopNOTUs <- names(sort(taxa_sums(top.taxa), TRUE)[1:20])
    BottumNOTUs <- names(taxa_sums(top.taxa))[which(!names(taxa_sums(top.taxa)) %in% TopNOTUs)]
    merged_physeq = merge_taxa(top.taxa, BottumNOTUs, 2)

    mdf = psmelt(merged_physeq); names(mdf)[names(mdf) == level] <- "level"
    mdf$OTU[which(is.na(mdf$level))] <- "aaaOther"
    mdf$level[which(is.na(mdf$level))] <- "aaaOther"
    aggr_mdf <- aggregate(Abundance ~ sample + run + method + level, data = mdf, sum)

    labs = aggr_mdf$level; labs[labs=="aaaOther"] <- "Other"
    cols = scales::hue_pal()(length(unique(labs))); cols[unique(labs) == "Other"] <- "#CCCCCC"

    p = ggplot(aggr_mdf, aes_string(x = "method", y = "Abundance", fill = "level"))
    p = p + scale_fill_manual(name = level, labels = unique(labs), values = cols)
    p = p + facet_grid(paste0("", aggr_mdf$run) ~ paste0("", aggr_mdf$sample))
    p = p + geom_bar(stat = "identity", position = "stack",  color = "black", size = 0.1)
    p = p + guides(fill=guide_legend(ncol=1))
    p = p + labs(x = "Method", y = "Absolute abundance")
    p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0, size = 5))

    ggsave(paste0("plots/aabund-", level, "-by-sample.pdf"), plot=p, device="pdf")

    p = ggplot(aggr_mdf, aes_string(x = "sample", y = "Abundance", fill = "level")) 
    p = p + scale_fill_manual(name = level, labels = unique(labs), values = cols)
    p = p + facet_grid(paste0("", aggr_mdf$run) ~ paste0("", aggr_mdf$method))
    p = p + geom_bar(stat = "identity", position = "stack",  color = "black", size = 0.1)
    p =	p + guides(fill=guide_legend(ncol=1))  
    p = p + labs(x = "Sample", y = "Absolute abundance")
    p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0, size = 5))

    ggsave(paste0("plots/aabund-", level, "-by-method.pdf"), plot=p, device="pdf")
}


for (level in colnames(taxmat)){
    top.taxa <- tax_glom(pphyseq, level)
    TopNOTUs <- names(sort(taxa_sums(top.taxa), TRUE)[1:20])
    BottumNOTUs <- names(taxa_sums(top.taxa))[which(!names(taxa_sums(top.taxa)) %in% TopNOTUs)]
    merged_physeq = merge_taxa(top.taxa, BottumNOTUs, 2)

    mdf = psmelt(merged_physeq); names(mdf)[names(mdf) == level] <- "level"
    mdf$OTU[which(is.na(mdf$level))] <- "aaaOther"
    mdf$level[which(is.na(mdf$level))] <- "aaaOther"
    aggr_mdf <- aggregate(Abundance ~ sample + run + method + level, data = mdf, sum)

    labs = aggr_mdf$level; labs[labs=="aaaOther"] <- "Other"
    cols = scales::hue_pal()(length(unique(labs))); cols[unique(labs) == "Other"] <- "#CCCCCC"

    p = ggplot(aggr_mdf, aes_string(x = "method", y = "Abundance", fill = "level"))
    p = p + scale_fill_manual(name = level, labels = unique(labs), values = cols)
    p = p + facet_grid(paste0("", aggr_mdf$run) ~ paste0("", aggr_mdf$sample))
    p = p + geom_bar(stat = "identity", position = "stack",  color = "black", size = 0.1)
    p =	p + guides(fill=guide_legend(ncol=1))  
    p = p + labs(x = "Method", y = "Absolute abundance")
    p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0, size = 5))

    ggsave(paste0("plots/rabund-", level, "-by-sample.pdf"), plot=p, device="pdf")

    p = ggplot(aggr_mdf, aes_string(x = "sample", y = "Abundance", fill = "level"))
    p = p + scale_fill_manual(name = level, labels = unique(labs), values = cols)
    p = p + facet_grid(paste0("", aggr_mdf$run) ~ paste0("", aggr_mdf$method))
    p = p + geom_bar(stat = "identity", position = "stack",  color = "black", size = 0.1)
    p =	p + guides(fill=guide_legend(ncol=1))  
    p = p + labs(x = "Sample", y = "Relative abundance")
    p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0, size = 5))

    ggsave(paste0("plots/rabund-", level, "-by-method.pdf"), plot=p, device="pdf")
}
