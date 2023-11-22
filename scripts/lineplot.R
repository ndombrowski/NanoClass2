#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

args = commandArgs(trailingOnly=TRUE)

dl = list()
for (i in 1:length(args)){
  # Get metadata
  lab = strsplit(args[i], "/")[[1]][4]
  run = strsplit(args[i], "/")[[1]][2]
  sample = strsplit(lab, "[.]")[[1]][1]
  method = strsplit(lab, "[.]")[[1]][2]
  # Read data and compute total prop precision
  tbl <- read.table(args[i], header = T, sep = '\t', 
                    comment = "", row.names = 1)
  dt <- as.data.frame(colSums(tbl, na.rm=T)/
                      colSums(!is.na(tbl)))
  names(dt) <- "precision"
  # Add meta data columns
  dt$level <- row.names(dt)
  dt$rank <- factor(dt$level, ordered=T, 
      levels=c("Domain","Phylum","Class",
               "Order","Family","Genus"))
  dt$run <- rep(run, nrow(dt))
  dt$sample <- rep(sample, nrow(dt))
  dt$method <- rep(method, nrow(dt))
  # Add to list
  dl[[i]] <- dt
}

df <- do.call(rbind, dl)

theme_set(theme_bw())
p = ggplot(df, aes(x = rank, y = precision, 
           group = method, col = method)) + 
    geom_line() + geom_point() +
    labs(x = "", y = "Precision") + ylim(0,1) +
    theme(axis.text.x = element_text(
        angle = 90, vjust = 0.5, hjust=1))
p = p + facet_grid(paste0("Run: ",run) ~ 
                   paste0("",sample))
ggsave("plots/precision.pdf", plot=p, device="pdf") 


