#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

args = commandArgs(trailingOnly=TRUE)

dl = list()
for (i in 1:length(args)){
  # Get metadata
  lab = strsplit(args[i], "/")[[1]][3]
  run = strsplit(args[i], "/")[[1]][2]
  sample = strsplit(strsplit(lab, "[_]")[[1]][3], "[.]")[[1]][1]
  method = strsplit(lab, "[_]")[[1]][1]
  # Read data
  dt <- as.data.frame(read.table(args[i], header = T, sep = '\t',
                    comment = "")[,1])
  # Add meta data columns
  dt$run <- rep(run, nrow(dt))
  dt$sample <- rep(sample, nrow(dt))
  dt$method <- rep(method, nrow(dt))
  # Add to list
  dl[[i]] <- dt
}

df <- do.call(rbind, dl)
names(df)[1] <- "s"

theme_set(theme_bw())

p1 = ggplot(df, aes(x = method, y = s, fill=method)) +
    geom_bar(stat="identity", color="black") +
    labs(x = "", y = "Runtime (s)") +
    theme(axis.text.x = element_text(
        angle = 90, vjust = 0.5, hjust=1))
p1 = p1 + facet_grid(paste0("Run: ",run) ~
                   paste0("Sample: ",sample))
p1 = p1 + theme(legend.position="none")
ggsave("plots/runtime-by-sample.pdf", plot=p1, device="pdf")

q1 = ggplot(df, aes(x = method, y = s, fill=method)) +
    geom_bar(stat="identity", color="black") +
    scale_y_continuous(trans = "log10") +
    labs(x = "", y = "Runtime (s)") +
    theme(axis.text.x = element_text(
        angle = 90, vjust = 0.5, hjust=1))
q1 = q1 + facet_grid(paste0("Run: ",run) ~
                   paste0("Sample: ",sample))
q1 = q1 + theme(legend.position="none")
ggsave("plots/runtime_log-by-sample.pdf", plot=q1, device="pdf")


p2 = ggplot(df, aes(x = sample, y = s, fill=sample)) +
    geom_bar(stat="identity", color="black") +
    labs(x = "", y = "Runtime (s)") +
    theme(axis.text.x = element_text(
        angle = 90, vjust = 0.5, hjust=1))
p2 = p2 + facet_grid(paste0("Run: ",run) ~
                   paste0("Method: ",method))
p2 = p2 + theme(legend.position="none")
ggsave("plots/runtime-by-method.pdf", plot=p2, device="pdf")

q2 = ggplot(df, aes(x = sample, y = s, fill=sample)) +
    geom_bar(stat="identity", color="black") +
    scale_y_continuous(trans = "log10") +
    labs(x = "", y = "Runtime (s)") +
    theme(axis.text.x = element_text(
        angle = 90, vjust = 0.5, hjust=1))
q2 = q2 + facet_grid(paste0("Run: ",run) ~
                   paste0("Method: ",method))
q2 = q2 + theme(legend.position="none")
ggsave("plots/runtime_log-by-method.pdf", plot=q2, device="pdf")
