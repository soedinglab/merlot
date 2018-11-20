#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(merlot))
library("optparse")

option_list <- list(
    make_option(c("-i", "--in"), 
                type = "character",
                dest = "input",
                help = "The elastic tree calculated by the CalculateElasticTree function with a reduced input from ScaffoldTree (RDS object)",
                metavar = "/input/file1"),
    make_option(c("-f", "--full"), 
                type = "character",
                dest = "manifold",
                help = "The manifold coordinates of the cells (cells x manifold components) (RDS object)",
                metavar = "/input/file2"),
    make_option(c("-o", "--out"), 
                type = "character",
                help = "Output file.",
                metavar = "/output/file")
);

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

# print(opt)

tree <- merlot::inflate_elastic_tree(readRDS(opt$input),
                                     readRDS(opt$manifold))

saveRDS(tree, file = opt$out)
