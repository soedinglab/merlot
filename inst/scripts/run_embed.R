#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(merlot))
library("optparse")

option_list <- list(
    make_option(c("-i", "--in"), 
                type = "character",
                dest = "input",
                help = "Input file (elastic tree in RDS format).",
                metavar = "/input/file1"),
    make_option(c("-f", "--full"), 
                type = "character",
                dest = "genes",
                help = "Gene expression matrix on which the tree will be embedded (cells x genes)",
                metavar = "/input/file2"),
    make_option(c("-o", "--out"), 
                type = "character",
                help = "Output file.",
                metavar = "/output/file"),
    make_option(c("-l", "--lambda_0"),
                type = "double",
                dest = "lambda",
                default = 2.03e-09,
                help = "Principal elastic tree energy function parameter. (default: %default)."),
    make_option(c("-m", "--mu_0"),
                type = "double",
                dest = "mu",
                default = 0.00625,
                help = "Principal elastic tree energy function parameter (default: %default)."),
    make_option(c("-n", "--increaseFactor_mu"),
                type = "double",
                dest = "increase_mu",
                default = 20,
                help = "Factor by which mu_0 will be increased for the embedding (default: %default)."),
    make_option(c("-k", "--increaseFactor_lambda"),
                type = "double",
                dest = "increase_lambda",
                default = 20,
                help = "Factor by which lambda_0 will be increased for the embedding (default: %default)."),
    make_option(c("-c", "--NCores"),
                type = "integer",
                dest = "cores",
                default = 1,
                help = "Number of CPU cores to be used for the calculation (default: %default).")
);

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

# print(opt)

tree <- merlot::GenesSpaceEmbedding(readRDS(opt$genes),
                                    readRDS(opt$input),
                                    lambda_0 = opt$lambda,
                                    mu_0 = opt$mu,
                                    increaseFactor_mu = opt$increase_mu,
                                    increaseFactor_lambda = opt$increase_lambda,
                                    NCores = opt$cores)

saveRDS(tree, file = opt$out)
