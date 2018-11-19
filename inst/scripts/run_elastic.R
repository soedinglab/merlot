#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(merlot))
library("optparse")

option_list <- list(
    make_option(c("-i", "--in"), 
                type = "character",
                dest = "input",
                help = "Input file (scaffold tree in RDS format).",
                metavar = "/input/file"),
    make_option(c("-o", "--out"), 
                type = "character",
                help = "Output file.",
                metavar = "/output/file"),
    make_option(c("-n", "--N_yk"),
                type = "integer",
                dest = "nodes",
                default = 100,
                help = "Number of nodes for the elastic principal tree (default: %default)."),
    make_option(c("-l", "--lambda_0"),
                type = "double",
                dest = "lambda",
                default = 8e-10,
                help = "Principal elastic tree energy function parameter. (default: %default)."),
    make_option(c("-m", "--mu_0"),
                type = "double",
                dest = "mu",
                default = 0.0025,
                help = "Principal elastic tree energy function parameter (default: %default)."),
    make_option(c("-f", "--FixEndpoints"),
                type = "logical",
                dest = "fixed",
                default = FALSE,
                action = "store_true",
                help = "Whether or not the end points coordinates are fixed (default: %default)."),
    make_option("--NoBranchScaffoldNodes",
                type = "logical",
                dest = "scafnodes",
                default = TRUE,
                action = "store_false",
                help = "Whether or not to add middle scaffold nodes (default: %default)."),
    make_option(c("-c", "--NCores"),
                type = "integer",
                dest = "cores",
                default = 1,
                help = "Number of CPU cores to be used for the calculation (default: %default).")

);

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

# print(opt)

tree <- merlot::CalculateElasticTree(readRDS(opt$input),
                                     N_yk = opt$nodes,
                                     lambda_0 = opt$lambda,
                                     mu_0 = opt$mu,
                                     FixEndpoints = opt$fixed,
                                     NBranchScaffoldNodes = opt$scafnodes,
                                     NCores = opt$cores,
                                     plot = FALSE)

saveRDS(tree, file = opt$out)
