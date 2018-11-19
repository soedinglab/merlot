#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(merlot))
library("optparse")

option_list <- list(
    make_option(c("-i", "--in"), 
                type = "character",
                dest = "input",
                help = "Input file (significant dimensions on reduced manifold).",
                metavar = "/input/file"),
    make_option(c("-o", "--out"), 
                type = "character",
                help = "Output file.",
                metavar = "/output/file"),
    make_option(c("-n", "--NEndpoints"),
                type = "integer", 
                dest = "endpoints",
                default = -1,
                help = "Users can specify how many endpoints they want the algorithm to find. In case this variable is not defined all branches producing branches longer than sqrt(N/2) will be added to the tree structure."),
    make_option(c("-m", "--BranchMinLength"),
                type = "integer", 
                dest = "minlength",
                default = -1,
                help = "Minimum number of nodes a branch has to contain in order to be included in the tree structure. By default this value is set to sqrt(N/2) with N being the total number of cells in the dataset."),
    make_option(c("-s", "--BranchMinLengthSensitive"),
                type = "integer",
                dest = "minlength_sens",
                help = "Minimum length for a branch to be included in the tree. It reconstructs the topology of the tree and maps cells to the potential new branch to decide if the branch will be added or not. Suggested value: sqrt(N) with N being the number of cells in the dataset.",
                default = -1),
    make_option(c("-r", "--reduced"),
                type = "integer", 
                dest = "reduced",
                default = 0,
                help = "The number of clusters to group cells in. If set to 0, no clustering will be performed and the scaffold tree will be calculated on all cells (default).")
);

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

NEndpoints <- opt$endpoints
if (opt$endpoints == -1) NEndpoints = NULL

coordinates <- readRDS(opt$input)

scaffold <- merlot::CalculateScaffoldTree(coordinates,
                                          NEndpoints = NEndpoints,
                                          BranchMinLength = opt$minlength,
                                          BranchMinLengthSensitive = opt$minlength_sens,
                                          reduced = opt$reduced,
                                          python_location = "python3")

saveRDS(object = scaffold, file = opt$out)
