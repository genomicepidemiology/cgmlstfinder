# create_tree.R

# Description: R script infers a phylogeny based on an input matrix created by
#              cgMLSTFinder.

library(cluster)
library(ape)

# Get arguments
args <- commandArgs(trailingOnly = TRUE)
matrix_file = args[0]
output = args[1]

df <- read.csv(matrix_file, sep = "\t", row.names = 1, colClasses = "factor")

# plot(hclust(daisy(df, metric = "gower")), hang = -1)

tree <- as.phylo(hclust(daisy(df, metric = "gower")))
write.tree(phy = tree, file = output)
