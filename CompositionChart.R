### Run the parent file first to generate variables: DataPrep.R
### This will generate all of the composition results including 
### the top 1% charts and the 0.1% table

library(ggplot2)
library(microbiome)
library(latticeExtra)
library(gridExtra)
library(MASS)
library(colorspace)
library(plyr)
library(Hmisc)
library(scales)
library(reshape2)
library(data.table)

### Data Upload and Manipulation
qiime.bac.data = read.table("Data/QIIME/bac_01_tissue.csv", header = TRUE, sep = ",")
qiime.fun.data = read.table("Data/QIIME/fun_01_tissue.csv", header = TRUE, sep = ",")
dada2.new.bac.rare.tiss = merge_samples(dada2.new.bac.rare, "Tissue", fun = "sum")
dada2.new.fun.rare.tiss = merge_samples(dada2.new.fun.rare, "Tissue", fun = "sum")
dada2.new.bac.esv = matrix(nrow = nrow(tax_table(dada2.new.bac.01.tiss)), ncol = 4)
dada2.new.fun.esv = matrix(nrow = nrow(tax_table(dada2.new.fun.01.tiss)), ncol = 4)
rownames(dada2.new.bac.esv) = tax_table(dada2.new.bac.01.tiss)[,6]
colnames(dada2.new.bac.esv) = rownames(otu_table(dada2.new.bac.01.tiss))
rownames(dada2.new.fun.esv) = tax_table(dada2.new.fun.01.tiss)[,6]
colnames(dada2.new.fun.esv) = rownames(otu_table(dada2.new.fun.01.tiss))

# Assigning values to the ESV table, very messy
# Some lines vary because not all top taxa went to the Genera level

dada2.new.bac.esv[1,1] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[1]])) & (otu_table(dada2.new.bac.rare.tiss)[1] >= 1), na.rm = TRUE)
dada2.new.bac.esv[2,1] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[2]])) & (otu_table(dada2.new.bac.rare.tiss)[1] >= 1), na.rm = TRUE)
dada2.new.bac.esv[3,1] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[3]])) & (otu_table(dada2.new.bac.rare.tiss)[1] >= 1), na.rm = TRUE)
dada2.new.bac.esv[4,1] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[4]])) & (otu_table(dada2.new.bac.rare.tiss)[1] >= 1), na.rm = TRUE)
dada2.new.bac.esv[5,1] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[5]])) & (otu_table(dada2.new.bac.rare.tiss)[1] >= 1), na.rm = TRUE)
dada2.new.bac.esv[6,1] = sum(t(paste(tax_table(dada2.new.bac.rare)[,5], tax_table(dada2.new.bac.rare)[,6], sep = "") == paste(tax_table(dada2.new.bac.01.tiss)[,5][[6]], tax_table(dada2.new.bac.01.tiss)[,6][[6]], sep = "")) & (otu_table(dada2.new.bac.rare.tiss)[1] >= 1), na.rm = TRUE)
dada2.new.bac.esv[7,1] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[7]])) & (otu_table(dada2.new.bac.rare.tiss)[1] >= 1), na.rm = TRUE)
dada2.new.bac.esv[8,1] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[8]])) & (otu_table(dada2.new.bac.rare.tiss)[1] >= 1), na.rm = TRUE)
dada2.new.bac.esv[9,1] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[9]])) & (otu_table(dada2.new.bac.rare.tiss)[1] >= 1), na.rm = TRUE)

dada2.new.bac.esv[1,2] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[1]])) & (otu_table(dada2.new.bac.rare.tiss)[2] >= 1), na.rm = TRUE)
dada2.new.bac.esv[2,2] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[2]])) & (otu_table(dada2.new.bac.rare.tiss)[2] >= 1), na.rm = TRUE)
dada2.new.bac.esv[3,2] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[3]])) & (otu_table(dada2.new.bac.rare.tiss)[2] >= 1), na.rm = TRUE)
dada2.new.bac.esv[4,2] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[4]])) & (otu_table(dada2.new.bac.rare.tiss)[2] >= 1), na.rm = TRUE)
dada2.new.bac.esv[5,2] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[5]])) & (otu_table(dada2.new.bac.rare.tiss)[2] >= 1), na.rm = TRUE)
dada2.new.bac.esv[6,2] = sum(t(paste(tax_table(dada2.new.bac.rare)[,5], tax_table(dada2.new.bac.rare)[,6], sep = "") == paste(tax_table(dada2.new.bac.01.tiss)[,5][[6]], tax_table(dada2.new.bac.01.tiss)[,6][[6]], sep = "")) & (otu_table(dada2.new.bac.rare.tiss)[2] >= 1), na.rm = TRUE)
dada2.new.bac.esv[7,2] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[7]])) & (otu_table(dada2.new.bac.rare.tiss)[2] >= 1), na.rm = TRUE)
dada2.new.bac.esv[8,2] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[8]])) & (otu_table(dada2.new.bac.rare.tiss)[2] >= 1), na.rm = TRUE)
dada2.new.bac.esv[9,2] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[9]])) & (otu_table(dada2.new.bac.rare.tiss)[2] >= 1), na.rm = TRUE)

dada2.new.bac.esv[1,3] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[1]])) & (otu_table(dada2.new.bac.rare.tiss)[3] >= 1), na.rm = TRUE)
dada2.new.bac.esv[2,3] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[2]])) & (otu_table(dada2.new.bac.rare.tiss)[3] >= 1), na.rm = TRUE)
dada2.new.bac.esv[3,3] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[3]])) & (otu_table(dada2.new.bac.rare.tiss)[3] >= 1), na.rm = TRUE)
dada2.new.bac.esv[4,3] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[4]])) & (otu_table(dada2.new.bac.rare.tiss)[3] >= 1), na.rm = TRUE)
dada2.new.bac.esv[5,3] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[5]])) & (otu_table(dada2.new.bac.rare.tiss)[3] >= 1), na.rm = TRUE)
dada2.new.bac.esv[6,3] = sum(t(paste(tax_table(dada2.new.bac.rare)[,5], tax_table(dada2.new.bac.rare)[,6], sep = "") == paste(tax_table(dada2.new.bac.01.tiss)[,5][[6]], tax_table(dada2.new.bac.01.tiss)[,6][[6]], sep = "")) & (otu_table(dada2.new.bac.rare.tiss)[3] >= 1), na.rm = TRUE)
dada2.new.bac.esv[7,3] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[7]])) & (otu_table(dada2.new.bac.rare.tiss)[3] >= 1), na.rm = TRUE)
dada2.new.bac.esv[8,3] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[8]])) & (otu_table(dada2.new.bac.rare.tiss)[3] >= 1), na.rm = TRUE)
dada2.new.bac.esv[9,3] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[9]])) & (otu_table(dada2.new.bac.rare.tiss)[3] >= 1), na.rm = TRUE)

dada2.new.bac.esv[1,4] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[1]])) & (otu_table(dada2.new.bac.rare.tiss)[4] >= 1), na.rm = TRUE)
dada2.new.bac.esv[2,4] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[2]])) & (otu_table(dada2.new.bac.rare.tiss)[4] >= 1), na.rm = TRUE)
dada2.new.bac.esv[3,4] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[3]])) & (otu_table(dada2.new.bac.rare.tiss)[4] >= 1), na.rm = TRUE)
dada2.new.bac.esv[4,4] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[4]])) & (otu_table(dada2.new.bac.rare.tiss)[4] >= 1), na.rm = TRUE)
dada2.new.bac.esv[5,4] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[5]])) & (otu_table(dada2.new.bac.rare.tiss)[4] >= 1), na.rm = TRUE)
dada2.new.bac.esv[6,4] = sum(t(paste(tax_table(dada2.new.bac.rare)[,5], tax_table(dada2.new.bac.rare)[,6], sep = "") == paste(tax_table(dada2.new.bac.01.tiss)[,5][[6]], tax_table(dada2.new.bac.01.tiss)[,6][[6]], sep = "")) & (otu_table(dada2.new.bac.rare.tiss)[4] >= 1), na.rm = TRUE)
dada2.new.bac.esv[7,4] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[7]])) & (otu_table(dada2.new.bac.rare.tiss)[4] >= 1), na.rm = TRUE)
dada2.new.bac.esv[8,4] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[8]])) & (otu_table(dada2.new.bac.rare.tiss)[4] >= 1), na.rm = TRUE)
dada2.new.bac.esv[9,4] = sum(t((tax_table(dada2.new.bac.rare.tiss)[,6] == tax_table(dada2.new.bac.01.tiss)[,6][[9]])) & (otu_table(dada2.new.bac.rare.tiss)[4] >= 1), na.rm = TRUE)

# Convert count data to relative values
dada2.new.bac.esv[,1] = dada2.new.bac.esv[,1] / colSums(dada2.new.bac.esv)[1]
dada2.new.bac.esv[,2] = dada2.new.bac.esv[,2] / colSums(dada2.new.bac.esv)[2]
dada2.new.bac.esv[,3] = dada2.new.bac.esv[,3] / colSums(dada2.new.bac.esv)[3]
dada2.new.bac.esv[,4] = dada2.new.bac.esv[,4] / colSums(dada2.new.bac.esv)[4]

dada2.new.fun.esv[1,1] = sum(paste(tax_table(dada2.new.fun.rare)[,4], tax_table(dada2.new.fun.rare)[,5], sep = "") == paste(tax_table(dada2.new.fun.rare.tiss)[,4][[1]], tax_table(dada2.new.fun.rare.tiss)[,5][[1]], sep = "") & (otu_table(dada2.new.fun.rare.tiss)[1] >= 1), na.rm = TRUE)
dada2.new.fun.esv[2,1] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[2]]) & (otu_table(dada2.new.fun.rare.tiss)[1] >= 1), na.rm = TRUE)
dada2.new.fun.esv[3,1] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[3]]) & (otu_table(dada2.new.fun.rare.tiss)[1] >= 1), na.rm = TRUE)
dada2.new.fun.esv[4,1] = sum(paste(tax_table(dada2.new.fun.rare)[,5], tax_table(dada2.new.fun.rare)[,6], sep = "") == paste(tax_table(dada2.new.fun.rare.tiss)[,5][[4]], tax_table(dada2.new.fun.rare.tiss)[,6][[4]], sep = "") & (otu_table(dada2.new.fun.rare.tiss)[1] >= 1), na.rm = TRUE)
dada2.new.fun.esv[5,1] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[5]]) & (otu_table(dada2.new.fun.rare.tiss)[1] >= 1), na.rm = TRUE)
dada2.new.fun.esv[6,1] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[6]]) & (otu_table(dada2.new.fun.rare.tiss)[1] >= 1), na.rm = TRUE)
dada2.new.fun.esv[7,1] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[7]]) & (otu_table(dada2.new.fun.rare.tiss)[1] >= 1), na.rm = TRUE)
dada2.new.fun.esv[8,1] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[8]]) & (otu_table(dada2.new.fun.rare.tiss)[1] >= 1), na.rm = TRUE)
dada2.new.fun.esv[9,1] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[9]]) & (otu_table(dada2.new.fun.rare.tiss)[1] >= 1), na.rm = TRUE)
dada2.new.fun.esv[10,1] = sum(paste(tax_table(dada2.new.fun.rare)[,4], tax_table(dada2.new.fun.rare)[,5], sep = "") == paste(tax_table(dada2.new.fun.rare.tiss)[,4][[10]], tax_table(dada2.new.fun.rare.tiss)[,5][[10]], sep = "") & (otu_table(dada2.new.fun.rare.tiss)[1] >= 1), na.rm = TRUE)
dada2.new.fun.esv[11,1] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[11]]) & (otu_table(dada2.new.fun.rare.tiss)[1] >= 1), na.rm = TRUE)
dada2.new.fun.esv[12,1] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[12]]) & (otu_table(dada2.new.fun.rare.tiss)[1] >= 1), na.rm = TRUE)
dada2.new.fun.esv[13,1] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[13]]) & (otu_table(dada2.new.fun.rare.tiss)[1] >= 1), na.rm = TRUE)

dada2.new.fun.esv[1,2] = sum(paste(tax_table(dada2.new.fun.rare)[,4], tax_table(dada2.new.fun.rare)[,5], sep = "") == paste(tax_table(dada2.new.fun.rare.tiss)[,4][[1]], tax_table(dada2.new.fun.rare.tiss)[,5][[1]], sep = "") & (otu_table(dada2.new.fun.rare.tiss)[2] >= 1), na.rm = TRUE)
dada2.new.fun.esv[2,2] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[2]]) & (otu_table(dada2.new.fun.rare.tiss)[2] >= 1), na.rm = TRUE)
dada2.new.fun.esv[3,2] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[3]]) & (otu_table(dada2.new.fun.rare.tiss)[2] >= 1), na.rm = TRUE)
dada2.new.fun.esv[4,2] = sum(paste(tax_table(dada2.new.fun.rare)[,5], tax_table(dada2.new.fun.rare)[,6], sep = "") == paste(tax_table(dada2.new.fun.rare.tiss)[,5][[4]], tax_table(dada2.new.fun.rare.tiss)[,6][[4]], sep = "") & (otu_table(dada2.new.fun.rare.tiss)[2] >= 1), na.rm = TRUE)
dada2.new.fun.esv[5,2] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[5]]) & (otu_table(dada2.new.fun.rare.tiss)[2] >= 1), na.rm = TRUE)
dada2.new.fun.esv[6,2] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[6]]) & (otu_table(dada2.new.fun.rare.tiss)[2] >= 1), na.rm = TRUE)
dada2.new.fun.esv[7,2] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[7]]) & (otu_table(dada2.new.fun.rare.tiss)[2] >= 1), na.rm = TRUE)
dada2.new.fun.esv[8,2] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[8]]) & (otu_table(dada2.new.fun.rare.tiss)[2] >= 1), na.rm = TRUE)
dada2.new.fun.esv[9,2] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[9]]) & (otu_table(dada2.new.fun.rare.tiss)[2] >= 1), na.rm = TRUE)
dada2.new.fun.esv[10,2] = sum(paste(tax_table(dada2.new.fun.rare)[,4], tax_table(dada2.new.fun.rare)[,5], sep = "") == paste(tax_table(dada2.new.fun.rare.tiss)[,4][[10]], tax_table(dada2.new.fun.rare.tiss)[,5][[10]], sep = "") & (otu_table(dada2.new.fun.rare.tiss)[2] >= 1), na.rm = TRUE)
dada2.new.fun.esv[11,2] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[11]]) & (otu_table(dada2.new.fun.rare.tiss)[2] >= 1), na.rm = TRUE)
dada2.new.fun.esv[12,2] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[12]]) & (otu_table(dada2.new.fun.rare.tiss)[2] >= 1), na.rm = TRUE)
dada2.new.fun.esv[13,2] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[13]]) & (otu_table(dada2.new.fun.rare.tiss)[2] >= 1), na.rm = TRUE)

dada2.new.fun.esv[1,3] = sum(paste(tax_table(dada2.new.fun.rare)[,4], tax_table(dada2.new.fun.rare)[,5], sep = "") == paste(tax_table(dada2.new.fun.rare.tiss)[,4][[1]], tax_table(dada2.new.fun.rare.tiss)[,5][[1]], sep = "") & (otu_table(dada2.new.fun.rare.tiss)[3] >= 1), na.rm = TRUE)
dada2.new.fun.esv[2,3] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[2]]) & (otu_table(dada2.new.fun.rare.tiss)[3] >= 1), na.rm = TRUE)
dada2.new.fun.esv[3,3] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[3]]) & (otu_table(dada2.new.fun.rare.tiss)[3] >= 1), na.rm = TRUE)
dada2.new.fun.esv[4,3] = sum(paste(tax_table(dada2.new.fun.rare)[,5], tax_table(dada2.new.fun.rare)[,6], sep = "") == paste(tax_table(dada2.new.fun.rare.tiss)[,5][[4]], tax_table(dada2.new.fun.rare.tiss)[,6][[4]], sep = "") & (otu_table(dada2.new.fun.rare.tiss)[3] >= 1), na.rm = TRUE)
dada2.new.fun.esv[5,3] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[5]]) & (otu_table(dada2.new.fun.rare.tiss)[3] >= 1), na.rm = TRUE)
dada2.new.fun.esv[6,3] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[6]]) & (otu_table(dada2.new.fun.rare.tiss)[3] >= 1), na.rm = TRUE)
dada2.new.fun.esv[7,3] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[7]]) & (otu_table(dada2.new.fun.rare.tiss)[3] >= 1), na.rm = TRUE)
dada2.new.fun.esv[8,3] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[8]]) & (otu_table(dada2.new.fun.rare.tiss)[3] >= 1), na.rm = TRUE)
dada2.new.fun.esv[9,3] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[9]]) & (otu_table(dada2.new.fun.rare.tiss)[3] >= 1), na.rm = TRUE)
dada2.new.fun.esv[10,3] = sum(paste(tax_table(dada2.new.fun.rare)[,4], tax_table(dada2.new.fun.rare)[,5], sep = "") == paste(tax_table(dada2.new.fun.rare.tiss)[,4][[10]], tax_table(dada2.new.fun.rare.tiss)[,5][[10]], sep = "") & (otu_table(dada2.new.fun.rare.tiss)[3] >= 1), na.rm = TRUE)
dada2.new.fun.esv[11,3] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[11]]) & (otu_table(dada2.new.fun.rare.tiss)[3] >= 1), na.rm = TRUE)
dada2.new.fun.esv[12,3] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[12]]) & (otu_table(dada2.new.fun.rare.tiss)[3] >= 1), na.rm = TRUE)
dada2.new.fun.esv[13,3] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[13]]) & (otu_table(dada2.new.fun.rare.tiss)[3] >= 1), na.rm = TRUE)

dada2.new.fun.esv[1,4] = sum(paste(tax_table(dada2.new.fun.rare)[,4], tax_table(dada2.new.fun.rare)[,5], sep = "") == paste(tax_table(dada2.new.fun.rare.tiss)[,4][[1]], tax_table(dada2.new.fun.rare.tiss)[,5][[1]], sep = "") & (otu_table(dada2.new.fun.rare.tiss)[4] >= 1), na.rm = TRUE)
dada2.new.fun.esv[2,4] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[2]]) & (otu_table(dada2.new.fun.rare.tiss)[4] >= 1), na.rm = TRUE)
dada2.new.fun.esv[3,4] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[3]]) & (otu_table(dada2.new.fun.rare.tiss)[4] >= 1), na.rm = TRUE)
dada2.new.fun.esv[4,4] = sum(paste(tax_table(dada2.new.fun.rare)[,5], tax_table(dada2.new.fun.rare)[,6], sep = "") == paste(tax_table(dada2.new.fun.rare.tiss)[,5][[4]], tax_table(dada2.new.fun.rare.tiss)[,6][[4]], sep = "") & (otu_table(dada2.new.fun.rare.tiss)[4] >= 1), na.rm = TRUE)
dada2.new.fun.esv[5,4] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[5]]) & (otu_table(dada2.new.fun.rare.tiss)[4] >= 1), na.rm = TRUE)
dada2.new.fun.esv[6,4] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[6]]) & (otu_table(dada2.new.fun.rare.tiss)[4] >= 1), na.rm = TRUE)
dada2.new.fun.esv[7,4] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[7]]) & (otu_table(dada2.new.fun.rare.tiss)[4] >= 1), na.rm = TRUE)
dada2.new.fun.esv[8,4] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[8]]) & (otu_table(dada2.new.fun.rare.tiss)[4] >= 1), na.rm = TRUE)
dada2.new.fun.esv[9,4] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[9]]) & (otu_table(dada2.new.fun.rare.tiss)[4] >= 1), na.rm = TRUE)
dada2.new.fun.esv[10,4] = sum(paste(tax_table(dada2.new.fun.rare)[,4], tax_table(dada2.new.fun.rare)[,5], sep = "") == paste(tax_table(dada2.new.fun.rare.tiss)[,4][[10]], tax_table(dada2.new.fun.rare.tiss)[,5][[10]], sep = "") & (otu_table(dada2.new.fun.rare.tiss)[4] >= 1), na.rm = TRUE)
dada2.new.fun.esv[11,4] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[11]]) & (otu_table(dada2.new.fun.rare.tiss)[4] >= 1), na.rm = TRUE)
dada2.new.fun.esv[12,4] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[12]]) & (otu_table(dada2.new.fun.rare.tiss)[4] >= 1), na.rm = TRUE)
dada2.new.fun.esv[13,4] = sum(t(tax_table(dada2.new.fun.rare)[,6] == tax_table(dada2.new.fun.rare.tiss)[,6][[13]]) & (otu_table(dada2.new.fun.rare.tiss)[4] >= 1), na.rm = TRUE)

dada2.new.fun.esv[,1] = dada2.new.fun.esv[,1] / colSums(dada2.new.fun.esv)[1]
dada2.new.fun.esv[,2] = dada2.new.fun.esv[,2] / colSums(dada2.new.fun.esv)[2]
dada2.new.fun.esv[,3] = dada2.new.fun.esv[,3] / colSums(dada2.new.fun.esv)[3]
dada2.new.fun.esv[,4] = dada2.new.fun.esv[,4] / colSums(dada2.new.fun.esv)[4]

# Chart Builder
### Top 1% Chart
#dada2.old.bac.tiss.chart = abundances(dada2.old.bac.01.tiss, "compositional")[,c("Grain", "Outer Grain", "Husk", "Outer Husk")]
#row.names(dada2.old.bac.tiss.chart) = tax_table(dada2.old.bac.01.tiss)[,"Genus"]
#row.names(dada2.old.bac.tiss.chart)[1] = "Microbacteriaceae"
#dada2.old.bac.tiss.chart = dada2.old.bac.tiss.chart[match(c("Microbacteriaceae", "Kineococcus", "Microbacterium", "Aureimonas", "Methylobacterium", "Rhizobium", "Sphingomonas", "Pantoea", "Pluralibacter"), rownames(dada2.old.bac.tiss.chart)),]
#dada2.old.bac.tiss.chart = melt(dada2.old.bac.tiss.chart)
#dada2.old.bac.ggplot = ggplot(dada2.old.bac.tiss.chart, aes(x = Var2, y = value, fill = Var1))
#dada2.old.bac.ggplot + geom_col(colour = "black", width = 0.6) + scale_fill_manual(values = c(sequential_hcl(3, h = 10, c = 50, power = 1.5), sequential_hcl(4, h = 130, c = 50, power = 1.5), sequential_hcl(2, h = 250, c = 50, power = 1.5))) + labs(x = "Tissue", y = "Relative Abundance", fill = "Taxa") + theme(panel.background = element_blank())

#dada2.old.fun.tiss.chart = abundances(dada2.old.fun.01.tiss, "compositional")[,c("Grain", "Outer Grain", "Husk", "Outer Husk")]
#row.names(dada2.old.fun.tiss.chart) = tax_table(dada2.old.fun.01.tiss)[,"Genus"]
#dada2.old.fun.tiss.chart = dada2.old.fun.tiss.chart[match(c("Occultifur", "Bulleribasidium", "Cryptococcus", "Hannaella", "Papiliotrema", "Alternaria", "Cladosporium", "Didymella", "Exserohilum","Phaeosphaeria", "Sphaerulina", "Gibberella", "Nigrospora"), rownames(dada2.old.fun.tiss.chart)),]
#dada2.old.fun.tiss.chart = melt(dada2.old.fun.tiss.chart)
#dada2.old.fun.ggplot = ggplot(dada2.old.fun.tiss.chart, aes(x = Var2, y = value, fill = Var1))
#dada2.old.fun.ggplot + geom_col(colour = "black", width = 0.6) + scale_fill_manual(values = c(sequential_hcl(5, h = 130, c = 50, power = 1.5), sequential_hcl(8, h = 250, c = 50, power = 1.5))) + labs(x = "Tissue", y = "Relative Abundance", fill = "Taxa") + theme(panel.background = element_blank())

### ## New Bacteria

# Upload Data, reorganize
dada2.new.bac.tiss.chart = abundances(dada2.new.bac.01.tiss, "compositional")[,c("Grain", "Outer Grain", "Husk", "Outer Husk")]
dada2.new.bac.esv = dada2.new.bac.esv[,c("Grain", "Outer Grain", "Husk", "Outer Husk")]
# Set rownames as some genera names were messy
row.names(dada2.new.bac.tiss.chart) = tax_table(dada2.new.bac.01.tiss)[,"Genus"]
rownames(dada2.new.bac.tiss.chart) = c(rownames(dada2.new.bac.tiss.chart)[1], "Rhizobium", rownames(dada2.new.bac.tiss.chart)[c(3:5)], "Enterobacteriaceae", rownames(dada2.new.bac.tiss.chart)[c(7:9)])
rownames(dada2.new.bac.esv) = rownames(dada2.new.bac.tiss.chart)
# Sort both datasets
dada2.new.bac.tiss.chart = dada2.new.bac.tiss.chart[match(c("Curtobacterium", "Microbacterium", "Aureimonas", "Methylobacterium", "Rhizobium", "Sphingomonas", "Enterobacteriaceae", "Atlantibacter", "Franconibacter"), rownames(dada2.new.bac.tiss.chart)),]
dada2.new.bac.esv = dada2.new.bac.esv[match(c("Curtobacterium", "Microbacterium", "Aureimonas", "Methylobacterium", "Rhizobium", "Sphingomonas", "Enterobacteriaceae", "Atlantibacter", "Franconibacter"), rownames(dada2.new.bac.tiss.chart)),]
# "Melt" the data into something ggplot2 can use
dada2.new.bac.tiss.chart = melt(dada2.new.bac.tiss.chart)
dada2.new.bac.esv = melt(dada2.new.bac.esv)
# Construction of the plots, sequential_hcl was used to generate colors
dada2.new.bac.ggplot = ggplot(dada2.new.bac.tiss.chart, aes(x = Var2, y = value, fill = Var1))
dada2.new.bac.ggplot + geom_col(colour = "black", width = 0.6) + scale_fill_manual(values = c(sequential_hcl(3, h = 10, c = 50, power = 1.5), sequential_hcl(4, h = 130, c = 50, power = 1.5), sequential_hcl(2, h = 250, c = 50, power = 1.5))) + labs(x = "Tissue", y = "Relative Abundance", fill = "Taxa") + theme(panel.background = element_blank())
dada2.new.bac.esv.ggplot = ggplot(dada2.new.bac.esv, aes(x = Var2, y = value, fill = Var1))
dada2.new.bac.esv.ggplot + geom_col(colour = "black", width = 0.6) + scale_fill_manual(values = c(sequential_hcl(2, h = 10, c = 50, power = 1.5), sequential_hcl(4, h = 130, c = 50, power = 1.5), sequential_hcl(3, h = 250, c = 50, power = 1.5))) + labs(x = "Tissue", y = "Relative Abundance", fill = "Taxa") + theme(panel.background = element_blank())

tiff(filename = "BacComp.tif", units = "px", width = 2000, height = 1188, res = 360)
dada2.new.bac.final.ggplot = grid.arrange(dada2.new.bac.ggplot + theme(legend.position="none") + geom_col(colour = "black", width = 0.6) + scale_fill_manual(values = c(sequential_hcl(3, h = 130, c = 65, power = 1)[1:2], sequential_hcl(3, h = 10, c = 65, power = 1)[1:2], sequential_hcl(3, h = 45, c = 65, power = 1)[1:2], sequential_hcl(4, h = 250, c = 65, power = 1)[1:3])) + labs(x = "Read Abundance", y = "Relative Abundance") + theme(panel.background = element_blank()), dada2.new.bac.esv.ggplot + theme(legend.position="none") + geom_col(colour = "black", width = 0.6) + scale_fill_manual(values = c(sequential_hcl(3, h = 130, c = 65, power = 1)[1:2], sequential_hcl(3, h = 10, c = 65, power = 1)[1:2], sequential_hcl(3, h = 45, c = 65, power = 1)[1:2], sequential_hcl(4, h = 250, c = 65, power = 1)[1:3])) + labs(x = "ASV Abundance", y = " ") + theme(panel.background = element_blank()), nrow = 1)
dev.off()

### ## New Fungi

dada2.new.fun.tiss.chart = abundances(dada2.new.fun.01.tiss, "compositional")[,c("Grain", "Outer Grain", "Husk", "Outer Husk")]
dada2.new.fun.esv = dada2.new.fun.esv[,c("Grain", "Outer Grain", "Husk", "Outer Husk")]
row.names(dada2.new.fun.tiss.chart) = tax_table(dada2.new.fun.01.tiss)[,"Genus"]
row.names(dada2.new.fun.tiss.chart) = c("Pleosporales", rownames(dada2.new.fun.tiss.chart)[c(2:3)], "Didymellaceae", rownames(dada2.new.fun.tiss.chart)[c(5:9)], "Tremellales", rownames(dada2.new.fun.tiss.chart)[c(11:13)])
rownames(dada2.new.fun.esv) = rownames(dada2.new.fun.tiss.chart)
dada2.new.fun.tiss.chart = dada2.new.fun.tiss.chart[match(c("Pleosporales", "Cladosporium", "Didymellaceae", "Cercospora", "Alternaria", "Curvularia", "Gibberella", "Nigrospora", "Tremellales", "Bullera", "Hannaella", "Papiliotrema", "Saitozyma"), rownames(dada2.new.fun.tiss.chart)),]
dada2.new.fun.esv = dada2.new.fun.esv[match(c("Pleosporales", "Cladosporium", "Didymellaceae", "Cercospora", "Alternaria", "Curvularia", "Gibberella", "Nigrospora", "Tremellales", "Bullera", "Hannaella", "Papiliotrema", "Saitozyma"), rownames(dada2.new.fun.esv)),]
dada2.new.fun.tiss.chart = melt(dada2.new.fun.tiss.chart)
dada2.new.fun.esv = melt(dada2.new.fun.esv)
dada2.new.fun.ggplot = ggplot(dada2.new.fun.tiss.chart, aes(x = Var2, y = value, fill = Var1))
dada2.new.fun.ggplot + geom_col(colour = "black", width = 0.6) + scale_fill_manual(values = c(sequential_hcl(3, h = 10, c = 50, power = 1), sequential_hcl(3, h = 45, c = 50, power = 1), sequential_hcl(2, h = 130, c = 50, power = 1), sequential_hcl(1, h = 240, c = 50, power = 1), sequential_hcl(4, h = 280, c = 50, power = 1))) + labs(x = "Tissue", y = "Relative Abundance", fill = "Taxa") + theme(panel.background = element_blank())
dada2.new.fun.esv.ggplot = ggplot(dada2.new.fun.esv, aes(x = Var2, y = value, fill = Var1))
dada2.new.fun.esv.ggplot + geom_col(colour = "black", width = 0.6) + scale_fill_manual(values = c(sequential_hcl(3, h = 10, c = 50, power = 1), sequential_hcl(3, h = 45, c = 50, power = 1), sequential_hcl(2, h = 130, c = 50, power = 1), sequential_hcl(5, h = 250, c = 50, power = 1))) + labs(x = "Tissue", y = "Relative Abundance", fill = "Taxa") + theme(panel.background = element_blank())

tiff(filename = "FunComp.tif", units = "px", width = 2000, height = 1188, res = 360)
dada2.new.fun.final.ggplot = grid.arrange(dada2.new.fun.ggplot + theme(legend.position="none") + geom_col(colour = "black", width = 0.6) + scale_fill_manual(values = c(sequential_hcl(4, h = 10, c = 65, power = 1)[1:3], sequential_hcl(4, h = 45, c = 65, power = 1)[1:3], sequential_hcl(3, h = 130, c = 65, power = 1)[1:2], sequential_hcl(1, h = 240, c = 65, power = 1), sequential_hcl(5, h = 280, c = 50, power = 1)[1:4])) + labs(x = "Read Abundance", y = "Relative Abundance") + theme(panel.background = element_blank()), dada2.new.fun.esv.ggplot + theme(legend.position="none") + geom_col(colour = "black", width = 0.6) + scale_fill_manual(values = c(sequential_hcl(4, h = 10, c = 65, power = 1)[1:3], sequential_hcl(4, h = 45, c = 65, power = 1)[1:3], sequential_hcl(3, h = 130, c = 65, power = 1)[1:2], sequential_hcl(1, h = 240, c = 65, power = 1), sequential_hcl(5, h = 280, c = 65, power = 1)[1:4])) + labs(x = "ASV Abundance", y = " ") + theme(panel.background = element_blank()), nrow = 1)
dev.off()

###

colnames(qiime.bac.data) = c("Taxa", "Grain", "Outer Grain", "Husk", "Outer Husk")
qiime.bac.data$Taxa = factor(qiime.bac.data$Taxa, levels = qiime.bac.data$Taxa)
qiime.old.bac.tiss.chart = melt(qiime.bac.data)
qiime.old.bac.ggplot = ggplot(qiime.old.bac.tiss.chart, aes(x = variable, y = value, fill = Taxa))
qiime.old.bac.ggplot + geom_col(colour = "black", width = 0.6) + scale_fill_manual(values = c(sequential_hcl(1, h = 10, c = 50, power = 1.5), sequential_hcl(4, h = 130, c = 50, power = 1.5), sequential_hcl(5, h = 250, c = 50, power = 1.5))) + labs(x = "Tissue", y = "Relative Abundance", fill = "Taxa") + theme(panel.background = element_blank())

colnames(qiime.fun.data) = c("Taxa", "Grain", "Outer Grain", "Husk", "Outer Husk")
qiime.fun.data = qiime.fun.data[c(13,9,10,11,12,8,7,6,5,4,3,1,2),]
qiime.fun.data$Taxa = factor(qiime.fun.data$Taxa, levels = qiime.fun.data$Taxa)
qiime.old.fun.tiss.chart = melt(qiime.fun.data)
qiime.old.fun.ggplot = ggplot(qiime.old.fun.tiss.chart, aes(x = variable, y = value, fill = Taxa))
qiime.old.fun.ggplot + geom_col(colour = "black", width = 0.6) + scale_fill_manual(values = c(sequential_hcl(1, h = 10, c = 50, power = 1.5), sequential_hcl(5, h = 130, c = 50, power = 1.5), sequential_hcl(7, h = 250, c = 50, power = 1.5))) + labs(x = "Tissue", y = "Relative Abundance", fill = "Taxa") + theme(panel.background = element_blank())



#####
dada2.new.bac.01.total = prune_taxa(taxa_sums(dada2.new.bac.total) / sum(taxa_sums(dada2.new.bac.total)) >= 0.01, dada2.new.bac.total)
dada2.new.bac.chart = abundances(dada2.new.bac.01.total, "compositional")[,c("S16", "S20", "S2", "S6", "S10", "S14", "S18", "S22", "S3", "S7", "S11", "S15", "S19", "S23", "S1", "S5", "S9", "S13", "S17", "S21")]
row.names(dada2.new.bac.chart) = tax_table(dada2.new.bac.01.total)[,"Genus"]
row.names(dada2.new.bac.chart) = c(rownames(dada2.new.bac.chart)[1], "Rhizobium", rownames(dada2.new.bac.chart)[c(3:5)], "Enterobacteriaceae", rownames(dada2.new.bac.chart)[c(7:9)])
dada2.new.bac.chart = dada2.new.bac.chart[match(c("Curtobacterium", "Microbacterium", "Aureimonas", "Methylobacterium", "Rhizobium", "Sphingomonas", "Enterobacteriaceae", "Atlantibacter", "Franconibacter"), rownames(dada2.new.bac.chart)),]
colnames(dada2.new.bac.chart) = c("D-Gr", "E-Gr", "A-OGr", "B-OGr", "C-OGr", "D-OGr", "E-OGr", "F-OGr", "A-H", "B-H", "C-H", "D-H", "E-H", "F-H", "A-OH", "B-OH", "C-OH", "D-OH", "E-OH", "F-OH")
dada2.new.bac.chart = melt(dada2.new.bac.chart)
dada2.new.bac2.ggplot = ggplot(dada2.new.bac.chart, aes(x = Var2, y = value, fill = Var1))
dada2.new.bac2.ggplot + geom_col(colour = "black", width = 0.6) + scale_fill_manual(values = c(sequential_hcl(3, h = 10, c = 50, power = 1.5), sequential_hcl(4, h = 130, c = 50, power = 1.5), sequential_hcl(2, h = 250, c = 50, power = 1.5))) + labs(x = "Samples", y = "Relative Abundance", fill = "Taxa") + theme(panel.background = element_blank())

dada2.new.fun.01.total = prune_taxa(taxa_sums(dada2.new.fun.total) / sum(taxa_sums(dada2.new.fun.total)) >= 0.01, dada2.new.fun.total)
dada2.new.fun.chart = abundances(dada2.new.fun.01.total, "compositional")[,c("S4", "S8", "S12", "S16", "S20", "S24", "S2", "S6", "S10", "S14", "S18", "S22", "S3", "S7", "S11", "S15", "S19", "S23", "S1", "S5", "S9", "S13", "S17", "S21")]
row.names(dada2.new.fun.chart) = tax_table(dada2.new.fun.01.total)[,"Genus"]
row.names(dada2.new.fun.chart) = c("Pleosporales", rownames(dada2.new.fun.chart)[c(2:3)], "Didymellaceae", rownames(dada2.new.fun.chart)[c(5:9)], "Tremellales", rownames(dada2.new.fun.chart)[c(11:13)])
dada2.new.fun.chart = dada2.new.fun.chart[match(c("Pleosporales", "Cladosporium", "Didymellaceae", "Cercospora", "Alternaria", "Curvularia", "Gibberella", "Nigrospora", "Tremellales", "Bullera", "Hannaella", "Papiliotrema", "Saitozyma"), rownames(dada2.new.fun.chart)),]
colnames(dada2.new.fun.chart) = c("A-Gr", "B-Gr", "C-Gr", "D-Gr", "E-Gr", "F-Gr", "A-OGr", "B-OGr", "C-OGr", "D-OGr", "E-OGr", "F-OGr", "A-H", "B-H", "C-H", "D-H", "E-H", "F-H", "A-OH", "B-OH", "C-OH", "D-OH", "E-OH", "F-OH")
dada2.new.fun.chart = melt(dada2.new.fun.chart)
dada2.new.fun2.ggplot = ggplot(dada2.new.fun.chart, aes(x = Var2, y = value, fill = Var1))
dada2.new.fun2.ggplot + geom_col(colour = "black", width = 0.6) + scale_fill_manual(values = c(sequential_hcl(3, h = 10, c = 50, power = 1), sequential_hcl(3, h = 45, c = 50, power = 1), sequential_hcl(2, h = 130, c = 50, power = 1), sequential_hcl(1, h = 240, c = 50, power = 1), sequential_hcl(4, h = 280, c = 50, power = 1))) + labs(x = "Samples", y = "Relative Abundance", fill = "Taxa") + theme(panel.background = element_blank())
``