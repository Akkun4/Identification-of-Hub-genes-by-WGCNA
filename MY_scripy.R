BiocManager::install("BiocVersion")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = 3.13)

library(WGCNA)
library(tidyverse)
library(DESeq2)
BiocManager::install("GEOquery")
library(GEOquery)
library(genefilter)
library(gridExtra)


rm(list = ls())

pancreatic_cancer <- getGEO("GSE28735", GSEMatrix = TRUE) # read the pancreatic_cancer data set from geo

fpancreatic_cancer <- fData(pancreatic_cancer[[1]]) # get the fdata of pancreatic_cancer which contains information on genesymbols, description etc....
genesymbol <- fpancreatic_cancer[,"gene_assignment"] # assign the gene symbol in the fdata to a variable

exprs_pancreatic_cancer <- exprs(pancreatic_cancer[[1]]) # get the gene expression data 
probeid <- row.names(exprs_pancreatic_cancer) # assign the probeids of the genes to a variable

collapserows <- collapseRows(exprs_pancreatic_cancer,genesymbol,probeid)[[2]]# perform collapse rows function to remove duplicate rows of genes
row.names(exprs_pancreatic_cancer)<- genesymbol # rename the probe ids as their corresponding gene symbols
exprs1_pancreatic_cancer <- exprs_pancreatic_cancer[collapserows[,1],] # filter the expression data using the collapsed rows output

gsg <- goodSamplesGenes(t(exprs_pancreatic_cancer)) # check whether the genes expression values are missing in large number of samples or not
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

exprs_pancreatic_cancer <- exprs_pancreatic_cancer[gsg$goodGenes == TRUE,] # filter out the bad genes if any

#OUTLIERS
htree <- hclust(dist(t(exprs_pancreatic_cancer)), method = "average")
plot(htree)


# PCA
prcomp(t(exprs_pancreatic_cancer))

# filter GENES
f1 <- pOverA(0.25, log2(100))
f2 <- function(x) (IQR(x) > 0.5)
ff <- filterfun(f1, f2)
selected <- genefilter(exprs1_pancreatic_cancer, ff)
sum(selected)
esetSub <- exprs1_pancreatic_cancer[selected, ]

#OUTLIERS
htree <- hclust(dist(t(esetSub)), method = "average")
plot(htree)



exp_pancreatic_cancer <- esetSub[,-c(54,56,2)] 
#exp_psoriasis <- exprs1_psoriasis[,-c(54,56,2)]

pheno_pancreatic_cancer <- pData(phenoData(pancreatic_cancer[[1]]))

pheno_pancreatic_cancer <- pheno_pancreatic_cancer[,c(2,34)]
pheno_pancreatic_cancer <- pheno_pancreatic_cancer[-c(54,56,2),]

condition <- pheno_pancreatic_cancer$`tissue:ch1`

pheno_pan_can <- data.frame(condition)
row.names(pheno_pan_can) <- row.names(pheno_pancreatic_cancer)

#SOFT POWER
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(exp_psoriasis1,
                         powerVector = power,
                         networkType = "signed")

