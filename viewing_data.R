setwd("C:\\Users\\esthe\\CSC87200_FinalProject\\scRNAseq_DeconvolutionaNeuralNetwork-1")

library <- ['GSE201473', 'GSE246792', ]

#### EXPLORE SCRNA GENE EXPRESSION OMNIBUS DATA ###

install.packages("BiocManager")
install.packages("forcats")
install.packages("stringr")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("readr")
install.packages("tidyr")
install.packages("survminer")
BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("pheatmap")
BiocManager::install("org.Hs.eg.db")

library(GEOquery)

# Define the regular expression pattern
pattern <- "fetal|embryo"

# Search GEO for datasets matching the pattern
datasets <- getGEO(search_terms = pattern)


# Extract IDs of the datasets
dataset_ids <- sapply(datasets, function(dataset) attr(dataset, "GEO")$series_id)







## change my_id to be the dataset that you want.
my_id <- "GSE246792"
gse <- getGEO(my_id)

## check how many platforms used
length(gse)

gse <- gse[[1]]
gse

## if more than one dataset is present, you can analyse the other dataset by changing the number inside the [[...]]
## e.g. gse <- gse[[2]]


pData(gse) ## print the sample information
fData(gse) ## print the gene annotation
exprs(gse) ## print the expression data



# CHECK THE NORMALIZATION AND SCALES USED

## exprs get the expression levels as a data frame and get the distribution, retrieve the expression values as a data frame; with one column per-sample and one row per-gene.
summary(exprs(gse))



# visually check if data has been normalized
exprs(gse) <- log2(exprs(gse))
boxplot(exprs(gse),outline=FALSE)


library(dplyr)

sampleInfo <- pData(gse)
sampleInfo


## source_name_ch1 and characteristics_ch1.1 seem to contain factors we might need for the analysis. Let's pick just those columns

sampleInfo <- select(sampleInfo, source_name_ch1,characteristics_ch1.1)

## Optionally, rename to more convenient column names
sampleInfo <- rename(sampleInfo,group = source_name_ch1, patient=characteristics_ch1.1)

sampleInfo


library(pheatmap)
## argument use="c" stops an error if there are any missing data points

corMatrix <- cor(exprs(gse),use="c")
pheatmap(corMatrix) 


############################### LOADING IN GSE TERMS

gse_list <- list(read.csv('GSE_list.csv'))


























