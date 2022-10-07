#######################################################################################################
###################### GENERATING TABLES BASED ON ANTHRACYCLINES AND/OR TAXANES + CARBOPLATIN TREATMENT ############################
###################################################################################################
# calling libraries
library(biomaRt)
library(GEOquery)
library(dplyr)

# load series and platform data from GEO
gset <- getGEO("GSE106977", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL17586", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# subsetting based on treatment type
sampleinfo <- as.data.frame(pData(gset))
gset_carb <- gset[,sampleinfo$characteristics_ch1.1 != "treatment: Anthracyclines and/or taxanes"]

# make proper column names to match toptable 
fvarLabels(gset_carb) <- make.names(fvarLabels(gset_carb))

# group membership for all samples
gsms <- paste0("11101001111110111111011110101001110011101001101111",
     "10000010001111110011111110101101101011011110100101",
   "0011111001100001000")
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
exC <- exprs(gset_carb)
qx <- as.numeric(quantile(exC, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { exC[which(exC <= 0)] <- NaN
  exprs(gset_carb) <- log2(exC) }
  
# Obtain counts table rounding them
exprs(gset_carb)
counts_c <- exprs(gset_carb)
counts_rounded_c <- round(counts_c)
counts_modified_c <- counts_rounded_c

# This is to remove the .1 that has been added to many of the ids for some reason.
row.names(counts_modified_c) <- gsub(".1$", "", row.names(counts_modified_c))

# Modify the gene tag
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
probe2entrez <- getBM(attributes = c("affy_hta_2_0", "entrezgene_id"), filters = "affy_hta_2_0", values=row.names(counts_modified_c), mart=ensembl)
#load("probe2entrez.RData")# probe2entrez <- probe2entrez[! duplicated(probe2entrez[,1]),]entrez_rownames <- probe2entrez[,2][match(row.names(counts_modified), probe2entrez[,1])]
entrez_rownames <- probe2entrez[,2][match(row.names(counts_modified_c), probe2entrez[,1])]

# To eliminate duplicated IDs
counts_modified_c <- counts_modified_c[! is.na(entrez_rownames), ]
entrez_rownames <- entrez_rownames[!is.na(entrez_rownames)]
counts_modified_c <- counts_modified_c[! duplicated(entrez_rownames), ]
row.names(counts_modified_c) <- entrez_rownames[! duplicated(entrez_rownames)]

# generate the counts table
write.table(counts_modified_c, file = "carb_counts.txt", sep="\t", quote = FALSE)

# obtain target table
target_carb <- as.data.frame(pData(gset_carb))
targetable_carb <- select(target_carb, "pathological complete response:ch1", "treatment:ch1" )

# create third column with sample names
targetable_carb$sample <- row.names(targetable_carb)

# modify name columns
colnames(targetable_carb) <- c("treat", "treat_received", "sample")
targetable_carb[targetable_carb == "yes"] <- "Treat"
targetable_carb[targetable_carb == "no"] <- "Ctrl"

# modify order of column
new_target_c <- targetable_carb[,c(3,1,2)]

# generate the target table
write.table(new_target_c, file = "carb_target.txt", sep="\t", quote = FALSE, row.names = FALSE)

#######################################################################################################
###################### GENERATING TABLES BASED ON ANTHRACYCLINES AND/OR TAXANES  TREATMENT ############################
###################################################################################################

# subsetting based on treatment type
sampleinfo <- as.data.frame(pData(gset))
gset_nocarb <- gset[,sampleinfo$characteristics_ch1.1 == "treatment: Anthracyclines and/or taxanes"]

# make proper column names to match toptable 
fvarLabels(gset_nocarb) <- make.names(fvarLabels(gset_nocarb))

# group membership for all samples
gsms <- paste0("11101001111110111111011110101001110011101001101111",
     "10000010001111110011111110101101101011011110100101",
   "0011111001100001000")
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
exNC <- exprs(gset_nocarb)
qx <- as.numeric(quantile(exNC, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { exNC[which(exNC <= 0)] <- NaN
  exprs(gset_nocarb) <- log2(exNC) }
  
exprs(gset_nocarb)
counts_nc <- exprs(gset_nocarb)
counts_rounded_nc <- round(counts_nc)
counts_modified_nc <- counts_rounded_nc


# This is to remove the .1 that has been added to many of the ids for some reason.
row.names(counts_modified_nc) <- gsub(".1$", "", row.names(counts_modified_nc))

# Modify the gene tag
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
probe2entrez <- getBM(attributes = c("affy_hta_2_0", "entrezgene_id"), filters = "affy_hta_2_0", values=row.names(counts_modified_nc), mart=ensembl)
#load("probe2entrez.RData")# probe2entrez <- probe2entrez[! duplicated(probe2entrez[,1]),]entrez_rownames <- probe2entrez[,2][match(row.names(counts_modified_nc), probe2entrez[,1])]
entrez_rownames <- probe2entrez[,2][match(row.names(counts_modified_nc), probe2entrez[,1])]

# To eliminate duplicated IDs
counts_modified_nc <- counts_modified_nc[! is.na(entrez_rownames), ]
entrez_rownames <- entrez_rownames[!is.na(entrez_rownames)]
counts_modified_nc <- counts_modified_nc[! duplicated(entrez_rownames), ]
row.names(counts_modified_nc) <- entrez_rownames[! duplicated(entrez_rownames)]

#Generate the counts table
write.table(counts_modified_nc, file = "nocarb_counts.txt", quote = FALSE, sep="\t")

# Obtain target table
target_nocarb <- as.data.frame(pData(gset_nocarb))
targetable_nocarb <- select(target_nocarb, "pathological complete response:ch1", "treatment:ch1" )

# create third column with sample names
targetable_nocarb$sample <- row.names(targetable_nocarb)

# Modify name columns
colnames(targetable_nocarb) <- c("treat", "treat_received", "sample")
targetable_nocarb[targetable_nocarb == "yes"] <- "Treat"
targetable_nocarb[targetable_nocarb == "no"] <- "Ctrl"

# Modify order of column
new_target_nc <- targetable_nocarb[,c(3,1,2)]

# Generate the target table
write.table(new_target_nc, file = "nocarb_target.txt", sep="\t", quote = FALSE, row.names = FALSE)

  

