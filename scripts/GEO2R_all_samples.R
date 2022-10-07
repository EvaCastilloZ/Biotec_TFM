######################################################################################################
############################################# GEO2R data ############################################
######################################################################################################

################# Analysis based in yes/no pCR (all samples) #####################################
# loading libraries
library(biomaRt)
library(dplyr)
library(GEOquery)

# load series and platform data from GEO
gset <- getGEO("GSE106977", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL17586", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("11101001111110111111011110101001110011101001101111",
     "10000010001111110011111110101101101011011110100101",
   "0011111001100001000")
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }
 
# Obtain counts table rounding them
exprs(gset)
counts <- exprs(gset)
counts_modified <- round(counts)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# This is to remove the .1 that has been added to many of the ids for some reason.

row.names(counts_modified) <- gsub(".1$", "", row.names(counts_modified))

probe2entrez <- getBM(attributes = c("affy_hta_2_0", "entrezgene_id"), filters = "affy_hta_2_0", values=row.names(counts_modified), mart=ensembl)
#load("probe2entrez.RData")# probe2entrez <- probe2entrez[! duplicated(probe2entrez[,1]),]entrez_rownames <- probe2entrez[,2][match(row.names(counts_modified), probe2entrez[,1])]
entrez_rownames <- probe2entrez[,2][match(row.names(counts_modified), probe2entrez[,1])]


counts_modified <- counts_modified[! is.na(entrez_rownames), ]
entrez_rownames <- entrez_rownames[!is.na(entrez_rownames)]
counts_modified <- counts_modified[! duplicated(entrez_rownames), ]
row.names(counts_modified) <- entrez_rownames[! duplicated(entrez_rownames)]

#Generate the counts table
write.table(counts_modified, sep="\t", quote = FALSE, file = "counts.txt")

# Obtain target table
target <- as.data.frame(pData(gset))
targetable <- select (target, "pathological complete response:ch1", "treatment:ch1" )

# create third column with sample names
targetable$sample <- row.names(targetable)

# Modify name columns
colnames(targetable) <- c("treat", "treat_received", "sample")
targetable[targetable == "yes"] <- "Treat"
targetable[targetable == "no"] <- "Ctrl"

# Modify order of column
new_target <- targetable[,c(3,1,2)]
# Generate the target table
write.table(new_target, file = "target.txt", quote = FALSE, row.names = FALSE, sep="\t")

###################################################################################
############### Analysis based in treatment (all samples) ########################################

# load series and platform data from GEO
gset <- getGEO("GSE106977", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL17586", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("11101001111110111111011110101001110011101001101111",
     "10000010001111110011111110101101101011011110100101",
   "0011111001100001000")
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }
             
# Obtain counts table rounding them
exprs(gset)
counts <- exprs(gset)
counts_modified <- round(counts)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# This is to remove the .1 that has been added to many of the ids for some reason.

row.names(counts_modified) <- gsub(".1$", "", row.names(counts_modified))

probe2entrez <- getBM(attributes = c("affy_hta_2_0", "entrezgene_id"), filters = "affy_hta_2_0", values=row.names(counts_modified), mart=ensembl)
#load("probe2entrez.RData")# probe2entrez <- probe2entrez[! duplicated(probe2entrez[,1]),]entrez_rownames <- probe2entrez[,2][match(row.names(counts_modified), probe2entrez[,1])]
entrez_rownames <- probe2entrez[,2][match(row.names(counts_modified), probe2entrez[,1])]

counts_modified <- counts_modified[! is.na(entrez_rownames), ]
entrez_rownames <- entrez_rownames[!is.na(entrez_rownames)]
counts_modified <- counts_modified[! duplicated(entrez_rownames), ]
row.names(counts_modified) <- entrez_rownames[! duplicated(entrez_rownames)]

#Generate the counts table
write.table(counts_modified, sep="\t", quote = FALSE, file = "counts_treatment.txt")

# Obtain target table
target <- as.data.frame(pData(gset))
targetable <- select (target, "pathological complete response:ch1", "treatment:ch1" )

# create third column with sample names
targetable$sample <- row.names(targetable)

# Modify name columns
colnames(targetable) <- c("pCR", "treat", "sample")
targetable[targetable == "Anthracyclines and/or taxanes"] <- "Treat"
targetable[targetable == "Anthracyclines and/or taxanes + carboplatin"] <- "Ctrl"

# Modify order of column
new_target <- targetable[,c(3,1,2)]
# Generate the target table
write.table(new_target, file = "target_treatment.txt", sep="\t", quote = FALSE, row.names = FALSE)

