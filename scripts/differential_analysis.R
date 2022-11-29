############### DIFFERENTIAL EXPRESSION ANALYSIS ALL SAMPLES BASED ON YES/NO PCR ##################

#   Differential expression analysis with limma
library(GEOquery)
library(biomaRt)
library(limma)
library(umap)
library(dplyr)

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
if (LogC) { exC[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("POSITIVE PCR", "NEGATIVE PCR"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","SPOT_ID"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
probes <- tT2[tT2$adj.P.Val<0.01 & abs(tT2$logFC)>1,c("ID","logFC","P.Value","adj.P.Val")]

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# This is to remove the .1 that has been added to many of the ids for some reason.

row.names(dT) <- gsub(".1$", "", row.names(dT))


probe2entrez <- getBM(attributes = c("affy_hta_2_0", "entrezgene_id"), filters = "affy_hta_2_0", values=row.names(dT), mart=ensembl)
#load("probe2entrez.RData")# probe2entrez <- probe2entrez[! duplicated(probe2entrez[,1]),]entrez_rownames <- probe2entrez[,2][match(row.names(counts_modified), probe2entrez[,1])]
entrez_rownames <- probe2entrez[,2][match(row.names(dT), probe2entrez[,1])]

dT <- dT[! is.na(entrez_rownames), ]
entrez_rownames <- entrez_rownames[!is.na(entrez_rownames)]
dT <- dT[! duplicated(entrez_rownames), ]
row.names(dT) <- entrez_rownames[! duplicated(entrez_rownames)]

# Venn diagram of results
vennDiagram(dT, circle.col=palette()) 

############### DIFFERENTIAL EXPRESSION ANALYSIS ALL SAMPLES BASED ON TREATMENT ##################

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
sampleinfo <- as.data.frame(pData(gset))
sml <- as.character(as.integer(sampleinfo$characteristics_ch1.1 != "treatment: Anthracyclines and/or taxanes + carboplatin"))

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { exC[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("ANTH/TAX + CARB","ANTH/TAX"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","SPOT_ID"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
probes <- tT2[tT2$adj.P.Val<0.01 & abs(tT2$logFC)>1,c("ID","logFC","P.Value","adj.P.Val")]
# To obtain table of genes differentially expressed
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
row.names(probes) <- gsub(".1$", "", row.names(probes))
probe2entrez <- getBM(attributes = c("affy_hta_2_0", "entrezgene_id"), filters = "affy_hta_2_0", values=row.names(probes), mart=ensembl)
#load("probe2entrez.RData")# probe2entrez <- probe2entrez[! duplicated(probe2entrez[,1]),]entrez_rownames <- probe2entrez[,2][match(row.names(counts_modified), probe2entrez[,1])]
entrez_rownames <- probe2entrez[,2][match(row.names(probes), probe2entrez[,1])]

probes <- probes[! is.na(entrez_rownames), ]
entrez_rownames <- entrez_rownames[!is.na(entrez_rownames)]
probes <- probes[! duplicated(entrez_rownames), ]
row.names(probes) <- entrez_rownames[! duplicated(entrez_rownames)]

write.table(probes, file="DEGs.txt", sep="\t", quote = FALSE, row.names = FALSE)

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


# This is to remove the .1 that has been added to many of the ids for some reason.

row.names(dT) <- gsub(".1$", "", row.names(dT))


probe2entrez <- getBM(attributes = c("affy_hta_2_0", "entrezgene_id"), filters = "affy_hta_2_0", values=row.names(dT), mart=ensembl)
#load("probe2entrez.RData")# probe2entrez <- probe2entrez[! duplicated(probe2entrez[,1]),]entrez_rownames <- probe2entrez[,2][match(row.names(counts_modified), probe2entrez[,1])]
entrez_rownames <- probe2entrez[,2][match(row.names(dT), probe2entrez[,1])]

dT <- dT[! is.na(entrez_rownames), ]
entrez_rownames <- entrez_rownames[!is.na(entrez_rownames)]
dT <- dT[! duplicated(entrez_rownames), ]
row.names(dT) <- entrez_rownames[! duplicated(entrez_rownames)]


# Venn diagram of results
vennDiagram(dT, circle.col=palette()) 

############### DIFFERENTIAL EXPRESSION ANALYSIS ALL SAMPLES BASED ON YES/NO PCR CONTROLLING TREATMENT ##################

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
sampleinfo <- as.data.frame(pData(gset))
sml_nocarb <- sml[sampleinfo$characteristics_ch1.1 != "treatment: Anthracyclines and/or taxanes"]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { exC[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }
    
# assign samples to groups and set up design matrix
gs_nocarb <- factor(sml_nocarb)
groups <- make.names(c("POSITIVE PCR","NEGATIVE PCR"))
levels(gs_nocarb) <- groups
gset$group <- gs_nocarb
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs_nocarb)

fit <- lmFit(gset_nocarb, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","SPOT_ID"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


# This is to remove the .1 that has been added to many of the ids for some reason.

row.names(dT) <- gsub(".1$", "", row.names(dT))


probe2entrez <- getBM(attributes = c("affy_hta_2_0", "entrezgene_id"), filters = "affy_hta_2_0", values=row.names(dT), mart=ensembl)
#load("probe2entrez.RData")# probe2entrez <- probe2entrez[! duplicated(probe2entrez[,1]),]entrez_rownames <- probe2entrez[,2][match(row.names(counts_modified), probe2entrez[,1])]
entrez_rownames <- probe2entrez[,2][match(row.names(dT), probe2entrez[,1])]

dT <- dT[! is.na(entrez_rownames), ]
entrez_rownames <- entrez_rownames[!is.na(entrez_rownames)]
dT <- dT[! duplicated(entrez_rownames), ]
row.names(dT) <- entrez_rownames[! duplicated(entrez_rownames)]

# Venn diagram of results
vennDiagram(dT, circle.col=palette()) 


####### DIFFERENTIAL EXPRESSION ANALYSIS NO CARBOPLATIN SAMPLES, BASED ON YES/NO PCR #############
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
sampleinfo <- as.data.frame(pData(gset))
sml_nocarb <- sml[sampleinfo$characteristics_ch1.1 == "treatment: Anthracyclines and/or taxanes"]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }
  
  gset_nocarb <- gset[,sampleinfo$characteristics_ch1.1 == "treatment: Anthracyclines and/or taxanes"]
  
# assign samples to groups and set up design matrix
gs_nocarb <- factor(sml_nocarb)
groups <- make.names(c("POSITIVE PCR","NEGATIVE PCR"))
levels(gs_nocarb) <- groups
gset_nocarb$group <- gs_nocarb
design <- model.matrix(~group + 0, gset_nocarb)
colnames(design) <- levels(gs_nocarb)

fit <- lmFit(gset_nocarb, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","SPOT_ID"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


# This is to remove the .1 that has been added to many of the ids for some reason.

row.names(dT) <- gsub(".1$", "", row.names(dT))


probe2entrez <- getBM(attributes = c("affy_hta_2_0", "entrezgene_id"), filters = "affy_hta_2_0", values=row.names(dT), mart=ensembl)
#load("probe2entrez.RData")# probe2entrez <- probe2entrez[! duplicated(probe2entrez[,1]),]entrez_rownames <- probe2entrez[,2][match(row.names(counts_modified), probe2entrez[,1])]
entrez_rownames <- probe2entrez[,2][match(row.names(dT), probe2entrez[,1])]

dT <- dT[! is.na(entrez_rownames), ]
entrez_rownames <- entrez_rownames[!is.na(entrez_rownames)]
dT <- dT[! duplicated(entrez_rownames), ]
row.names(dT) <- entrez_rownames[! duplicated(entrez_rownames)]

# Venn diagram of results
vennDiagram(dT, circle.col=palette()) 

#######################################################################################
####### DIFFERENTIAL EXPRESSION ANALYSIS CARBOPLATIN SAMPLES, BASED ON YES/NO PCR #############
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
sampleinfo <- as.data.frame(pData(gset))
sml_carb <- sml[sampleinfo$characteristics_ch1.1 != "treatment: Anthracyclines and/or taxanes"]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }
  
  gset_carb <- gset[,sampleinfo$characteristics_ch1.1 != "treatment: Anthracyclines and/or taxanes"]
  
# assign samples to groups and set up design matrix
gs_carb <- factor(sml_carb)
groups <- make.names(c("POSITIVE PCR","NEGATIVE PCR"))
levels(gs_carb) <- groups
gset_carb$group <- gs_carb
design <- model.matrix(~group + 0, gset_carb)
colnames(design) <- levels(gs_carb)

fit <- lmFit(gset_carb, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","SPOT_ID"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


# This is to remove the .1 that has been added to many of the ids for some reason.

row.names(dT) <- gsub(".1$", "", row.names(dT))


probe2entrez <- getBM(attributes = c("affy_hta_2_0", "entrezgene_id"), filters = "affy_hta_2_0", values=row.names(dT), mart=ensembl)
#load("probe2entrez.RData")# probe2entrez <- probe2entrez[! duplicated(probe2entrez[,1]),]entrez_rownames <- probe2entrez[,2][match(row.names(counts_modified), probe2entrez[,1])]
entrez_rownames <- probe2entrez[,2][match(row.names(dT), probe2entrez[,1])]

dT <- dT[! is.na(entrez_rownames), ]
entrez_rownames <- entrez_rownames[!is.na(entrez_rownames)]
dT <- dT[! duplicated(entrez_rownames), ]
row.names(dT) <- entrez_rownames[! duplicated(entrez_rownames)]

# Venn diagram of results
vennDiagram(dT, circle.col=palette()) 



