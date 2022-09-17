######################################################################################################
############################################# GEO2R data ############################################
######################################################################################################

# loading libraries
library(biomaRt)
library(dplyr)
library(GEOquery)
library(limma)

# load series and platform data from GEO
gset <- getGEO("GSE106977", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL17586", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

gset_medians <- rowMedians(Biobase::exprs(gset))

hist_res <- hist(gset_medians, 100, col = "cornsilk1", freq = FALSE, 
            main = "Histogram of the median intensities", 
            border = "antiquewhite4",
            xlab = "Median intensities")

#gset <- gset[rowMedians(Biobase::exprs(gset)) > 5,]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("11101001111110111111011110101001110011101001101111",
     "10000010001111110011111110101101101011011110100101",
   "0011111001100001000")
sml <- strsplit(gsms, split="")[[1]]
#sampleinfo <- as.data.frame(pData(gset))
#sml_nocarb <- sml[sampleinfo$characteristics_ch1.1 == "treatment: Anthracyclines and/or taxanes"]
#sml_carb <- sml[sampleinfo$characteristics_ch1.1 != "treatment: Anthracyclines and/or taxanes"]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }
  
  # assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("yes pcr","no pcr"))
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
             

# Obtain counts table rounding them
exprs(gset)
counts <- exprs(gset)
counts_rounded <- round(counts)
counts_modified <- counts_rounded

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
write.table(counts_modified, sep="\t", quote = FALSE, file = "counts_modified.txt")

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
write.table(new_target, file = "target.txt", sep="\t", quote = FALSE, row.names = FALSE)

