# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Differential expression analysis with limma
setwd ("/home/eva/Escritorio/Biotec_TFM")
library(GEOquery)
library(limma)
library(umap)
library("dplyr")

# load series and platform data from GEO
# gset <- getGEO("GSE106977", GSEMatrix =TRUE, AnnotGPL=FALSE)
# if (length(gset) > 1) idx <- grep("GPL17586", attr(gset, "names")) else idx <- 1
# gset <- gset[[idx]]

#save(gset, file="gset_full.RData")
load("gset_full.RData")
gset <- gset[rowMedians(Biobase::exprs(gset)) > 5,]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
#gsms <- paste0("11101001111110111111011110101001110011101001101111",
   #   "10000010001111110011111110101101101011011110100101",
   # "0011111001100001000")
#sml <- strsplit(gsms, split="")[[1]]
sml_nocarb <- sml[sampleinfo$characteristics_ch1.1 == "treatment: Anthracyclines and/or taxanes"]
sml_carb <- sml[sampleinfo$characteristics_ch1.1 != "treatment: Anthracyclines and/or taxanes"]

sampleinfo <- as.data.frame(pData(gset))

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }
  
gset_medians <- rowMedians(Biobase::exprs(gset))

hist_res <- hist(gset_medians, 100, col = "cornsilk1", freq = FALSE, 
            main = "Histogram of the median intensities", 
            border = "antiquewhite4",
            xlab = "Median intensities")
            
# Divide gset based on treatment
gset_nocarb <- gset[,sampleinfo$characteristics_ch1.1 == "treatment: Anthracyclines and/or taxanes"]
gset_carb <- gset[,sampleinfo$characteristics_ch1.1 != "treatment: Anthracyclines and/or taxanes"]

# assign samples to groups and set up design matrix (NOCARB)
gs_nocarb <- factor(sml_nocarb)
groups <- make.names(c("POSITIVE PCR","NEGATIVE PCR"))
levels(gs_nocarb) <- groups
gset_nocarb$group <- gs_nocarb
design <- model.matrix(~group + 0, gset_nocarb)
colnames(design) <- levels(gs_nocarb)

fit <- lmFit(gset_nocarb, design)  # fit linear model

# assign samples to groups and set up design matrix (CARB)
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
# Genes for adj.P.Val<0.01 y P.Value>1 in absolute value
tT2[tT2$adj.P.Val<0.01 & abs(tT2$logFC)>1,c("ID","logFC","P.Value","adj.P.Val")]


hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
  ylab = "Number of genes", main = "P-adj value distribution")
  

# DATAFRAME CREATED BY TREATMENT AND PCR
CARB <- as.data.frame(pData(gset_carb))
NOCARB <- as.data.frame(pData(gset_nocarb))
CARBTABLE <- select (CARB, "pathological complete response:ch1")
NOCARBTABLE <- select (NOCARB, "pathological complete response:ch1")
write.table(CARBTABLE)
write.table(NOCARBTABLE)

#save(CARBTABLE, file= "CARBTABLE.RData")
load("CARBTABLE.RData")

#save(NOCARBTABLE, file = "NOCARBTABLE.RData")
load("NOCARBTABLE.RData")

# TARGET TABLE
target <- as.data.frame(pData(gset))
targetable <- select (target, "pathological complete response:ch1", "treatment:ch1" )
write.table(targetable, file = "completetargettable.R", sep="\t", quote = FALSE, row.names = FALSE)

# create third column with sample names
targetable$sample <- row.names(targetable)

# modify name columns
colnames(targetable) <- c("treat", "treat_received", "sample")
targetable[targetable == "yes"] <- "Treat"
targetable[targetable == "no"] <- "Ctrl"

#modify order of column
TARGET <- targetable[,c(3,1,2)]
write.table(TARGET, file = "target.txt", sep="\t", quote = FALSE, row.names = FALSE)


#COUNTS TABLE
exprs(gset)
counts <- exprs(gset)
Countsrounded <- round(counts)
write.table(Countsrounded, sep="\t", quote = FALSE, file = "counts.txt")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)

# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
  highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

################################################################
# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
dev.new(width=3+ncol(gset)/6, height=5)
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE106977", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")
dev.off()

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE106977", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
col=1:nlevels(gs), title="Group", pt.cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

# mean-variance trend, helps to see if precision weights are needed
plotSA(fit2, main="Mean variance trend, GSE106977")





write.csv(tT2, "/home/eva/Dropbox/EvaCastillo/tablapCR")

TABLEYESNO <- read.csv("tabla.txt", sep=",", row.names=1)

# para modificar un objeto o vector. Donde ponía yes ahora pondrá treat y donde pone no ahora ctrl.
> CARBTABLE[CARBTABLE == "yes"] <- "Treat"
> CARBTABLE[CARBTABLE == "no"] <- "Ctrl"
> CARBTABLE
write.table(CARBTABLE, file = "CARBTABLE.R", sep="\t", quote = FALSE)

NOCARBTABLE[NOCARBTABLE == "yes"] <- "Treat"
NOCARBTABLE[NOCARBTABLE == "no"] <- "Ctrl"
NOCARBTABLE
write.table(NOCARBTABLE, file = "NOCARBTABLE.R", sep="\t", quote = FALSE)



