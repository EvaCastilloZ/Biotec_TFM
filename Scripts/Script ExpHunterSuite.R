cd ~/software/Biotec_TFM/Data
R
library("biomaRt")
library(GEOquery)
library(ExpHunterSuite)
library(limma)
library(umap)
library("dplyr")

load("gset_full.RData")
gset <- gset[rowMedians(Biobase::exprs(gset)) > 5,]

exprs(gset)
counts <- exprs(gset)
counts_rounded <- round(counts)
counts_modified <- counts_rounded

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


# This is to remove the .1 that has been added to many of the ids for some reason.

row.names(counts_modified) <- gsub(".1$", "", row.names(counts_modified))


probe2entrez <- getBM(attributes = c("affy_hta_2_0", "entrezgene_id"), filters = "affy_hta_2_0", values=row.names(counts_modified), mart=ensembl); save(probe2entrez, file="probe2entrez.RData")
#load("probe2entrez.RData")# probe2entrez <- probe2entrez[! duplicated(probe2entrez[,1]),]entrez_rownames <- probe2entrez[,2][match(row.names(counts_modified), probe2entrez[,1])]
entrez_rownames <- probe2entrez[,2][match(row.names(counts_modified), probe2entrez[,1])]


counts_modified <- counts_modified[! is.na(entrez_rownames), ]
entrez_rownames <- entrez_rownames[!is.na(entrez_rownames)]
counts_modified <- counts_modified[! duplicated(entrez_rownames), ]
row.names(counts_modified) <- entrez_rownames[! duplicated(entrez_rownames)]

write.table(counts_modified, sep="\t", quote = FALSE, file = "counts_modified.txt")

 # Code to load target file
table_of_target <- read.csv("../Data/target.txt", sep="\t")
# Code to load counts file
counts_no_rownames <- read.table("../Data/counts_modified.txt", sep="\t", row.names=1, header=TRUE)

# To run coexpression analysis (with WGCNA)
hunter_results <- main_degenes_Hunter(raw=counts_no_rownames,
                                        target=table_of_target, 
                                        modules="DW", 
                                        output_files = "RESULTS")
                                        
write_expression_report(exp_results=hunter_results) 
#save(dea_coexpression, file="../Data/hunter_results.RData")
load("hunter_results.RData")
 
hunter_results$DE_all_genes[hunter_results$DE_all_genes == "NOT_DEG"] <- "PREVALENT_DEG"
#save(hunter_results, file="../Data/hunter_results_genetagmodified.RData")
load("hunter_results_genetagmodified.RData")

# Code to run functional analysis
#functional_analysis <- main_functional_hunter(dea_4_modules, 'Human', enrich_dbs = c("MF", "BP","Reactome"), enrich_methods = "ORA")
func_results <- main_functional_hunter(hunter_results, 'Human', 
                                                           enrich_dbs = c("MF","BP","CC"), 
                                                           enrich_methods = "ORA",
                                                           input_gene_id = "ENTREZID" )
write_enrich_files(func_results)
write_functional_report(hunter_results=hunter_results,
                        func_results=func_results)


