# Loading library
library(ExpHunterSuite)
 
 # Code to load target file
table_of_target <- read.csv("target.txt", sep="\t")
# Code to load counts file
counts_no_rownames <- read.table("counts_modified.txt", sep="\t", row.names=1, header=TRUE)

# To run coexpression analysis (with WGCNA)
hunter_results_coexp <- main_degenes_Hunter(raw=counts_no_rownames,
                                        target=table_of_target, 
                                        modules="DELNW", 
                                        output_files = "RESULTS")
                                        
# we recommend use the save function because its size                 
                                        
# Code to obtain expression report                          
write_expression_report(exp_results=hunter_results_coexp)
 
 
hunter_results_coexp$DE_all_genes[hunter_results_coexp$DE_all_genes == "NOT_DEG"] <- "PREVALENT_DEG"

# Code to run functional analysis
func_results <- main_functional_hunter(hunter_results_coexp, 'Human', 
                                                           enrich_dbs = c("MF","BP","CC"), 
                                                           enrich_methods = "ORA",
                                                           input_gene_id = "ENTREZID" )
                        
# Code to obtain functional report                                   
write_functional_report(hunter_results=hunter_results_coexp,
                        func_results=func_results,
                        report = "fci")


