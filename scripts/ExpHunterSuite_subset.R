#######################################################################################################
###################### RUNNING ANTHRACYCLINES AND/OR TAXANES + CARBOPLATIN ANALYSIS ############################
###################################################################################################
# loading library
library(ExpHunterSuite)
 
 # code to load target file
table_of_target <- read.csv("carb_target.txt", sep="\t")
# Code to load counts file
counts_no_rownames <- read.table("carb_counts.txt", header=TRUE, row.names=1, sep="\t")

# to run diferential and coexpression analysis
hunter_results_coexp <- main_degenes_Hunter(modules="DW", 
					    raw=counts_no_rownames,
                                            target=table_of_target)
                                                                           
                                        
# we recommend use the save function because the following function often fails               
                                        
# code to obtain expression report                          
write_expression_report(exp_results=hunter_results_coexp)
 
 
hunter_results_coexp$DE_all_genes[hunter_results_coexp$DE_all_genes == "NOT_DEG"] <- "PREVALENT_DEG"

# Code to run functional analysis
func_results <- main_functional_hunter(hunter_results_coexp, 'Human', 
                                                           enrich_dbs = c("MF","BP","CC"), 
                                                           enrich_methods = "ORA",
                                                           input_gene_id = "ENTREZID" )
                        
# Code to obtain functional report                                   
write_functional_report(func_results=func_results, 
			hunter_results=hunter_results_coexp,
                        report = "fci")

## NOTE: Before using the next analysis, note you should modify the report names manually because their function doesn't let modifying the file name. Also it's possible to specify the path to save the report with output_files in both functions.

#######################################################################################################
###################### RUNNING ANTHRACYCLINES AND/OR TAXANES ANALYSIS ############################
###################################################################################################

# code to load target file
table_of_target <- read.csv("nocarb_target.txt", sep="\t")
# Code to load counts file
counts_no_rownames <- read.table("nocarb_counts.txt", header=TRUE, row.names=1, sep="\t")

# to run diferential and coexpression analysis
hunter_results_coexp <- main_degenes_Hunter(modules="DW", 
					    raw=counts_no_rownames,
                                            target=table_of_target)
                                        
# we recommend use the save function because the following function often fails               
                                        
# to obtain expression report                          
write_expression_report(exp_results=hunter_results_coexp)
 
hunter_results_coexp$DE_all_genes[hunter_results_coexp$DE_all_genes == "NOT_DEG"] <- "PREVALENT_DEG"

# to run functional analysis
func_results <- main_functional_hunter(hunter_results_coexp, 'Human', 
                                                           enrich_dbs = c("MF","BP","CC"), 
                                                           enrich_methods = "ORA",
                                                           input_gene_id = "ENTREZID" )
                        
# Code to obtain functional report                                   
write_functional_report(hunter_results=hunter_results_coexp,
                        func_results=func_results,
                        report = "fci")
