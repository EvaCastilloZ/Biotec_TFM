library(ExpHunterSuite)
 # Code to load target file
table_of_target <- read.csv("../Data/target.txt", sep="\t")
# Code to load counts file
counts_no_rownames <- read.table("../Data/counts.txt", sep="\t", row.names=1, header=TRUE)
# Code to run expression analysis
#dea_4_modules <- main_degenes_Hunter(raw=counts_no_rownames, 
                                     # target=tableoftarget,
                                     # modules="DNLE",
                                     # output_files = "RESULTS")

# To run coexpression analysis (with WGCNA)
dea_coexpression <- main_degenes_Hunter(raw=counts_no_rownames,
                                        target=table_of_target, 
                                        modules="DW", 
                                        output_files = "RESULTS")
save(dea_coexpression, file="../Data/dea_coexpression.RData")
load("dea_coexpression.RData")
# Code to make a report on the expression analysis
#write_expression_report(exp_results=dea_4_modules) 
write_expression_report(exp_results=dea_coexpression) 
# Code to run functional analysis
#functional_analysis <- main_functional_hunter(dea_4_modules, 'Human', enrich_dbs = c("MF", "BP","Reactome"), enrich_methods = "ORA")
coexpression_functional_analysis <- main_functional_hunter(dea_coexpression,
                                                           'Human', 
                                                           enrich_dbs = c("MF","BP","Reactome"), 
                                                           enrich_methods = "ORA")
# Code to make a report on the functional analysis
print(wd())
write_enrich_files(func_results=functional_analysis)
write_enrich_files(func_results=coexpression_functional_analysis)
write_functional_report(hunter_results=dea_4_modules,
                        func_results=functional_analysis)
write_functional_report(hunter_results=dea_coexpression,
                        func_results=coexpression_functional_analysis)
