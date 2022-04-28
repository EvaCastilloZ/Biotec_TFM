# Specify the work directory
 setwd ("/(where target and counts file are)")
# Code to load target file
 read.csv("target.txt")
tableoftarget <- read.csv("target.txt", sep="\t")
# Code to load counts file
read.table("counts.txt")
tableofcountswrn <- read.table("counts.txt", sep="\t", row.names=1, header=TRUE)
# Code to run expression analysis (and coexpression)
library(ExpHunterSuite)
DEAcomplete <- main_degenes_Hunter(raw=tableofcountswrn, target=tableoftarget, modules="DNLE", output_files = "RESULTS")
DEAcompleteCOexp <- main_degenes_Hunter(raw=tableofcountswrn, target=tableoftarget, modules="DNLEW", output_files = "RESULTS")
# Code to make a report on the expression analysis
write_expression_report(exp_results=DEAcompleteCOexp) 
# Code to run functional analysis
Functionalanalysis <- functional_hunter(DEAcomplete, 'Mouse', func_annot_db = "gKR", GO_subont = "BMC", analysis_type= "o")
FunctionalanalysisCOexp <- functional_hunter(DEAcompleteCOexp, 'Mouse', func_annot_db = "gKR", GO_subont = "BMC", analysis_type= "o")
# Code to make a report on the functional analysis
print(wd())
write_enrich_files(func_results=Functionalanalysis)
write_functional_report(hunter_results=DEAcomplete, func_results=Functionalanalysis)
# Code to make a report on the functional analysis (for coexpression analysis)
> write_functional_report(hunter_results=DEAcompleteCOexp, func_results=FunctionalanalysisCOexp)
