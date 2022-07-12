 setwd ("/home/eva/Escritorio/")
 # Code to load target file
 read.csv("TARGET.txt")
tableoftarget <- read.csv("TARGET.txt", sep="\t")
# Code to load counts file
read.table("COUNTS.txt")
tableofcountswrn <- read.table("COUNTS.txt", sep="\t", row.names=1, header=TRUE)
# Code to run expression analysis
library(ExpHunterSuite)
#DEAcomplete <- main_degenes_Hunter(raw=tableofcountswrn, target=tableoftarget, modules="DNLE", output_files = "RESULTS")

# TO RUN WGCNA too
DEAcompleteCOexp <- main_degenes_Hunter(raw=tableofcountswrn, target=tableoftarget, modules="DW", output_files = "RESULTS")
# Code to make a report on the expression analysis
#write_expression_report(exp_results=DEAcomplete) 
write_expression_report(exp_results=DEAcompleteCOexp) 
# Code to run functional analysis
#Functionalanalysis <- main_functional_hunter(DEAcomplete, 'Human')
FunctionalanalysisCOexp <- main_functional_hunter(DEAcompleteCOexp, 'Human')
# Code to make a report on the functional analysis
print(wd())
write_enrich_files(func_results=Functionalanalysis)
write_functional_report(hunter_results=DEAcomplete, func_results=Functionalanalysis)


> write_functional_report(hunter_results=DEAcompleteCOexp, func_results=FunctionalanalysisCOexp)
