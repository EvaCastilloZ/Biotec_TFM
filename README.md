Script to run expression and functional analysis for a file with packages from ExpHunterSuite (Bioconductor)

# REQUERIMENTS
- R version 4.2 (or later)
- Bioconductor v. 3.15

# Installation
Install latest R-Version
Go to page https://cloud.r-project.org/ and install the latest R version on your computer. Install also the latest Biodoncudctor version in http://bioconductor.org/install
# Install ExpHunterSuite and required R-packages
To download ExpHunterSuite: https://bioconductor.org/packages/release/workflows/html/ExpHunterSuite.html

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ExpHunterSuite")

# Run ExpHunterSuite
Expression and functional analysis with the run_ehs.R script
