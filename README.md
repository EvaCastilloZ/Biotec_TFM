# DOCUMENTATION - README
> Master's end project that tries to find genes with differential expression between many patients that has developed triple negative breast cancer (TNBC).

## REQUERIMENTS
R version 4.2 (or later)

## Installation

### Install latest R-Version 
Go to page https://cloud.r-project.org/ and install the latest R version on your computer.

### Install Bioconductor
To install the latest Bioconductor version, go to page http://bioconductor.org/install/

### Install ExpHunterSuite and required R-packages
To install ExpHunterSuite from console, we need to use the devtools to install R packages from GitHub. It can be done with these commands:

``` bash
devtools::install_github("seoanezonjic/ExpHunterSuite")
```

To install the latest versions of all R-packages required to run ExpHunterSuite, use the necessary_libraries.R located in Scripts folder.

## Obtaining the counts and target table
To obtain these dataframes, we need to run the GEO2R_script.R located in the Scripts folder.

## Run ExpHunterSuite
Expression and functional analysis with the "script_ExpHunterSuite.R" script and obtain functional and expression reports (and/or coexpression adding module W to main_degenes_hunter function. 


