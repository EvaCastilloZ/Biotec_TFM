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

To install required R-packages, try to install them via BiocManager (http://bioconductor.org/install/) calling BiocManager::install() or by CRAN. For CRAN, you can type the following:
```bash
install.packages(c("biomaRt", "dplyr", "GEOquery", "limma"))
```

## Obtaining the counts and target table from GEO2R
To obtain these dataframes, we need to run the GEO2R_script.R located in the Scripts folder.

## Run ExpHunterSuite and obtaining reports
Run expression and functional analysis with the "script_ExpHunterSuite.R" script and obtain functional and expression reports (and/or coexpression adding module W to main_degenes_hunter function). 

## CITATIONS

Perkins J, Seoane Zonjic P, Moreno Jabato F, Córdoba Caballero J,
  Rojano Rivera E, Bautista Moreno R, Claros M, Gonzalez Gayte I,
  García Ranea J (2022). _ExpHunterSuite: Package For The Comprehensive
  Analysis Of Transcriptomic Data_. R package version 0.99.13.

