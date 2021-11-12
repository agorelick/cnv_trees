
## set this to the location where the github repo was downloaded
#setwd('GITHUB-REPO-LOCATION')
setwd('~/lab_repos/cnv_trees_tmp')

## install any missing required R packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!"Biobase" %in% installed.packages()[,"Package"]) BiocManager::install("Biobase")
if(!"ACE" %in% installed.packages()[,"Package"]) BiocManager::install("ACE")
if(!"data.table" %in% installed.packages()[,"Package"]) install.packages('data.table')
if(!"ggplot2" %in% installed.packages()[,"Package"]) install.packages('ggplot2')

## load required R packages
library(Biobase)
library(data.table)
library(ACE)
library(ggplot2)

## source required functions
source('R/func.R')

## create expected output directories
if(!dir.exists('output')) dir.create('output')

## process raw data from ACE to create CNV segment matrices for each subject
source('R/process_C146.R')
source('R/process_C154.R')
source('R/process_C157.R')
source('R/process_E4.R')



