
## set this to the location where the github repo was downloaded
#setwd('GITHUB-REPO-LOCATION')
setwd('~/Dropbox/Naxerova lab/Alex/cnv_trees')

## install any missing required R packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!"Biobase" %in% installed.packages()[,"Package"]) BiocManager::install("Biobase")
if(!"ACE" %in% installed.packages()[,"Package"]) BiocManager::install("ACE")
if(!"data.table" %in% installed.packages()[,"Package"]) install.packages('data.table')
if(!"ggplot2" %in% installed.packages()[,"Package"]) install.packages('ggplot2')

## source required functions
source(here::here('R/func.R'))

## create expected output directories
if(!dir.exists(here('output/C146/fits'))) dir.create(here('output/C146/fits'),recursive=T)
if(!dir.exists(here('output/C154/fits'))) dir.create(here('output/C154/fits'),recursive=T)
if(!dir.exists(here('output/C157/fits'))) dir.create(here('output/C157/fits'),recursive=T)
if(!dir.exists(here('output/C159/fits'))) dir.create(here('output/C159/fits'),recursive=T)
if(!dir.exists(here('output/C161/fits'))) dir.create(here('output/C161/fits'),recursive=T)
if(!dir.exists(here('output/C186/fits'))) dir.create(here('output/C186/fits'),recursive=T)
if(!dir.exists(here('output/E4/fits'))) dir.create(here('output/E4/fits'),recursive=T)


## process raw data from ACE to create CNV segment matrices for each subject
source(here('R/process_C146.R'))
source(here('R/process_C154.R'))
source(here('R/process_C157.R'))
source(here('R/process_C159.R'))
source(here('R/process_C161.R'))
source(here('R/process_C186.R'))
source(here('R/process_E4.R'))



