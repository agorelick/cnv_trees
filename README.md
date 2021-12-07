## Regenerate somatic copy number alteration (SCNA) analyses for [peritoneal metastases paper]

Use this repo to process shallow WGS-based SCNA data. This will automatically do the following:
- Fit the purity/ploidy of each sample (starting from 1Mb binned counts from ACE R-package).
- Generate copy number segments using purity/ploidy corrected binned counts.
- Generate SCNA euclidean distance matrix/trees and compare to pre-generated poly-G angular distance data for each patient.
- Generate SCNA heatmap for each patient.

## Instructions

1. In a mac/linux terminal, clone the repo:
`git clone https://github.com/agorelick/cnv_trees`

2. Open R and install the accompanying *polyG* R-package:
```r
install.packages('cnv_trees/rpkg',type='src',repos=NULL)
```

3. 'cd' to the cloned repo, run the `run.R` script. That's it!

`
cd cnv_trees
Rscript run.R
`

