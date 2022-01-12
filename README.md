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
install.packages('cnv_trees/pkg',type='src',repos=NULL)
```

3. 'cd' to the cloned repo, run the `run.R` script. That's it!

`
cd cnv_trees
Rscript run.R
`

## The 'polyG' R-package

The polyG package contained in the `pkg/` directory has two main purposes. 

First, it contains many 'convenience' functions that are used in various analyses in this project. 

Second, it provides functions to generate 'annotated_phylo' objects, which is a custom object built around standard 'phylo' type objects, but containing additional metadata for each leaf/sample on the tree. This main additional information is the "group" to which each leaf on the phylogenetic tree belongs. Groups are typically the sample-type, such as: Normal, Primary, Locoregional, Peritoneum, Lung, Liver, and Distant (other). Group information is used for sample-type-specific analyses, as well as for quickly plotting standardized trees (labels colored by sample type and rooted at the Normal sample).

Some basic usage of "annotated_phylo" objects is as follows (assuming we have a distance matrix "dm"). NB: the row/column names in the distance matrix should be the standard "new" Naxerova lab naming scheme for samples: [Site][lesion][sample][autopsy]. For example, Liv1a has type='Liv', lesion=1, sample='a', autopsy=FALSE. Lun2-A is type='Lun', lesion=2, sample=NA, autopsy=TRUE.

```r
## use the group_samples() function to parse the sample names and define groups as desired
## here, lung samples are a distinct group from other distant-mets, but liver/peritoneum are combined with other distant mets.
groups <- group_samples(dm,color=T,lun=T,liv=F,per=F)

## generate the annotated_phylo object
tree <- annotated_phylo(dm, groups)

## print tree to see the information it contains (data from the angular distance matrix for C100):
R> tree

Annotated phylogenetic tree with 17 tips and 16 internal nodes.

Rooted; includes branch lengths.

    barcode order type lesion sample autopsy           group
 1:   Liv1a     1  Liv      1      a   FALSE Distant (other)
 2:   Liv1b     2  Liv      1      b   FALSE Distant (other)
 3:   Liv1c     3  Liv      1      c   FALSE Distant (other)
 4:   Liv2a     4  Liv      2      a   FALSE Distant (other)
 5:   Liv3a     5  Liv      3      a   FALSE Distant (other)
 6:   Liv3b     6  Liv      3      b   FALSE Distant (other)
 7:   Liv3c     7  Liv      3      c   FALSE Distant (other)
 8:   Liv3d     8  Liv      3      d   FALSE Distant (other)
 9:      L1     9    L      1          FALSE    Locoregional
10:      L2    10    L      2          FALSE    Locoregional
11:      L3    11    L      3          FALSE    Locoregional
12:      L4    12    L      4          FALSE    Locoregional
13:      P1    13    P      1          FALSE         Primary
14:      P2    14    P      2          FALSE         Primary
15:      P3    15    P      3          FALSE         Primary
16:      P4    16    P      4          FALSE         Primary
17:      N1    17    N      1          FALSE          Normal


## plot a standardized tree (uses 'ggtree' R-package)
plot(tree)
```






