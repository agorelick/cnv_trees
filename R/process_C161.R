source(here::here('R/func.R'))

subject <- 'C161'
sex <- 'female'

## because hacked to include X,Y, this is now as a list of 'obj', where each 'obj' has 1 sample
obj_list <- readRDS(here(paste0('original_data/',subject,'_1000kbp_withXYMT.rds'))) 
obj_list <- rename_samples(obj_list, subject)

## remove results from previous run
purity_file=here(paste0('output/',subject,'/fits/purity_ploidy.txt'))
if(file.exists(purity_file)) file.remove(purity_file)

samples <- names(obj_list)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# looking through the fits, we choose ploidy=4
# target homdel at 8p, 2 copies of 3 and 4q, 3 copies of 4p
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

samples <- names(obj_list)

## A1
refit(obj_list[[1]], samplename=samples[1], sex=sex, ploidy=4, purity=0.28, save=T, output_dir=here(paste0('output/',subject)))

## A2, odd, looks more like the rest with 0 8p and 2 4q, is this diff time from other lun?
refit(obj_list[[2]], samplename=samples[2], sex=sex, ploidy=4, purity=0.18, save=T, output_dir=here(paste0('output/',subject)))

## A5
refit(obj_list[[3]], samplename=samples[3], sex=sex, ploidy=4, purity=0.18, save=T, output_dir=here(paste0('output/',subject)))

## A6, odd, looks more like the rest with 0 8p and 2 4q, is this diff time from other lun?
refit(obj_list[[4]], samplename=samples[4], sex=sex, ploidy=4, purity=0.34, save=T, output_dir=here(paste0('output/',subject)))

## B1 
refit(obj_list[[5]], samplename=samples[5], sex=sex, ploidy=4, purity=0.21, save=T, output_dir=here(paste0('output/',subject)))

## B4 
refit(obj_list[[6]], samplename=samples[6], sex=sex, ploidy=4, purity=0.22, save=T, output_dir=here(paste0('output/',subject)))

## B6 
refit(obj_list[[7]], samplename=samples[7], sex=sex, ploidy=4, purity=0.17, save=T, output_dir=here(paste0('output/',subject)))

## B7
refit(obj_list[[8]], samplename=samples[8], sex=sex, ploidy=4, purity=0.16, save=T, output_dir=here(paste0('output/',subject)))

## H1
refit(obj_list[[9]], samplename=samples[9], sex=sex, ploidy=4, purity=0.24, save=T, output_dir=here(paste0('output/',subject)))

## L2
refit(obj_list[[10]], samplename=samples[10], sex=sex, ploidy=4, purity=0.44, save=T, output_dir=here(paste0('output/',subject)))

## L3
refit(obj_list[[11]], samplename=samples[11], sex=sex, ploidy=4, purity=0.23, save=T, output_dir=here(paste0('output/',subject)))

## ? L4
refit(obj_list[[12]], samplename=samples[12], sex=sex, ploidy=4, purity=0.30, save=T, output_dir=here(paste0('output/',subject)))

## Ld1 
refit(obj_list[[13]], samplename=samples[13], sex=sex, ploidy=4, purity=0.27, save=T, output_dir=here(paste0('output/',subject)))

## ? Ld2 
refit(obj_list[[14]], samplename=samples[14], sex=sex, ploidy=4, purity=0.17, save=T, output_dir=here(paste0('output/',subject)))

## N1
refit(obj_list[[15]], samplename=samples[15], sex=sex, ploidy=2, purity=1, save=T, output_dir=here(paste0('output/',subject)))

## O1
refit(obj_list[[16]], samplename=samples[16], sex=sex, ploidy=4, purity=0.23, save=T, output_dir=here(paste0('output/',subject)))

## O2
refit(obj_list[[17]], samplename=samples[17], sex=sex, ploidy=4, purity=0.28, save=T, output_dir=here(paste0('output/',subject)))

## P2
refit(obj_list[[18]], samplename=samples[18], sex=sex, ploidy=4, purity=0.32, save=T, output_dir=here(paste0('output/',subject)))

## P3
refit(obj_list[[19]], samplename=samples[19], sex=sex, ploidy=4, purity=0.33, save=T, output_dir=here(paste0('output/',subject)))

## P4
refit(obj_list[[20]], samplename=samples[20], sex=sex, ploidy=4, purity=0.41, save=T, output_dir=here(paste0('output/',subject)))

## ? P5: chr3 should be 2-3 copies
refit(obj_list[[21]], samplename=samples[21], sex=sex, ploidy=4, purity=0.28, save=T, output_dir=here(paste0('output/',subject)))

## P6
refit(obj_list[[22]], samplename=samples[22], sex=sex, ploidy=4, purity=0.37, save=T, output_dir=here(paste0('output/',subject)))

## P7
refit(obj_list[[23]], samplename=samples[23], sex=sex, ploidy=4, purity=0.32, save=T, output_dir=here(paste0('output/',subject)))

## P8
refit(obj_list[[24]], samplename=samples[24], sex=sex, ploidy=4, purity=0.38, save=T, output_dir=here(paste0('output/',subject)))

## P9
refit(obj_list[[25]], samplename=samples[25], sex=sex, ploidy=4, purity=0.43, save=T, output_dir=here(paste0('output/',subject)))

## ? Pa1 (either 8p is negative or 3, 4q, 10q all subclonal ~2.3 copies)
refit(obj_list[[26]], samplename=samples[26], sex=sex, ploidy=4, purity=0.22, save=T, output_dir=here(paste0('output/',subject)))

## S1
refit(obj_list[[27]], samplename=samples[27], sex=sex, ploidy=4, purity=0.30, save=T, output_dir=here(paste0('output/',subject)))

## SB1
refit(obj_list[[28]], samplename=samples[28], sex=sex, ploidy=4, purity=0.35, save=T, output_dir=here(paste0('output/',subject)))

## SB2
refit(obj_list[[29]], samplename=samples[29], sex=sex, ploidy=4, purity=0.28, save=T, output_dir=here(paste0('output/',subject)))

## St1
refit(obj_list[[30]], samplename=samples[30], sex=sex, ploidy=4, purity=0.21, save=T, output_dir=here(paste0('output/',subject)))

## St2
refit(obj_list[[31]], samplename=samples[31], sex=sex, ploidy=4, purity=0.33, save=T, output_dir=here(paste0('output/',subject)))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# process cnv data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

info <- get_bins_and_segments(obj_list, purity_file, sex)
samples <- names(obj_list)
info <- process_SCNA_data(samples, info, subject)

## make CNV segment heatmap
p <- cnv_heatmap(info$mat, info$seg, info$distance_matrix, this.subject=subject)
ggsave(here(paste0('output/',subject,'/',subject,'_cnv_segment_heatmap.pdf')),width=11,height=8)

## save plots comparing CNV 5+ Mb segment phylogeny to poly-G angular distance phylogeny
set.seed(42)
p2 <- compare_matrices(info$distance_matrix,subject, R=1e4)
ggsave(here(paste0('output/',subject,'/',subject,'_segment_euclidean_matrix_comparison.pdf')),width=7,height=6)
tree2.1 <- compare_trees(info$distance_matrix,subject, tree_method='nj')$plot
ggsave(here(paste0('output/',subject,'/',subject,'_segment_euclidean_nj_tree_comparison.pdf')),width=10,height=8)


## generate bootstrapped trees
set.seed(42)
chr_trees <- get_bootstrapped_trees(subject, type='chromosome')
set.seed(42)
seg_trees <- get_bootstrapped_trees(subject, type='segment')
set.seed(42)
polyg_trees <- get_bootstrapped_trees(subject, type='polyg')
out <- list(chr=chr_trees, seg=seg_trees, polyg=polyg_trees)
saveRDS(out,file=here(paste0('output/',subject,'/',subject,'_bootstrapped_trees.rds')))


## save plots comparing CNV (binned) phylogeny to poly-G angular distance phylogeny
bins <- info$bins[,c('sample','bin','copies'),with=F]
bins <- dcast(sample ~ bin, value.var='copies', data=bins)
bins <- d2m(bins)
bins <- as.matrix(dist(bins, method='euclidean'))

set.seed(42)
p2 <- compare_matrices(bins,subject, R=1e4)
ggsave(here(paste0('output/',subject,'/',subject,'_binned_euclidean_matrix_comparison.pdf')),width=7,height=6)
tree2.1 <- compare_trees(bins,subject, tree_method='nj')$plot
ggsave(here(paste0('output/',subject,'/',subject,'_binned_euclidean_nj_tree_comparison.pdf')),width=10,height=8)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# test tree similarity
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message('Testing tree similarity ...')
mat1 <- read_distance_matrix(here(paste0('output/',subject,'/',subject,'_cnv_distance_matrix.txt')))
mat2 <- read_distance_matrix(here(paste0('original_data/polyG/ad_matrix_oldnames/',subject,'_angular_dist_matrix_w_root_usedmarkers_repreReplicate_oldnames.txt')))
mat2 <- update_distance_matrix_samples(mat2)


## unnormalized
set.seed(42)
p1 <- test_tree_similarity(mat1,mat2,title=paste(subject,'SCNA segment vs Poly-G tree'),nperm=10000)
ggsave(here(paste0('output/',subject,'/',subject,'_cnv_segment_polyG_tree_similarity_GRF.pdf')),width=7,height=4.5)






