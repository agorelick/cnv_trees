source(here::here('R/func.R'))

subject <- 'C154'
sex <- 'female'

## because hacked to include X,Y, this is now as a list of 'obj', where each 'obj' has 1 sample
obj_list <- readRDS(here(paste0('original_data/',subject,'_1000kbp_withXYMT.rds'))) 
obj_list <- rename_samples(obj_list, subject)

## remove results from previous run
purity_file=here(paste0('output/',subject,'/fits/purity_ploidy.txt'))
if(file.exists(purity_file)) file.remove(purity_file)

samples <- names(obj_list)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# fit purity and ploidy for each sample
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## B1
refit(obj_list[[1]], samplename=samples[1], sex=sex, ploidy=2, purity=0.58, save=T, output_dir=here(paste0('output/',subject)))

## B3: strange, no longer X del, maybe this is from a separate part of the primary?
refit(obj_list[[2]], samplename=samples[2], sex=sex, ploidy=2, purity=0.5, save=T, output_dir=here(paste0('output/',subject)))

## B4:
refit(obj_list[[3]], samplename=samples[3], sex=sex, ploidy=2, purity=0.42, save=T, output_dir=here(paste0('output/',subject)))

## B5:
refit(obj_list[[4]], samplename=samples[4], sex=sex, ploidy=2, purity=0.51, save=T, output_dir=here(paste0('output/',subject)))

## B7:
refit(obj_list[[5]], samplename=samples[5], sex=sex, ploidy=2, purity=0.45, save=T, output_dir=here(paste0('output/',subject)))

## B8: 4q del is not clear
refit(obj_list[[6]], samplename=samples[6], sex=sex, ploidy=2, purity=0.29, save=T, output_dir=here(paste0('output/',subject)))

## B9:
refit(obj_list[[7]], samplename=samples[7], sex=sex, ploidy=2, purity=0.31, save=T, output_dir=here(paste0('output/',subject)))

## L11:
refit(obj_list[[8]], samplename=samples[8], sex=sex, ploidy=2, purity=0.65, save=T, output_dir=here(paste0('output/',subject)))

## L13:
refit(obj_list[[9]], samplename=samples[9], sex=sex, ploidy=2, purity=0.52, save=T, output_dir=here(paste0('output/',subject)))

## L14:
refit(obj_list[[10]], samplename=samples[10], sex=sex, ploidy=2, purity=0.57, save=T, output_dir=here(paste0('output/',subject)))

## L2a: X diploid, weird that L2b is different!
refit(obj_list[[11]], samplename=samples[11], sex=sex, ploidy=2, purity=0.66, save=T, output_dir=here(paste0('output/',subject)))

## L2b: subclonal X del
refit(obj_list[[12]], samplename=samples[12], sex=sex, ploidy=2, purity=0.30, save=T, output_dir=here(paste0('output/',subject)))

## L5
refit(obj_list[[13]], samplename=samples[13], sex=sex, ploidy=2, purity=0.65, save=T, output_dir=here(paste0('output/',subject)))

## N1: X is oddly low, could there be tumor contamination?
refit(obj_list[[14]], samplename=samples[14], sex=sex, ploidy=2, purity=1, save=T, output_dir=here(paste0('output/',subject)))

## P11: clonal X del, 8q dup, 4q del
refit(obj_list[[15]], samplename=samples[15], sex=sex, ploidy=2, purity=0.54, save=T, output_dir=here(paste0('output/',subject)))

## P12 clonal X del, 3 copies of 8q; 4q del
refit(obj_list[[16]], samplename=samples[16], sex=sex, ploidy=2, purity=0.55, save=T, output_dir=here(paste0('output/',subject)))

## P2
refit(obj_list[[17]], samplename=samples[17], sex=sex, ploidy=2, purity=0.25, save=T, output_dir=here(paste0('output/',subject)))

## P3: difficult to fit. P=0.48: X del clonal, but other events subclonal; P=0.4: fits not perfect but rounded give delX, del 4q, dup8q
refit(obj_list[[18]], samplename=samples[18], sex=sex, ploidy=2, purity=0.40, save=T, output_dir=here(paste0('output/',subject)))

## P4
refit(obj_list[[19]], samplename=samples[19], sex=sex, ploidy=2, purity=0.38, save=T, output_dir=here(paste0('output/',subject)))

## P7
refit(obj_list[[20]], samplename=samples[20], sex=sex, ploidy=2, purity=0.40, save=T, output_dir=here(paste0('output/',subject)))

## P8: 4q del, 8q dup, X subclonal deletion
refit(obj_list[[21]], samplename=samples[21], sex=sex, ploidy=2, purity=0.45, save=T, output_dir=here(paste0('output/',subject)))

## P9
refit(obj_list[[22]], samplename=samples[22], sex=sex, ploidy=2, purity=0.32, save=T, output_dir=here(paste0('output/',subject)))

## TD1
refit(obj_list[[23]], samplename=samples[23], sex=sex, ploidy=2, purity=0.64, save=T, output_dir=here(paste0('output/',subject)))



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
p1 <- test_tree_similarity(mat1,mat2,title=paste(subject,'SCNA segment vs Poly-G tree'),method='grf',nperm=10000)
ggsave(here(paste0('output/',subject,'/',subject,'_cnv_segment_polyG_tree_similarity_GRF.pdf')),width=7,height=4.5)



