source(here::here('R/func.R'))

subject <- 'C186'
sex <- 'male'

## because hacked to include X,Y, this is now as a list of 'obj', where each 'obj' has 1 sample
obj_list <- readRDS(here(paste0('original_data/',subject,'_1000kbp_withXYMT.rds'))) 
obj_list <- rename_samples(obj_list, subject)

## remove results from previous run
purity_file=here(paste0('output/',subject,'/fits/purity_ploidy.txt'))
if(file.exists(purity_file)) file.remove(purity_file)

samples <- names(obj_list)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# looking through the fits, we choose ploidy=4
# 1 copy deletion in 1q and 12
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## A1
refit(obj_list[[1]], samplename=samples[1], sex=sex, ploidy=4, purity=0.33, save=T, output_dir=here(paste0('output/',subject)))

## A2
refit(obj_list[[2]], samplename=samples[2], sex=sex, ploidy=4, purity=0.35, save=T, output_dir=here(paste0('output/',subject)))

## A3
refit(obj_list[[3]], samplename=samples[3], sex=sex, ploidy=4, purity=0.36, save=T, output_dir=here(paste0('output/',subject)))

## A4
refit(obj_list[[4]], samplename=samples[4], sex=sex, ploidy=4, purity=0.41, save=T, output_dir=here(paste0('output/',subject)))

## A5
refit(obj_list[[5]], samplename=samples[5], sex=sex, ploidy=4, purity=0.36, save=T, output_dir=here(paste0('output/',subject)))

## CA1 
refit(obj_list[[6]], samplename=samples[6], sex=sex, ploidy=4, purity=0.58, save=T, output_dir=here(paste0('output/',subject)))

## CA2
refit(obj_list[[7]], samplename=samples[7], sex=sex, ploidy=4, purity=0.42, save=T, output_dir=here(paste0('output/',subject)))

## M10
refit(obj_list[[8]], samplename=samples[8], sex=sex, ploidy=4, purity=0.75, save=T, output_dir=here(paste0('output/',subject)))

## M3
refit(obj_list[[9]], samplename=samples[9], sex=sex, ploidy=4, purity=0.50, save=T, output_dir=here(paste0('output/',subject)))

## M5
refit(obj_list[[10]], samplename=samples[10], sex=sex, ploidy=4, purity=0.65, save=T, output_dir=here(paste0('output/',subject)))

## M6
refit(obj_list[[11]], samplename=samples[11], sex=sex, ploidy=4, purity=0.71, save=T, output_dir=here(paste0('output/',subject)))

## M7
refit(obj_list[[12]], samplename=samples[12], sex=sex, ploidy=4, purity=0.41, save=T, output_dir=here(paste0('output/',subject)))

## N3
refit(obj_list[[13]], samplename=samples[13], sex=sex, ploidy=2, purity=1, save=T, output_dir=here(paste0('output/',subject)))

## P1
refit(obj_list[[14]], samplename=samples[14], sex=sex, ploidy=4, purity=0.53, save=T, output_dir=here(paste0('output/',subject)))

## P2
refit(obj_list[[15]], samplename=samples[15], sex=sex, ploidy=4, purity=0.55, save=T, output_dir=here(paste0('output/',subject)))

## P3
refit(obj_list[[16]], samplename=samples[16], sex=sex, ploidy=4, purity=0.60, save=T, output_dir=here(paste0('output/',subject)))

## P7 (chr12 at 3 copies)
refit(obj_list[[17]], samplename=samples[17], sex=sex, ploidy=4, purity=0.30, save=T, output_dir=here(paste0('output/',subject)))


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

## bootstrap SCNA tree
p <- bootstrap_cnv_tree(info$mat, B=1000, this.subject=subject,collapse_threshold=0)
ggsave(here(paste0('output/',subject,'/',subject,'_cnv_segment_tree_bootstrapped.pdf')),width=9,height=8)


## save the bootstrap values for to a table
bs <- p$data$bootstrap
bs <- bs[!is.na(bs)]
bs_vals <- data.table(subject=subject,bs=bs)
write_tsv(bs_vals, here(paste0('output/',subject,'/',subject,'_segment_euclidean_nj_tree_bootstrap_values.txt')))

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
mat2 <- read_distance_matrix(here(paste0('original_data/polyG/',subject,'_ad_matrix.txt')))

## unnormalized
set.seed(42)
p1 <- test_tree_similarity(mat1,mat2,title=paste(subject,'SCNA segment vs Poly-G tree'),nperm=10000)
ggsave(here(paste0('output/',subject,'/',subject,'_cnv_segment_polyG_tree_similarity.pdf')),width=7,height=4.5)

## normalized
set.seed(42)
p2 <- test_tree_similarity(mat1,mat2,title=paste(subject,'SCNA segment vs Poly-G tree'),nperm=10000,normalize=T)
ggsave(here(paste0('output/',subject,'/',subject,'_cnv_segment_polyG_tree_similarity_normalized.pdf')),width=7,height=4.5)


