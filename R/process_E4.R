rm(list=ls())
source(here::here('R/func.R'))

subject <- 'E4'
sex <- 'female'

## because hacked to include X,Y, this is now as a list of 'obj', where each 'obj' has 1 sample
obj_list <- readRDS(here(paste0('original_data/',subject,'_1000kbp_withXYMT.rds'))) 
obj_list <- rename_samples(obj_list, subject)

## remove results from previous run
purity_file=here(paste0('output/',subject,'/fits/purity_ploidy.txt'))
if(file.exists(purity_file)) file.remove(purity_file)

samples <- names(obj_list)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# use refit function to fit the purity and ploidy for each sample, generating PDFs
# and saving the purity/ploidy to a table as you go
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Per8
refit(obj_list[[1]], samplename=samples[1], sex=sex, ploidy=2, purity=0.34, save=T, output_dir=here(paste0('output/',subject)))

## Per11
refit(obj_list[[2]], samplename=samples[2], sex=sex, ploidy=2, purity=0.28, save=T, output_dir=here(paste0('output/',subject)))

## Per4
refit(obj_list[[3]], samplename=samples[3], sex=sex, ploidy=2, purity=0.26, save=T, output_dir=here(paste0('output/',subject)))

## Liv1
refit(obj_list[[4]], samplename=samples[4], sex=sex, ploidy=2, purity=0.67, save=T, output_dir=here(paste0('output/',subject)))

## TD1
refit(obj_list[[5]], samplename=samples[5], sex=sex, ploidy=2, purity=0.62, save=T, output_dir=here(paste0('output/',subject)))

## N1 
refit(obj_list[[6]], samplename=samples[6], sex=sex, ploidy=2, purity=1, save=T, output_dir=here(paste0('output/',subject)))

## P7 
refit(obj_list[[7]], samplename=samples[7], sex=sex, ploidy=2, purity=0.32, save=T, output_dir=here(paste0('output/',subject)))

## P8 
refit(obj_list[[8]], samplename=samples[8], sex=sex, ploidy=2, purity=0.53, save=T, output_dir=here(paste0('output/',subject)))

## P3 
refit(obj_list[[9]], samplename=samples[9], sex=sex, ploidy=2, purity=0.61, save=T, output_dir=here(paste0('output/',subject)))

## Per9 
refit(obj_list[[10]], samplename=samples[10], sex=sex, ploidy=2, purity=0.13, save=T, output_dir=here(paste0('output/',subject)))

## Per18 
refit(obj_list[[11]], samplename=samples[11], sex=sex, ploidy=2, purity=0.15, save=T, output_dir=here(paste0('output/',subject)))

## Per17
refit(obj_list[[12]], samplename=samples[12], sex=sex, ploidy=2, purity=0.19, save=T, output_dir=here(paste0('output/',subject)))

## Per10
refit(obj_list[[13]], samplename=samples[13], sex=sex, ploidy=2, purity=0.18, save=T, output_dir=here(paste0('output/',subject)))

## Per21
refit(obj_list[[14]], samplename=samples[14], sex=sex, ploidy=2, purity=0.10, save=T, output_dir=here(paste0('output/',subject)))

## Per12
refit(obj_list[[15]], samplename=samples[15], sex=sex, ploidy=2, purity=0.22, save=T, output_dir=here(paste0('output/',subject)))

## Per6
refit(obj_list[[16]], samplename=samples[16], sex=sex, ploidy=2, purity=0.20, save=T, output_dir=here(paste0('output/',subject)))

## Per14
refit(obj_list[[17]], samplename=samples[17], sex=sex, ploidy=2, purity=0.19, save=T, output_dir=here(paste0('output/',subject)))

## Per1
refit(obj_list[[18]], samplename=samples[18], sex=sex, ploidy=2, purity=0.24, save=T, output_dir=here(paste0('output/',subject)))

## Per2
refit(obj_list[[19]], samplename=samples[19], sex=sex, ploidy=2, purity=0.24, save=T, output_dir=here(paste0('output/',subject)))

## Per3
refit(obj_list[[20]], samplename=samples[20], sex=sex, ploidy=2, purity=0.19, save=T, output_dir=here(paste0('output/',subject)))

## Per15
refit(obj_list[[21]], samplename=samples[21], sex=sex, ploidy=2, purity=0.21, save=T, output_dir=here(paste0('output/',subject)))

## Per16
refit(obj_list[[22]], samplename=samples[22], sex=sex, ploidy=2, purity=0.25, save=T, output_dir=here(paste0('output/',subject)))

## Per19
refit(obj_list[[23]], samplename=samples[23], sex=sex, ploidy=2, purity=0.18, save=T, output_dir=here(paste0('output/',subject)))

## L5
refit(obj_list[[24]], samplename=samples[24], sex=sex, ploidy=2, purity=0.23, save=T, output_dir=here(paste0('output/',subject)))

## TD3
refit(obj_list[[25]], samplename=samples[25], sex=sex, ploidy=2, purity=0.13, save=T, output_dir=here(paste0('output/',subject)))

## L6
refit(obj_list[[26]], samplename=samples[26], sex=sex, ploidy=2, purity=0.12, save=T, output_dir=here(paste0('output/',subject)))

## L1
refit(obj_list[[27]], samplename=samples[27], sex=sex, ploidy=2, purity=0.43, save=T, output_dir=here(paste0('output/',subject)))

## L2
refit(obj_list[[28]], samplename=samples[28], sex=sex, ploidy=2, purity=0.18, save=T, output_dir=here(paste0('output/',subject)))

## TD4
refit(obj_list[[29]], samplename=samples[29], sex=sex, ploidy=2, purity=0.11, save=T, output_dir=here(paste0('output/',subject)))

## L3
refit(obj_list[[30]], samplename=samples[30], sex=sex, ploidy=2, purity=0.35, save=T, output_dir=here(paste0('output/',subject)))

## L4
refit(obj_list[[31]], samplename=samples[31], sex=sex, ploidy=2, purity=0.26, save=T, output_dir=here(paste0('output/',subject)))

## P1
refit(obj_list[[32]], samplename=samples[32], sex=sex, ploidy=2, purity=0.35, save=T, output_dir=here(paste0('output/',subject)))

## P9
refit(obj_list[[33]], samplename=samples[33], sex=sex, ploidy=2, purity=0.11, save=T, output_dir=here(paste0('output/',subject)))

## P4
refit(obj_list[[34]], samplename=samples[34], sex=sex, ploidy=2, purity=0.40, save=T, output_dir=here(paste0('output/',subject)))

## P5
refit(obj_list[[35]], samplename=samples[35], sex=sex, ploidy=2, purity=0.42, save=T, output_dir=here(paste0('output/',subject)))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# process cnv data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## remove any problematic samples based on manual curation
bad_samples <- c("Per9","Per21","TD3","L6","L2","TD4")
good_samples <- names(obj_list)[!names(obj_list) %in% bad_samples]
fit_file=here(paste0('output/',subject,'/fits/purity_ploidy.txt'))
fits <- fread(fit_file)
fits <- fits[barcode %in% good_samples,]
write_tsv(fits, fit_file)
obj_list <- obj_list[fits$barcode] ## names must match fits$barcode


info <- get_bins_and_segments(obj_list, fit_file, sex)

manual_override_chrX <- c('Per17','Per10','Per12','Per16','Per15','Per2','P7','Per6','P1','L3','L5','L1','Per19','Per18','P9','Per14','Per1')

override_X_bins <- function(tmp) {
    tmp$copynumbers[tmp$chr=='X'] <- 2
    tmp$segments[tmp$chr=='X'] <- 2
    tmp$adjustedcopynumbers[tmp$chr=='X'] <- 2
    tmp$intcopy[tmp$chr=='X'] <- 2
    tmp$meancopy[tmp$chr=='X'] <- 2
    tmp
}
override_X_segs <- function(tmp) {
    tmp$Copies[tmp$Chromosome=='X'] <- 2
    tmp$Segment_Mean[tmp$Chromosome=='X'] <- 2
    tmp$Segment_Mean2[tmp$Chromosome=='X'] <- 2
    tmp$Segment_SE[tmp$Chromosome=='X'] <- NA
    tmp$P_log10[tmp$Chromosome=='X'] <- NA
    tmp
}

for(this.sample in manual_override_chrX) {
    message(this.sample)
    info$bins[[this.sample]] <- override_X_bins(info$bins[[this.sample]])
    info$segments[[this.sample]] <- override_X_segs(info$segments[[this.sample]])
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# process cnv data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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



## save the bootstrap values for to a table
p <- bootstrap_cnv_tree(info$mat, B=1000, this.subject=subject,collapse_threshold=0)
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




## bootstrap SCNA tree
set.seed(42)
for(bs in seq(0,95,by=5)) {
    message(bs)
    if(bs < 10) {
        bslab <- paste0('0',bs)
    } else {
        bslab <- as.character(bs)
    }
    p <- bootstrap_cnv_tree(info$mat, B=1000, this.subject=subject,collapse_threshold=bs)
    ggsave(here(paste0('output/',subject,'/collapsed_trees/',subject,'_cnv_segment_tree_bootstrapped_',bslab,'.pdf')),width=9,height=8)
}


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


