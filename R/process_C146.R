source(here::here('R/func.R'))

subject <- 'C146'
sex <- 'male'

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

## B1b very difficult because < 1 copies of X and negative Y
refit(obj_list[[1]], samplename=samples[1], sex=sex, ploidy=2, purity=0.54, save=T, output_dir=here(paste0('output/',subject)))

## LN1: 0.5 gives Y=0, X=1, 3p and 4 exactly at 1.
refit(obj_list[[2]], samplename=samples[2], sex=sex, ploidy=2, purity=0.5, save=T, output_dir=here(paste0('output/',subject)))

## LN5: 
refit(obj_list[[3]], samplename=samples[3], sex=sex, ploidy=2, purity=0.55, save=T, output_dir=here(paste0('output/',subject)))

## N1
refit(obj_list[[4]], samplename=samples[4], sex=sex, ploidy=2, purity=1, save=T, output_dir=here(paste0('output/',subject)))

## P1
refit(obj_list[[5]], samplename=samples[5], sex=sex, ploidy=2, purity=0.57, save=T, output_dir=here(paste0('output/',subject)))

## P2
refit(obj_list[[6]], samplename=samples[6], sex=sex, ploidy=2, purity=0.4, save=T, output_dir=here(paste0('output/',subject)))

## P5
refit(obj_list[[7]], samplename=samples[7], sex=sex, ploidy=2, purity=0.37, save=T, output_dir=here(paste0('output/',subject)))

## P7
refit(obj_list[[8]], samplename=samples[8], sex=sex, ploidy=2, purity=0.47, save=T, output_dir=here(paste0('output/',subject)))

## P8
refit(obj_list[[9]], samplename=samples[9], sex=sex, ploidy=2, purity=0.49, save=T, output_dir=here(paste0('output/',subject)))

## TD2
refit(obj_list[[10]], samplename=samples[10], sex=sex, ploidy=2, purity=0.57, save=T, output_dir=here(paste0('output/',subject)))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# process cnv data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

info <- process_copynumber_data(obj_list, fit_file=here(paste0('output/',subject,'/fits/purity_ploidy.txt')), sex=sex, this.subject=subject, min_segment_bins=5,field='meancopy', R=10000, ncpus=4)

## make CNV segment heatmap
p <- cnv_heatmap(info$mat, info$seg, info$distance_matrix, this.subject=subject)
ggsave(here(paste0('output/',subject,'/',subject,'_cnv_segment_heatmap.pdf')),width=11,height=9)


