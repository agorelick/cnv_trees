source(here::here('R/func.R'))

subject <- 'E4'
sex <- 'female'

## because hacked to include X,Y, this is now as a list of 'obj', where each 'obj' has 1 sample
obj_list <- readRDS(here(paste0('original_data/',subject,'_1000kbp_withXYMT.rds'))) 
obj_list <- rename_samples(obj_list, subject)

## remove results from previous run
purity_file=here(paste0('output/',subject,'/purity_ploidy.txt'))
if(file.exists(purity_file)) file.remove(purity_file)

samples <- names(obj_list)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# use refit function to fit the purity and ploidy for each sample, generating PDFs
# and saving the purity/ploidy to a table as you go
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# B11
refit(obj_list[[1]], samplename=samples[1], sex=sex, ploidy=2, purity=0.34, save=T, output_dir=here(paste0('output/',subject)))

## B19
refit(obj_list[[2]], samplename=samples[2], sex=sex, ploidy=2, purity=0.28, save=T, output_dir=here(paste0('output/',subject)))

## B6
refit(obj_list[[3]], samplename=samples[3], sex=sex, ploidy=2, purity=0.26, save=T, output_dir=here(paste0('output/',subject)))

## H1 
refit(obj_list[[4]], samplename=samples[4], sex=sex, ploidy=2, purity=0.67, save=T, output_dir=here(paste0('output/',subject)))

## L3 
refit(obj_list[[5]], samplename=samples[5], sex=sex, ploidy=2, purity=0.62, save=T, output_dir=here(paste0('output/',subject)))

## N1 
refit(obj_list[[6]], samplename=samples[6], sex=sex, ploidy=2, purity=1, save=T, output_dir=here(paste0('output/',subject)))

## P11 
refit(obj_list[[7]], samplename=samples[7], sex=sex, ploidy=2, purity=0.32, save=T, output_dir=here(paste0('output/',subject)))

## P14 
refit(obj_list[[8]], samplename=samples[8], sex=sex, ploidy=2, purity=0.53, save=T, output_dir=here(paste0('output/',subject)))

## P14 
refit(obj_list[[9]], samplename=samples[9], sex=sex, ploidy=2, purity=0.61, save=T, output_dir=here(paste0('output/',subject)))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# process cnv data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

info <- process_copynumber_data(obj_list, fit_file=here(paste0('output/',subject,'/purity_ploidy.txt')), sex=sex, this.subject=subject, min_segment_bins=5,field='meancopy', R=10000, ncpus=4)

## make CNV segment heatmap
p <- cnv_heatmap(info$mat, info$seg, info$distance_matrix, this.subject=subject)
ggsave(here(paste0('output/',subject,'/',subject,'_cnv_segment_heatmap.pdf')),width=11,height=9)


