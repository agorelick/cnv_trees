subject <- 'C157'
sex <- 'male'

obj_list <- readRDS(paste0('original_data/',subject,'_1000kbp_withXYMT.rds')) ## because hacked to include X,Y, this is now as a list of 'obj', where each 'obj' has 1 sample
names(obj_list) <- gsub(subject,'',gsub('_aligned','',names(obj_list)))

## remove results from previous run
if(file.exists(paste0('output/',subject,'/purity_ploidy.txt'))) file.remove(paste0('output/',subject,'/purity_ploidy.txt'))

samples <- names(obj_list)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# fit purity and ploidy for each sample
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## B1
refit(obj_list[[1]], samplename=samples[1], sex='male', ploidy=2, purity=0.34, save=T, output_dir=paste0('output/',subject))

## B2
refit(obj_list[[2]], samplename=samples[2], sex='male', ploidy=2, purity=0.38, save=T, output_dir=paste0('output/',subject))

## Di1
refit(obj_list[[3]], samplename=samples[3], sex='male', ploidy=2, purity=0.25, save=T, output_dir=paste0('output/',subject))

## Di2
refit(obj_list[[4]], samplename=samples[4], sex='male', ploidy=2, purity=0.24, save=T, output_dir=paste0('output/',subject))

## LN4
refit(obj_list[[5]], samplename=samples[5], sex='male', ploidy=2, purity=0.43, save=T, output_dir=paste0('output/',subject))

## N1
refit(obj_list[[6]], samplename=samples[6], sex='male', ploidy=2, purity=1, save=T, output_dir=paste0('output/',subject))

## P10
refit(obj_list[[7]], samplename=samples[7], sex='male', ploidy=2, purity=0.38, save=T, output_dir=paste0('output/',subject))

## P1
refit(obj_list[[8]], samplename=samples[8], sex='male', ploidy=2, purity=0.52, save=T, output_dir=paste0('output/',subject))

## P2
refit(obj_list[[9]], samplename=samples[9], sex='male', ploidy=2, purity=0.53, save=T, output_dir=paste0('output/',subject))

## P3 
refit(obj_list[[10]], samplename=samples[10], sex='male', ploidy=2, purity=0.48, save=T, output_dir=paste0('output/',subject))

## P4
refit(obj_list[[11]], samplename=samples[11], sex='male', ploidy=2, purity=0.60, save=T, output_dir=paste0('output/',subject))

## P6
refit(obj_list[[12]], samplename=samples[12], sex='male', ploidy=2, purity=0.49, save=T, output_dir=paste0('output/',subject))

## P9
refit(obj_list[[13]], samplename=samples[13], sex='male', ploidy=2, purity=0.24, save=T, output_dir=paste0('output/',subject))

## TD1
refit(obj_list[[14]], samplename=samples[14], sex='male', ploidy=2, purity=0.50, save=T, output_dir=paste0('output/',subject))

## TD2a
refit(obj_list[[15]], samplename=samples[15], sex='male', ploidy=2, purity=0.53, save=T, output_dir=paste0('output/',subject))

## TD2b
refit(obj_list[[16]], samplename=samples[16], sex='male', ploidy=2, purity=0.51, save=T, output_dir=paste0('output/',subject))

## TD3
refit(obj_list[[17]], samplename=samples[17], sex='male', ploidy=2, purity=0.76, save=T, output_dir=paste0('output/',subject))

## TD7
refit(obj_list[[18]], samplename=samples[18], sex='male', ploidy=2, purity=0.47, save=T, output_dir=paste0('output/',subject))




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# process cnv data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

info <- get_cnv_segments(obj_list, fit_file=paste0('output/',subject,'/purity_ploidy.txt'), sex=sex, min_segment_bins=5,field='meancopy')

mat <- info$mat
mat <- cbind(sample=rownames(mat), as.data.table(mat))
write_tsv(mat,paste0('output/',subject,'_cnv_matrix.txt'))

segments <- info$segs
write_tsv(segments,paste0('output/',subject,'_cnv_segments.txt'))

bins <- info$bins
write_tsv(bins,paste0('output/',subject,'_cnv_bins.txt'))




