subject <- 'E4'
sex <- 'female'

obj_list <- readRDS(paste0('original_data/',subject,'_1000kbp_withXYMT.rds')) ## because hacked to include X,Y, this is now as a list of 'obj', where each 'obj' has 1 sample
names(obj_list) <- gsub(subject,'',gsub('_aligned','',names(obj_list)))

## remove results from previous run
if(file.exists('refits/purity_ploidy.txt')) file.remove('refits/purity_ploidy.txt')

samples <- names(obj_list)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# use refit function to fit the purity and ploidy for each sample, generating PDFs
# and saving the purity/ploidy to a table as you go
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## B11
refit(obj_list[[1]], samplename=samples[1], sex='female', ploidy=2, purity=0.34, save=T, output_dir=paste0('output/',subject))

## B19
refit(obj_list[[2]], samplename=samples[2], sex='female', ploidy=2, purity=0.28, save=T, output_dir=paste0('output/',subject))

## B6
refit(obj_list[[3]], samplename=samples[3], sex='female', ploidy=2, purity=0.26, save=T, output_dir=paste0('output/',subject))

## H1 
refit(obj_list[[4]], samplename=samples[4], sex='female', ploidy=2, purity=0.67, save=T, output_dir=paste0('output/',subject))

## L3 
refit(obj_list[[5]], samplename=samples[5], sex='female', ploidy=2, purity=0.62, save=T, output_dir=paste0('output/',subject))

## N1 
refit(obj_list[[6]], samplename=samples[6], sex='female', ploidy=2, purity=1, save=T, output_dir=paste0('output/',subject))

## P11 
refit(obj_list[[7]], samplename=samples[7], sex='female', ploidy=2, purity=0.32, save=T, output_dir=paste0('output/',subject))

## P14 
refit(obj_list[[8]], samplename=samples[8], sex='female', ploidy=2, purity=0.53, save=T, output_dir=paste0('output/',subject))

## P14 
refit(obj_list[[9]], samplename=samples[9], sex='female', ploidy=2, purity=0.61, save=T, output_dir=paste0('output/',subject))



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



