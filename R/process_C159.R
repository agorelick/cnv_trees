source(here::here('R/func.R'))

subject <- 'C159'
sex <- 'male'

## because hacked to include X,Y, this is now as a list of 'obj', where each 'obj' has 1 sample
obj_list <- readRDS(here(paste0('original_data/',subject,'_1000kbp_withXYMT.rds'))) 
obj_list <- rename_samples(obj_list, subject)

## remove results from previous run
purity_file=here(paste0('output/',subject,'/fits/purity_ploidy.txt'))
if(file.exists(purity_file)) file.remove(purity_file)

samples <- names(obj_list)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# looking through the fits, we choose ploidy=2 as most likely for this sample.
# The following events seem clonal (in every sample). We will use them to decide each sample's purity.
# 1q: 1 copy duplication
# 6q (before the end of chromosome): 1 copy deletion
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# clonal events based on primary tumors:
# - 4: del
# clonal in many but not all primary tumors:
# - 14q: del
# - 18: del (18q deep del subclonal in P1w)


## A1
refit(obj_list[[1]], samplename=samples[1], sex=sex, ploidy=3, purity=0.38, save=T, output_dir=here(paste0('output/',subject)))

## A1w
refit(obj_list[[2]], samplename=samples[2], sex=sex, ploidy=3, purity=0.30, save=T, output_dir=here(paste0('output/',subject)))

## A2a
refit(obj_list[[3]], samplename=samples[3], sex=sex, ploidy=3, purity=0.35, save=T, output_dir=here(paste0('output/',subject)))

## A2b
refit(obj_list[[4]], samplename=samples[4], sex=sex, ploidy=3, purity=0.17, save=T, output_dir=here(paste0('output/',subject)))

## A2w
refit(obj_list[[5]], samplename=samples[5], sex=sex, ploidy=3, purity=0.06, save=T, output_dir=here(paste0('output/',subject)))

## Ab1
refit(obj_list[[6]], samplename=samples[6], sex=sex, ploidy=3, purity=0.69, save=T, output_dir=here(paste0('output/',subject)))

## Ad1
refit(obj_list[[7]], samplename=samples[7], sex=sex, ploidy=3, purity=0.50, save=T, output_dir=here(paste0('output/',subject)))

## Ad2
refit(obj_list[[8]], samplename=samples[8], sex=sex, ploidy=3, purity=0.21, save=T, output_dir=here(paste0('output/',subject)))

## AR1w
refit(obj_list[[9]], samplename=samples[9], sex=sex, ploidy=3, purity=0.54, save=T, output_dir=here(paste0('output/',subject)))

## AR3w (0% tumor)
refit(obj_list[[10]], samplename=samples[10], sex=sex, ploidy=2, purity=1, save=T, output_dir=here(paste0('output/',subject)))

## B1
refit(obj_list[[11]], samplename=samples[11], sex=sex, ploidy=3, purity=0.37, save=T, output_dir=here(paste0('output/',subject)))

## B2w
refit(obj_list[[12]], samplename=samples[12], sex=sex, ploidy=3, purity=0.59, save=T, output_dir=here(paste0('output/',subject)))

## B3
refit(obj_list[[13]], samplename=samples[13], sex=sex, ploidy=3, purity=0.5, save=T, output_dir=here(paste0('output/',subject)))

## B3w
refit(obj_list[[14]], samplename=samples[14], sex=sex, ploidy=3, purity=0.45, save=T, output_dir=here(paste0('output/',subject)))

## B4
refit(obj_list[[15]], samplename=samples[15], sex=sex, ploidy=3, purity=0.58, save=T, output_dir=here(paste0('output/',subject)))

## D1w
refit(obj_list[[16]], samplename=samples[16], sex=sex, ploidy=3, purity=0.58, save=T, output_dir=here(paste0('output/',subject)))

## Di1
refit(obj_list[[17]], samplename=samples[17], sex=sex, ploidy=3, purity=0.75, save=T, output_dir=here(paste0('output/',subject)))

## Di2
refit(obj_list[[18]], samplename=samples[18], sex=sex, ploidy=3, purity=0.56, save=T, output_dir=here(paste0('output/',subject)))

## Di3
refit(obj_list[[19]], samplename=samples[19], sex=sex, ploidy=3, purity=0.29, save=T, output_dir=here(paste0('output/',subject)))

## H1w
refit(obj_list[[20]], samplename=samples[20], sex=sex, ploidy=3, purity=0.70, save=T, output_dir=here(paste0('output/',subject)))

## H3
refit(obj_list[[21]], samplename=samples[21], sex=sex, ploidy=3, purity=0.49, save=T, output_dir=here(paste0('output/',subject)))

## H6b
refit(obj_list[[22]], samplename=samples[22], sex=sex, ploidy=3, purity=0.59, save=T, output_dir=here(paste0('output/',subject)))

## H7a
refit(obj_list[[23]], samplename=samples[23], sex=sex, ploidy=3, purity=0.61, save=T, output_dir=here(paste0('output/',subject)))

## H7b
refit(obj_list[[24]], samplename=samples[24], sex=sex, ploidy=3, purity=0.69, save=T, output_dir=here(paste0('output/',subject)))

## H8
refit(obj_list[[25]], samplename=samples[25], sex=sex, ploidy=3, purity=0.43, save=T, output_dir=here(paste0('output/',subject)))

## LD1
refit(obj_list[[26]], samplename=samples[26], sex=sex, ploidy=3, purity=0.61, save=T, output_dir=here(paste0('output/',subject)))

## LD1w
refit(obj_list[[27]], samplename=samples[27], sex=sex, ploidy=3, purity=0.32, save=T, output_dir=here(paste0('output/',subject)))

## LD2
refit(obj_list[[28]], samplename=samples[28], sex=sex, ploidy=3, purity=0.54, save=T, output_dir=here(paste0('output/',subject)))

## LD2w
refit(obj_list[[29]], samplename=samples[29], sex=sex, ploidy=3, purity=0.13, save=T, output_dir=here(paste0('output/',subject)))

## LD3
refit(obj_list[[30]], samplename=samples[30], sex=sex, ploidy=3, purity=0.44, save=T, output_dir=here(paste0('output/',subject)))

## LD3w
refit(obj_list[[31]], samplename=samples[31], sex=sex, ploidy=3, purity=0.37, save=T, output_dir=here(paste0('output/',subject)))

## LD4
refit(obj_list[[32]], samplename=samples[32], sex=sex, ploidy=3, purity=0.47, save=T, output_dir=here(paste0('output/',subject)))

## LD5
refit(obj_list[[33]], samplename=samples[33], sex=sex, ploidy=3, purity=0.40, save=T, output_dir=here(paste0('output/',subject)))

## LN13
refit(obj_list[[34]], samplename=samples[34], sex=sex, ploidy=3, purity=0.12, save=T, output_dir=here(paste0('output/',subject)))

## LN1
refit(obj_list[[35]], samplename=samples[35], sex=sex, ploidy=3, purity=0.52, save=T, output_dir=here(paste0('output/',subject)))

## LN2
refit(obj_list[[36]], samplename=samples[36], sex=sex, ploidy=3, purity=0.44, save=T, output_dir=here(paste0('output/',subject)))

## LN4
refit(obj_list[[37]], samplename=samples[37], sex=sex, ploidy=3, purity=0.22, save=T, output_dir=here(paste0('output/',subject)))

## LN5
refit(obj_list[[38]], samplename=samples[38], sex=sex, ploidy=3, purity=0.73, save=T, output_dir=here(paste0('output/',subject)))

## LN6
refit(obj_list[[39]], samplename=samples[39], sex=sex, ploidy=3, purity=0.15, save=T, output_dir=here(paste0('output/',subject)))

## Me1
refit(obj_list[[40]], samplename=samples[40], sex=sex, ploidy=3, purity=0.70, save=T, output_dir=here(paste0('output/',subject)))

## N1
refit(obj_list[[41]], samplename=samples[41], sex=sex, ploidy=2, purity=1, save=T, output_dir=here(paste0('output/',subject)))

## P10
refit(obj_list[[42]], samplename=samples[42], sex=sex, ploidy=3, purity=0.79, save=T, output_dir=here(paste0('output/',subject)))

## P1w
refit(obj_list[[43]], samplename=samples[43], sex=sex, ploidy=3, purity=0.30, save=T, output_dir=here(paste0('output/',subject)))

## P2
refit(obj_list[[44]], samplename=samples[44], sex=sex, ploidy=3, purity=0.60, save=T, output_dir=here(paste0('output/',subject)))

## P2w
refit(obj_list[[45]], samplename=samples[45], sex=sex, ploidy=3, purity=0.62, save=T, output_dir=here(paste0('output/',subject)))

## P3
refit(obj_list[[46]], samplename=samples[46], sex=sex, ploidy=3, purity=0.54, save=T, output_dir=here(paste0('output/',subject)))

## P3w
refit(obj_list[[47]], samplename=samples[47], sex=sex, ploidy=3, purity=0.65, save=T, output_dir=here(paste0('output/',subject)))

## P4
refit(obj_list[[48]], samplename=samples[48], sex=sex, ploidy=3, purity=0.65, save=T, output_dir=here(paste0('output/',subject)))

## P4w
refit(obj_list[[49]], samplename=samples[49], sex=sex, ploidy=3, purity=0.82, save=T, output_dir=here(paste0('output/',subject)))

## P5
refit(obj_list[[50]], samplename=samples[50], sex=sex, ploidy=3, purity=0.55, save=T, output_dir=here(paste0('output/',subject)))

## P6
refit(obj_list[[51]], samplename=samples[51], sex=sex, ploidy=3, purity=0.55, save=T, output_dir=here(paste0('output/',subject)))

## P7
refit(obj_list[[52]], samplename=samples[52], sex=sex, ploidy=3, purity=0.45, save=T, output_dir=here(paste0('output/',subject)))

## P8
refit(obj_list[[53]], samplename=samples[53], sex=sex, ploidy=3, purity=0.57, save=T, output_dir=here(paste0('output/',subject)))

## P9
refit(obj_list[[54]], samplename=samples[54], sex=sex, ploidy=3, purity=0.41, save=T, output_dir=here(paste0('output/',subject)))

## Pa1
refit(obj_list[[55]], samplename=samples[55], sex=sex, ploidy=3, purity=0.60, save=T, output_dir=here(paste0('output/',subject)))

## Pg1
refit(obj_list[[56]], samplename=samples[56], sex=sex, ploidy=3, purity=0.38, save=T, output_dir=here(paste0('output/',subject)))

## S1w
refit(obj_list[[57]], samplename=samples[57], sex=sex, ploidy=3, purity=0.44, save=T, output_dir=here(paste0('output/',subject)))

## SB1
refit(obj_list[[58]], samplename=samples[58], sex=sex, ploidy=3, purity=0.37, save=T, output_dir=here(paste0('output/',subject)))

## Sp1
refit(obj_list[[59]], samplename=samples[59], sex=sex, ploidy=3, purity=0.38, save=T, output_dir=here(paste0('output/',subject)))

## TD1
refit(obj_list[[60]], samplename=samples[60], sex=sex, ploidy=3, purity=0.22, save=T, output_dir=here(paste0('output/',subject)))

## X1a
refit(obj_list[[61]], samplename=samples[61], sex=sex, ploidy=3, purity=0.38, save=T, output_dir=here(paste0('output/',subject)))

## X1b
refit(obj_list[[62]], samplename=samples[62], sex=sex, ploidy=3, purity=0.37, save=T, output_dir=here(paste0('output/',subject)))

## X1c
refit(obj_list[[63]], samplename=samples[63], sex=sex, ploidy=3, purity=0.55, save=T, output_dir=here(paste0('output/',subject)))

## X2w
refit(obj_list[[64]], samplename=samples[64], sex=sex, ploidy=3, purity=0.59, save=T, output_dir=here(paste0('output/',subject)))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# process cnv data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

info <- process_copynumber_data(obj_list, fit_file=here(paste0('output/',subject,'/fits/purity_ploidy.txt')), sex=sex, this.subject=subject, min_segment_bins=5,field='meancopy', R=10000, ncpus=4)

## make CNV segment heatmap
p <- cnv_heatmap(info$mat, info$seg, info$distance_matrix, this.subject=subject)
ggsave(here(paste0('output/',subject,'/',subject,'_cnv_segment_heatmap.pdf')),width=11,height=9)




