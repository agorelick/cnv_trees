source(here::here('R/func.R'))

subject <- 'E4'
sex <- 'female'

## because hacked to include X,Y, this is now as a list of 'obj', where each 'obj' has 1 sample
obj_list <- readRDS(here(paste0('original_data/',subject,'_1000kbp_withXYMT.rds'))) 
obj_list <- rename_samples(obj_list, subject)

## remove results from previous run
purity_file=here(paste0('output/',subject,'/fits/purity_ploidy.txt'))
#if(file.exists(purity_file)) file.remove(purity_file)

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

#fits <- fread(here(paste0('output/',subject,'/fits/purity_ploidy.txt')))
#fits <- fits[1:9,]
#fits <- fits[purity >= 0.2,]
#obj_list <- obj_list[fits$barcode]
#write_tsv(fits,here(paste0('output/',subject,'/fits/purity_ploidy.txt')))

info <- process_copynumber_data(obj_list, fit_file=here(paste0('output/',subject,'/fits/purity_ploidy.txt')), sex=sex, this.subject=subject, min_segment_bins=5,field='meancopy', R=10000, ncpus=4)

## make CNV segment heatmap
p <- cnv_heatmap(info$mat, info$seg, info$distance_matrix, this.subject=subject)
ggsave(here(paste0('output/',subject,'/',subject,'_cnv_segment_heatmap.pdf')),width=11,height=8)

#x <- read_distance_matrix(here('output/E4/results_all_samples/E4_cnv_matrix.txt'))

## bootstrap SCNA tree
p <- bootstrap_cnv_tree(info$mat, B=1000, this.subject=subject)
ggsave(here(paste0('output/',subject,'/',subject,'_cnv_segment_tree_bootstrapped.pdf')),width=9,height=8)






run_scrap=F

if(run_scrap==T) {

    pdf('~/Desktop/E4_allsamples_bootstrapped.pdf',width=8,height=10)
    plot(tree, main = "E4 SCNA tree (all sWGS samples) with 1000 bootstrap values\nno subclonal SCNAs, no chrX")
    drawSupportOnEdges(boot)
    dev.off()



    #fits <- fread(here(paste0('output/',subject,'/fits/purity_ploidy.txt')))
    #fits <- fits[purity >= 0.2,]
    #x <- x[fits$barcode,]
    x <- round(x)
    rows <- rownames(x)
    rng <- range(x)
    values <- rng[1]:rng[2]
    x <- apply(x, 2, as.character)
    rownames(x) <- rows
    segments <- colnames(x)
    segments <- segments[!grepl('X',segments)]
    x <- x[,segments]

    ## get maximum parsimony tree
    set.seed(42)
    pd <- as.phyDat(x,levels=as.character(values),type='USER',ambiguity='-')
    tree <- pratchet(pd, trace=0)
    tree <- acctran(tree, pd)
    tree <- phytools::reroot(tree, node.number=grep('N1',tree$tip.label))

    #dm <- as.matrix(distTips(tree,method='nNodes'))
    dm <- as.matrix(distTips(tree,method='patristic'))
    ad <- read_distance_matrix(here('original_data/polyG/E4_ad_matrix.txt'))
    common_samples <- intersect(rownames(dm),rownames(ad))
    common_samples <- common_samples[!grepl('^N',common_samples)]
    sample_levels <- sort(common_samples, decreasing=T)
    dm_subset <- dm[common_samples,common_samples]
    ad_subset <- ad[common_samples,common_samples]
    tst <- dist_similarity(test_dist=dm_subset, ref_dist=ad_subset, nperm=1000, cpus=4, return_only_pval=F, method='spearman')

    this.subject <- 'E4'
    ## heatmap
    d_subset <- as.data.table(reshape2::melt(tst$test_dist))
    d_subset$data <- 'CNV'
    d_subset$Var1 <- factor(d_subset$Var1, levels=sample_levels)
    d_subset$Var2 <- factor(d_subset$Var2, levels=sample_levels)
    g_subset <- as.data.table(reshape2::melt(tst$ref_dist))
    g_subset$data <- 'poly-G'
    g_subset$Var1 <- factor(g_subset$Var1, levels=sample_levels)
    g_subset$Var2 <- factor(g_subset$Var2, levels=sample_levels)

    p1 <- ggplot(d_subset, aes(x=Var1,y=Var2)) +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0)) +
        theme_ang(base_size=10) +
        geom_tile(aes(fill=value)) +
        theme(
              axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(), 
              legend.position='right') +
scale_fill_gradient(low='white',high='steelblue',name='Euclidean\ndistance') +
labs(x=NULL,y=NULL,subtitle=paste(this.subject,'SCNA distance matrix'))

    p2 <- ggplot(g_subset, aes(x=Var1,y=Var2)) +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0)) +
        theme_ang(base_size=10) +
        geom_tile(aes(fill=value)) +
        theme(
              axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(), 
              legend.position='right') +
scale_fill_gradient(low='white',high='steelblue',name='Angular\ndistance') +
labs(x=NULL,y=NULL,subtitle=paste(this.subject,'Polyguanine distance matrix'))
    p_top <- plot_grid(p1, p2, ncol=2)
p_space <- ggplot() + theme_nothing()
p_bottom <- plot_grid(p_space, tst$plot, p_space, nrow=1, rel_widths=c(1,3,1))
p <- plot_grid(p_top, p_bottom, ncol=1, rel_heights=c(1,1.2))
p


common_samples <- intersect(rownames(dm),rownames(ad))
dm_subset <- dm[common_samples,common_samples]
ad_subset <- ad[common_samples,common_samples]

## define groups for plot
groups <- group_samples(dm_subset,color=T,lun=F,liv=F,per=T)

## CNV tree
tree1 <- annotated_phylo(dm_subset, groups)
tree$colors <- tree1$colors
tree$tip.annotations <- tree1$tip.annotations
tree$distance.matrix <- tree1$distance.matrix
class(tree) <- 'annotated_phylo'
tree$edge.length <- NULL

p1 <- plot(tree, legend.position='none', angle=T) + labs(subtitle='SCNA maximum parsimony tree')

## poly-G tree
tree2 <- annotated_phylo(ad_subset, groups)
p2 <- plot(tree2, legend.position='none', angle=T) + labs(subtitle='Polyguanine angular distance')

p <- plot_grid(p1,p2,ncol=2)
ggsave('~/Desktop/E4_maxpars_vs_ad_tree_min20purity.pdf',width=10,height=8)

}



