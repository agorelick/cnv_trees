rm(list=ls())

source(here::here('R/func.R'))


## define helper functions
plot_tree <- function(m) {
    set.seed(42)
    pd <- as.phyDat(m,levels=as.character(0:10),type='USER',ambiguity='-')
    dm <- dist.hamming(pd,exclude='pairwise')
    x <- as.matrix(dm)
    hammingtree <- NJ(dm)
    tree <- pratchet(pd, start=hammingtree)
    tree <- acctran(tree, pd)
    tree
}

plot_heatmap <- function(mat2, chr, xpct=0.025, title, ylab, na.color='black') {
    chr2 <- copy(chr)
    chr2 <- chr2[chr %in% c(1:18,'X','Y'),c('chr','midpoint'),with=F]
    chr2 <- rbind(chr2, data.table(midpoint=mean(chr$midpoint[chr$chr %in% c(19:22)]),chr='19-22'))
    chr2 <- chr2[order(midpoint),]
    sample_levels <- levels(mat2$cell_IDs)
    nrows <- length(sample_levels)
    xpos <- -1*xpct*nrows
    dels <- c('#2d89e5','#99ccff')
    names(dels) <- c('0','1')
    neutral <- 'white'
    names(neutral) <- '2'
    amps <- c(brewer.pal(8,'YlOrRd'))
    names(amps) <- as.character(3:10)
    cols <- c(dels,neutral,amps)

    p_heatmap <- ggplot(mat2) + 
        scale_x_continuous(breaks=c(0,chr$chr_end), label=c('',as.character(chr$chr)), expand=c(0,0)) + 
        scale_y_continuous(breaks=1:length(sample_levels),
                           labels=sample_levels,expand=c(0,0),position='right') + 
        geom_rect(aes(xmin=start,xmax=end,
                      ymin=samplenum-0.5,ymax=samplenum+0.5,fill=copy.number),size=0.125) + 
        geom_text(data=chr2,aes(x=midpoint,label=chr),y=xpos,size=3) +
        scale_fill_manual(values=cols,name='Copies',na.value=na.color) + 
        labs(title=title,
             x='\nGenomic position',y=ylab) +
        theme_ang(base_size=12) +
        theme(axis.text.x=element_blank(),axis.line.x=element_blank(),axis.line.y=element_blank(),
              axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank()) + 
    geom_vline(xintercept=c(0,chr$chr_end),size=0.5) +
    coord_cartesian(clip = "off")

    p_heatmap
}

collapse_chr <- function(chr) {
    ## for each chr (or chr-arm), get the start/end positions
    chr_start <- min(chr$start)
    chr_end <- max(chr$end)
    midpoint <- (chr_start + chr_end)/2
    list(chr_start=chr_start,chr_end=chr_end,midpoint=midpoint)
}

parse_segment <- function(pos) {
    s <- strsplit(pos,'[:]')[[1]]
    chr <- s[1]
    s <- strsplit(s[2],'[-]')[[1]]
    start <- as.integer(s[1])
    end <- as.integer(s[2])
    list(chr=chr,start=start,end=end) 
}

get_character_matrix <- function(dat) {
    ## convert it to a matrix of character copy-number [0-10], '-' means NA; rows are cells, columns are segments
    ## this will be used to construct maximum parsimony tree
    dat <- dat[order(chr,start,end),]
    dat$segment <- factor(dat$segment, levels=unique(dat$segment))
    m <- dat[,c('sample','cell_IDs','segment','copy.number'),with=F]
    rng <- range(m$copy.number,na.rm=T)
    values <- as.character(rng[1]:rng[2])
    m$copy.number <- as.character(m$copy.number)
    m <- d2m(dcast(cell_IDs ~ segment, value.var='copy.number', data=m))
    m[is.na(m)] <- '-'
    m
}


get_chr <- function(cytobands_or_available_bins) {
    ## for plotting, we use the cumulative sum of each chromosome's position along the complete genome
    ## if available_bins is provided, then we need to reduce the chromosome lengths to the total retained length (this is necessary for the position of the chromosome labels and borders to be correct in the heatmap)

    chr <- cytobands_or_available_bins[,collapse_chr(.SD),by=c('chr')]
    chr$chr <- factor(chr$chr, levels=c(1:22,'X')) 
    chr <- chr[order(chr),]
    chr$chr_start <- chr$chr_start / 1e6
    chr$chr_end <- chr$chr_end / 1e6
    chr$midpoint <- chr$midpoint / 1e6
    chr$chr_end <- cumsum(chr$chr_end) 
    for(i in 2:nrow(chr)) {
        chr$chr_start[i] <- chr$chr_end[i-1]
        chr$midpoint[i] <- chr$midpoint[i] + chr$chr_end[i-1]
    }
    chr
}


prep_data_for_heatmap <- function(m, cytobands, shrink_to_available_bins=F) { 
    ## create a data.table object 'dat' which will be used for plotting the heatmap
    dat <- t(m)
    l <- lapply(as.character(rownames(dat)),parse_segment)
    toadd2 <- rbindlist(l)
    dat <- cbind(toadd2, dat)

    ## get the chromosome start/end positons
    if(shrink_to_available_bins==T) {
        chr <- get_chr(dat[,c('chr','start','end'),with=F])
    } else {
        chr <- get_chr(cytobands)   
    }

    ## add chromosome arm to the bin data
    arms <- cytobands[,collapse_chr(.SD),by=c('chr','arm')]
    arms[,midpoint:=NULL]
    arms[,chr_start:=chr_start / 1e6]
    arms[,chr_end:=chr_end / 1e6]
    dat$start <- (dat$start / 1e6)
    dat$end <- (dat$end / 1e6)
    setkey(dat,'chr','start','end')
    setkey(arms,'chr','chr_start','chr_end')
    dat <- foverlaps(dat, arms, type='within')
    dat[,c('chr_start','chr_end'):=NULL]
    dat$chr <- factor(dat$chr, levels=c(1:22,'X'))

    ## add a new segment field for the Mb-scaled coordinates
    dat[,segment:=paste0(chr,':',round(start),'-',round(end))]
    dat <- dat[order(chr,start,end),]

    ## add the cumulative genomic positions and chromosome-arm to the copy number data
    dat <- merge(dat, chr[,c('chr','chr_start'),with=F], by='chr', all.x=T)
    dat[,start:=start + chr_start]
    dat[,end:=end + chr_start]
    dat[,chr_start:=NULL]
    dat <- melt(dat, id.vars=c('segment','chr','arm','start','end'))
    setnames(dat,c('variable','value'),c('cell_IDs','copy.number'))
   
    ## return both the bin data and the chr object 
    list(dat=dat,chr=chr)
}


# ~~~~~~~~

## load the unfiltered aneufinder output, add segment ID
ew035 <- fread(here('original_data/scDNA/EW035_bins_allsamples_unfiltered.txt.gz'))
setnames(ew035,'seqnames','chr')
ew035$chr <- factor(ew035$chr, levels=c(1:22,'X'))
ew035[,segment:=paste0(chr,':',start,'-',end)]
cell_counts1 <- as.data.frame(table(ew035$sample[!duplicated(ew035$cell_IDs)]))

## map of fresh-frozen sample names to FFPE old/new sample names
sample_info <- fread(here('original_data/scDNA/E4_FF_FFPE_conversion.txt'))

## get the start/end positions of each chromosome (we use the hg19 cytobands to also get the start/end positions of the individual chromosome arms)
cytobands <- fread(here('original_data/cytoBand.txt'))
names(cytobands) <- c('chr','start','end','band','info')
cytobands$arm <- strtrim(cytobands$band,1)
cytobands[,chr:=gsub('chr','',chr)]
cytobands$chr <- factor(cytobands$chr, levels=c(1:22,'X'))
cytobands <- cytobands[!is.na(chr),]
cytobands <- cytobands[order(chr),]

## generate matrix of cell vs bin with character copy number '0',...'10',('-'=NA)
m <- get_character_matrix(ew035) 

## format SCNA data for the heatmap
heatmap_data <- prep_data_for_heatmap(m, cytobands, shrink_to_available_bins=F)
dat <- heatmap_data$dat
chr <- heatmap_data$chr

## plot the maximum parsimony tree (takes ~1 min) using all cells/segments and the accompanying heatmap
tree <- plot_tree(m)
p_tree <- ggtree(tree, layout='rect')

## plot the heatmap with cells aligned to the tree
plotdat <- as.data.table(p_tree$data)
plotdat <- plotdat[isTip==T]
plotdat <- plotdat[order(y,decreasing=T),]
dat$cell_IDs <- factor(dat$cell_IDs, levels=rev(plotdat$label))
dat$samplenum <- as.integer(dat$cell_IDs)
p_heatmap <- plot_heatmap(dat, chr, title='E4 single cell SCNA matrix and maximum parsimony tree',ylab=paste0('Cells (N=',nrow(m),')'),na.color='#bfbfbf')
p <- plot_grid(p_tree, p_heatmap, align='h', rel_widths=c(1,6))
ggsave(here('output/E4_singlecell/heatmap_unfiltered.pdf'),width=10,height=7)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# based on the above heatmap, define clonal CNVs and
# then subset for cells with at least one of these
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

check_for_clonal_events <- function(d) {
    dels <- c('0','1')
    dups <- as.character(c(3:10))
    dup1q <- any(d$copy.number %in% dups & d$charm=='1q',na.rm=T)
    dup6q <- any(d$copy.number %in% dups & d$charm=='6q',na.rm=T)
    del6q <- any(d$copy.number %in% dels & d$charm=='6q',na.rm=T)
    del12q <- any(d$copy.number %in% dels & d$charm=='12q',na.rm=T)
    dup15 <- any(d$copy.number %in% dups & d$chr=='15',na.rm=T)
    bins_avail <- sum(d$copy.number!='-')
    frac_cna_bins <- sum(d$copy.number %in% c(dels,dups)) / bins_avail
    frac_NA_bins <- mean(d$copy.number=='-')
    list(dup1q=dup1q,dup6q=dup6q,del6q=del6q,del12q=del12q,dup15=dup15,
         frac_cna_bins=frac_cna_bins,frac_NA_bins=frac_NA_bins)
}
d <- copy(dat)
d[,charm:=paste0(chr,arm)]
res <- d[,check_for_clonal_events(.SD),by=c('cell_IDs')]
res$clonal_events <- rowSums(res[,(2:6),with=F])

## add tissue-type and sample name annotations
extract_sample <- function(cell_ID) {
    str <- strsplit(cell_ID,'_')[[1]]
    str[2]
}
res$sample <- sapply(as.character(res$cell_IDs), extract_sample)
res <- merge(res, sample_info[,c('sample','tissue'),with=F], by='sample', all.x=T)
res[sample=='R',tissue:='primary']
res[sample=='Q',tissue:='ovarian']

## subset for cancer cells based on those with at least 3/5 of the above clonal SCNAs
res <- res[clonal_events >= 3,]

## exclude extremely noisy cells according to if they either have (1) an outlier fraction of bins with SCNAs, or (2) an outlier fraction of bins with missing copy.number data. Detect outliers based on modified Z-score as described by B. Iglewicz and D.C. Hoaglin, How to Detect and Handle Outliers (American Society for Quality Control, Milwaukee, WI, 1993).
modified_zscore <- 0.6745*(res$frac_cna_bins - median(res$frac_cna_bins)) / mad(res$frac_cna_bins)
res$outlier_frac_SCNA <- F
res[abs(modified_zscore) > 3.5,outlier_frac_SCNA:=T]
modified_zscore <- 0.6745*(res$frac_NA_bins - median(res$frac_NA_bins)) / mad(res$frac_NA_bins)
res$outlier_frac_NA <- F
res[abs(modified_zscore) > 3.5,outlier_frac_NA:=T]
res[outlier_frac_SCNA==F & outlier_frac_NA==F,status:='Retain']
res[outlier_frac_SCNA==T & outlier_frac_NA==F,status:='% SCNA bin outlier']
res[outlier_frac_SCNA==F & outlier_frac_NA==T,status:='% NA bin outlier']
res[outlier_frac_SCNA==T & outlier_frac_NA==T,status:='Both % SCNA and NA bin outlier']
tbl <- table_freq(res$status)
res <- merge(res, tbl, by.x='status', by.y='value', all.x=T)
res[,status:=paste0(status,' (n=',N,')')]
p <- ggplot(res,aes(x=frac_NA_bins,y=frac_cna_bins)) +
    geom_point(aes(fill=status),pch=21,color='black',stroke=0.25) +
    theme_ang(base_size=12) +
    scale_fill_brewer(palette='Accent',name='Cell status') +
    labs(x='Frac of bins with missing copy number',y='Frac bins with SCNAs')
ggsave(here('output/E4_singlecell/cell_outlier_status_scatterplot.pdf'),width=7,height=3.5)

res <- res[outlier_frac_SCNA==F & outlier_frac_NA==F,]
res$cell_IDs <- as.character(res$cell_IDs)
xtabs(~clonal_events+tissue,data=res)

## add number of filtered cells for each sample to the cell_counts table and write the table
tbl <- table_freq(res$sample)
cell_counts2 <- merge(cell_counts1, tbl, by.x='Var1', by.y='value', all.x=T)
names(cell_counts2) <- c('sample','n_unfiltered','n_filtered')
cell_counts2 <- as.data.table(cell_counts2)
cell_counts2[is.na(n_filtered),n_filtered:=0]
cell_counts2 <- cell_counts2[order(sample),]
out <- merge(sample_info, cell_counts2, by='sample', all.x=T)
write_tsv(out,here('output/E4_singlecell/cell_counts.txt'))

## get an updated heatmap and parsimony tree after first round of filtering
tree2 <- plot_tree(m[rownames(m) %in% res$cell_IDs,])
p_tree2 <- ggtree(tree2, layout='rect')
plotdat <- as.data.table(p_tree2$data)
plotdat <- plotdat[isTip==T]
plotdat <- plotdat[order(y,decreasing=T),]
dat_filtered <- dat[cell_IDs %in% plotdat$label,]
dat_filtered$cell_IDs <- factor(dat_filtered$cell_IDs, levels=rev(plotdat$label))
dat_filtered$samplenum <- as.integer(dat_filtered$cell_IDs)
p_heatmap2 <- plot_heatmap(dat_filtered, chr, title='Filtered E4 single cell SCNA matrix and maximum parsimony tree',ylab=paste0('Cells (N=',nrow(res),')'), na.color='#bfbfbf')
p2 <- plot_grid(p_tree2, p_heatmap2, align='h', rel_widths=c(1,6))
ggsave(here('output/E4_singlecell/heatmap_filtered_partialexample.pdf'),width=10,height=7)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# using the filtered cells, get the euclidean distance trees
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dat_filtered <- m[rownames(m) %in% res$cell_IDs,]
rows <- rownames(dat_filtered)
dat_filtered <- apply(dat_filtered,2,as.integer)
rownames(dat_filtered) <- rows
dat_filtered <- reshape2::melt(dat_filtered)
dat_filtered <- as.data.table(dat_filtered)
names(dat_filtered) <- c('cell_IDs','segment','copy.number')
parse_id <- function(s) {
    s <- strsplit(s,'_')[[1]]
    list(subject=s[1],sample=s[2],cell=s[3])
}
l <- rbindlist(lapply(as.character(dat_filtered$cell_IDs), parse_id))
dat_filtered <- cbind(l, dat_filtered)

## add tissue-type and sample name annotations
dat_filtered <- merge(dat_filtered, sample_info, by='sample', all.x=T)
dat_filtered[sample=='R',sample_name:='primary tumor']
dat_filtered[sample=='R',oldname:='E4P-R']
dat_filtered[sample=='R',newname:='E4P-R']
dat_filtered[sample=='R',tissue:='primary']
dat_filtered[sample=='R',color:='#008c45']

## get map of cell to sample (before O1,2,3 -> O)
cells_samples <- dat_filtered[!duplicated(cell_IDs),c('cell_IDs','sample'),with=F]

## combine O1,O2,O3 together (because same sample) for imputation. This will not exist after this section and shall be redone.
dat_filtered[sample %in% c('O1','O2','O3'),sample:='O']

## first, remove any frequently-missing segments (NA in 20%+ of cells)
prop_with_NA <- function(d) {
    propNA=mean(is.na(d$copy.number))
    list(propNA=propNA)
}
bins <- dat_filtered[,prop_with_NA(.SD),by=c('segment')]
bins_with_recurrent_NAs <- as.character(bins$segment[bins$propNA >= 0.2])
dat_filtered <- dat_filtered[!segment %in% bins_with_recurrent_NAs,]
dat_filtered$segment <- factor(dat_filtered$segment, levels=unique(dat_filtered$segment))

## now within each sample, we impute NA values based on the nearest neighbor (within same tissue) where it is available.
impute_within_type <- function(this.tissue, data) {
    message(this.tissue)

    ## first subset the data for cells from this tissue type
    d <- data[sample %in% this.tissue]
    d <- data.table::dcast(data=d, segment ~ cell_IDs, value.var='copy.number')
    d <- d2m(d)

    ## get distance between each pair of samples (allow NAs) 
    nn <- (as.matrix(dist(t(d),diag=F,method='euclidean',upper=T)))
    nn <- cbind(cell=rownames(nn), as.data.table(nn))

    ## extract the cell/segment combos we need to impute
    d_noNA <- d[!is.na(rowSums(d)),]
    toimpute <- as.data.table(reshape2::melt(d[is.na(rowSums(d)),]))
    names(toimpute) <- c('segment','cell','cn')

    ## annotate the cell/segments to be imputed with the cells' distance to all other cells
    toimpute <- merge(toimpute, nn, by='cell', all.x=T)
    toimpute <- melt(toimpute, id.vars=c('cell','segment','cn'))
    setnames(toimpute,c('variable','value'),c('neighbor','distance'))
    toimpute <- toimpute[cell!=neighbor,]

    ## for each neighboring cell, merge in its CN value for the given segment
    toadd <- melt(cbind(segment=rownames(d),as.data.table(d)),id.vars='segment')
    toimpute <- merge(toimpute, toadd, by.x=c('neighbor','segment'), by.y=c('variable','segment'), all.x=T)

    ## finally, for any cell/segment combos with missing values, impute them with the value for the closest neighbor for which that value is available.
    impute <- function(tmp) {
        tmp <- tmp[order(distance,decreasing=F),]
        if(any(!is.na(tmp$cn))) {
            new_value <- unique(tmp$cn)
        } else {
            tmp <- tmp[!is.na(value),]
            if(nrow(tmp) > 0) {
                new_value <- tmp$value[tmp$distance==min(tmp$distance)]
                if(length(new_value) > 1) new_value=sample(new_value,1)
            } else {
                new_value <- as.numeric(NA)
            }
        }
        list(new_value=as.numeric(new_value))
    }
    imputed <- toimpute[,impute(.SD),by=c('cell','segment')]

    ## merge the cell/segments that did not require imputation with the imputed portion
    imputed <- d2m(dcast(segment ~ cell, value.var='new_value', data=imputed))
    d_imputed <- rbind(d_noNA, imputed)
    d_imputed <- d_imputed[levels(data$segment),]
    d_imputed <- cbind(tissue=this.tissue, segment=rownames(d_imputed), as.data.table(d_imputed))
    d_imputed <- melt(d_imputed, id.vars=c('tissue','segment'))
}
tbl <- table_freq(dat_filtered$sample[!duplicated(dat_filtered$cell_IDs)])
imputed_tissues <- tbl$value[tbl$N > 1]
unimputed_tissues <- tbl$value[tbl$N == 1]
l <- lapply(imputed_tissues, impute_within_type,dat_filtered)
imputed <- rbindlist(l)
setnames(imputed,c('variable','value'),c('cell_IDs','copy.number'))
setnames(imputed,'tissue','sample')
unimputed <- dat_filtered[sample %in% unimputed_tissues,c('sample','segment','cell_IDs','copy.number'),with=F]
combined <- rbind(imputed, unimputed)

## finally, remove any segments that were unimputable
segments <- combined[,prop_with_NA(.SD),by=c('segment')]
bins_unimputable <- as.character(segments$segment[segments$propNA > 0])
combined <- combined[!segment %in% bins_unimputable,]

## add in any lost bins as totally missing
bad_bins <- unique(c(bins_unimputable,bins_with_recurrent_NAs))
to_expand <- combined[!duplicated(cell_IDs),]
to_expand[,segment:=NULL]
expand_bin <- function(bad_bin, to_expand) {
    to_expand$segment <- bad_bin
    to_expand$copy.number <- as.integer(NA)
    to_expand
}
bad_bins_expanded <- rbindlist(lapply(bad_bins, expand_bin, to_expand))
combined <- rbind(combined, bad_bins_expanded)

#imputed <- copy(dat_filtered)
combined[,sample:=NULL]
## parse the segments and reformat as a wide matrix
l <- rbindlist(lapply(as.character(combined$segment), parse_segment))
combined <- cbind(combined, l)
combined <- merge(combined, cells_samples, by='cell_IDs', all.x=T)
m_filtered_imputed <- get_character_matrix(combined)

## format SCNA data for the heatmap
filtered_imputed_heatmap_data <- prep_data_for_heatmap(m_filtered_imputed, cytobands, shrink_to_available_bins=F)
filtered_imputed_dat <- filtered_imputed_heatmap_data$dat
filtered_imputed_chr <- filtered_imputed_heatmap_data$chr

## plot the maximum parsimony tree (takes ~1 min) using all cells/segments and the accompanying heatmap
filtered_imputed_tree <- plot_tree(m_filtered_imputed)
p_filtered_imputed_tree <- ggtree(filtered_imputed_tree, layout='rect')

## plot the heatmap with cells aligned to the tree
plotdat <- as.data.table(p_filtered_imputed_tree$data)
plotdat <- plotdat[isTip==T]
plotdat <- plotdat[order(y,decreasing=T),]
filtered_imputed_dat$cell_IDs <- factor(filtered_imputed_dat$cell_IDs, levels=rev(plotdat$label))
filtered_imputed_dat$samplenum <- as.integer(filtered_imputed_dat$cell_IDs)
p_filtered_imputed_heatmap <- plot_heatmap(filtered_imputed_dat, chr, title='Filtered, imputed E4 single cell SCNA matrix and maximum parsimony tree',ylab=paste0('Cells (N=',nrow(m_filtered_imputed),')'),na.color='#bfbfbf')

p <- plot_grid(p_filtered_imputed_tree, p_filtered_imputed_heatmap, align='h', rel_widths=c(1,6))
ggsave(here('output/E4_singlecell/heatmap_filtered_imputed.pdf'),width=10,height=7)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get the euclidean distance tree and root it at the normal (diploid pseudo-cell)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

treedat <- copy(combined)
treedat <- treedat[order(cell_IDs,chr,start,end),]
treedat$segment <- factor(treedat$segment, levels=unique(treedat$segment))
x <- dcast(cell_IDs ~ segment, value.var='copy.number', data=treedat)
x <- d2m(x)
root <- rep(2,ncol(x))
x <- rbind(x, root)
filtered_output <- cbind(cell_IDs=rownames(x), as.data.table(x))
write_tsv(filtered_output,here('output/E4_singlecell/EW035_bin_matrix_filtered_imputed.txt'))

dm <- dist(x, method='euclidean')
dm <- as.matrix(dm)
finaltree <- nj(dm)
finaltree <- phytools::reroot(finaltree, node.number=grep('root',finaltree$tip.label))

## info is data.table containing tissue types for plotting with ggtree
info <- data.table(label=finaltree$tip.label)
info$pos <- 1:nrow(info)
info <- merge(info, dat_filtered[,c('cell_IDs','tissue','sample','sample_name','oldname','newname','color'),with=F], by.x='label', by.y='cell_IDs', all.x=T)
info[label=='root',tissue:='normal']
info[label=='root',sample:='Normal']
info[label=='root',sample_name:='Normal']
info[label=='root',oldname:='Normal']
info[label=='root',newname:='Normal']
info[,sample_label:=paste0(oldname,' (',newname,') ',sample_name)]
info[label=='root',sample_label:='Normal']
info[label=='root',color:='black']

## plot euclidean distance rooted tree highlighting cells from each sample individually
plot_tree_for_sample <- function(this.oldname, info, finaltree) {
    plot_title <- unique(info$sample_label[info$oldname==this.oldname & !is.na(info$sample_label)])
    dot_color <- unique(info$color[info$oldname==this.oldname & !is.na(info$sample_label)])
    p <- ggtree(finaltree,layout='ape',color='grey',size=0.5)
    p <- p %<+% info
    tmp <- as.data.frame(p$data)
    tmp <- tmp[tmp$oldname==this.oldname,]
    n_cells <- length(unique(tmp$label[!is.na(tmp$label)]))
    p <- p + geom_rootpoint(color="black", size=2) 
    p <- p + geom_tippoint(size=2, color=dot_color, data=tmp) + labs(title=paste0(plot_title,'; N=',n_cells,' cells'))
    pdf_file <- here(paste0('output/E4_singlecell/finaltree_',gsub('[/]','-',this.oldname),'.pdf'))
    ggsave(pdf_file,width=7,height=6)
    message(pdf_file)
}
samples <- unique(info$oldname)
samples <- samples[samples!='Normal']
trash <- lapply(samples, plot_tree_for_sample, info, finaltree)

## color based on the samples' tissue type
tissue_colors <- info[!duplicated(tissue),]
cols <- tissue_colors$color
names(cols) <- tissue_colors$tissue
p <- ggtree(finaltree,layout='ape',color='grey',size=0.5)
p <- p %<+% info
p <- p + geom_rootpoint(color="black", size=2) 
p <- p + geom_tippoint(size=2, pch=16, aes(color=tissue))
p <- p + scale_color_manual(values=cols,name='Tissue type')
p <- p + labs(title='E4 filtered cells (colored by sample tissue type)')
ggsave(here('output/E4_singlecell/finaltree_allsamples.pdf'),width=10,height=6)


## color based on individual sample
#labels <- unique(info$sample_label)
#per_samples <- grep('Per',labels,value=T)
#other_samples <- labels[!grepl('Per',labels)]
#cols_peritoneal <- c(brewer.pal(8,'Accent'),'white')
#names(cols_peritoneal) <- per_samples
#cols_other <- c(brewer.pal(7,'Set1'),'black')
#names(cols_other) <- other_samples
#
#tmp_peritoneal <- tmp[tmp$sample_label %in% per_samples,]
#tmp_other <- tmp[!tmp$sample_label %in% per_samples,]
#p <- p + geom_tippoint(data=tmp_other, size=2, pch=19, aes(color=sample_label)) 
#p <- p + geom_tippoint(data=tmp_peritoneal, size=2, pch=21, aes(fill=sample_label)) 
#p <- p + scale_fill_manual(values=cols_peritoneal,name='Peritoneal samples') 
#p <- p + scale_color_manual(values=cols_other,name='Other samples') 
#p <- p + labs(title='E4 filtered cells (all samples)')
#ggsave(here('output/E4_singlecell/finaltree_allsamples.pdf'),width=10,height=6)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generate pseudo-bulk CNV tree based on the mean/median copy number per bin in
# each sample, compare this to the poly-G tree. NB: the samples won't perfectly overlap
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d <- fread(here('output/E4_singlecell/EW035_bin_matrix_filtered_imputed.txt'))
get_sample <- function(s) strsplit(s,'_')[[1]][2]
d$sample <- sapply(d$cell_IDs, get_sample)
d[cell_IDs=='root',sample:='Normal']
d[,cell_IDs:=NULL]
d <- merge(d, sample_info[,c('sample','oldname'),with=F], by='sample', all.x=T)
d[sample=='Normal',oldname:='E4N1']
d <- d[!is.na(oldname),]

L1L2 <- d[oldname=='E4L1/L2',]
L1 <- copy(L1L2)
L2 <- copy(L1L2)
L1$oldname <- 'E4L1'
L2$oldname <- 'E4L2'

L3L4 <- d[oldname=='E4L3/L4',]
L3 <- copy(L3L4)
L4 <- copy(L3L4)
L3$oldname <- 'E4L3'
L4$oldname <- 'E4L4'

P1P2 <- d[oldname=='E4P1/P2',]
P1 <- copy(P1P2)
P2 <- copy(P1P2)
P1$oldname <- 'E4P1'
P2$oldname <- 'E4P2'

P3P14 <- d[oldname=='E4P3/P14',]
P3 <- copy(P3P14)
P14 <- copy(P3P14)
P3$oldname <- 'E4P3'
P14$oldname <- 'E4P14'

d <- d[!oldname %in% c('E4L1/L2','E4L3/L4','E4P1/P2','E4P3/P14')]
d <- rbind(d,L1,L2,L3,L4,P1,P2,P3,P14)

tbl <- table_freq(d$oldname)
samples_with_enough_cells <- tbl$value[tbl$N >= 1 | tbl$value=='E4N1']
d <- d[oldname %in% samples_with_enough_cells]



## load the angular distance matrix
ad <- fread(here('original_data/polyG/E4_ad_matrix.txt'))
ad <- melt(ad,id.vars='V1')
names(ad) <- c('s1','s2','ad')
info <- fread(here('original_data/sample_info.txt'))
info <- info[subject=='E4',c('oldname','barcode'),with=F]
ad <- merge(ad, info, by.x='s1', by.y='barcode', all.x=T)
setnames(ad,'oldname','oldname1')
ad <- merge(ad, info, by.x='s2', by.y='barcode', all.x=T)
setnames(ad,'oldname','oldname2')


## only consider samples common to both the poly-G and scDNA data
valid_samples <- c(intersect(d$oldname, ad$oldname1))


## get the angular distance matrix with the oldnames
ad <- ad[oldname1 %in% valid_samples & oldname2 %in% valid_samples,]
ad_dm <- dcast(oldname1 ~ oldname2, value.var='ad', data=ad)
ad_dm <- d2m(ad_dm)
ad_dm <- ad_dm[valid_samples,valid_samples]


## convert the sc data to euclidean distance matrix also with oldnames
d <- d[oldname %in% c(valid_samples),]
d[,sample:=NULL]


get_mean <- function(d) {
    #x <- apply(d,2,median,na.rm=T)
    #x <- as.numeric(x)
    x <- colMeans(d)    
    x <- as.list(x)
    names(x) <- names(d)
    x
}
d <- d[,get_mean(.SD),by=oldname]
sc_mat <- copy(d)
rows <- sc_mat$oldname
sc_mat[,c('oldname'):=NULL]
sc_mat <- as.matrix(sc_mat)
rownames(sc_mat) <- rows
sc_dm <- as.matrix(dist(sc_mat,method='euclidean'))
sc_dm <- sc_dm[valid_samples,valid_samples]


## load the angular distance matrix
info <- fread(here('original_data/sample_info.txt'))
info <- info[subject=='E4',c('oldname','barcode'),with=F]
tmp <- data.table(oldname=valid_samples)
tmp$pos <- 1:nrow(tmp)
tmp <- merge(tmp, info, by='oldname', all.x=T)
tmp <- tmp[order(pos),]
rownames(sc_dm) <- tmp$barcode
colnames(sc_dm) <- tmp$barcode
rownames(ad_dm) <- tmp$barcode
colnames(ad_dm) <- tmp$barcode

sc_groups <- group_samples(rownames(sc_dm),color=T)
sc_tree <- annotated_phylo(sc_dm, sc_groups)
ad_groups <- group_samples(rownames(ad_dm),color=T)
ad_tree <- annotated_phylo(ad_dm, ad_groups)
sc_plot <- plot(sc_tree)
ad_plot <- plot(ad_tree)
p <- plot_grid(sc_plot, ad_plot, ncol=1)
ggsave(here('output/E4_singlecell/pseudobulk_cna_vs_polyG_tree.pdf'),width=8,height=10)












## the single cell data includes merges of lesions individually profiled by poly-G, so let's split them up into identical lesions to better compare them in the trees.
single_samples <- d[!grepl('[/]',oldname),]
single_samples[,sample:=NULL]
segments <- names(d)[!names(d) %in% c('sample','oldname')]

## we will average together O1/O2/O3, which correspond to L3/L4:
L3L4 <- d[sample %in% c('O1','O2','O3'),c(segments),with=F]
L3L4 <- as.data.table(get_mean(L3L4))
L3 <- copy(L3L4)
L3$oldname <- 'E4L3'
L4 <- copy(L3L4)
L4$oldname <- 'E4L4'

## we will duplicated sample S, which correspond to P1/P2:
P1P2 <- d[sample %in% 'S',c(segments),with=F]
P1 <- copy(P1P2)
P1$oldname <- 'E4P1'
P2 <- copy(P1P2)
P2$oldname <- 'E4P2'

## we will duplicated sample T, which correspond to P3/P14:
P3P14 <- d[sample %in% 'T',c(segments),with=F]
P3 <- copy(P3P14)
P3$oldname <- 'E4P3'
P14 <- copy(P3P14)
P14$oldname <- 'E4P14'

## merge the new duplicated samples with the single-sample data
sc <- rbind(single_samples, L3, L4, P1, P2, P3, P14)

## remove samples with too few cells
tbl <- table_freq(sc$oldname)
valid_samples <- tbl$value[tbl$N >= 1 | tbl$value=='E4N1']
sc <- sc[oldname %in% valid_samples,]

get_mean <- function(d) {
    x <- apply(d,2,median,na.rm=T)
    x <- as.numeric(x)
    #x <- colMeans(d)    
    x <- as.list(x)
    names(x) <- names(d)
    x
}
res <- sc[,get_mean(.SD),by=oldname]

## load the angular distance matrix
ad <- fread(here('original_data/polyG/E4_ad_matrix.txt'))
ad <- melt(ad,id.vars='V1')
names(ad) <- c('s1','s2','ad')
info <- fread(here('original_data/sample_info.txt'))
info <- info[subject=='E4',c('oldname','barcode'),with=F]
ad <- merge(ad, info, by.x='s1', by.y='barcode', all.x=T)
setnames(ad,'oldname','oldname1')
ad <- merge(ad, info, by.x='s2', by.y='barcode', all.x=T)
setnames(ad,'oldname','oldname2')



## convert the sc data to euclidean distance matrix
sc_mat <- copy(sc)
rows <- sc_mat$oldname
sc_mat[,oldname:=NULL]
sc_mat <- as.matrix(sc_mat)
rownames(sc_mat) <- rows
sc_dm <- as.matrix(dist(sc_mat,method='euclidean'))

## get the angular distance matrix with the oldnames
ad_dm <- dcast(oldname1 ~ oldname2, value.var='ad', data=ad)
ad_dm <- d2m(ad_dm)

library(polyG)
valid_samples <- intersect(rownames(sc_dm),rownames(ad_dm))
sc_dm <- sc_dm[valid_samples,valid_samples]
ad_dm <- ad_dm[valid_samples,valid_samples]

tmp <- data.table(oldname=valid_samples)
tmp$pos <- 1:nrow(tmp)
tmp <- merge(tmp, info, by='oldname', all.x=T)
tmp <- tmp[order(pos),]
rownames(sc_dm) <- tmp$barcode
colnames(sc_dm) <- tmp$barcode
rownames(ad_dm) <- tmp$barcode
colnames(ad_dm) <- tmp$barcode

tst <- dist_similarity(test_dist=sc_dm, ref_dist=ad_dm, nperm=1000, cpus=4, return_only_pval=F, method='spearman')
groups <- group_samples(rownames(sc_dm),color=T)
sc_tree <- annotated_phylo(tst$test_dist, groups)
ad_tree <- annotated_phylo(tst$ref_dist, groups)
sc_plot <- plot(sc_tree)
ad_plot <- plot(ad_tree)
p <- plot_grid(sc_plot, ad_plot, ncol=2)




