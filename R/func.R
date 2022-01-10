## load required R packages
library(polyG)
library(RColorBrewer)
library(cowplot)


rename_samples <- function(obj_list,this.subject) {
    ## rename samples in the plots
    names(obj_list) <- gsub(subject,'',gsub('_aligned','',names(obj_list)))
    samples <- data.table(oldname=names(obj_list))
    samples$pos <- 1:nrow(samples)
    map <- fread(here('original_data/sample_info.txt'),select=c('subject','oldname','barcode'))
    map <- map[subject==this.subject]
    map[,oldname:=gsub(this.subject,'',oldname)]
    if(this.subject=='C186') {
        map$oldname <- paste0('LM1',map$oldname)
        map[oldname=='LM1N3',barcode:='N1']
    }
    samples <- merge(samples, map, by='oldname', all.x=T)
    samples <- samples[order(pos)]
    samples <- samples$barcode
    names(obj_list) <- samples
    obj_list
}

process_copynumber_data <- function(obj_list, fit_file, sex, this.subject, min_segment_bins=5, field='intcopy', map_file=here('original_data/sample_info.txt'), R, ncpus) {
    #browser()
    options(scipen = 99)
    segments <- list()
    bins <- list()
    fits <- fread(fit_file)
    if(!all(names(obj_list)==fits$barcode)) stop('Names of the obj_list should match the *barcode* column of the fit file!')

    samples <- names(obj_list)    
    chr_info <- sex_chr_info(sex)
    included_samples <- intersect(names(obj_list),samples)
    obj_list <- obj_list[included_samples]
    fits <- fits[barcode %in% included_samples,]

    ## loop through data object and collect relevant bin-level information
    for (i in 1:length(included_samples)){
        sample <- included_samples[i]
        object <- obj_list[[i]]
        purity <- fits$purity[fits$barcode==sample]
        ploidy <- fits$ploidy[fits$barcode==sample]
        message(sample,'; purity=',purity,'; ploidy=',ploidy)
        template <- objectsampletotemplate(object, index=1)
        include_chr <- unique(template$chr)[chr_info$plot_chr_included]
        template <- template[template$chr %in% include_chr,]
        segmentdf <- getadjustedsegments(template, ploidy=ploidy, cellularity = purity, sgc=chr_info$sgc) 
        bindata <- template
        bindata <- bindata[!is.na(bindata$copynumbers),]
        bindata$chr <- as.character(bindata$chr)
        bindata$intcopy <- NA
        bindata$meancopy <- NA

        for (j in 1:dim(segmentdf)[1]){
            chr <- as.character(segmentdf$Chromosome[j])
            start <- segmentdf$Start[j]
            end <- segmentdf$End[j]
            bindata$intcopy[which(bindata$chr==chr & bindata$start >= start & bindata$end <= end)] <- segmentdf$Copies[j]
            bindata$meancopy[which(bindata$chr==chr & bindata$start >= start & bindata$end <= end)] <- segmentdf$Segment_Mean[j]
        }

        segments[[i]] <- segmentdf 
        bins[[i]] <- bindata
    }
    
    # first we create a wide table of int-copies for each bin, merged across all samples
    merge_field_across_samples <- function(bins,included_samples,field) {
        b <- bins[[1]][,c('bin','chr','start','end')]
        extract_value <- function(bins_for_subject) data.frame(value=bins_for_subject[[field]])
        bin_list <- lapply(bins[1:length(bins)], extract_value)
        bin_data <- do.call(cbind,bin_list)
        names(bin_data) <- included_samples
        out <- cbind(b, bin_data)
        out
    }
    b <- (merge_field_across_samples(bins, included_samples, field=field))
    b <- as.data.table(b)

    ## next in each sample we check where any breakpoints occurred, then we format this data into a matrix with only the copy number changes in each sample
    check_breakpoints_per_sample <- function(sample, b) {
        out <- b[,c('bin','chr','start','end')]
        out$sample <- sample
        out$value <- b[[sample]]
        out$change <- c(0,diff(b[[sample]]))
        out
    }
    l0 <- rbindlist(lapply(included_samples, check_breakpoints_per_sample, b))
    l <- data.table::dcast(data=l0, bin + chr + start + end ~ sample, value.var='change')
    l[,id:=paste0(chr,':',start,'-',end)]

    ## now we get the unique breakpoints, based on changes in any sample's copy number
    m <- as.matrix(l[,(included_samples),with=F])
    l$breakpoint <- rowSums(m != 0) > 0
    l$segment <- cumsum(l$breakpoint) + 1

    ## get the regions for each breakpoint
    summarize_segments <- function(l) {
        seg_start <- min(l$start)
        seg_end <- max(l$end)
        n_bins <- nrow(l)
        list(seg_start=seg_start,seg_end=seg_end,n_bins=n_bins)
    }
    segs <- l[,summarize_segments(.SD),by=c('segment','chr')]
    segs <- segs[n_bins >= min_segment_bins,] ## exclude segments less than X bins long

    ## reshape the int-copy number data so that we can get each sample's int-copy number within each segment
    setkey(segs,'chr','seg_start','seg_end')
    setkey(l0,'chr','start','end')
    segs <- foverlaps(l0, segs, type='within')
    segs[,c('start','end','bin','change'):=NULL]
    segs$segid <- paste0(segs$chr,':',segs$seg_start,'-',segs$seg_end)
    segs <- segs[!is.na(segment),]

    ## summarize the unique integer copy number per segment per sample
    summarize_segments_again <- function(segs) {
        value <- unique(segs$value)
        list(value=value, n=length(value))
    }
    dat <- segs[,summarize_segments_again(.SD),by=c('segment','chr','seg_start','seg_end','n_bins','segid','sample')]
    dat <- data.table::dcast(segment + segid + chr + seg_start + seg_end + n_bins ~ sample, value.var='value', data=dat)
    dat <- dat[order(segment),]
    m <- as.matrix(dat[,(7:ncol(dat)),with=F])
    rownames(m) <- dat$segid
    m <- t(m)

    ## get adjusted bins for each sample
    run_for_sample <- function(this.sample, obj_list, fits, sex) {
        message(this.sample)
        chrs <- sex_chr_info(sex)
        ind <- which(fits$barcode==this.sample)
        object <- obj_list[[which(names(obj_list)==this.sample)]]
        cellularity=fits$purity[ind]
        ploidy=fits$ploidy[ind]
        bins <- as.data.table(objectsampletotemplate(object, index=1))
        bins <- bins[!is.na(copynumbers),]
        bins$gc <- rep(2, nrow(bins))
        bins[chr %in% chrs$sgc,gc:=1]
        bins[chr=='MT',gc:=NA]
        bins <- bins[!chr %in% chrs$model_chr_excluded]
        standard <- median(bins$segments,na.rm=T)
        bins$adjustedcopynumbers <- bins$copynumbers * (ploidy + 2/cellularity - 2)/standard - bins$gc/cellularity + bins$gc
        setkey(bins,'chr','start','end')
        template <- objectsampletotemplate(object, index = 1)
        segmentdf <- getadjustedsegments(template, cellularity = cellularity, ploidy=ploidy, sgc=chrs$sgc)
        seg <- as.data.table(segmentdf)
        seg <- seg[,c('Chromosome','Start','End','Num_Bins','Copies'),with=F]
        setkey(seg,'Chromosome','Start','End')
        bins <- foverlaps(bins, seg, type='within')
        bins[,segment:=paste0(chr,':',Start,'-',End)]
        bins[is.na(Num_Bins),segment:=NA]
        bins$sample <- this.sample
        bins
    }

    message('Processing bins ...')
    bin_list <- lapply(samples, run_for_sample, obj_list, fits, sex)
    bins <- rbindlist(bin_list)
   
    ## test for subclonal CNV segments in each sample 
    get_subclonal_segments <- function(seg, bins, threshold, padding) {
        ## load the segments and get the number of samples with detectable CNVs at each segment
        seg <- as.data.table(reshape2::melt(seg))
        parse_segment <- function(s) {
            str <- strsplit(s,'[:]')[[1]]
            chr <- str[1]
            str <- strsplit(str[2],'[-]')[[1]]
            start <- str[1]
            end <- str[2]
            list(chr=chr,start=start,end=end)
        }
        l <- rbindlist(lapply(as.character(seg$Var2), parse_segment))
        alleles <- cbind(seg, l)
        alleles <- alleles[chr %in% 1:22,]
        alleles$copies_expected <- 2
        alleles$copies_observed <- round(alleles$value)
        alleles[copies_observed < copies_expected, state:='del']
        alleles[copies_observed > copies_expected, state:='dup']
        alleles[copies_observed == copies_expected, state:='bal']
        alleles$state <- factor(alleles$state, levels=c('del','bal','dup'))
        alleles <- dcast(Var2 ~ state, fun.aggregate=length, data=alleles)
        alleles$ambiguous <- F
        alleles[dup > 0 & del > 0, ambiguous:=T]
        alleles <- alleles[ambiguous==F,]
        alleles <- alleles[del > 0 | dup > 0,]
        alleles[,ambiguous:=NULL]
        names(alleles)[1] <- 'segment'

        ## annotate each bin with the segment data
        parse_segment <- function(s) {
            str <- strsplit(s,'[:]')[[1]]
            chr <- str[1]
            str <- strsplit(str[2],'[-]')[[1]]
            start <- str[1]
            end <- str[2]
            list(chr=chr,start=start,end=end)
        }
        l <- lapply(as.character(seg$Var2), parse_segment)
        l <- rbindlist(l)
        l$end <- gsub('M','',l$end)
        l$start <- as.integer(l$start)
        l$end <- as.integer(l$end)
        seg2 <- as.data.table(cbind(seg, l))
        setnames(seg2,'Var2','segment')
        setnames(seg2,'Var1','sample')
        bins <- bins[,c('sample','chr','start','end','adjustedcopynumbers'),with=F]
        bins$id <- paste(bins$sample,bins$chr)
        seg2$id <- paste(seg2$sample,seg2$chr)
        setkey(bins,'id','start','end')
        setkey(seg2,'id','start','end')
        bins <- foverlaps(bins, seg2, type='within')
        setnames(bins,c('i.start','i.end'),c('bin.start','bin.end'))
        bins <- bins[!is.na(chr),]
        bins$adjustedintegercopies <- round(bins$value)
        bins <- bins[chr %in% 1:22,]
        bins <- merge(bins, alleles, by='segment', all.x=T)
        bins[is.na(del),del:=F]
        bins[is.na(bal),bal:=F]
        bins[is.na(dup),dup:=F]
        bins$chr <- factor(bins$chr, levels=c(1:22,'X','Y'))
        bins <- bins[order(sample,chr,bin.start,bin.end),]

        ## here we use a t.test to check if each segment has greater/lower mean than the expected integer number of copies, while also constraining the test to segments where the mean would be an (overall) deletion or duplication appropriately
        test_segment <- function(bins) {
            current_mean <- unique(bins$value)
            expected_copies <- unique(bins$adjustedintegercopies)

            ## check that this segment has the expected direction, otherwise we skip it.
            if(unique(bins$del) > 0 & unique(bins$dup) == 0) {
                expected_direction <- 'deletion'
            } else if(unique(bins$dup) > 0 & unique(bins$del) == 0) {
                expected_direction <- 'duplication'
            }    

            if(expected_direction=='deletion' & current_mean > 2) {
                ## this segment was observed as a deletion in some sample, but here the segmean is > 2, 
                ## so this can't be a deletion.
                p <- as.numeric(NA)

            } else if(expected_direction=='duplication' & current_mean < 2) {
                ## this segment was observed as a duplication in some sample, but here the segmean is < 2, 
                ## so this can't be a duplication.
                p <- as.numeric(NA)

            } else if(current_mean < expected_copies - padding) {
                ## the segment has LOWER mean than its nearest integer
                ## therefore we will check if it its subclonally deleted
                p=t.test(bins$adjustedcopynumbers, mu=expected_copies, alternative='less')$p.value

            } else if(current_mean > expected_copies + padding) {
                ## the segment has LOWER mean than its nearest integer
                ## therefore we will check if it its subclonally deleted
                p=t.test(bins$adjustedcopynumbers, mu=expected_copies, alternative='greater')$p.value

            } else {
                p <- as.numeric(NA)

            } 
            list(p=p,mean=current_mean,expected_copies=expected_copies,expected_direction=expected_direction,padding=padding) 
        }
        test <- bins[segment %in% alleles$segment,]
        test <- test[,test_segment(.SD),by=c('sample','segment')]

        ## correct each segment for multiple hypothesis across samples
        test_across_samples <- function(bins) {
            p.adj <- p.adjust(bins$p,method='BH')
            bins$p.adj=p.adj
            bins
        }
        test <- test[,test_across_samples(.SD),by=c('segment')]

        ## format output
        out <- bins[,c('sample','segment','chr','start','end','value','adjustedintegercopies'),with=F]
        out <- merge(out, test[,c('sample','segment','p.adj','expected_copies','expected_direction','padding'),with=F], by=c('sample','segment'), all.x=T)
        out$subclonal <- out$p.adj < threshold
        out[is.na(subclonal),subclonal:=F]
        out$id <- paste(out$sample, out$segment)
        out <- out[!duplicated(id),]
        out[,id:=NULL]
        out <- out[order(segment,sample),]
        setnames(out,c('p.adj','expected_copies','expected_direction','subclonal'),c('sc.p.adj','sc.copies.exp','sc.direction.exp','sc.call'))
        out$sc.threshold <- threshold
        out[sc.call==F,copies:=adjustedintegercopies]
        out[sc.call==T,copies:=round(value,1)]
        out[copies < 0, copies:=0]
        out
    }
    segs <- get_subclonal_segments(m, bins, threshold=0.01, padding=0.1)

    ## load the map file to recode the sample names
    #map <- fread(map_file,select=c('subject','oldname','barcode'))
    #map <- map[subject==this.subject,] 
    #map$oldname <- gsub(subject,'',map$oldname)

    ## update sample names in the 'mat' object
    #m <- cbind(sample=rownames(m), as.data.table(m))
    #tmp <- data.table(oldname=m$sample)
    #tmp$pos <- 1:nrow(tmp)
    #tmp <- merge(tmp, map, by='oldname', all.x=T)
    ##tmp[is.na(barcode),barcode:=paste0(oldname,'?')]
    #tmp <- tmp[order(pos),]
    #m$sample <- tmp$barcode

    ## update sample names in the 'segs' object
    #segs$sample <- as.character(segs$sample)
    #tmp <- data.table(oldname=segs$sample)
    #tmp$pos <- 1:nrow(tmp)
    #tmp <- merge(tmp, map, by='oldname', all.x=T)
    ##tmp[is.na(barcode),barcode:=paste0(oldname,'?')]
    #tmp <- tmp[order(pos),]
    #segs$sample <- tmp$barcode

    ## update sample names in the 'bins' object
    #bins$sample <- as.character(bins$sample)
    #tmp <- data.table(oldname=bins$sample)
    #tmp$pos <- 1:nrow(tmp)
    #tmp <- merge(tmp, map, by='oldname', all.x=T)
    ##tmp[is.na(barcode),barcode:=paste0(oldname,'?')]
    #tmp <- tmp[order(pos),]
    #bins$sample <- tmp$barcode

    ## create a distance matrix from the binned segments
    seg_matrix <- copy(m)
    distance_matrix <- dist(seg_matrix,method='euclidean')
    distance_matrix <- as.matrix(distance_matrix)
   
    ## save output from CNV data
    m_out <- cbind(barcode=rownames(m),as.data.table(m))
    write_tsv(m_out,here(paste0('output/',this.subject,'/',this.subject,'_cnv_matrix.txt')))
    write_tsv(segs,here(paste0('output/',this.subject,'/',this.subject,'_cnv_segments.txt')))
    write_tsv(bins,here(paste0('output/',this.subject,'/',this.subject,'_cnv_bins.txt')))
    write_distance_matrix(dm_df=distance_matrix,filepath=here(paste0('output/',this.subject,'/',this.subject,'_cnv_distance_matrix.txt')))

    ## save plots comparing CNV 5+ Mb segment phylogeny to poly-G angular distance phylogeny
    set.seed(42)
    p2 <- compare_matrices(distance_matrix,this.subject, R, ncpus)
    ggsave(here(paste0('output/',this.subject,'/',this.subject,'_segment_euclidean_matrix_comparison.pdf')),width=7,height=6)
    info2.1 <- compare_trees(distance_matrix,this.subject, tree_method='nj')
    tree2.1 <- info2.1$plot
    ggsave(here(paste0('output/',this.subject,'/',this.subject,'_segment_euclidean_nj_tree_comparison.pdf')),width=10,height=8)

    list(mat=m, bins=bins, segs=segs, distance_matrix=distance_matrix)
}



cnv_heatmap <- function(mat, seg, distance_matrix, this.subject) {

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # heatmap summarizing clonal/subclonal CNV segments
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ## matrix of CNV segments with seg-means
    #mat <- fread(mat_file)
    #mat <- d2m(mat)

    mat <- as.data.table(reshape2::melt(mat))
    names(mat) <- c('sample','segment','seg.mean')

    ## annotated CNV calls
    #seg <- fread(seg_file)
    seg <- seg[,c('sample','segment','sc.call'),with=F]
    mat <- merge(mat, seg, by=c('sample','segment'), all.x=T)
    mat[is.na(sc.call),sc.call:=F] ## only affects chrX
    mat[seg.mean < 0,sc.call:=F]
    mat[seg.mean < 0,seg.mean:=0]
    mat$copies <- ''
    mat[sc.call==F & round(seg.mean) %in% 0:5,copies:=as.character(round(seg.mean))]
    mat[copies=='5' | seg.mean >= 5,copies:='5+']
    mat[copies=='' & sc.call==T,copies:=paste0('(',floor(seg.mean),'-',ceiling(seg.mean),')')]
    mat$copies <- factor(mat$copies, levels=c('0','(0-1)','1','(1-2)','2','(2-3)','3','(3-4)','4','(4-5)','5+'))

    ## merge data
    mat$sample <- factor(mat$sample, levels=unique(mat$sample))
    mat[copies %in% as.character(0:4),subclonal:=F]
    mat[!copies %in% as.character(0:4),subclonal:=T]
    mat[is.na(copies),subclonal:=F]
    mat[copies=='5+',subclonal:=F]

    ## parse genomic position
    parse <- function(pos) {
        s <- strsplit(pos,'[:]')[[1]]
        chr <- s[1]
        s <- strsplit(s[2],'[-]')[[1]]
        start <- as.integer(s[1])
        end <- as.integer(s[2])
        list(chr=chr,start=start,end=end) 
    }
    l <- lapply(as.character(mat$segment),parse)
    toadd2 <- rbindlist(l)
    mat2 <- cbind(mat, toadd2)
    mat2$start <- round(mat2$start / 1e6)
    mat2$end <- round(mat2$end / 1e6)
    mat2[,segment:=paste0(chr,':',start,'-',end)]
    mat2$chr <- factor(mat2$chr, levels=c(1:22,'X','Y'))
    mat2 <- mat2[order(sample,chr,start,end),]
    valid_chr <- as.character(unique(mat2$chr))
    mat2$chr <- factor(mat2$chr, levels=valid_chr)
    segments <- mat2[!duplicated(segment),c('segment','chr','start','end'),with=F]
    mat2$sample <- as.character(mat2$sample)

    ## use the copy-number euclidean distance tree for ordering the samples
    tree <- nj(distance_matrix)
    tree <- phytools::reroot(tree, node.number=grep('^N[0-9]',tree$tip.label))
    p_tree <- ggtree(tree,layout='rect',size=0.5)
    dat <- as.data.frame(p_tree$data)
    dat <- dat[dat$isTip==T,]
    dat <- dat[order(dat$y,decreasing=F),]
    mat2$sample <- factor(mat2$sample, levels=dat$label)

    ## define color scheme
    cols <- c(rev(brewer.pal(9,'Blues')[c(2,4,6,8)]),'white','#FEE8DE',brewer.pal(9,'Reds')[c(3,5,7,9)],'#460000')
    #cols <- c(rev(brewer.pal(9,'Blues')[c(2,2,4,4)]),'white','white',brewer.pal(9,'Reds')[c(3,3,6,6,9)])
    names(cols) <- levels(mat$copies)

    ## expand each segment to to non-NA length
    expand_segment <- function(mat3) {
        start <- min(mat3$start)
        end <- max(mat3$end)
        delta <- start
        start <- start - delta
        end <- end - delta
        list(start=start,end=end)
    }
    tmp <- mat2[,expand_segment(.SD),by=c('sample','segment','chr','seg.mean','sc.call','copies','subclonal')]

    ## concatenate the segments over any gaps
    concat_segments_per_sample <- function(tmp) {
        for(i in 2:nrow(tmp)) {
            tmp$start[i] = tmp$start[i] + tmp$end[i-1]
            tmp$end[i] = tmp$start[i] + tmp$end[i]         
        }
        tmp 
    }
    tmp2 <- tmp[,concat_segments_per_sample(.SD),by=sample]

    tmp2$chr <- factor(tmp2$chr, levels=valid_chr)
    tmp2 <- tmp2[order(sample,chr,start,end),]
    tmp2$samplenum <- as.integer(tmp2$sample)
    tmp2$midpoint <- (tmp2$start + tmp2$end) / 2

    ## get new chr positions
    collapse_chr <- function(tmp2) {
        chr_start <- min(tmp2$start)
        chr_end <- max(tmp2$end)
        midpoint <- (chr_start + chr_end)/2
        list(chr_start=chr_start,chr_end=chr_end,midpoint=midpoint)
    }
    chr2 <- tmp2[,collapse_chr(.SD),by=c('chr')]
    chr3 <- copy(chr2)
    chr3 <- chr3[chr %in% c(1:18,'X','Y'),c('chr','midpoint'),with=F]
    chr3 <- rbind(chr3, data.table(midpoint=mean(chr2$midpoint[chr2$chr %in% c(19:22)]),chr='19-22'))
    chr3 <- chr3[order(midpoint),]
    sample_levels <- levels(tmp2$sample)

    p_heatmap <- ggplot(tmp2) + 
        scale_x_continuous(breaks=c(0,chr2$chr_end), label=c('',as.character(chr2$chr)), expand=c(0,0)) + 
        scale_y_continuous(breaks=1:length(sample_levels),labels=sample_levels,expand=c(0,0),position='right') + 
        geom_rect(aes(xmin=start,xmax=end,ymin=samplenum-0.5,ymax=samplenum+0.5,fill=copies),size=0.125) + 
        geom_point(data=tmp2[subclonal==T],aes(x=midpoint,y=samplenum),pch=16,size=0.3)  +
        geom_text(data=chr3,aes(x=midpoint,label=chr),y=0.25,size=3.5) +
        scale_fill_manual(values=cols,name='Copies',na.value='black') + 
        labs(title=paste(this.subject,'SCNA matrix and euclidean distance tree'),x='\nGenomic position',y='Sample') +
        theme_ang(base_size=12) +
        theme(axis.text.x=element_blank(),axis.line.x=element_blank(),axis.line.y=element_blank(),
              axis.ticks.x=element_blank(),axis.ticks.y=element_blank()) + 
        geom_hline(yintercept=0.5+(0:length(sample_levels)),size=0.25) + 
        geom_vline(xintercept=c(0,chr2$chr_end),size=0.25) +
        coord_cartesian(clip = "off")
    p <- plot_grid(p_tree, p_heatmap, align='h', rel_widths=c(1,6))
    p
}












## side-by-side comparison of distance matrix heatmaps from CNVs and polyG
compare_matrices <- function(cnv_distance_matrix, this.subject, R, ncpus) {
    set.seed(42)

    g <- read_distance_matrix(here(paste0('original_data/polyG/',this.subject,'_ad_matrix.txt')))
    common_samples <- intersect(rownames(cnv_distance_matrix),rownames(g))
    common_samples <- common_samples[!grepl('^N',common_samples)]
    sample_levels <- sort(common_samples, decreasing=T)
    d_subset <- cnv_distance_matrix[common_samples,common_samples]
    g_subset <- g[common_samples,common_samples]
    tst <- dist_similarity(test_dist=d_subset, ref_dist=g_subset, nperm=R, cpus=ncpus, return_only_pval=F, method='spearman')

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
}


compare_trees <- function(cnv_distance_matrix, this.subject, tree_method) {
    ## align CNV and polyG distance matrices
    g <- read_distance_matrix(here(paste0('original_data/polyG/',this.subject,'_ad_matrix.txt')))
    common_samples <- intersect(rownames(cnv_distance_matrix),rownames(g))
    sample_levels <- sort(common_samples, decreasing=T)
    d_subset <- cnv_distance_matrix[common_samples,common_samples]
    g_subset <- g[common_samples,common_samples]

    ## test matrix similarity
    tst <- dist_similarity(test_dist=d_subset, ref_dist=g_subset, nperm=1000, cpus=4, return_only_pval=F)

    ## define groups for plot
    groups <- group_samples(tst$test_dist,color=T,lun=T,liv=F,per=F)

    ## CNV tree
    tree1 <- annotated_phylo(tst$test_dist, groups, method=tree_method)
    p1 <- plot(tree1, legend.position='none', angle=F) + labs(subtitle=paste('     ',this.subject,'SCNA euclidean distance'))

    ## poly-G tree
    tree2 <- annotated_phylo(tst$ref_dist, groups, method=tree_method)
    p2 <- plot(tree2, legend.position='none', angle=F) + labs(subtitle=paste('     ',this.subject,'Polyguanine angular distance'))

    p <- plot_grid(p1,p2,ncol=2)

    list(plot=p, d_subset=d_subset, g_subset=g_subset)
}


sample_info <- function(query.subject = NA) {
    info <- fread(file=here(paste0("original_data/sample_info.txt")))
    #info[is.na(newname), `:=`(newname, "n/a")]
    if (!is.na(query.subject)) {
        info[subject == query.subject]
    }
    else {
        info
    }
}


sex_chr_info <- function(sex) {
    get_sgc <- function(sex) {
        if (sex == "male") {
            c("X", "Y")
        }
        else if (sex == "female") {
            c()
        }
    }
    sgc <- get_sgc(sex)
    if ("Y" %in% sgc) {
        chrsubset = 1:24
        excluded_chr = c("X", "Y", "MT")
    }
    else {
        chrsubset = 1:23
        excluded_chr = c("Y", "MT")
    }
    list(sgc=sgc, plot_chr_included=chrsubset, model_chr_excluded=excluded_chr)
}


refit <- function(object, sex, purity=NA, ploidy, samplename=NA, sampleindex=1, save=F, bottom=-1, cap=10, output_dir='.') {
    get_sgc <- function(sex) {
        if(sex=='male') {
            c('X','Y')
        } else if(sex=='female') {
            c()
        }
    }
    sgc <- get_sgc(sex)
    if('Y' %in% sgc) {
        ## 1X, 1Y
        chrsubset=1:24
        excluded_chr=c('X','Y','MT')
    } else {
        ## 2X, 0Y
        chrsubset=1:23
        excluded_chr=c('Y','MT')
    }

    if(is.na(samplename)) {
        pd <- Biobase::pData(object)
        samplename <- gsub("_.*", "",pd$name)
    }

    if(is.na(purity)) {
        model <- singlemodel(object, QDNAseqobjectsample=sampleindex, ploidy=ploidy, exclude=excluded_chr); model$minima
        return(model$minima)
    }

    title = paste0(samplename,': ploidy=',ploidy,', purity=',purity)
    model <- singlemodel(object, QDNAseqobjectsample=sampleindex, ploidy=ploidy, exclude=excluded_chr)
    p <- singleplot(object, QDNAseqobjectsample=sampleindex, cellularity=purity, ploidy=ploidy, bottom=bottom, cap=cap, sgc=sgc, onlyautosomes=F, chrsubset=chrsubset, title=title)

    if(save==F) {
        ## return the results
        p
    } else {
        refits_dir=file.path(output_dir,'fits')
        if(!dir.exists(refits_dir)) dir.create(refits_dir,recursive=T)

        ## save the plot for this fit
        plot_file <- paste0(refits_dir,'/',samplename,'_N=',ploidy,'_cellularity=',purity,'.pdf')
        message('Saving plot: ',plot_file)
        ggsave(p,dev=cairo_pdf,filename=plot_file,width=8,height=5)

        ## save data (append it to existing file along with the date of the entry)
        fit_file <- paste0(refits_dir,'/purity_ploidy.txt')
        if(!file.exists(fit_file)) cat('barcode\tpurity\tploidy\tdate\n',file=fit_file,append=F) 
        when <- as.character(format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
        message('Saving data: ',fit_file)
        cat(paste(samplename,purity,ploidy,when,sep='\t'),'\n',file=fit_file,append=T)
    }
}
