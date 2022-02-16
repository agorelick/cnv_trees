## load required R packages
require(RColorBrewer)
require(cowplot)
require(TreeTools)
require(here)
require(data.table)
require(ACE)
require(ggplot2)
require(RColorBrewer)
require(ggtree)
require(phytools)
require(TreeDist)

## first extract the adjusted copy number bins and segments for each sample
get_bins_and_segments <- function(obj_list, fit_file, sex) {
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

    ## loop through each sample, generating bins and segments with purity/ploidy-adjusted values
    for (i in 1:length(included_samples)){
        sample <- included_samples[i]
        object <- obj_list[[i]]
        purity <- fits$purity[fits$barcode==sample]
        ploidy <- fits$ploidy[fits$barcode==sample]
        message(sample,'; purity=',purity,'; ploidy=',ploidy)
    
        ## get adjusted segments
        template <- objectsampletotemplate(object, index=1)
        include_chr <- unique(template$chr)[chr_info$plot_chr_included]
        template <- template[template$chr %in% include_chr,]
        segmentdf <- getadjustedsegments(template, ploidy=ploidy, cellularity = purity, sgc=chr_info$sgc) 

        ## adjust the bin's copy number with the same formula as for getting adjustedsegments
        bindata <- template
        bindata$chr <- as.character(bindata$chr)
        standard <- median(template$segments,na.rm=T)
        gc <- rep(2, nrow(template))
        gc[template$chr %in% chr_info$sgc] <- 1
        bindata$adjustedcopynumbers <- bindata$copynumbers * (ploidy + 2/purity - 2)/standard - gc/purity + gc
        bindata <- bindata[!is.na(bindata$copynumbers),]

        ## add in the segment's adjusted intcopy and meancopy to the bins
        bindata$intcopy <- NA
        bindata$meancopy <- NA
        for (j in 1:dim(segmentdf)[1]){
            chr <- as.character(segmentdf$Chromosome[j])
            start <- segmentdf$Start[j]
            end <- segmentdf$End[j]
            bindata$intcopy[which(bindata$chr==chr & bindata$start >= start & bindata$end <= end)] <- segmentdf$Copies[j]
            bindata$meancopy[which(bindata$chr==chr & bindata$start >= start & bindata$end <= end)] <- segmentdf$Segment_Mean[j]
        }

        ## add the bins and segments to our running lists for each sample
        segments[[i]] <- segmentdf 
        bins[[i]] <- bindata
    }
    names(segments) <- included_samples
    names(bins) <- included_samples
    list(segments=segments, bins=bins)
}


test_subclonal_segments_and_annotate_bins <- function(this.sample, info) {
    ## for each sample, test whether each segment has a significantly subclonal copy number alteration based on: 
    ## 1. segmean is more than 0.1 away from the nearest integer copy number
    ## 2. adjust p-value from a one-sided t-test in the expected direction is significant
    ## 3. the copy number is less than 5

    message('Testing for subclonal segments in sample ', this.sample)
    segs <- as.data.table(info$segments[[this.sample]])
    segs[,segment:=paste0(Chromosome,':',Start,'-',End)]
    bins <- as.data.table(info$bins[[this.sample]])
    bins[,binID:=paste0(chr,':',start,'-',end)]

    setkey(bins,'chr','start','end')
    setkey(segs,'Chromosome','Start','End')
    tmp <- foverlaps(bins, segs, type='any')
    bin_to_segment <- tmp[!duplicated(binID),c('binID','segment'),with=F]
    bins <- merge(bins, bin_to_segment, by='binID', all.x=T)
    bins[,binID:=NULL]
    test_subclonal_segment <- function(tmp) {
        options(scipen=6)
        this_cn <- as.numeric(unique(tmp$Segment_Mean))
        diff_from_clonal <- as.numeric(abs(this_cn - round(this_cn)))
        int_cn <- as.integer(unique(tmp$Copies))
        p <- tryCatch({
            t.test(tmp$adjustedcopynumbers, mu=int_cn, alternative='two.sided')$p.value
        },error=function(e) {
            as.numeric(NA)
        })
        list(this_cn=this_cn,diff_from_clonal=diff_from_clonal,int_cn=int_cn,p=p)
    }
    res <- tmp[,test_subclonal_segment(.SD),by=c('chr','Start','End','segment','Num_Bins','P_log10')]
    res$sc.p.adj <- p.adjust(res$p,method='BH')
    res[,p:=NULL]
    res[diff_from_clonal < 0.1 | int_cn >= 5 | chr %in% c('X','Y','MT'),c('sc.p.adj'):=list(NA)]
    res[sc.p.adj >= 0.01 | is.na(sc.p.adj),sc.call:=F]
    res[sc.p.adj < 0.01,sc.call:=T]
    res[sc.call==F,copies:=as.numeric(int_cn)]
    res[sc.call==T,copies:=as.numeric(round(this_cn,1))]
    res[copies < 0,copies:=as.numeric(0)]

    ## merge the subclonality results to the bins
    bins <- merge(bins, res[,c('segment','sc.p.adj','sc.call','copies'),with=F], by='segment', all.x=T)
    bins$chr <- factor(bins$chr, levels=c(1:22,'X','Y'))
    bins <- bins[order(chr,start,end),]                   
    bins$sample <- this.sample
    bins
}


## now we use the floating-point copy number to get a new list of segments
get_new_segments_from_bins <- function(bins_list, min_segment_bins=5) {

    # first we create a wide table of copies for each bin, merged across all samples
    # we will use the 'copies' field, containing either integer copies, or the floating-point copy number in the case of subclonal SCNAs
    merge_field_across_samples <- function(bins,included_samples,field) {
        b <- bins[[1]][,c('bin','chr','start','end')]
        extract_value <- function(bins_for_subject) data.frame(value=bins_for_subject[[field]])
        bin_list <- lapply(bins[1:length(bins)], extract_value)
        bin_data <- do.call(cbind,bin_list)
        names(bin_data) <- included_samples
        out <- cbind(b, bin_data)
        out
    }
    b <- (merge_field_across_samples(bins_list, samples, field='copies'))
    b <- as.data.table(b)

    ## next in each sample we check where any breakpoints occurred, then we format this data into a matrix with only the copy number changes in each sample
    check_breakpoints_per_sample <- function(sample, b) {
        out <- b[,c('bin','chr','start','end')]
        out$sample <- sample
        out$value <- b[[sample]]
        out$change <- c(0,diff(b[[sample]]))
        out
    }
    l0 <- rbindlist(lapply(samples, check_breakpoints_per_sample, b))
    l <- data.table::dcast(data=l0, bin + chr + start + end ~ sample, value.var='change')
    l[,id:=paste0(chr,':',start,'-',end)]

    ## now we get the unique breakpoints, based on changes in any sample's copy number
    m <- as.matrix(l[,(samples),with=F])
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

    ## remove duplicate segments from merging
    segs[,uniqueID:=paste(sample,segid)]
    segs <- segs[!duplicated(uniqueID),]
    segs[,uniqueID:=NULL]
}

get_copynumber_matrix_from_segs <- function(segs,field) { 
    ## summarize the unique floating-point copy number per segment per sample
    summarize_segments_again <- function(segs,field) {
        myvalue <- unique(segs[[field]])
        list(myvalue=myvalue, n=length(myvalue))
    }
    dat <- segs[,summarize_segments_again(.SD,field),by=c('segment','chr','seg_start','seg_end','n_bins','segid','sample')]
    dat <- data.table::dcast(segment + segid + chr + seg_start + seg_end + n_bins ~ sample, value.var='myvalue', data=dat)
    dat <- dat[order(segment),]
    m <- as.matrix(dat[,(7:ncol(dat)),with=F])
    rownames(m) <- dat$segid
    m <- t(m)
    m
}


get_new_segments_and_annotate_with_subclonality <- function(bins_list) {
    message('Generate new segment based on the integer (and subclonal) copy number bins ...')
    segs <- get_new_segments_from_bins(bins_list)

    ## annotate the new segments with subclonality info extracted from the bins
    bins <- rbindlist(bins_list)
    tmp <- bins[,c('sample','bin','chr','start','end','sc.p.adj','sc.call','copies'),with=F]
    segs[,sample_chr:=paste(sample,chr)]
    tmp[,sample_chr:=paste(sample,chr)]
    setkey(segs,'sample_chr','seg_start','seg_end')
    setkey(tmp,'sample_chr','start','end')
    tmp <- foverlaps(tmp, segs, type='any')
    tmp[,uniqueID:=paste(sample,chr,start,end)]

    collapse_segment_subclonality <- function(tmp) {
        tmp[!duplicated(copies),] 
    }
    tmp <- tmp[,collapse_segment_subclonality(.SD),by=c('sample','segment','chr','seg_start','seg_end')]
    tmp <- tmp[!is.na(segid) & !is.na(sample),]
    tmp[,c('i.chr','i.sample','start','end','uniqueID'):=NULL]
    tmp <- tmp[order(chr,seg_start,seg_end,sample),]
    tmp$segid <- factor(tmp$segid, levels=unique(tmp$segid))
    tmp$segment <- as.integer(tmp$segid)
    tmp[,sample_chr:=NULL]
    tmp
}


update_bins_with_new_segments <- function(bins_list, segs) {
    ## update the 'segment' field in bins to reflect the new segments

    bins <- rbindlist(bins_list)
    bins[,segment:=NULL]
    front_fields <- names(bins)
    bins[,sample_chr:=paste(sample,chr)]
    setkey(bins,'sample_chr','start','end')

    toadd <- segs[,c('sample','chr','seg_start','seg_end','segid'),with=F]
    toadd[,sample_chr:=paste(sample,chr)]
    toadd[,c('sample','chr'):=NULL]
    setkey(toadd,'sample_chr','seg_start','seg_end')
    bins <- foverlaps(bins, toadd, type='any')
    bins[,sample_chr:=NULL]
    back_fields <- names(bins)[!names(bins) %in% front_fields]
    fields <- c(front_fields, back_fields)
    fields <- fields[fields %in% names(bins)]
    bins <- bins[,(fields),with=F]
    bins
}


process_SCNA_data <- function(samples, info, this.subject) {
    ## get bins with purity/ploidy-corrected copy number and annotations for subclonality
    bins_list <- lapply(samples, test_subclonal_segments_and_annotate_bins, info)
    names(bins_list) <- samples

    ## get new copy number segments based on the bins' integer or subclonal copy numbers
    segs <- get_new_segments_and_annotate_with_subclonality(bins_list)

    ## get the sample/segment copy number matrix
    mat <- get_copynumber_matrix_from_segs(segs,'copies')

    ## update the 'segment' field in bins to reflect the new segments
    bins <- update_bins_with_new_segments(bins_list, segs)

    ## create a distance matrix from the binned segments
    seg_matrix <- copy(mat)
    distance_matrix <- dist(seg_matrix,method='euclidean')
    distance_matrix <- as.matrix(distance_matrix)

    ## save output from CNV data
    mat_out <- cbind(barcode=rownames(mat),as.data.table(mat))
    write_tsv(mat_out,here(paste0('output/',this.subject,'/',this.subject,'_cnv_matrix.txt')))
    write_tsv(segs,here(paste0('output/',this.subject,'/',this.subject,'_cnv_segments.txt')))
    write_tsv(bins,here(paste0('output/',this.subject,'/',this.subject,'_cnv_bins.txt')))
    write_distance_matrix(dm_df=distance_matrix,filepath=here(paste0('output/',this.subject,'/',this.subject,'_cnv_distance_matrix.txt')))

    list(mat=mat, bins=bins, segs=segs, distance_matrix=distance_matrix)
}


rename_samples <- function(obj_list,this.subject) {
    ## rename samples in the plots
    names(obj_list) <- gsub('_aligned','',names(obj_list))
    samples <- data.table(Sample_ID=names(obj_list))
    samples$pos <- 1:nrow(samples)
    map <- fread(here('original_data/sample_info.txt'),select=c('Subject','Sample_ID','barcode'))
    map <- map[Subject==this.subject]
    samples <- merge(samples, map, by='Sample_ID', all.x=T)
    samples <- samples[order(pos)]
    names(obj_list) <- samples$barcode
    obj_list
}


process_copynumber_data <- function(obj_list, fit_file, sex, this.subject, min_segment_bins=5, field='intcopy', map_file=here('original_data/sample_info.txt'), R) {
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

    list(mat=m, bins=bins, segs=segs, distance_matrix=distance_matrix)
}



test_tree_similarity <- function(mat1,mat2,method,nperm=1000,title=NULL,return_plot=T) {
    require(TreeDist)
    require(Quartet)

    ## subset both matrices to the common samples 
    common_samples <- intersect(rownames(mat1),rownames(mat2))
    mat1 <- mat1[common_samples,common_samples]
    mat2 <- mat2[common_samples,common_samples]
    tree1 <- nj(mat1)
    tree2 <- nj(mat2)

    permute_test_and_get_shared_info <- function(i, mat1, tree2) {
        original_samples <- rownames(mat1)
        if(i==0) {
            permuted_samples <- copy(original_samples)
        } else {
            permuted_samples <- sample(original_samples,replace=F)
        }
        perm_mat1 <- mat1[permuted_samples,permuted_samples]    
        rownames(perm_mat1) <- original_samples; colnames(perm_mat1) <- original_samples
        perm_tree1 <- nj(perm_mat1)        

        if(method=='grf') {
            tree_grf_similarity <- function(test_tree, ref_tree) {
                test_info <- SplitwiseInfo(test_tree)
                ref_info <- SplitwiseInfo(ref_tree) 
                shared_info <- SharedPhylogeneticInfo(test_tree, ref_tree)
                data.table(test_info=test_info, ref_info=ref_info, shared_info=shared_info)
            }
            info <- tree_grf_similarity(test_tree=perm_tree1, ref_tree=tree2)
        } else {
            info <- as.data.frame(QuartetStatus(perm_tree1, cf=tree2))
        }
        info$perm <- i
        info
    }
    l <- lapply(0:nperm, permute_test_and_get_shared_info, mat1, tree2)
    res <- rbindlist(l)

    if(method=='quartet') {
        ## for Quartet
        res$similarity <- res$s / res$Q
    } else {
        ## for GRF
        m <- as.matrix(res[,c('test_info','ref_info'),with=F])
        res$total_info <- apply(m,1,min)
        #res$similarity <- res$shared_info / res$ref_info # normalize to the reference
        res$similarity <- res$shared_info / res$total_info
    }
    obs <- res[perm==0,]
    exp <- res[perm > 0,]
    numerator <- sum(exp$similarity >= obs$similarity) + 1
    denominator <- nrow(exp) + 1    
    pval <- numerator / denominator
    pval <- prettyNum(pval,digits=3)
    mean_exp <- prettyNum(mean(exp$similarity),digits=3)
    confint_exp <- prettyNum(quantile(exp$similarity,c(0.025,0.975)),digits=3)

    label <- paste0('Expected: ',mean_exp,', 95%CI: [',confint_exp[1],'-',confint_exp[2],']; Observed: ',prettyNum(obs$similarity,digits=3),'; P-value: ',pval,'; ',nperm,' permutations')

    p <- ggplot(exp,aes(x=similarity)) +
        geom_histogram(bins=80,fill='#a6a6a6',size=1,color='white') +
        theme_ang(base_size=12) +
        geom_segment(x=obs$similarity,xend=obs$similarity,y=0,yend=Inf,color='red')

    p <- p + labs(x='Normalized tree similarity (shared phylogenetic info)',y='N random permutations',title=title,subtitle=label)   
    p <- p + scale_x_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25),expand=c(0,0))
    p <- p + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))

    if(return_plot==T) {
        p
    } else {
        res
    }
}


get_bootstrapped_trees <- function(this.subject, type) {
    ## load the SCNA table and generate both segment and chromosome-level bootstrap values.
    ## we will compare these to the poly-G bootstrap values to show that chromosome-level are most consistent across samples (like poly-G) as opposed to segment-level.
    set.seed(42)

    if(type!='polyg') {
        message(paste0('Bootstrapping SCNA tree (',type,' resampling) ...'))
    } else {
        message('Bootstrapping poly-G tree (marker resampling) ...')    
    }

    get_valid_samples <- function(this.subject) {
        info <- fread(here('original_data/sample_info.txt'))
        info <- info[Subject==this.subject]
        list(polyg_samples=info$Sample_ID, scna_samples=info$barcode)
    }
    samples <- get_valid_samples(this.subject)

    if(type=='chromosome') {
        ## load the original SCNA sample segment matrix
        scna_segments <- fread(here(paste0('output/',this.subject,'/',this.subject,'_cnv_segments.txt')))
        valid_samples <- samples$scna_samples
        scna_segments <- scna_segments[sample %in% valid_samples,]

        get_trees_bootstrapped_by_chromosome <- function(i, mat) { 
            ## keep the segments for X and Y labeled separately but group them under XY for the sampling
            mat[chr %in% 'X',segid:=paste0('XY.',segid)]
            mat[chr %in% 'Y',segid:=paste0('XY.',segid)]
            mat[chr %in% c('X','Y'),chr:='XY']
            mat$chr <- factor(mat$chr,levels=c(1:22,'XY'))
            chr_with_replacement <- c(1:22,'XY')
            if(i > 0) chr_with_replacement <- sample(c(1:22,'XY'),replace=T)
            data_per_chromosome <- function(this.chr, mat) {
                mat <- mat[chr %in% this.chr] 
                mat  
            }
            list_resampled_mat <- lapply(chr_with_replacement, data_per_chromosome, mat)
            for(i in 1:length(list_resampled_mat)) list_resampled_mat[[i]]$newchr <- i 
            resampled_mat <- rbindlist(list_resampled_mat)
            resampled_mat[,segid:=paste0(newchr,'.',segid)]
            resampled_mat[,c('chr','newchr'):=NULL]
            resampled_mat <- d2m(resampled_mat)
            distance_matrix <- dist(t(resampled_mat),method='euclidean')
            tree <- nj(distance_matrix)
        }
        mat <- dcast(chr + segid ~ sample, value.var='copies', data=scna_segments)
        bstrees <- lapply(1:1000, get_trees_bootstrapped_by_chromosome, mat)
        tree <- get_trees_bootstrapped_by_chromosome(0, mat)
        bstrees <- TreeTools::as.multiPhylo(bstrees)

        ## root tree before adding the boostrap values!
        refind <- grep(paste0('^Normal'),tree$tip.label)
        tree <- root(tree,outgroup=refind,resolve.root=TRUE)
        tree <- addConfidences(tree, bstrees) 
        tree$node.label <- round(100*tree$node.label)


    } else if(type=='segment') { 
        ## load the original SCNA sample segment matrix
        scna_segments <- fread(here(paste0('output/',this.subject,'/',this.subject,'_cnv_segments.txt')))
        valid_samples <- samples$scna_samples
        scna_segments <- scna_segments[sample %in% valid_samples,]

        get_trees_bootstrapped_by_segment <- function(i, mat2) { 
            if(i > 0) mat2 <- mat2[,sample(colnames(mat2),replace=T)]
            distance_matrix <- dist(mat2,method='euclidean')
            tree <- nj(distance_matrix)
        }
        mat2 <- dcast(sample ~ segid, value.var='copies', data=scna_segments)
        mat2 <- d2m(mat2)
        bstrees <- lapply(1:1000, get_trees_bootstrapped_by_segment, mat2)
        tree <- get_trees_bootstrapped_by_segment(0, mat2)
        bstrees <- TreeTools::as.multiPhylo(bstrees)

        ## root tree before adding the boostrap values!
        refind <- grep(paste0('^Normal'),tree$tip.label)
        tree <- root(tree,outgroup=refind,resolve.root=TRUE)
        tree <- addConfidences(tree, bstrees) 
        tree$node.label <- round(100*tree$node.label)

    } else if(type=='polyg') {
        orig_matrix <- read_distance_matrix(here(paste0('original_data/polyG/ad_matrix_oldnames/',this.subject,'_angular_dist_matrix_w_root_usedmarkers_repreReplicate_oldnames.txt')))
        tree <- nj(orig_matrix)

        ## root tree before adding the boostrap values!
        refind <- grep(paste0('^',this.subject,'N[0-9]$'),tree$tip.label)
        tree <- root(tree,outgroup=refind,resolve.root=TRUE)

        ## now add bootstrap values
        load(here(paste0('original_data/polyG/bstrees/',this.subject,'_bstrees_usedmarkers_repreReplicate_NJ.RData'))) ## load bstrees
        bstrees <- TreeTools::as.multiPhylo(bstrees)
        tree <- addConfidences(tree, bstrees) 
        tree$node.label <- round(100*tree$node.label)

        ## update sample names after bootstrapping
        update_tree_labels <- function(tree) { 
            tmp <- data.table(Sample_ID=tree$tip.label)
            tmp$pos <- 1:nrow(tmp)
            map <- data.table(Sample_ID=samples$polyg_samples, newname=samples$scna_samples)
            tmp <- merge(tmp, map, by='Sample_ID', all.x=T)
            tmp <- tmp[order(pos),]
            tree$tip.label <- tmp$newname
            tree
        }
        tree <- update_tree_labels(tree)
        for(i in 1:length(bstrees)) bstrees[[i]] <- update_tree_labels(bstrees[[i]])
    }
    list(tree=tree, bstrees=bstrees)
}



cnv_heatmap <- function(mat, seg, distance_matrix, this.subject) {
    #browser()
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # heatmap summarizing clonal/subclonal CNV segments
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ## matrix of CNV segments with seg-means
    #mat <- fread(mat_file)
    #mat <- d2m(mat)

    mat <- as.data.table(reshape2::melt(mat))
    names(mat) <- c('sample','segment','seg.mean')

    ## annotated CNV calls
    seg <- seg[,c('sample','segid','sc.call'),with=F]
    mat <- merge(mat, seg, by.x=c('sample','segment'), by.y=c('sample','segid'), all.x=T)

    mat$copies <- as.character(mat$seg.mean)
    mat[seg.mean >= 5,copies:='5+']
    mat[sc.call==T,copies:=paste0('(',floor(seg.mean),'-',ceiling(seg.mean),')')]
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
    tree <- phytools::reroot(tree, node.number=grep('^Normal[0-9]',tree$tip.label))
    p_tree <- ggtree(tree,layout='rect',size=0.5)
    dat <- as.data.frame(p_tree$data)
    dat <- dat[dat$isTip==T,]
    dat <- dat[order(dat$y,decreasing=F),]
    mat2$sample <- factor(mat2$sample, levels=dat$label)

    ## define color scheme
    cols <- c(rev(brewer.pal(9,'Blues')[c(2,4,6,8)]),'white','#FEE8DE',brewer.pal(9,'Reds')[c(3,5,7,9)],'#52000a')
    #cols <- c("#08519C","#4292C6","#9ECAE1","#DEEBF7",'white','#FEE0D2','#FC9272','#EF3B2C','#A50F15','#840c10','#63090c')
    #cols <- c(rev(brewer.pal(9,'Blues')[c(2,4,6,8)]),'white',brewer.pal(9,'OrRd')[c(2,4,6,8,9)])
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
        geom_point(data=tmp2[subclonal==T],aes(x=midpoint,y=samplenum),pch=16,color='black',size=0.3)  +
        geom_point(data=tmp2[subclonal==T & copies=='(4-5)'],aes(x=midpoint,y=samplenum),pch=16,color='white',size=0.3)  +
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


update_distance_matrix_samples <- function(g) {
    map <- fread(here('original_data/sample_info.txt'),select=c('Subject','Sample_ID','barcode'))
    tmp <- data.table(Sample_ID=rownames(g))
    tmp$pos <- 1:nrow(tmp)
    tmp <- merge(tmp, map[,c('Sample_ID','barcode'),with=F], by='Sample_ID', all.x=T)
    tmp <- tmp[order(pos),]
    rownames(g) <- tmp$barcode; colnames(g) <- tmp$barcode
    g
}



## side-by-side comparison of distance matrix heatmaps from CNVs and polyG
compare_matrices <- function(cnv_distance_matrix, this.subject, R) {
    set.seed(42)

    g <- read_distance_matrix(here(paste0('original_data/polyG/ad_matrix_oldnames/',this.subject,'_angular_dist_matrix_w_root_usedmarkers_repreReplicate_oldnames.txt')))
    g <- update_distance_matrix_samples(g)

    common_samples <- intersect(rownames(cnv_distance_matrix),rownames(g))
    common_samples <- common_samples[!grepl('^Normal',common_samples)]
    sample_levels <- sort(common_samples, decreasing=T)
    d_subset <- cnv_distance_matrix[common_samples,common_samples]
    g_subset <- g[common_samples,common_samples]
    tst <- dist_similarity(test_dist=d_subset, ref_dist=g_subset, nperm=R, return_only_pval=F, method='spearman')

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
    g <- read_distance_matrix(here(paste0('original_data/polyG/ad_matrix_oldnames/',this.subject,'_angular_dist_matrix_w_root_usedmarkers_repreReplicate_oldnames.txt')))
    g <- update_distance_matrix_samples(g)

    common_samples <- intersect(rownames(cnv_distance_matrix),rownames(g))
    sample_levels <- sort(common_samples, decreasing=T)
    d_subset <- cnv_distance_matrix[common_samples,common_samples]
    g_subset <- g[common_samples,common_samples]

    ## test matrix similarity
    tst <- dist_similarity(test_dist=d_subset, ref_dist=g_subset, nperm=1000, return_only_pval=F)

    ## define groups for plot
    groups <- group_samples(tst$test_dist,color=T,lun=F,liv=F,per=T)

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


d2m <- function(dt) {
    ## assumes first column should be rownames of a matrix made from columns 2:ncol
    rows <- dt[[1]]
    if('data.table' %in% class(dt)) {
        dt <- dt[,c(2:ncol(dt)),with=F]
    } else if(class(dt)=='data.frame') {
        ## assume this is a data.frame
        dt <- dt[,c(2:ncol(dt))]
    } else {
        stop('enter a data.table or a data.frame')
    }
    m <- as.matrix(dt)
    rownames(m) <- rows
    m
}



distance_matrix_correlation <- function(test_dist_full, ref_dist_full, method, return_only_r=F) {
    ref_dist <- ref_dist_full[order(rownames(ref_dist_full)),order(colnames(ref_dist_full))]    
    test_dist <- test_dist_full[order(rownames(test_dist_full)),order(colnames(test_dist_full))]    

    ## subset the distances for common samples 
    common_samples <- sort(intersect(rownames(test_dist), rownames(ref_dist)))
    test_dist <- test_dist[common_samples,common_samples]
    ref_dist <- ref_dist[common_samples,common_samples]
    
    ## get a long table of sample:sample distances with columns for test and ref data
    ref_dist_long <- ref_dist; test_dist_long <- test_dist; 
    ref_dist_long[upper.tri(ref_dist_long)] <- NA
    test_dist_long[upper.tri(test_dist_long)] <- NA
    test_dist_long <- as.data.table(reshape2::melt(test_dist_long))
    ref_dist_long <- as.data.table(reshape2::melt(ref_dist_long))
    merged <- merge(test_dist_long,ref_dist_long,by=c('Var1','Var2'))
    merged <- merged[Var1!=Var2 & !is.na(value.x) & !is.na(value.y)] 
    merged[,c('Var1','Var2'):=NULL]
    names(merged) <- c('test_dist','ref_dist')

    ## get the correlation
    r <- cor(merged$test_dist,merged$ref_dist,method=method)    

    if(return_only_r==T) {
        out <- r
    } else {
        out <- list(r=r, merged=merged, test_dist=test_dist, ref_dist=ref_dist, test_dist_full=test_dist_full, ref_dist_full=ref_dist_full)
    }
    out
}


dist_similarity <- function(test_dist, ref_dist, nperm=100, return_only_pval=T, method='spearman') {

    shuffled_correlation <- function(i, test_dist, ref_dist, return_only_r, method) { 
        if(i > 0) {
            ## shuffle among all the samples with data (not only the common ones)
            orig_samples <- rownames(test_dist)
            shuffled_samples <- sample(orig_samples,replace=F)
            test_dist <- test_dist[shuffled_samples,shuffled_samples]
            rownames(test_dist) <- orig_samples
            colnames(test_dist) <- orig_samples
        }
        info <- distance_matrix_correlation(test_dist, ref_dist, return_only_r=F, method=method)
        r <- info$r
        merged <- info$merged
        test_dist <- info$test_dist
        ref_dist <- info$ref_dist
        test_dist_full <- info$test_dist_full
        ref_dist_full <- info$ref_dist_full
        ## return either just the distance or the complete data
        if(return_only_r) { 
            r
        } else {
            list(test_dist=test_dist, ref_dist=ref_dist, test_dist_full=test_dist_full, ref_dist_full=ref_dist_full, merged=merged, method=method, r=r)
        }
    }

    observed_info <- shuffled_correlation(0, test_dist, ref_dist, return_only_r=F, method=method)
    ref_dist <- observed_info$ref_dist
    test_dist <- observed_info$test_dist
    observed <- observed_info$r
    merged <- observed_info$merged

    permuted <- unlist(lapply(1:nperm, shuffled_correlation, test_dist, ref_dist, return_only_r=T, method=method))
    perms_as_correlated <- sum(permuted >= observed)
    p.value <- (perms_as_correlated + 1) / (nperm + 1) ## add pseudo-count to avoid p=0

    if(return_only_pval==T) {
        p.value
    } else {
        tmp <- data.table(permuted=permuted)
        plot <- ggplot(tmp, aes(x=permuted)) +
            geom_histogram(fill='black',color='white',bins=30) +
            theme_ang(base_size=12)  +
            geom_vline(xintercept=observed,linetype='dashed',color='red') +
            labs(x='Correlation of sample-sample distances',y='N permutations',
                 subtitle=paste0('Observed R=',round(observed,3),'; P-value = ',round(p.value,log10(nperm)))
            )
        list(p.value=p.value, observed=observed, permuted=permuted, merged=merged, ref_dist=ref_dist, test_dist=test_dist, plot=plot)
    }
}


read_distance_matrix <- function(file,return.as.matrix=T) {
    ## read txt file with a saved distance matrix (i.e. table with named rows and cols and numeric distances as cells)
    distance_matrix <- fread(file)
    rows <- distance_matrix[[1]]
    distance_matrix <- distance_matrix[,(2:ncol(distance_matrix)),with=F]
    m <- as.matrix(distance_matrix)
    rownames(m) <- rows
    if(return.as.matrix==F) {
        as.dist(m,diag=T)
    } else {
        m
    }
}


write_tsv <- function(d,file,sep='\t',quote=F,row.names=F,...) write.table(d,file=file,sep=sep,quote=quote,row.names=row.names,...)


write_distance_matrix <- function(filepath,dm_df) {
    write.table(dm_df,file=filepath,sep="\t",quote=FALSE,col.names=NA)
}


theme_ang <- function (base_size = 11, base_line_size = base_size/22, base_rect_size = base_size/22) {
    require(ggplot2)
    theme_bw(base_size = base_size, base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace%
    theme(
          line = element_line(colour = "black", size = base_line_size, linetype = 1, lineend = "round"), 
          text = element_text(colour = "black", size = base_size, lineheight = 0.9, hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(), debug = F), 
          axis.text = element_text(colour = "black", size = rel(0.8)), 
          axis.ticks = element_line(colour = "black", size = rel(1)), 
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black", size = rel(1)), 
          legend.key = element_blank(),
          strip.background = element_blank())
}


break_axis <- function(y, maxlower, minupper=NA, lowerticksize, upperticksize, ratio_lower_to_upper) {
    if(is.na(minupper)) {
        breakpos <- maxlower
        lowerticklabels <- seq(0,breakpos,by=lowerticksize); lowerticklabels
        upperticklabels <- seq(breakpos+upperticksize,max(y)+upperticksize,by=upperticksize); upperticklabels
        ticklabels <- c(lowerticklabels, upperticklabels); ticklabels
        lowertickpos <- lowerticklabels
        uppertickspacing <- ratio_lower_to_upper * lowerticksize
        uppertickpos <- breakpos + ((1:length(upperticklabels))*uppertickspacing)
        tickpos <- c(lowertickpos, uppertickpos)
        newy <- as.numeric(y)
        ind <- newy > breakpos
        newy[ind] <- breakpos + uppertickspacing*((newy[ind]-breakpos) / upperticksize)
        list(newy=newy, breaks=tickpos, labels=ticklabels, limits=range(tickpos))
    } else {
        lowerticklabels <- seq(0,maxlower,by=lowerticksize); lowerticklabels
        upperticklabels <- seq(minupper,max(y)+upperticksize,by=upperticksize); upperticklabels
        ticklabels <- c(lowerticklabels, upperticklabels); ticklabels
        lowertickpos <- lowerticklabels
        uppertickspacing <- ratio_lower_to_upper * lowerticksize
        uppertickpos <- maxlower + 0.5*lowerticksize + ((1:length(upperticklabels))*uppertickspacing)
        tickpos <- c(lowertickpos, uppertickpos)
        newy <- as.numeric(y)
        ind <- newy > maxlower
        newy[ind] <- maxlower + 0.5*lowerticksize + 1*uppertickspacing + uppertickspacing*((newy[ind]-minupper) / upperticksize)
        list(newy=newy, breaks=tickpos, labels=ticklabels, limits=range(tickpos))
    }
}


extract_gglegend <- function(p){
    require(ggplot2)
    require(cowplot)

    ## extract the legend from a ggplot object
    tmp <- ggplot_gtable(ggplot_build(p))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    if(length(leg) > 0) leg <- tmp$grobs[[leg]]
    else leg <- NULL
    leg

    ## return the legend as a ggplot object
    legend <- cowplot::ggdraw() + cowplot::draw_grob(grid::grobTree(leg))
    plot <- p + theme(legend.position='none')
    list(plot=plot,legend=legend)
}


table_freq <- function(value) {
    if(is.null(value) | length(value)==0) {
        tbl <- data.table(value=NA,N=NA)
    } else {
        tbl <- as.data.table(table(value))
        tbl <- tbl[order(tbl$N,decreasing=T),]
    }
    tbl
}


annotated_phylo <- function(distance.matrix,groups,method='nj',root_at_normal=T) { 
    ## annotated_phylo: a class based on ape's nj() with extra annotations about the type of each leaf
    require(data.table)
    require(ape)
    require(phangorn)
    require(phytools)
    require(ggtree)

    if(method=='nj') {
        tree <- ape::nj(distance.matrix)
    } else if(method=='upgma') {
        tree <- phangorn::upgma(distance.matrix)
    } else {
        stop('Method must be either nj or upgma')
    }
    ## constructor, takes a phylo class object and adds group information.
    ## groups must be a named character vector, where the character-elements are the groups, and the names are samples (should match tree$tip.label).
    stopifnot(class(tree)=='phylo')
    info <- parse_barcode(tree$tip.label)
    info$order <- 1:nrow(info)
    info <- merge(info, groups, by='barcode', all.x=T)
    info <- info[order(order),]
    colors <- info[!duplicated(group)]
    usual_groups <- c('Normal','Primary','Locoregional','Peritoneum','Lung','Liver','Distant (other)')
    additional_groups <- unique(info$group[!info$group %in% usual_groups])
    ordered_groups <- c(usual_groups, additional_groups)
    ordered_groups <- unique(ordered_groups[ordered_groups %in% info$group])
    colors$group <- factor(colors$group, levels=ordered_groups)
    colors <- colors[order(group)]
    color_scheme <- colors$color
    names(color_scheme) <- as.character(colors$group)
    fields <- c('barcode','order','type','lesion','sample','autopsy','group')
    fields <- fields[fields %in% names(info)]
    info <- info[,(fields),with=F]

    ## attempt to root the tree
    if(root_at_normal==T) { 
        root.barcode=grep('^N',info$barcode,value=T)
        node.number <- which(tree$tip.label==root.barcode)
        if(length(node.number)!=1) {
            n_roots <- length(node.number)
            warning(n_roots,' root samples found! Need exactly one, proceeding unrooted.')
        } else {
            tree <- phytools::reroot(tree, node.number=node.number)
        }
    }

    ## add the information extracted from the barcodes to the tree
    tree$tip.annotations <- info
    tree$distance.matrix <- distance.matrix
    tree$colors <- color_scheme

    ## extend the class of the new object
    class(tree) <- c('annotated_phylo',class(tree)) 
    tree
}


print.annotated_phylo <- function(x, printlen=30) {
    nb.tip <- length(x$tip.label)
    nb.node <- x$Nnode
    cat(paste("\nAnnotated phylogenetic tree with", nb.tip, "tips and",
              nb.node, "internal nodes.\n"))
    rlab <- if (is.rooted(x))
        "Rooted"
    else "Unrooted"
    cat("\n", rlab, "; ", sep = "")
    blen <- if (is.null(x$edge.length))
        "no branch lengths."
    else "includes branch lengths."
    cat(blen, "\n\n", sep = "")
    print(head(x$tip.annotations, printlen))
}


group_samples <- function(input,lun=T,liv=T,per=T,primary_autopsy_is_distant=T,highlight_peritoneum_when_multi=T,color=F) {
    if(any(c('data.frame','matrix') %in% class(input))) {
        barcodes <- rownames(input)
    } else if('character' %in% class(input)) {
        barcodes <- input
    }
    info <- parse_barcode(barcodes)    
    info[grepl('^N[0-9]',barcode) | grepl('^Normal[0-9]',barcode),group:='Normal']
    info[grepl('^P[0-9]',barcode) | grepl('^PT[0-9]',barcode),group:='Primary']
    info[grepl('^L[0-9]',barcode) | grepl('^LN[0-9]',barcode) | grepl('^TD[0-9]',barcode),group:='Locoregional']
    info[grepl('^Lun[0-9]',barcode),group:='Lung']
    info[grepl('^Liv[0-9]',barcode),group:='Liver']
    #info[grepl('^Ad[0-9]',barcode) | grepl('^AD[0-9]',barcode),group:='Adenoma']
    info[grepl('^Per[0-9]',barcode) | grepl('^Di[0-9]',barcode) | grepl('^Om[0-9]',barcode) | grepl('^PerOv[0-9]',barcode),group:='Peritoneum']
    info[is.na(group),group:='Distant (other)']

    if(lun==F) info[group=='Lung',group:='Distant (other)']
    if(liv==F) info[group=='Liver',group:='Distant (other)']
    if(per==F) info[group=='Peritoneum',group:='Distant (other)']
    #if(tumor_deposit_with_lymph==F) info[grepl('^TD[0-9]',barcode),group:='Tumor deposit']
    if(primary_autopsy_is_distant==T) info[group=='Primary' & autopsy==T,group:='Distant (other)']
    out <- info[,c('barcode','group'),with=F]
   
    if(highlight_peritoneum_when_multi) {
        ## default coloring showing all major types, but highlighting peritoneum
        out[group=='Normal',color:='black']
        out[group=='Primary',color:='#008c45']
        out[group=='Locoregional',color:='#eb5b2b']
        out[group=='Liver',color:='#4c86c6']
        out[group=='Lung',color:='#ea6a8c']
        out[group=='Peritoneum',color:='#fab31d']
        out[group=='Distant (other)',color:='#534797']    
    } else {
        ## default coloring showing all major types, but highlighting lung
        out[group=='Normal',color:='black']
        out[group=='Primary',color:='#008c45']
        out[group=='Locoregional',color:='#eb5b2b']
        out[group=='Liver',color:='#4c86c6']
        out[group=='Peritoneum',color:='#ea6a8c'] ## find alternative?
        out[group=='Lung',color:='#fab31d']
        out[group=='Distant (other)',color:='#534797']    
    }

    ## depending on which were included in input argumens, highlight single type
    if(lun==T & liv==F & per==F) {
        out[group=='Lung',color:='#fab31d']
        out[group %in% 'Distant (other)',color:='#4c86c6']

    } else if(lun==F & liv==T & per==F) {
        out[group=='Liver',color:='#fab31d']
        out[group %in% 'Distant (other)',color:='#4c86c6']

    } else if(lun==F & liv==F & per==T) {
        out[group=='Peritoneum',color:='#fab31d']
        out[group %in% 'Distant (other)',color:='#4c86c6']
    }
    if(color==F) out[,color:=NULL]
    out 
}


plot.annotated_phylo <- function(tree,cex=2.5,angle=F,layout='ape',legend.position='none',suppress_tip_labs=F, xpad=0.1, fontface=1, ...){ 
    stopifnot('annotated_phylo' %in% class(tree))
    info <- data.table(label=tree$tip.label,group=tree$tip.annotations$group)
    color_scheme <- tree$colors
    tree$group <- NULL
    class(tree) <- 'phylo'
    p <- ggtree(tree, layout=layout, ...) %<+% info
    if(suppress_tip_labs==F & angle==F) p <- p + geom_tiplab(aes(color=group),fontface=fontface,size=cex,angle=F)
    if(suppress_tip_labs==F & angle==T) p <- p + geom_tiplab(aes(color=group,angle=angle),fontface=fontface,size=cex)
    p <- p + scale_color_manual(values=color_scheme,name=NULL) + theme(legend.position=legend.position)
    xr <- range(p$data$x)
    xwidth <- xr[2]-xr[1]
    xmin <- xr[1] - xpad*xwidth
    xmax <- xr[2] + xpad*xwidth
    p <- p + xlim(c(xmin,xmax))
    if('title' %in% names(tree)) {
        tree_title <- tree$title
        p <- p + labs(title=tree_title)
    }
    p 
}


expand_tree <- function(tree) {
    barcode_order <- tree$tip.label
    x <- as.data.table(reshape2::melt(tree$distance.matrix))
    names(x) <- c('barcode.1','barcode.2','distance')
    anno <- tree$tip.annotations
    if('colors' %in% names(tree)) {
        cols <- data.table(group=names(tree$colors),color=tree$colors)
    } else {
        cols <- data.table(group=unique(anno$group))
        cols$color=rainbow(n=nrow(cols))
    }
    anno <- merge(anno, cols, by='group', all.x=T)
    x <- merge(x, anno, by.x='barcode.1', by.y='barcode', all.x=T)
    setnames(x,c('type','lesion','sample','autopsy','group','color'),c('type.1','lesion.1','sample.1','autopsy.1','group.1','color.1'))
    x <- merge(x, anno, by.x='barcode.2', by.y='barcode', all.x=T)
    setnames(x,c('type','lesion','sample','autopsy','group','color'),c('type.2','lesion.2','sample.2','autopsy.2','group.2','color.2'))
    x <- x[,c('barcode.1','barcode.2','distance','group.1','type.1','lesion.1','sample.1','autopsy.1','color.1','group.2','type.2','lesion.2','sample.2','autopsy.2','color.2'),with=F]
    x$barcode.1 <- factor(x$barcode.1, levels=barcode_order)
    x$barcode.2 <- factor(x$barcode.2, levels=barcode_order)
    x <- x[order(barcode.1,barcode.2),]
    x$barcode.1 <- as.character(x$barcode.1)
    x$barcode.2 <- as.character(x$barcode.2)
    x$order <- 1:nrow(x)
    x
}


contract_tree <- function(x) {
    ## remake the distance matrix
    barcode_order <- unique(x$barcode.1)
    dm <- data.table::dcast(barcode.1 ~ barcode.2, value.var='distance', data=x)
    rows <- dm$barcode.1
    dm[,barcode.1:=NULL]
    dm <- as.matrix(dm)
    rownames(dm) <- rows
    dm <- dm[barcode_order,barcode_order]

    ## remake the groups
    groups <- x[,c('barcode.1','group.1','color.1'),with=F]
    names(groups) <- gsub('[.]1','',names(groups))
    groups <- groups[!duplicated(barcode),]

    ## remake the annotated_phylo object
    tree <- annotated_phylo(dm, groups)
    tree
}


parse_barcode <- function(barcodes) {
    ## extract the main tissue type, lesion, and sample number from each sample's barcode
    .parse_barcode <- function(barcode) {
        str <- strsplit(gsub("([A-Za-z]*)([0-9]*)([A-Za-z]*)", "\\1 \\2 \\3", barcode), " ")[[1]]
        list(barcode=barcode,type=str[1],lesion=str[2],sample=str[3])
    }
    s <- rbindlist(lapply(barcodes, .parse_barcode))
    s[is.na(sample),sample:='']
    s[grepl('-A$',barcode),autopsy:=T]
    s[!grepl('-A$',barcode),autopsy:=F]
    s$sample <- gsub('-A$','',s$sample)
    s
}


collapse_tree_reiter <- function(tree,break_tie_method='least_similar',iter.max=500) {
    # Collapse full phylogenies to the one-sample-per-lesion trees as used in Reiter et al, Nature Genetics 2020.

    select_one <- function(candidates,info,method) {
        get_avg_distance <- function(info) {
            mu <- mean(info$distance,na.rm=T)
            data.table(mu=mu)
        }
        tmp <- info[barcode.1 %in% candidates,get_avg_distance(.SD),by=barcode.1]
        tmp <- tmp[order(mu,decreasing=T),]

        if(method=='least_similar') {
            ## for each candidate sample, get its average distance to all other samples. Keep the candidate that is most far away from all other samples. Remove the other candidates from the data and return this updated data.
            least_similar <- tmp$barcode.1[1]
            drop_tmp <- tmp$barcode.1[!tmp$barcode.1 %in% least_similar]
            if(length(drop_tmp) > 0) info <- info[!barcode.1 %in% drop_tmp & !barcode.2 %in% drop_tmp,]

        } else if(method=='random') {
            random_choice <- sample(tmp$barcode.1,1)
            drop_tmp <- tmp$barcode.1[!tmp$barcode.1 %in% random_choice]
            if(length(drop_tmp) > 0) info <- info[!barcode.1 %in% drop_tmp & !barcode.2 %in% drop_tmp,]    
        }
        info
    }
   
    get_comparison_id <- function(tmp) {
        b1 <- tmp$barcode.1
        b2 <- tmp$barcode.2
        b <- sort(c(b1,b2))
        id <- paste(b[1],b[2],sep=':')
        tmp$comparison <- id
        tmp$intra_lesion <- tmp$main_tissue.1==tmp$main_tissue.2 & tmp$lesion.1==tmp$lesion.2 & tmp$sample.1!=tmp$sample.2
        tmp
    }

    get_clades_for_lesion <- function(lesion, info) {
        ## get the number of nodes separating every pair of samples on the nj tree
        mat <- contract_distance_matrix(info)
        tree <- nj(mat) 
        all_barcodes_for_lesion <- info[lesion.1==lesion & lesion.2==lesion,c('barcode.1','barcode.2'),with=F]
        n <- nrow(all_barcodes_for_lesion)

        if(n < 2) {
            list(largest_clades=NULL, samples_not_in_largest_clades=NA)
        } else {
            check_same_clade <- function(i,all_barcodes_for_lesion,tree,info) {
                out <- all_barcodes_for_lesion[i,]
                sample.1 <- out$barcode.1
                sample.2 <- out$barcode.2
                if(sample.1==sample.2) {
                    in_clade <- T
                } else {
                    sample_lesion <- info$lesion.1[info$barcode.1==sample.1 & info$barcode.2==sample.2]
                    node <- fastMRCA(sp1=sample.1,sp2=sample.2,tree=tree)
                    clade <- extract.clade(tree,node)
                    samples_in_clade <- clade$tip.label
                    lesions_in_clade <- unique(info$lesion.1[info$barcode.1 %in% samples_in_clade & info$barcode.2 %in% samples_in_clade])
                    in_clade <- ifelse(all(lesions_in_clade %in% sample_lesion),T,F)
                }
                out$in_clade <- in_clade
                out
            }
            clades_dat <- rbindlist(lapply(1:n, check_same_clade, all_barcodes_for_lesion, tree, info))
            clades_dat$distance <- as.integer(clades_dat$in_clade==F)
            clades_dat <- contract_distance_matrix(clades_dat) ## 1 means separate, 0 means together
            f=function(s) paste(s,collapse='')
            distances <- apply(clades_dat, 1, f)
            clade_size <- table(distances)
            max_clade_size <- max(clade_size)
            if(max_clade_size > 1) {
                largest_clades <- names(clade_size[clade_size==max_clade_size]) ## there may be multiple clades of the same maximum size
                get_samples_in_clade <- function(largest_clade, distances) {
                    names(distances[distances==largest_clade])            
                }
                largest_clades <- lapply(largest_clades, get_samples_in_clade, distances) 
            } else {
                ## in case multiple samples are all the same distance apart
                largest_clades <- list(names(distances))
            }
            samples_in_largest_clades <- unlist(largest_clades)
            all_samples <- rownames(clades_dat)
            samples_not_in_largest_clades <- all_samples[!all_samples %in% samples_in_largest_clades]
            list(largest_clades=largest_clades, samples_not_in_largest_clades=samples_not_in_largest_clades)
        }
    }

    contract_distance_matrix <- function(info) {
        ## recover the distance matrix from the expanded (long) table
        dm <- data.table::dcast( barcode.1 ~ barcode.2, value.var='distance', data=info)
        rows <- dm$barcode.1
        dm[,barcode.1:=NULL]
        dm <- as.matrix(dm)
        rownames(dm) <- rows    
        dm
    }

    ## repeat this process until it cannot be collapsed any further (force stop at 500 iterations)
    size <- length(tree$tip.label)
    run <- T 
    r <- 0
    while(r < iter.max & run==T) { 
        ## expand the distance matrix
        info <- expand_tree(tree)
        info$lesion.1 <- paste0(info$type.1,info$lesion.1)
        info$lesion.2 <- paste0(info$type.2,info$lesion.2)
        info$groupl.1 <- tolower(info$group.1)
        info$groupl.2 <- tolower(info$group.2)
        if(sum(info$groupl.1=='normal')==0) warning('No normal samples found!')
        if(sum(info$groupl.1=='primary')==0) warning('No primary samples found!')

        ## any samples not organized into the specified groups will be grouped according to the type from their barcodes
        info[groupl.1=='',groupl.1:=type.1]
        info[groupl.2=='',groupl.2:=type.2]

        ## remove extra germline normals
        normals <- sort(unique(info$barcode.1[info$groupl.1=='normal']))
        if(length(normals) > 1) info <- select_one(normals,info,method=break_tie_method)

        ## remove non-cancer samples
        info <- info[!groupl.1 %in% 'adenoma' & !groupl.2 %in% 'adenoma',]

        ## prune non-representative samples for lesions. in the case of lesions with multiple samples where the majority are in a single clade, remove the lesion's other samples which aren't on this clade.
        multiple_samples <- which(info$groupl.1==info$groupl.2 & info$lesion.1==info$lesion.2 & info$sample.1!=info$sample.2 & !info$groupl.1 %in% c('normal','primary'))
        lesions_with_multiple_samples <- sort(unique(info$lesion.1[multiple_samples]))

        for(l in lesions_with_multiple_samples) { 
            clade_info <- get_clades_for_lesion(l,info)
            to_drop <- clade_info$samples_not_in_largest_clades
            to_drop <- to_drop[!is.na(to_drop)]
            clades_to_collapse <- clade_info$largest_clades
            num_clades_to_collapse <- length(clades_to_collapse)
            if(length(to_drop) > 0) { 
                message(' -> Pruning: ',paste(to_drop,collapse=', '))
                info <- info[!barcode.1 %in% to_drop & !barcode.2 %in% to_drop,]
            }

            for(i in 1:num_clades_to_collapse) { 
                candidates <- clades_to_collapse[[i]]
                if(length(candidates) > 1) {
                    message('Collapsing clade: ',paste(candidates,collapse=', '))
                    info <- select_one(candidates,info,method=break_tie_method)
                } 
            }
        }
        tree <- contract_tree(info)
        newsize <- length(tree$tip.label)
        if(newsize==size) {
            run <- F
        } else {
            size <- newsize
            r <- r+1
        }
    }
    if(r >= iter.max) warning('Maximum iterations reached!')
    tree
}


collapse_tree_all <- function(tree) {
    # Take a distance matrix which may contain multiple distinct samples per lesion and return a list of all possible matrices where each lesion is represented by a single sample.
    # param tree A annotated_phylo class object.
    # return a list of all unique possible trees where each lesion is represented by a single sample

    ## get expanded data
    x <- expand_tree(tree)
    x$lesion.1 <- paste0(x$type.1,x$lesion.1)
    x$lesion.2 <- paste0(x$type.2,x$lesion.2)
    x$groupl.1 <- tolower(x$group.1)
    x$groupl.2 <- tolower(x$group.2)

    ## split tip annotations into Primary and Non-primary samples
    orig_info <- tree$tip.annotations
    orig_matrix <- tree$distance.matrix
    info <- tree$tip.annotations[group!='Primary',]
    info$lesion <- paste0(info$type,info$lesion)
    infoP <- tree$tip.annotations[group=='Primary',]
    infoP$lesion <- paste0(infoP$type,infoP$lesion)
  
    ## get a list of samples per each lesion
    lesions <- unique(info$lesion)
    f <- function(lesion, info) info$barcode[info$lesion==lesion]
    l <- lapply(unique(lesions), f, info)
    names(l) <- lesions

    ## process the primary tumors the same way, but here we force it to treat every primary sample as its own lesion
    barcodes <- unique(infoP$barcode)
    fP <- function(barcode, infoP) infoP$barcode[infoP$barcode==barcode]
    lP <- lapply(unique(barcodes), fP, infoP)
    names(lP) <- infoP$barcode
    l <- append(l, lP)

    ## get all combinations of samples where each lesion is represented by a single sample
    grid <- as.matrix(expand.grid(l))
    n_subsets <- nrow(grid)
    
    ## get the subset distance matrix for samples in each combination
    subset_tree <- function(i, grid, orig_matrix) {
        sample_subset <- as.character(grid[i,])

        ## subset distance matrix
        subset_matrix <- orig_matrix[sample_subset,sample_subset] 
        subset_groups <- group_samples(rownames(subset_matrix),liv=T,lun=T,per=T,primary_autopsy_is_distant=T,color=T)
        subset_tree <- annotated_phylo(subset_matrix,subset_groups)
        subset_tree
    }
    subset_trees <- lapply(1:n_subsets, subset_tree, grid, orig_matrix)
    subset_trees
}



