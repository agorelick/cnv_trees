
sample_type7 <- c("#55823b","#bf2026","#2f75b6","#df812b","#bf61a6","#6f3996","#395829","#000000","#bfbfbf")
names(sample_type7) <- c("Normal","Primary","Locoregional","Liver","Lung","Peritoneum","Distant (other)","Other","n/a")

sample_type4 <- c("#55823b","#bf2026","#2f75b6","#395829","#000000","#bfbfb")
names(sample_type4) <- c("Normal","Primary","Locoregional","Distant","Other","n/a")

colors <- list(sample_type7=sample_type7, sample_type4=sample_type4)




get_cnv_segments <- function(obj_list, fit_file, sex, min_segment_bins=10, field='intcopy') {
    options(scipen = 99)
    segments <- list()
    bins <- list()

    fits <- fread(fit_file)
    if(!all(names(obj_list)==fits$sample)) stop('Names of the obj_list should match the *sample* column of the fit file!')

    samples <- names(obj_list)    
    chr_info <- sex_chr_info(sex)
    included_samples <- intersect(names(obj_list),samples)
    obj_list <- obj_list[included_samples]
    fits <- fits[sample %in% included_samples,]

    ## loop through data object and collect relevant bin-level information
    for (i in 1:length(included_samples)){
        sample <- included_samples[i]
        object <- obj_list[[i]]
        purity <- fits$purity[fits$sample==sample]
        ploidy <- fits$ploidy[fits$sample==sample]
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
        ind <- which(fits$sample==this.sample)
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
        #browser()
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

    list(mat=m, bins=bins, segs=segs)
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



write_tsv <- function(d,file,sep='\t',quote=F,row.names=F,...) write.table(d,file=file,sep=sep,quote=quote,row.names=row.names,...)





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
        message(title)
        p
    } else {
        refits_dir=file.path(output_dir)
        if(!dir.exists(refits_dir)) dir.create(refits_dir,recursive=T)

        ## save the plot for this fit
        plot_file <- paste0(refits_dir,'/',samplename,'_N=',ploidy,'_cellularity=',purity,'.pdf')
        message('Saving plot: ',plot_file)
        ggsave(p,dev=cairo_pdf,filename=plot_file,width=8,height=5)

        ## save data (append it to existing file along with the date of the entry)
        fit_file <- paste0(refits_dir,'/purity_ploidy.txt')
        if(!file.exists(fit_file)) cat('sample\tpurity\tploidy\tdate\n',file=fit_file,append=F) 
        when <- as.character(format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
        message('Saving data: ',fit_file)
        cat(paste(samplename,purity,ploidy,when,sep='\t'),'\n',file=fit_file,append=T)
    }
}
