##' d2m
##' @export
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



##' distance_matrix_correlation
##' @export
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



##' dist_similarity
##' @export
dist_similarity <- function(test_dist, ref_dist, nperm=100, cpus=1, return_only_pval=T, method='spearman') {

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

    if(cpus==1) {
        permuted <- unlist(lapply(1:nperm, shuffled_correlation, test_dist, ref_dist, return_only_r=T, method=method))
    } else {
        require(parallel)
        permuted <- unlist(mclapply(1:nperm, shuffled_correlation, test_dist, ref_dist, return_only_r=T, mc.cores=cpus, method=method))
    }
    
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


##' read_distance_matrix
##' @export
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


##' write_tsv
##' @export
write_tsv <- function(d,file,sep='\t',quote=F,row.names=F,...) write.table(d,file=file,sep=sep,quote=quote,row.names=row.names,...)


##' write_distance_matrix
##' @export
write_distance_matrix <- function(filepath,dm_df) {
    write.table(dm_df,file=filepath,sep="\t",quote=FALSE,col.names=NA)
}


##' theme_ang
##' @export
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


##' break_axis
##' @export
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


##' extract_gglegend
##' @export
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



##' table_freq
##' @export
table_freq <- function(value) {
    if(is.null(value) | length(value)==0) {
        tbl <- data.table(value=NA,N=NA)
    } else {
        tbl <- as.data.table(table(value))
        tbl <- tbl[order(tbl$N,decreasing=T),]
    }
    tbl
}



