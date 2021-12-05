##' annotated_phylo
##' @export
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
    colors$group <- factor(colors$group, levels=c('Normal','Primary','Locoregional','Peritoneum','Lung','Liver','Distant (other)'))
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

##' print.annotated_phylo
##' @export
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


##' group_samples
##' @export
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


##' plot.annotated_phylo
##' @export
plot.annotated_phylo <- function(tree,cex=2.5,angle=F,layout='ape',legend.position='none',suppress_tip_labs=F, xpad=0.1, ...){ 
    stopifnot('annotated_phylo' %in% class(tree))
    info <- data.table(label=tree$tip.label,group=tree$tip.annotations$group)
    color_scheme <- tree$colors
    tree$group <- NULL
    class(tree) <- 'phylo'
    p <- ggtree(tree, layout=layout, ...) %<+% info
    if(suppress_tip_labs==F & angle==F) p <- p + geom_tiplab(aes(color=group),fontface=1,size=cex,angle=F)
    if(suppress_tip_labs==F & angle==T) p <- p + geom_tiplab(aes(color=group,angle=angle),fontface=1,size=cex)
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

##' expand_tree
##' @export
expand_tree <- function(tree) {
    barcode_order <- tree$tip.label
    x <- as.data.table(reshape2::melt(tree$distance.matrix))
    names(x) <- c('barcode.1','barcode.2','distance')
    anno <- tree$tip.annotations
    x <- merge(x, anno, by.x='barcode.1', by.y='barcode', all.x=T)
    setnames(x,c('type','lesion','sample','autopsy','group'),c('type.1','lesion.1','sample.1','autopsy.1','group.1'))
    x <- merge(x, anno, by.x='barcode.2', by.y='barcode', all.x=T)
    setnames(x,c('type','lesion','sample','autopsy','group'),c('type.2','lesion.2','sample.2','autopsy.2','group.2'))
    x <- x[,c('barcode.1','barcode.2','distance','group.1','type.1','lesion.1','sample.1','autopsy.1','group.2','type.2','lesion.2','sample.2','autopsy.2'),with=F]
    x$barcode.1 <- factor(x$barcode.1, levels=barcode_order)
    x$barcode.2 <- factor(x$barcode.2, levels=barcode_order)
    x <- x[order(barcode.1,barcode.2),]
    x$barcode.1 <- as.character(x$barcode.1)
    x$barcode.2 <- as.character(x$barcode.2)
    x$order <- 1:nrow(x)
    x
}

##' contract_tree
##' @export
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
    tmp <- x[,c('barcode.1','type.1','lesion.1','sample.1','autopsy.1','group.1'),with=F]
    names(tmp) <- gsub('[.]1','',names(tmp))
    groups <- tmp$barcode
    names(groups) <- tmp$group
    
    ## remake the annotated_phylo object
    tree <- annotated_phylo(dm, groups)
    tree
}

##' parse_barcode
##' @export
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

##' collapse_tree_reiter
##' @export
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

##' collapse_tree_all
##'
##' Take a distance matrix which may contain multiple distinct samples per lesion and return a list of all possible matrices where each lesion is represented by a single sample.
##' @param tree A annotated_phylo class object.
##' @export
collapse_tree_all <- function(tree) {
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

        ## subset annotations
        subset_info <- orig_info[barcode %in% sample_subset,]
        subset_info <- subset_info[!duplicated(barcode),]

        ## subset annotated_phylo
        subset_groups <- orig_info$barcode
        names(subset_groups) <- orig_info$group
        subset_tree <- annotated_phylo(subset_matrix,subset_groups)
        subset_tree
    }
    subset_trees <- lapply(1:n_subsets, subset_tree, grid, orig_matrix)
    subset_trees
}



