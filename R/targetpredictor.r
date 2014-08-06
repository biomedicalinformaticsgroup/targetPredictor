#######
#  module integrating mirna target predictions for human, mouse and rat
#######

require(stringr)    #string operations
require(plyr)       #dataframe operations
require(targetPredictor.db) #data

options(gsubfn.engine = "R")


#' @docType package
#' @name targetPredictor
#' @title Integrating miRNA target predictions from multiple sources
#' @details It is a package with tools to facilitate implementation of workflows requiring miRNA prediction through access to multiple prediction results (DIANA, Targetscan, PicTar and Miranda) and their aggregation. 
#' Three aggregation methods are available: minimum, maximum and geometric mean, additional parameters provide further tuning of the results. 
#' Predictions are available for Homo sapiens, Mus musculus and Rattus norvegicus (the last one through homology translation).
#' @import plyr stringr targetPredictor.db
#' @author Maciej Pajak \email{m.pajak@@sms.ed.ac.uk}, Ian Simpson 
#' @examples
#' #direct targets in mouse aggregated from all sources:
#' targets_mouse <- getPredictedTargets('let-7a',species='mmu', method='geom') 
#' #homology-translated targets in rat aggregated from all sources
#' targets_rat <- getPredictedTargets('let-7a',species='mmu', method='geom') 
NULL


#####################
#  MAIN FUNCTION
#####################

#' @title Get aggregated ordered list of predicted targets for miRNA
#'
#' @details This method performs aggregation of target lists from multiple sources.
#' Aggregated list is more accurate than any list from a single source.
#' Multiple aggregation methods are available.
#' Direct target data from four sources for Human and Mouse is supplied through targetPredictor.db package, for Rat targets are derived through homology translations whereever direct ones are not available.
#'
#' @param mirna miRNA in a standard format
#' @param sources a list of sources to use for aggregation, default c('pictar','diana','targetscan','miranda')
#' @param species species in a standard three-letter acronym, 'mmu' and 'hsa' available as direct targets, 'rno' as homology translations, default 'mmu'
#' @param min_src minimum number of sources required for a target to be considered, default 2
#' @param method method of aggregation - choose from 'min', 'max', 'geom', default 'min' is a minimum of ranks, 'max' is a maximum of ranks, and 'geom' is based on geometric mean of the ranks, it proves to be the most accurate.
#' @param promote add weights to improve accuracy of the method, default TRUE
#' @param ... any optional arguments
#' @return a data.frame object where row names are entrez IDs of target genes, ranks from individual sources and aggregated rank are shown in columns
#' @export
#' @author Maciej Pajak \email{m.pajak@@sms.ed.ac.uk}
#' @references
#' Friedman, R. C., Farh, K. K.-H., Burge, C. B., and Bartel, D. P. (2009). Most mammalian mRNAs are conserved targets of microRNAs. Genome research, 19(1):92-105.
#' @references
#' Griffiths-Jones, S., Saini, H. K., van Dongen, S., and Enright, A. J. (2008). miRBase: tools for microRNA genomics. Nucleic acids research, 36(Database issue):D154-8.
#' @references
#' Lall, S., Grun, D., Krek, A., Chen, K., Wang, Y.-L., Dewey, C. N., ... Rajewsky, N. (2006). A genome-wide map of conserved microRNA targets in C. elegans. Current biology : CB, 16(5):460-71.
#' @references
#'  Maragkakis, M., Vergoulis, T., Alexiou, P., Reczko, M., Plomaritou, K., Gousis, M., ... Hatzigeorgiou, A. G. (2011). DIANA-microT Web server upgrade supports Fly and Worm miRNA target prediction and bibliographic miRNA to disease association. Nucleic Acids Research, 39(Web Server issue), W145-8.
#' @examples
#' ##Not run:
#' targets <- getPredictedTargets('let-7a',species='mmu') 
#' head(targets) #top of the list with minimum aggregation
#' targets2 <- getPredictedTargets('let-7a',species='mmu', method='geom') 
#' head(targets2) #top of the list with geometric mean aggregation
getPredictedTargets <- function(mirna, sources=c('pictar','diana','targetscan','miranda'), species = 'mmu', min_src = 2, method = 'min', promote = TRUE, ...) {
  
  if (!(species %in% c('mmu','hsa','rno','dme'))) {
    warning(paste('species ',species,' not supported in the current version',sep=''))
    return(NULL)
  }
   
  #conL <- targetPredictor.db:::.getTpConnection()
  n_sources <- length(sources)
  
  if (n_sources<min_src) {
    warning('min_src > n_sources!')
    return(NULL)
  }
  
  if (min_src<=0) {
    warning('min_src changed to 1')
    min_src <- 1
  }
  
  if (nchar(mirna)<3) {
    warning('unrecognised miRNA')
    return(NULL)
  }

  
  l_outputs <- list()
  
  cols <- c('GeneID','score')
  
  if (substr(mirna,1,3)==species) {
    mirna=substr(mirna,5,15)
  }

  for (src in sources) {
    l_outputs[[src]] <- getTargetsFromSource(mirna, species, source = src)
  }
  
  #it creates non-unique colnames, hence warning surpression
  merged.scores <- suppressWarnings(Reduce(function(...) .mergeRename(..., by='GeneID', all=TRUE), l_outputs))
  
  if (dim(merged.scores)[1]<1) { # 
    warning(paste('no targets found for mirna ', mirna, sep = ''))
    return(NULL)
  }
  
  valid_srcs <- (1:(n_sources+1))[colSums(merged.scores,na.rm = TRUE)>0]
  
  n_valid_srcs <- length(valid_srcs)-1
  
  if (n_valid_srcs<min_src) {
    warning(paste('algos which returned a target list<min_src, min_src reduced to ',n_valid_srcs,sep=''))
    min_src <- n_valid_srcs
  }  
  
  merged.scores <- merged.scores[,valid_srcs]
  
  #take over by ranking function.
  result <- aggregateRanks(merged.scores, n_valid_srcs, min_src, method = method, promote=promote)
  
  return(result)
}


####################
#  Ranking function
####################
#' @title Auxiliary targetPredictor functions
#'
#' @details This function performs aggregation phase of target prediction in getPredictedTargets()
#' @param ranks data.frame with ordered scored
#' @param n_valid_srcs number of valid sources in the dataset
#' @param min_srcs minimum acceptable number fo sources
#' @param method 'min','max', or 'geom'
#' @param promote add weights to improve accuracy of the method, default TRUE
#' @return data.frame object with aggregate ranks
#' @author Maciej Pajak \email{m.pajak@@sms.ed.ac.uk} 
aggregateRanks <- function(ranks, n_valid_srcs, min_src, method = 'min', promote=TRUE) {

    if (n_valid_srcs>1) {
    target_found <- rowSums(!is.na(ranks[,2:(n_valid_srcs+1)])) #n_sources
  } else {
    target_found <- 1*(!is.na(ranks[,2:(n_valid_srcs+1)])) #n_sources
  }
  
  ranks <- ranks[target_found>=min_src,]
  ranks <- ranks[!duplicated(ranks$GeneID),]
  ranks <- ranks[!is.na(ranks$GeneID),]
  
  row.names(ranks) <- ranks$GeneID
  
  data1 <- ranks[,2:(n_valid_srcs+1)] #n_sources
  data1 <- as.data.frame(data1)
  num.gene <- dim(data1)[1]
  
  if (num.gene<1) {
    warning('not enough targets overlapping between sources, reduce min_srcs parameter or add sources')
    return(NULL)
  } else if (num.gene==1) {
    warning('only one target found')
    result <- data.frame(matrix(1,nrow=1,ncol=n_valid_srcs+2))
    colnames(result) <- c(paste('source_',1:n_valid_srcs,sep=''),'rank_product','rank_final')
    row.names(result) <- row.names(data1)
    return(result)
  }
  
  #we need to reverse rank them
  rank_ind <- apply(data1, 2, function(x) (num.gene+1) - rank(x, ties.method='average', na.last = FALSE))
  
  
  if (method=='geom') {
        result = .aggregateGeom(data1, rank_ind, promote)
    }  else if (method=='max') {
        result = .aggregateMinMax(data1, rank_ind, minmax=method)
    } else  { #even if it's something invalid, still use min
        result = .aggregateMinMax(data1, rank_ind, minmax=method)
    } #more methods possible to add, so far geom suffices and scores well
  
    colnames(result) <- c(paste('source_',1:n_valid_srcs,sep=''),'rank_product','rank_final') #coating
    return(result)

}



#####################
#  Auxiliary
#####################
 

.aggregateGeom <- function(data1, rank_ind, promote=TRUE) {
    rank_ind[is.na(data1)] <- 1
    num.rank <- apply(is.na(data1) == FALSE, 1, sum)
    if (promote==TRUE) {
      rank_mean <- (1/num.rank)*(apply(rank_ind, 1, prod))^(1/num.rank)
    } else {
      rank_mean <- (apply(rank_ind, 1, prod))^(1/num.rank)
    }
    rank_mean[num.rank == 0] <- NA
    rank_final <- rank(rank_mean)
    rank_ind[is.na(data1)] <- NA
    result <- cbind(rank_ind,rank_mean,rank_final)  
    result <- result[order(result[,'rank_final']),]
    return(result)
}

.aggregateMinMax <- function(data1, rank_ind, minmax='min') {
    rank_ind[is.na(data1)] <- NA
    if (minmax=='max') {
        rank_mean <- apply(rank_ind, 1, max, na.rm=TRUE)
    } else {
        rank_mean <- apply(rank_ind, 1, min, na.rm=TRUE)
    }
    rank_final <- rank(rank_mean)
    result <- cbind(rank_ind,rank_mean,rank_final)
    result <- result[order(result[,'rank_final']),]
}


.mergeRename <- function(...) {
  merged <- merge(...)
  colnames(merged) <- c(colnames(merged)[1],paste('V',1:(length(names(merged))-1),sep=''))
  return(merged)
}

