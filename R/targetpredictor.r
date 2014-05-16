#######
#	To be a module mining multiple sources of target predictions
#	to start with - mouse - Mus musculus or human - Homo Sapiens
#
#######

require(stringr)		#string operations, essential
require(plyr) 			#dataframe operations, do i actually use it?
#require(RSQLite)		#local sqlite db, actually I'm not using it here

require(targetPredictor.db) #DATA


options(sqldf.driver = "SQLite") # as per FAQ #7 force SQLite
options(gsubfn.engine = "R")


#####################
#####################
#	MAIN FUNCTION
#####################
#####################



#' Get aggregated ordered list of predicted targets for miRNA
#'
#' This method performs aggregation of target lists from multiple sources.
#' Aggregated list is more accurate than any list from a single source.
#' Multiple aggregation methods are available.
#' 
#'
#' @param mirna mirna in a standard format
#' @param sources a list of sources to use for aggregation, default c('pictar','diana','targetscan','miranda')
#' @param species species in a standard three-letter acronym, default 'mmu'
#' @param min_src minimum number of sources required for a target to be considered, default 2
#' @param method method of aggregation - choose from c('min', 'max', 'geom'), default 'min' is a minimum of ranks, 'max' is a maximum of ranks, and 'geom' is based on geometric mean of the ranks, it proves to be the most accurate.
#' @param promote add weights to improve accuracy of the method, default TRUE
#' @param ... any optional arguments
#' @return getPredictedTargets a data.frame object where row names are entrez IDs of target genes
#' @export
#' @author Maciej Pajak \email{m.pajak@@sms.ed.ac.uk}
#' @references
#' Friedman, R. C., Farh, K. K.-H., Burge, C. B., and Bartel, D. P. (2009). Most mammalian mRNAs are conserved targets of microRNAs. Genome research, 19(1):92-105.
#' @references
#' Griffiths-Jones, S., Saini, H. K., van Dongen, S., and Enright, A. J. (2008). miRBase: tools for microRNA genomics. Nucleic acids research, 36(Database issue):D154-8.
#' @references
#' Lall, S., Grun, D., Krek, A., Chen, K., Wang, Y.-L., Dewey, C. N., ... Rajewsky, N. (2006). A genome-wide map of conserved microRNA targets in C. elegans. Current biology : CB, 16(5):460-71.
#' @references
#'	Maragkakis, M., Vergoulis, T., Alexiou, P., Reczko, M., Plomaritou, K., Gousis, M., ... Hatzigeorgiou, A. G. (2011). DIANA-microT Web server upgrade supports Fly and Worm miRNA target prediction and bibliographic miRNA to disease association. Nucleic Acids Research, 39(Web Server issue), W145-8.
#' @examples
#' ##Not run:
#' targets = getPredictedTargets('let-7a',species='mmu') 
#' head(targets) #top of the list
getPredictedTargets = function(mirna,sources=c('pictar','diana','targetscan','miranda'),species='mmu', min_src=2, method='min', promote=TRUE, ...) {
	
	if (!(species %in% c('mmu','hsa','rno','dme'))) {
		message(paste('species ',species,' not supported in the current version',sep=''))
		return(NULL)
	}
	 
	conL = get_connection()
	
	n_sources = length(sources)
	
	message(paste('sources: ',n_sources,sep=''))
	if (n_sources<min_src) {
		#!! either 'repair' it or inform the user that they're stupid
		# well in this case rather the latter, right?
		message('min_src > n_sources!')
		return(NULL)
	}
	
	if (min_src<=0) {
		message('min_src changed to 1')
		min_src = 1
	}
	
	if (nchar(mirna)<3) {
		message('unrecognised miRNA')
		return(NULL)
	}
	
	#h_outputs = hash()
	
	l_outputs = list()
	
	cols=c('GeneID','score')
	
	if (substr(mirna,1,3)==species) {
		mirna=substr(mirna,5,15)
	}
	
	#!!extract ifs below into a wrapper on getters DONE
	#this is the important step - obtaining all target lists
	for (src in sources) {
	#!! but what if theres just one source
		l_outputs[[src]]=getTargetsFromSource(mirna,species,conL,source=src)
	}
	
	
	#watch out - it creates non-unique colnames!
	merged.scores = suppressWarnings(Reduce(function(...) merge_rename(..., by='GeneID', all=TRUE), l_outputs))
	
	#print(merged.scores)
	#print(dim(merged.scores))
	
	if (dim(merged.scores)[1]<1) { # 
		message(paste('no targets found for mirna ',mirna,sep=''))
		return(NULL)
	}
	
	valid_srcs = (1:(n_sources+1))[colSums(merged.scores,na.rm = TRUE)>0]
	#print(valid_srcs)
	
	n_valid_srcs = length(valid_srcs)-1
	#print(n_valid_srcs)
	
	if (n_valid_srcs<min_src) {
		#!!technically that's not user's mistake, it should get repaired automatically
		#and some sort of warning should be given
		message(paste('algos which returned a target list<min_src, min_src reduced to ',n_valid_srcs,sep=''))
		min_src = n_valid_srcs
		#return(NULL)
	}	
	
	merged.scores = merged.scores[,valid_srcs]
	
	if (n_valid_srcs>1) {
		target_found = rowSums(!is.na(merged.scores[,2:(n_valid_srcs+1)])) #n_sources
	} else {
		target_found = 1*(!is.na(merged.scores[,2:(n_valid_srcs+1)])) #n_sources
	}
	#max_appears=max(target_found) #possibly? used it in the paper
	
	
	merged.scores2 = merged.scores[target_found>=min_src,]
	
	#duplicate removal needs to go at the level of a single result list
	#!! so it means that we don't need this line...
	merged.scores3 = merged.scores2[!duplicated(merged.scores2$GeneID),]
	
	merged.scores3 = merged.scores3[!is.na(merged.scores3$GeneID),]
  
	row.names(merged.scores3) = merged.scores3$GeneID
	
	#RELEVANT LINES COPIED FROM RP
	data1 = merged.scores3[,2:(n_valid_srcs+1)] #n_sources
	
	data1 = as.data.frame(data1)
	
	#data2 = data1[,NULL]
	
	num.gene = dim(data1)[1]
	
	#message(dim(data1))
	
	
	if (num.gene<1) {
		message('not enough targets overlapping between sources, reduce min_srcs parameter or add sources')
		return(NULL)
	} else if (num.gene==1) {
		message('only one target found')
		result = data.frame(matrix(1,nrow=1,ncol=n_valid_srcs+2))
		colnames(result)=c(paste('source_',1:n_valid_srcs,sep=''),'rank_product','rank_final')
		row.names(result) = row.names(data1)
		return(result)
	}
	
	#we need to reverse rank them
	rank.rep = apply(data1, 2, function(x) (num.gene+1) - rank(x, ties.method='average', na.last = FALSE))
	num.rank = apply(is.na(data1) == FALSE, 1, sum)
	
	
	
	if (method=='geom') {
		rank.rep[is.na(data1)] = 1
		
		if (promote==TRUE) {
			rank.prod = (1/num.rank)*(apply(rank.rep, 1, prod))^(1/num.rank)
		} else {
			rank.prod = (apply(rank.rep, 1, prod))^(1/num.rank)
		}
		rank.prod[num.rank == 0] = NA #icing - all NAs won't happen anyway, !!so delete?
		
		rank.ori <- rank(rank.prod)
		
		
		rank.rep[is.na(data1)] = NA
		result = cbind(rank.rep,rank.prod,rank.ori)  #,pval.downin2
 		result = result[order(result[,'rank.ori']),]
	
	}  else if (method=='max') {
	
		rank.rep[is.na(data1)] = NA
		rank.prod = apply(rank.rep, 1, max, na.rm=TRUE)
		rank.ori <- rank(rank.prod)
		result = cbind(rank.rep,rank.prod,rank.ori)
		result = result[order(result[,'rank.ori']),]
	
	} else  { #even if it's something invalid, still use min
	
		rank.rep[is.na(data1)] = NA
		rank.prod = apply(rank.rep, 1, min, na.rm=TRUE)
		rank.ori <- rank(rank.prod)
		result = cbind(rank.rep,rank.prod,rank.ori)
		result = result[order(result[,'rank.ori']),]
			
	} 
	
	
	dbDisconnect(conL)
	
	colnames(result)=c(paste('source_',1:n_valid_srcs,sep=''),'rank_product','rank_final') #coating
	return(result)
}




#####################
#	
#####################

#' Auxiliary internal TargetPredictor functions
#'
#' This function is not meant to be called directly by the user
#'
#' @param ... parameters to merge
#' @export
#' @author Maciej Pajak \email{m.pajak@@sms.ed.ac.uk}	
merge_rename = function(...) {
	merged = merge(...)
	colnames(merged) = c(colnames(merged)[1],paste('V',1:(length(names(merged))-1),sep=''))
	return(merged)
}

####################
#	useless stuff
####################
#bin?
# generate_csv = function(mir,species='mmu',...) {
	
	# ao = getPredictedTargets(mir, species=species,method='geom',min_src=2)
	# xo = as.data.frame(ao)
	# xo$entrez = row.names(ao)
	# xo = xo[,c('entrez','rank.prod','rank.ori','pval.downin2')]
	
	# filename = paste(species,'-',mir,'.csv',sep='')
	# write.csv(xo, filename, row.names=FALSE)
	
# }
