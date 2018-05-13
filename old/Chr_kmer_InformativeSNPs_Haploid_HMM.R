#### Packages ####
library(data.table)
library(reshape2)
library(gtools)
library(splitstackshape)
library(ggplot2)
library(useful)

############## Functions ############

## Create a 3d matrix with the conditional probabilities by SNP ##
SelectInformativeSNPs <- function(Chr_Strands, Prob_Pops, diff_quantile = 0.1){
	## Prob_pops is a one-dimensional vector with the probability of each population


	print("Finding Informative SNPs")
	Alleles <- data.frame(REF = "0", ALT = "1")


	## Initially filled with 1's. Number of rows equal to number of populations (with each row labeled as the name of a particular population),
	## number of columns equal to number of Alleles (with each column labeled as the name of an allele),
	## and number of ??? equal to number of SNPs (with each SNP labeled as the ID of the corresponding SNP)
	conditional_prob_matrices <-array(1,c(length(Prob_Pops),ncol(Alleles),nrow(Chr_Strands)),dimnames=list(Pop=colnames(ADMIX2),Allele=colnames(Alleles),SNP=Chr_Strands$ID))


	## Calculate population conditional probabilities of Alleles by SNP #
	for (i in 1:ncol(Alleles)){ ## Is it necessary to iterate through each Allele? Or do we only really use information from the REF allele


		## RowSums is a vector where the Nth value is equal to the number of columns in row N that have Allele i
		## Prob_Allele is a vector with the same length as RowSums, where each value is initialized to 0
		## (each item in RowSums and Prob_Allele corresponds to one gene)
		Allele_df <- data.frame(RowSums=rowSums(Chr_Strands[,10:ncol(Chr_Strands)] == as.character(Alleles[,i])),
													Prob_Allele=numeric(nrow(Chr_Strands)), row.names=Chr_Strands$ID)

		## Nth value in Prob_Allele is now equal to the Nth value of RowSums, divided by number of people?
		## This is equal to the probability of a random person having that allele for Nth gene
		Allele_df$Prob_Allele <- Allele_df$RowSums/(nrow(ADMIX2)) ## Prob(Allele[i])


		for (j in 1:length(Prob_Pops)){
			## sum_over_individuals(#Allele[i] per Individual * %Pop[j] per Individual)
			Pop_Sum <- rowSums(sweep(Chr_Strands[,10:ncol(Chr_Strands)] == as.character(Alleles[,i]),MARGIN=2,ADMIX2[,j],'*',check.margin=FALSE))
			Pop_WeightedMean <- ifelse(Pop_Sum==0 & Allele_df$RowSums == 0, 0,Pop_Sum/Allele_df$RowSums) ## Prob(Pop[j]|Allele[i])
			conditional_prob_matrices[j,i,] <- (Pop_WeightedMean * Allele_df$Prob_Allele) / (Prob_Pops[j]) ## Prob(Allele[i]|Pop1)
		}
	}

	## Select only top SNPs based on Allele freq differences between Pops $
	conditional_prob_differences <- abs(conditional_prob_matrices[1,1,] - conditional_prob_matrices[2,1,])
	diff_quantile_value <- quantile(conditional_prob_differences, ifelse(diff_quantile <= 0.5,1-diff_quantile,diff_quantile))
	top_conditional_prob_diff <- subset(conditional_prob_differences, abs(conditional_prob_differences) >= diff_quantile_value)
	Chr_Strands <- subset(Chr_Strands,Chr_Strands$ID %in% names(top_conditional_prob_diff))
	return(Chr_Strands)
}

SplitStringLogic <- function(text, k) {
	### split text string into substrings each of length k
	sst <- strsplit(text, "")[[1]]
	if (k == 2){
		paste0(sst[c(T, F)],sst[c(F, T)])
	} else if (k == 3){
		paste0(sst[c(T, F, F)],sst[c(F, T, F)], sst[c(F, F, T)])
	} else if (k == 4){
		paste0(sst[c(T, F, F, F)],sst[c(F, T, F, F)], sst[c(F, F, T, F)], sst[c(F, F, F, T)])
	} else if (k == 5){
		paste0(sst[c(T, F, F, F, F)],sst[c(F, T, F, F, F)], sst[c(F, F, T, F, F)], sst[c(F, F, F, T, F)], sst[c(F, F, F, F, T)])
	}
}

CreateKmerHaplotypes <- function(Chr_Strands, kmer){
	print(paste0("Creating ",kmer,"-mer haplotypes"))
	if (kmer != 1){
		remainder <- (nrow(Chr_Strands) %% kmer)
		if (remainder != 0) {
			### remove extra rows
			Chr_Strands <- Chr_Strands[-((nrow(Chr_Strands)-remainder+1):nrow(Chr_Strands)),]
			# print(nrow(Chr_Strands))
			}

		### get every kth row (starting with first row)
		Chr_Strand_Substrings <- Chr_Strands[seq(1,nrow(Chr_Strands),by=kmer),]
		colnames(Chr_Strand_Substrings)[2] <- "POS_Start"
		# print(nrow(Chr_Strand_Substrings))
		# print(length(Chr_Strands[seq(kmer,nrow(Chr_Strands),by=kmer),"POS"]))
		Chr_Strand_Substrings$POS_End <- Chr_Strands[seq(kmer,nrow(Chr_Strands),by=kmer),"POS"]

		### move last column to the third column position and shift everything else down
		### TODO: why?
		### col1 col2 coln col3 col4 ... col(n-1)
		Chr_Strand_Substrings <- Chr_Strand_Substrings[,c(1,2,ncol(Chr_Strand_Substrings),3:(ncol(Chr_Strand_Substrings)-1))]
		for(i in 11:ncol(Chr_Strand_Substrings)){
			Chr_Strands_AsStrings <- apply(Chr_Strands[i-1],2,paste,collapse="")
			Chr_Strand_Substrings[,i] <- SplitStringLogic(Chr_Strands_AsStrings, kmer)
		}
		rm(Chr_Strands_AsStrings)
	} else {
		Chr_Strand_Substrings <- Chr_Strands
		colnames(Chr_Strand_Substrings)[2] <- "POS_Start"
		Chr_Strand_Substrings$POS_End <- Chr_Strand_Substrings$POS_Start
		Chr_Strand_Substrings <- Chr_Strand_Substrings[,c(1,2,ncol(Chr_Strand_Substrings),3:(ncol(Chr_Strand_Substrings)-1))]
		}
	return(Chr_Strand_Substrings)
}

CreateEmissionMatrix <- function(Chr_Strand_Substrings,ADMIX2,kmer_haplotypes){
	print("Creating Emission Matrix")
	emissionmatrices <-array(1,c(ncol(ADMIX2),length(kmer_haplotypes),nrow(Chr_Strand_Substrings)),dimnames=list(Pop=colnames(ADMIX2),Allele=kmer_haplotypes,SNP=Chr_Strand_Substrings$ID))
	for (i in 1:length(kmer_haplotypes)){
		haplotype_df <- data.frame(RowSums=rowSums(Chr_Strand_Substrings[,11:ncol(Chr_Strand_Substrings)] == kmer_haplotypes[i]),
															 Prob_Haplotype=numeric(nrow(Chr_Strand_Substrings)), row.names=Chr_Strand_Substrings$ID)
		haplotype_df$Prob_Haplotype <- haplotype_df$RowSums/(nrow(ADMIX2)) ### Prob(haplotype[i])
		for (j in 1:ncol(ADMIX2)){
			Pop_Sum <- rowSums(sweep(Chr_Strand_Substrings[,11:ncol(Chr_Strand_Substrings)] == kmer_haplotypes[i],MARGIN=2,ADMIX2[,j],'*',check.margin=FALSE)) # sum_over_individuals(#haplotype[i] per Individual * %Pop[j] per Individual)
			Pop_WeightedMean <- ifelse(Pop_Sum==0 & haplotype_df$RowSums == 0, 0,Pop_Sum/haplotype_df$RowSums) # Prob(Pop[j]|haplotype[i])
			emissionmatrices[j,i,] <- (Pop_WeightedMean * haplotype_df$Prob_Haplotype) / (Prob_Pops[j]) # Prob(REFREF|Pop1)
		}
	}
	## Convert to log2 emission matrices ##
	log2_emissionmatrices <- log2(emissionmatrices)
	return(log2_emissionmatrices)
}

viterbi_init <- function(next_haplotype,init_prob,emissionmatrix, kmer_haplotypes){
	haplotype_index <- which(kmer_haplotypes == next_haplotype)
	next_prob <- init_prob
	for (pop_index in 1:length(next_prob)){
		next_prob[pop_index] <- emissionmatrix[pop_index,haplotype_index] + init_prob[pop_index]
	}
	return(next_prob)
}

viterbi <- function(next_haplotype,prev_prob,emissionmatrix,transitionmatrix, kmer_haplotypes){
	haplotype_index <- which(kmer_haplotypes == next_haplotype)
	next_prob <- prev_prob
	for (pop_index in 1:length(next_prob)){
		next_prob[pop_index] <- emissionmatrix[pop_index,haplotype_index] + max(prev_prob + transitionmatrix[pop_index,])
	}  
	return(next_prob)
}

ForwardHMM <- function(Chr_Strand_Substrings_Array, init_prob, log2_emissionmatrices, log2_transitionmatrix, kmer_haplotypes){
	Chr_States_Forward <- Chr_Strand_Substrings_Array
	next_prob <- viterbi_init(as.character(Chr_Strand_Substrings_Array[1]),log2(init_prob),log2_emissionmatrices[,,1], kmer_haplotypes) #FORWARDS
	Chr_States_Forward[1] <- names(next_prob)[which.max(next_prob)] #FORWARDS
	prev_prob <- next_prob
	for (j in 2:length(Chr_Strand_Substrings_Array)){ #FORWARD
		next_prob <- viterbi(as.character(Chr_Strand_Substrings_Array[j]),prev_prob,log2_emissionmatrices[,,j],log2_transitionmatrix, kmer_haplotypes)
		Chr_States_Forward[j] <- names(next_prob)[which.max(next_prob)]
		prev_prob <- next_prob
	}
	print("Finished Forward Direction")
	return(Chr_States_Forward)
}

BackwardHMM <- function(Chr_Strand_Substrings_Array, init_prob, log2_emissionmatrices, log2_transitionmatrix, kmer_haplotypes){
	Chr_States_Backward <- Chr_Strand_Substrings_Array
	next_prob <- viterbi_init(as.character(Chr_Strand_Substrings_Array[length(Chr_Strand_Substrings_Array)]),log2(init_prob),log2_emissionmatrices[,,length(Chr_Strand_Substrings_Array)], kmer_haplotypes)
	Chr_States_Backward[length(Chr_Strand_Substrings_Array)] <- names(next_prob)[which.max(next_prob)] #FORWARDS
	prev_prob <- next_prob
	for (j in (length(Chr_Strand_Substrings_Array)-1):1){
		next_prob <- viterbi(as.character(Chr_Strand_Substrings_Array[j]),prev_prob,log2_emissionmatrices[,,j],log2_transitionmatrix, kmer_haplotypes)
		Chr_States_Backward[j] <- names(next_prob)[which.max(next_prob)]
		prev_prob <- next_prob
	}
	print("Finished Backward Direction")
	return(Chr_States_Backward)
}

ForwardHMM_Prob <- function(Chr_Strand_Substrings_Array, init_prob, log2_emissionmatrices, log2_transitionmatrix, kmer_haplotypes){
	Chr_States_Forward <- Chr_Strand_Substrings_Array
	Chr_States_Forward_Prob <- Chr_Strand_Substrings_Array
	next_prob <- viterbi_init(as.character(Chr_Strand_Substrings_Array[1]),log2(init_prob),log2_emissionmatrices[,,1], kmer_haplotypes) #FORWARDS
	Chr_States_Forward_Prob[1] <- (next_prob[1]/next_prob[2])
	Chr_States_Forward[1] <- names(next_prob)[which.max(next_prob)] #FORWARDS
	prev_prob <- next_prob
	for (j in 2:length(Chr_States_Forward)){ #FORWARD
		next_prob <- viterbi(as.character(Chr_Strand_Substrings_Array[j]),prev_prob,log2_emissionmatrices[,,j],log2_transitionmatrix, kmer_haplotypes)
		Chr_States_Forward_Prob[j] <- (next_prob[1]/next_prob[2])
		Chr_States_Forward[j] <- names(next_prob)[which.max(next_prob)]
		prev_prob <- next_prob
	}
	print("Finished Forward Direction")
	return(list(Chr_States_Forward,Chr_States_Forward_Prob))
}

BackwardHMM_Prob <- function(Chr_Strand_Substrings_Array, init_prob, log2_emissionmatrices, log2_transitionmatrix, kmer_haplotypes){
	Chr_States_Backward <- Chr_Strand_Substrings_Array
	Chr_States_Backward_Prob <- Chr_Strand_Substrings_Array
	next_prob <- viterbi_init(as.character(Chr_Strand_Substrings_Array[length(Chr_Strand_Substrings_Array)]),log2(init_prob),log2_emissionmatrices[,,length(Chr_Strand_Substrings_Array)], kmer_haplotypes)
	Chr_States_Backward_Prob[length(Chr_Strand_Substrings_Array)] <- (next_prob[1]/next_prob[2])
	Chr_States_Backward[length(Chr_Strand_Substrings_Array)] <- names(next_prob)[which.max(next_prob)] 
	prev_prob <- next_prob
	for (j in (length(Chr_Strand_Substrings_Array)-1):1){
		next_prob <- viterbi(as.character(Chr_Strand_Substrings_Array[j]),prev_prob,log2_emissionmatrices[,,j],log2_transitionmatrix, kmer_haplotypes)
		Chr_States_Backward_Prob[j] <- (next_prob[1]/next_prob[2])
		Chr_States_Backward[j] <- names(next_prob)[which.max(next_prob)]
		prev_prob <- next_prob
	}
	print("Finished Backward Direction")
	return(list(Chr_States_Backward,Chr_States_Backward_Prob))
}


ForwardLocalAncestry <- function(Chr_Strand_Substrings, ADMIX2, Individual_IDs, log2_emissionmatrices, kmer, recomb_rate = 0.001){
	### "00000", "00001", "00011", etc. (length of each string == kmer)
	kmer_haplotypes <- apply(permutations(n = 2, r = kmer, v = c("0","1"),repeats=TRUE), 1, paste, collapse="")
	Chr_States_Forward <- Chr_Strand_Substrings[,c("POS_Start","POS_End",Individual_IDs)]

	### clear out all but first two columns of Chr_States_Forward
	Chr_States_Forward[,3:ncol(Chr_States_Forward)] <- rep(NA, (ncol(Chr_States_Forward)-2)*nrow(Chr_States_Forward))

	for (id in Individual_IDs) {
		print(paste0("Running Viterbi for ", id))
		### Initial State ###
		init_prob <- ADMIX2[id,]

		### Create transition matrix for each individual ###
		transitionmatrix <- array(data=rep(0,ncol(ADMIX2)^2),dim=c(ncol(ADMIX2),ncol(ADMIX2)),dimnames = list(Previous=colnames(ADMIX2),Next=colnames(ADMIX2)))
		for (pop_index in 1:ncol(ADMIX2)){
			transitionmatrix[-(pop_index),pop_index] <- as.numeric(recomb_rate * init_prob[pop_index]) 
		}
		diag(transitionmatrix) <- (1 - rowSums(transitionmatrix))
		log2_transitionmatrix <- log2(transitionmatrix)
		### Forward Direction HMM ###
		Chr_States_Forward[,id] <- ForwardHMM(Chr_Strand_Substrings[,id], init_prob, log2_emissionmatrices, log2_transitionmatrix, kmer_haplotypes)
	}
	return(Chr_States_Forward)
}

BackwardLocalAncestry <- function(Chr_Strand_Substrings, ADMIX2, Individual_IDs, log2_emissionmatrices, kmer, recomb_rate = 0.001){
	### "00000", "00001", "00011", etc. (length of each string == kmer)
	kmer_haplotypes <- apply(permutations(n = 2, r = kmer, v = c("0","1"),repeats=TRUE), 1, paste, collapse="")
	Chr_States_Backward <- Chr_Strand_Substrings[,c("POS_Start","POS_End",Individual_IDs)]

	### clear out all but first two columns of Chr_States_Forward
	Chr_States_Backward[,3:ncol(Chr_States_Backward)] <- rep(NA, (ncol(Chr_States_Backward)-2)*nrow(Chr_States_Backward))
	for (id in Individual_IDs) {
		print(paste0("Running Viterbi for ", id))
		### Initial State ###
		init_prob <- ADMIX2[id,]

		### Create transition matrix for each individual ###
		transitionmatrix <- array(data=rep(0,ncol(ADMIX2)^2),dim=c(ncol(ADMIX2),ncol(ADMIX2)),dimnames = list(Previous=colnames(ADMIX2),Next=colnames(ADMIX2)))
		for (pop_index in 1:ncol(ADMIX2)){
			transitionmatrix[-(pop_index),pop_index] <- as.numeric(recomb_rate * init_prob[pop_index]) 
		}
		diag(transitionmatrix) <- (1 - rowSums(transitionmatrix))
		log2_transitionmatrix <- log2(transitionmatrix)
		### Backward Direction HMM ###
		Chr_States_Backward[,id] <- BackwardHMM(Chr_Strand_Substrings[,id], init_prob, log2_emissionmatrices, log2_transitionmatrix, kmer_haplotypes)
	}
	return(Chr_States_Backward)
}

BothDirectionsLocalAncestry <- function(Chr_Strand_Substrings, ADMIX2, Individual_IDs, log2_emissionmatrices, kmer, recomb_rate = 0.001){
	kmer_haplotypes <- apply(permutations(n = 2, r = kmer, v = c("0","1"),repeats=TRUE), 1, paste, collapse="")
	Chr_States_Forward <- Chr_Strand_Substrings[,c("POS_Start","POS_End",Individual_IDs)]
	Chr_States_Forward[,3:ncol(Chr_States_Forward)] <- rep(NA, (ncol(Chr_States_Forward)-2)*nrow(Chr_States_Forward))
	Chr_States_Backward <- Chr_States_Forward
	Chr_States_Both_Directions <- Chr_States_Forward
	for (id in Individual_IDs) {
		print(paste0("Running Viterbi for ", id))
		### Initial State ###
		init_prob <- ADMIX2[id,]

		### Create transition matrix for each individual ###
		transitionmatrix <- array(data=rep(0,ncol(ADMIX2)^2),dim=c(ncol(ADMIX2),ncol(ADMIX2)),dimnames = list(Previous=colnames(ADMIX2),Next=colnames(ADMIX2)))
		for (pop_index in 1:ncol(ADMIX2)){
			transitionmatrix[-(pop_index),pop_index] <- as.numeric(recomb_rate * init_prob[pop_index]) 
		}
		diag(transitionmatrix) <- (1 - rowSums(transitionmatrix))
		log2_transitionmatrix <- log2(transitionmatrix)
		### Forward Direction HMM ###
		Chr_States_Forward[,id] <- ForwardHMM(Chr_Strand_Substrings[,id], init_prob, log2_emissionmatrices, log2_transitionmatrix, kmer_haplotypes)
		### Backward Direction HMM ###
		Chr_States_Backward[,id] <- BackwardHMM(Chr_Strand_Substrings[,id], init_prob, log2_emissionmatrices, log2_transitionmatrix, kmer_haplotypes)
		### Combine Both Directions ###
		Chr_States_Both_Directions[,id] <- ifelse(Chr_States_Forward[,id] == Chr_States_Backward[,id],
																							Chr_States_Backward[,id],"Unknown")
	}
	return(Chr_States_Both_Directions)
}

BothDirectionsLocalAncestry_Prob <- function(Chr_Strand_Substrings, ADMIX2, Individual_IDs, log2_emissionmatrices, kmer, recomb_rate = 0.001){
	kmer_haplotypes <- apply(permutations(n = 2, r = kmer, v = c("0","1"),repeats=TRUE), 1, paste, collapse="")
	Chr_States_Forward <- Chr_Strand_Substrings[,c("POS_Start","POS_End",Individual_IDs)]
	Chr_States_Forward[,3:ncol(Chr_States_Forward)] <- rep(NA, (ncol(Chr_States_Forward)-2)*nrow(Chr_States_Forward))
	Chr_States_Forward_Prob <- Chr_States_Forward
	Chr_States_Backward <- Chr_States_Forward
	Chr_States_Backward_Prob <- Chr_States_Forward
	Chr_States_Both_Directions <- Chr_States_Forward
	for (id in Individual_IDs) {
		print(paste0("Running Viterbi for ", id))
		### Initial State ###
		init_prob <- ADMIX2[id,]

		### Create transition matrix for each individual ###
		transitionmatrix <- array(data=rep(0,ncol(ADMIX2)^2),dim=c(ncol(ADMIX2),ncol(ADMIX2)),dimnames = list(Previous=colnames(ADMIX2),Next=colnames(ADMIX2)))
		for (pop_index in 1:ncol(ADMIX2)){
			transitionmatrix[-(pop_index),pop_index] <- as.numeric(recomb_rate * init_prob[pop_index]) 
		}
		diag(transitionmatrix) <- (1 - rowSums(transitionmatrix))
		log2_transitionmatrix <- log2(transitionmatrix)

		### Forward Direction HMM ###
		ForwardHMM_Results <- ForwardHMM_Prob(Chr_Strand_Substrings[,id], init_prob, log2_emissionmatrices, log2_transitionmatrix, kmer_haplotypes)
		Chr_States_Forward[,id] <- ForwardHMM_Results[[1]]
		Chr_States_Forward_Prob[,id] <- unlist(ForwardHMM_Results[[2]])
		### Backward Direction HMM ###
		BackwardHMM_Results <- BackwardHMM_Prob(Chr_Strand_Substrings[,id], init_prob, log2_emissionmatrices, log2_transitionmatrix, kmer_haplotypes)
		Chr_States_Backward[,id] <- BackwardHMM_Results[[1]]
		Chr_States_Backward_Prob[,id] <- unlist(BackwardHMM_Results[[2]])

		### Combine Both Directions ###
		Chr_States_Both_Directions[,id] <- ifelse(Chr_States_Forward[,id] == Chr_States_Backward[,id],
																							Chr_States_Backward[,id],"Unknown")
	}
	return(list(Chr_States_Both_Directions,Chr_States_Forward,Chr_States_Forward_Prob,Chr_States_Backward,Chr_States_Backward_Prob))
}


PlotLocalAncestry <- function(LocalAncestryDataFrame, kmer, imagename=NULL){
	LocalAncestryDataFrame_long <- melt(LocalAncestryDataFrame, id.vars=c("POS_Start","POS_End"),value.name="EstimatedAncestry")
	imagename <- ifelse(is.null(imagename), paste0("InformativeSNPs",kmer,"merHaplotype_Haploid.png"),imagename)
	plots<-ggplot(LocalAncestryDataFrame_long,aes(x=POS_Start,y=variable,color=EstimatedAncestry))+
		geom_point(size=0.5)+
		xlab("Chr Position")+
		ylab("")+
		scale_color_manual(values=c("#FF0000", "#0000FF", "#73008C"),                         
											 name="Estimated Ancestry",
											 breaks=c("Pop1", "Pop2", "Unknown"),
											 labels=c("Ancestry 1", "Ancestry 2", "Ambiguous"))+
		ggtitle(paste0("Informative SNPs, ", kmer, "-mer haplotypes"))+
		theme_bw(base_size = 5)+
		theme(legend.key.size = unit(0.1, "in"))
	ggsave(plots,filename=imagename,height = 1.5, width= 5, units="in")
}