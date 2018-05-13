#####








#### Don't look at this, look at Test_Local_Ancestry_Inference_ASW_from_3Pops.R in RScripts folder







####

source("../bin/Chr_kmer_InformativeSNPs_Haploid_HMM.R")

args <- commandArgs(trailingOnly=TRUE)
homologous_filename <- args[1]
reference_filename <- args[2]
all_admixture_filename <- args[3]
admixture_filename <- args[4]
out_filename <- args[5]
kmer <- as.numeric(args[6])

#### Read in Strands + Select only the ASW samples ####
print("Loading Data")
ASW_Chr_Strands <- read.table(homologous_filename, head=T, sep="\t")
ASW_Chr_Strands <- subset(ASW_Chr_Strands,!(duplicated(ASW_Chr_Strands$ID)|duplicated(ASW_Chr_Strands$ID,fromLast = TRUE))) ### why are all duplicated IDs completely removed?
Ref_Chr_Strands <- read.table(reference_filename, head=T, sep="\t")
Ref_Chr_Strands <- subset(Ref_Chr_Strands,!(duplicated(Ref_Chr_Strands$ID)|duplicated(Ref_Chr_Strands$ID,fromLast = TRUE)))
Merged_Chr_Strands <- merge(ASW_Chr_Strands,Ref_Chr_Strands,by=c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"))

###### Function to Create Test Chromosomes ######
Create_Test_Chromosomes <- function(Chr_Strands,Pops,Ref_IDs,Num_Switch){
	Test_Pos_Chromosome <- sample(Chr_Strands$POS, Num_Switch, replace=FALSE)
	Test_Pos_Chromosome <- sort(Test_Pos_Chromosome)
	Test_Pos_Chr_index_stop <- sapply(1:Num_Switch, function(x){which(Chr_Strands$POS == Test_Pos_Chromosome[x])})
	Test_Pos_Chr_index_start <- c(0,Test_Pos_Chr_index_stop)
	Test_Pos_Chr_index_stop <- c(Test_Pos_Chr_index_stop,nrow(Chr_Strands))
	Test_Chr <- numeric(nrow(Chr_Strands))
	True_Chr <- character(nrow(Chr_Strands))
	for (j in 1:length(Pops)){
		for (i in seq(j,length(Test_Pos_Chr_index_start),length(Pops))){
			Test_Chr[(Test_Pos_Chr_index_start[i]+1):Test_Pos_Chr_index_stop[i]] <- Chr_Strands[(Test_Pos_Chr_index_start[i]+1):Test_Pos_Chr_index_stop[i],Ref_IDs[j]]
			True_Chr[(Test_Pos_Chr_index_start[i]+1):Test_Pos_Chr_index_stop[i]] <- rep(Pops[j],length((Test_Pos_Chr_index_start[i]+1):Test_Pos_Chr_index_stop[i]))
		}
	}
	return(list(Test_Chr,True_Chr))
}
#### Test Strand 1 ####
print("Creating Test Strands")
Known_Sample1 <- Create_Test_Chromosomes(Merged_Chr_Strands,Pops = c("Pop2","Pop1"),Ref_IDs = c("NA06984_1","NA18486_1"),Num_Switch = 4)
Test_1 <- Known_Sample1[[1]]
True_1 <- Known_Sample1[[2]]

#### Test Strand 2 ####
Known_Sample2 <- Create_Test_Chromosomes(Merged_Chr_Strands,Pops = c("Pop2","Pop1"),Ref_IDs = c("NA06984_2","NA18486_2"),Num_Switch = 6)
Test_2 <- Known_Sample2[[1]]
True_2 <- Known_Sample2[[2]]

True_LA <- data.frame(POS_Start=ASW_Chr_Strands$POS,
											POS_End=ASW_Chr_Strands$POS,True_1,True_2)

#### Read in Admixture results for all chromosomes ####
All_ADMIX2<-read.table(all_admixture_filename,head=F)
num_pops <- ncol(All_ADMIX2)
All_ADMIX2_ordered <- All_ADMIX2[,order(colMeans(All_ADMIX2),decreasing=TRUE)]
colnames(All_ADMIX2_ordered)<-paste0("Pop",1:ncol(All_ADMIX2_ordered))

#### Read in Admixture (k=2) results for single chromosome ####
print("Loading Admixture")
ADMIX2_unordered <-read.table(admixture_filename,head=F)
Admix_chr_diff_all <- apply(ADMIX2_unordered, 2, function(x) {sum(All_ADMIX2_ordered$Pop1-x)})
ADMIX2 <- ADMIX2_unordered[,order(Admix_chr_diff_all)]
colnames(ADMIX2)<-paste0("Pop",1:ncol(ADMIX2))
rownames(ADMIX2)<-colnames(ASW_Chr_Strands)[10:ncol(ASW_Chr_Strands)]
Test_ADMIX <- rbind((table(True_LA$True_1)/nrow(True_LA)),(table(True_LA$True_2)/nrow(True_LA)))
Test_IDs <- c("Test_1","Test_2")
row.names(Test_ADMIX) <- Test_IDs
ADMIX2_with_Test <- rbind(ADMIX2,Test_ADMIX)
Prob_Pops <- (colSums(ADMIX2_with_Test)/nrow(ADMIX2_with_Test)) #Prob(Pop1) & Prob(Pop2), respectively

#### Update ASW_Chr_Strands with Test Individuals ####
ASW_Chr_Strands_with_Test <- data.frame(ASW_Chr_Strands,Test_1 = Test_1, Test_2 = Test_2)

#### Select the most informative SNPs ####
Chr_Strands_with_Test_top <- SelectInformativeSNPs(ASW_Chr_Strands_with_Test, Prob_Pops, diff_quantile = 0.1)

#### Create new data frame with k consecutive snps together (k-mers) ####
Chr_Strand_with_Test_Substrings <- CreateKmerHaplotypes(Chr_Strands_with_Test_top ,kmer)

#### Create Emission Matrices for each k-mer ####
kmer_haplotypes <- apply(permutations(n = 2, r = kmer, v = c("0","1"),repeats=TRUE), 1, paste, collapse="")
log2_emissionmatrices <- CreateEmissionMatrix(Chr_Strand_with_Test_Substrings,ADMIX2_with_Test,kmer_haplotypes)

#### Local Ancestry Inference ####
Train_with_Test_LA <- BothDirectionsLocalAncestry_Prob(Chr_Strand_with_Test_Substrings, ADMIX2_with_Test, Test_IDs, log2_emissionmatrices, kmer, recomb_rate = 0.002)
Train_with_Test_LA_Both <- Train_with_Test_LA[[1]]
True_Full <- data.frame(POS=ASW_Chr_Strands$POS, True_1, True_2)
True_Substrings <- merge(Chr_Strand_with_Test_Substrings[,1:3], True_Full, by.x="POS_Start", by.y="POS")
Test_LA_Reordered <- data.frame(Train_with_Test_LA_Both[,1:3], True_1=True_Substrings$True_1, Test_2=Train_with_Test_LA_Both[,4], True_2=True_Substrings$True_2)
write.table(Test_LA_Reordered, out_filename, quote=F, sep="\t", row.names = F)
out_png_filename <- gsub("txt", "png", out_filename)
PlotLocalAncestry(LocalAncestryDataFrame = Test_LA_Reordered, kmer, imagename=out_png_filename)
