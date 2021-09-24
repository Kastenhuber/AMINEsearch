#AMINEsearch
#Ted Kastenhuber, last updated 2021-09-17

#setup packages
#one-time only. installs the human genome (~850 MB download).
if (FALSE){source("http://www.bioconductor.org/biocLite.R")
biocLite("BSgenome")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
biocLite("Biostrings")
biocLite("biomaRt")
biocLite("hom.Hs.inp.db")
biocLite("org.Mm.eg.db")
biocLite("AnnotationDbi")
install.packages('readr')
install.packages('tidyverse')
install.packages('foreach')
install.packages('doParallel')
}

#load packages
library(BSgenome.Hsapiens.UCSC.hg19) 
library(Biostrings)
library(biomaRt)
library(hom.Hs.inp.db) 
library(org.Mm.eg.db) 
library(AnnotationDbi)
library(readr)
library(Biostrings)
library(tidyverse)
library(foreach)
library(doParallel)


###############################################################################################################
#input files - to be specified by end user. local or remote
setwd(folder) # user defined PATH

#load .maf file for mutations - can use any standard .maf, double check heading names are consistent for new inputs. 
#make sure that the following columns are present with these exact names: 
# c("Hugo_Symbol","Chromosome","Start_Position","End_Position","HGVSp_Short","Reference_Allele","Tumor_Seq_Allele2","Variant_Classification")	
#check variant classifications listed in format c("Nonsense_Mutation", "Missense_Mutation", "Splice_Site", "Nonstop_Mutation")
#maffile="20200616-IMPACT-data_mutations_extended_patched.maf.txt"
#read.delim(maffile,stringsAsFactors=FALSE)->maf

#mutation data source
maffile="20200616-IMPACT-data_mutations_extended_patched.maf.txt"   #name or describe mutation data source
maf=read.table(file="./input/20200616-IMPACT-data_mutations_extended_patched.maf.txt",stringsAsFactors=FALSE,sep="\t",header=TRUE)

#alterative data source:
#maf = read.table(url("https://wcm.box.com/shared/static/uxoyunnl2nzg0cd68lmxjr93fjdotwnb.txt"),stringsAsFactors=FALSE,sep="\t",header=TRUE)

#define Cas9 and editor space
cas9s=read.table(file="./input/cas9_properties_20210828.txt",stringsAsFactors=FALSE,sep="\t",header=TRUE)
edit=read.table(file="./input/editors_20210828.txt",stringsAsFactors=FALSE,sep="\t",header=TRUE)

#alterative data source:
#cas9s = read.table(url("https://wcm.box.com/shared/static/dcocqm6xwfat74tbrv18lnxztmrhoyoi.txt"),stringsAsFactors=FALSE,sep="\t",header=TRUE)
#edit = read.table(url("https://wcm.box.com/shared/static/hrplnq9l089jjzp7urr97orebdajwawj.txt"),stringsAsFactors=FALSE,sep="\t",header=TRUE)


#settings - to be specified by end user
minObservations=6
registerDoParallel(1) #how many cores to use for parallelization. To check how many cores available use 'detectCores()'



#worklist
humansearch=TRUE
mousesearch=TRUE
correction_lib=FALSE 


#load inputs specific to mousesearch
##load("REFSEQ_IMPACT.Rdata") #replace with full refseq for general case 
#RefSeq <- readRDS("REFSEQ_FULL_Hs2Mm.rds") #human to mouse orthology for mousesearch
if (mousesearch==TRUE){
	#load 4 data objects for mouse transcript selection and annotation: RefSeq, mf, mfanno, select_mouse_transcript_annotation
	RefSeq <- readRDS("./reference_files/REFSEQ_FULL_Hs2Mm.rds") #human to mouse orthology for mousesearch
	load("./reference_files/select_mouse_transcript_annotation.Rdata") #select_mouse_transcript_annotation
	load("./reference_files/mfanno.gencode.vM25.pc_transcripts.annotation_select.Rdata") #mfanno
	mf=readDNAStringSet("./reference_files/gencode.vM25.pc_transcripts_select.fa.gz")

	#alterative data source:
	#RefSeq <- readRDS(url("https://wcm.box.com/shared/static/ep5jn1c4l4qsuv46ehbdp14wxud35umd.rds")) #human to mouse orthology for mousesearch
	#temp <- tempfile(fileext = ".Rdata")
	#download.file("https://wcm.box.com/shared/static/q0kdo2qa3s7l94c6mg5p0du4s3k8vyjx.rdata", temp)
	#load(temp)
}



###############################################################################################################
#maf to candidates file
#make sure that the following columns are present with these exact names: 
#Hugo_Symbol	 Chromosome Start_Position End_Position	HGVSp_Short Reference_Allele Variant_Allele Variant_Classification	
if("Tumor_Seq_Allele2" %in% colnames(maf)){colnames(maf)[which(colnames(maf)=="Tumor_Seq_Allele2")]="Variant_Allele"}
req.cols=c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "HGVSp_Short", "Reference_Allele", "Variant_Allele", "Variant_Classification")
if ( all(req.cols %in% colnames(maf))==F ){stop(paste("\n","Error. Please ensure the following columns in maf are present:", req.cols[which(req.cols %in% colnames(maf)==F)]) )}

#check variant classifications listed in format c("Nonsense_Mutation", "Missense_Mutation", "Splice_Site", "Nonstop_Mutation")
req.var.class=c("Nonsense_Mutation", "Missense_Mutation", "Splice_Site", "Nonstop_Mutation")
error.var.class=c("Nonsense Mutation", "Missense Mutation", "Splice Site", "Nonstop Mutation","Nonsense", "Missense", "Splice", "Nonstop")
if ( any(req.var.class %in% maf$Variant_Classification)==F ){stop(paste("\n","Error. Please ensure the Variant_Classification in maf follows:", paste(req.var.class, collapse = ", ")) )}
if ( any(error.var.class %in% maf$Variant_Classification) ){stop(paste("\n","Error. Please ensure the Variant_Classification in maf follows:", paste(req.var.class, collapse = ", ")) )}


maf$ID=paste(maf$Hugo_Symbol,maf$HGVSp_Short,maf$Chromosome,maf$Start_Position,maf$change,sep=" ")
maf$ID_AA=paste(maf$Hugo_Symbol,maf$HGVSp_Short,sep=" ")
maf$recurrentAAsub=duplicated(maf$ID_AA)
maf$change=paste(maf$Reference_Allele,maf$Variant_Allele,sep="")
if (minObservations>1){
candidates0 =maf[which(maf$Variant_Classification %in% c("Nonsense_Mutation", "Missense_Mutation", "Splice_Site", "Nonstop_Mutation") & maf$ID_AA %in% maf$ID_AA[which(maf$recurrentAAsub==TRUE)] & maf$Reference_Allele %in% c("A","C","G","T") & maf$Variant_Allele %in% c("A","C","G","T")),]
}
if (minObservations==1){candidates0 =maf[which(maf$Variant_Classification %in% c("Nonsense_Mutation", "Missense_Mutation", "Splice_Site", "Nonstop_Mutation") &  maf$Reference_Allele %in% c("A","C","G","T") & maf$Variant_Allele %in% c("A","C","G","T")),]}



#condense to one entry per mutant at nt level (already selected for recurrent mutations at AA level)
candidates0$duplicated_nt_change=duplicated(candidates0$ID)   
candidates=candidates0[which(candidates0$duplicated_nt_change==FALSE),]	
for (c in 1:nrow(candidates)){
	candidates$numberMutationsObserved[c]=length(which(maf$ID==candidates$ID[c]))
	candidates$numberSubstitutionsObserved[c]=length(which(maf$ID_AA==candidates$ID_AA[c]))
	}
keepcols=c("Hugo_Symbol","Chromosome","Start_Position","End_Position","Strand","Reference_Allele","Variant_Allele","Variant_Classification","HGVSp_Short","ID","ID_AA","change", "numberMutationsObserved","numberSubstitutionsObserved","flank")
candidates=candidates[which(candidates$numberSubstitutionsObserved>=minObservations),which(colnames(candidates) %in% keepcols)]


# Function to pull flanking sequence. Defaults to +/- 30 bp
getflank <- function(chr, position,allele, offset=30) {
leftflank <- getSeq(Hsapiens,paste("chr",chr,sep=""),position-offset,position-1)
rightflank <- getSeq(Hsapiens,paste("chr",chr,sep=""),position+1,position+offset)
paste(leftflank,allele,rightflank,sep="")
}

#get flanking sequence for each candidate mutation
candidates$flank=NA
for (m in 1:nrow(candidates)){
	candidates$flank[m]=getflank(candidates$Chromosome[m],candidates$Start_Position[m],candidates$Reference_Allele[m])
	}


if (correction_lib==TRUE){
	WT=candidates$Reference_Allele
	MUT=candidates$Variant_Allele
	candidates$Reference_Allele=MUT
	candidates$Variant_Allele=WT
	candidates$change=paste0(MUT,WT)	
	substr(candidates$flank,31,31)<-MUT
	substr(candidates$mouse_flank,31,31)<-MUT
	candidates$correction_lib=TRUE
}


if ("maffile" %in% ls()){
saveRDS(candidates,file=paste0(substr(Sys.time(),0,10),"_from",unlist(strsplit(maffile,split="[.]"))[1],".candidates.rds"))} else saveRDS(candidates,file=paste0(substr(Sys.time(),0,10),".candidates.rds"))



#candidates=readRDS("2021-08-28_from20200616-IMPACT-data_mutations_extended_patched.candidates.rds")

###############################################################################################################
#human guide search

#load candidates if not already loaded
#load(file="2021-08-28_20200616-IMPACT-data_mutations_extended_patched.candidates.Rdata")

if (humansearch==TRUE){
	
	
	
ptm <- proc.time()
#make parallel analysis of individual mutations
run.idx <- 1:nrow(candidates)
output <- foreach(m = run.idx) %dopar% {
	
#procedure for each mutation	
targets=candidates[0,]
t=0	
#for (m in 1:nrow(candidates)){	# replaced by %dopar%
	for (c in 1:nrow(cas9s)){
		#define rules for this specific cas9
		cas=cas9s$Cas_ortholog[c] #which cas9 ortholog?
		PAM=unlist(strsplit(cas9s$PAM[c],split=",")) #cas9-specific PAM required
		GL=cas9s$gRNA_Length[c] #length of the guide
			for (e in 1:nrow(edit)){
				#define rules for this specific editor
				editor=edit$Editor[e] #pick editor
				PDmin=GL-edit$window_end[e]+1 #editor-specific window	
				if(edit$window_end[e]<0){PDmin =GL-edit$window_end[e]}	
				PDmax=GL-edit$window_start[e]+1
				if(edit$window_start[e]<0){PDmax=GL-edit$window_start[e]} #fixed for -1 editing
				BE=unlist(strsplit(edit$BE[e],split="[|]")) #which bases can be edited by this enzyme?			
				#does this mutation match this enzyme conversion?
				if (candidates$change[m] %in% BE){
				forward=BE[1]
				reverse=BE[2]
				#is the desired editing on the reverse strand? If so look at the reverse complement

				tmpflank=NA
				if (candidates$change[m]==forward){tmpflank=candidates$flank[m]} 
				if (candidates$change[m]==reverse){tmpflank=as.character(reverseComplement(DNAString(candidates$flank[m])))} 
					for (s in PDmin:PDmax){
					#scan through the possible guide spacing and check for a PAM
						for (p in 1:length(PAM)){ 
						PL=nchar(PAM[p]) #length of the PAM	
						if(countPattern(DNAString(PAM[p]),DNAString(substr(tmpflank,31+s,31+s+PL-1)),fixed=FALSE)>0){ #fixed=FALSE to include wildcards
						t=t+1
						targets[t,1:ncol(candidates)]=candidates[m,]
						targets$flank[t]=tmpflank
						targets$sgRNA[t]=substr(targets$flank[t],31+s-GL,30+s) 
						targets$sgRNA_G[t]=paste0("G",substr(targets$sgRNA[t],2,nchar(targets$sgRNA[t]))) #substitute first nt of guide and guide within sensor with 'G'					
						targets$sensorv4[t]=paste0("CATAGCGTACACGTCTCACACC", targets$sgRNA_G[t], "GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTT", substr(targets$flank[t],31+s-GL-10,30+s+10), "GAATTCTAGATCCGGTCGTCAACGGCAC")
						targets$window[t]=substr(targets$flank[t],30+s-GL+edit$window_start[e],30+s-GL+edit$window_end[e]) #version for edits in/after guide
						if (edit$window_start[e]<0){targets$window[t]=substr(targets$flank[t],31+s-GL+edit$window_start[e],30+s-GL+edit$window_end[e])} #fixed if edited outside of guide range						
						targets$PAM[t]=as.character(PAM[p])
						targets$cas9[t]= cas
						targets$editor[t]= editor
						targets$adjacentEdit[t]=countPattern(substr(BE[1],0,1),DNAString(targets$window[t]))-1 
						targets$targetBase_position[t]=GL-s+1 			
						if (targets $targetBase_position[t]<=0){targets $targetBase_position[t]= targets $targetBase_position[t]-1} #fixed if edited outside of guide range (5')
						targets$sgtop[t]=paste0("CACC",as.character(targets$sgRNA_G[t]))
						targets$sgbot[t]=paste0("AAAC", as.character(reverseComplement(DNAString(targets$sgRNA_G[t]))))
						}#end for if statement
						}#end for pam
					}#end of sliding window
				} #end of base change check 
			}#end for editor
	}#end for cas9 ortholog
	targets
}#end for candidate mutation site
#row.names(targets)=seq(from=1,to=nrow(targets))




### collect parallel output
names(output) <- candidates$ID[run.idx]
final.res.hu <- output %>%
  tibble(Name = names(.), content = .) %>%
  mutate(nrow = map_dbl(content, nrow)) %>%
  filter(nrow > 0) %>%
  dplyr::select(-nrow) %>%
  unnest(content) %>%
  as.data.frame

proc.time() - ptm


keepcols=c("HGVSp_Short","Variant_Classification", "Chromosome", "Start_Position", "End_Position","Hugo_Symbol", "ID", "change", "numberMutationsObserved", "flank", "sgRNA","sgRNA_G", "window", "PAM", "cas9","editor","adjacentEdit","targetBase_position","sensorv4","correction_lib","sgtop","sgbot")
final.res.hu = final.res.hu[,which(colnames(final.res.hu) %in% keepcols)]





#write.table(targets,file=paste0("IMPACT_",substr(Sys.time(),0,10),"AMINEsearch_human.txt"),sep="\t",row.names=FALSE) #output to txt option

if ("maffile" %in% ls()){
saveRDS(final.res.hu,file=paste0(substr(Sys.time(),0,10),"_from",unlist(strsplit(maffile,split="[.]"))[1],".AMINEsearch_human.targets.rds"))
} else saveRDS(final.res.hu,file=paste0(substr(Sys.time(),0,10),".AMINEsearch_human.targets.rds"))


}#end humansearch


###########################################################################################################
#mouse guide search
if (mousesearch==TRUE){
	


#pre-processing mouse orthology, alignments, and AA consequences
#load candidates if not already loaded

candidates=candidates[order(candidates$Hugo_Symbol),]

#refer to mouse ensembl to get mouse orthologs
#load Refseq (mouse>human orthologs)

candidates$human.strand= RefSeq$strand[match(candidates$Hugo_Symbol, RefSeq$external_gene_name)]
candidates$mouse.ensg=candidates$mouse.symbol= NA
candidates$mouse.symbol= RefSeq$mmusculus_homolog_associated_gene_name[match(candidates$Hugo_Symbol, RefSeq$external_gene_name)]
candidates$mouse.ensg = RefSeq$mmusculus_homolog_ensembl_gene[match(candidates$Hugo_Symbol, RefSeq$external_gene_name)]
candidates$mouse.chromosome=RefSeq$mmusculus_homolog_chromosome[match(candidates$Hugo_Symbol, RefSeq$external_gene_name)]
candidates$mouse.start=RefSeq$mmusculus_homolog_chrom_start[match(candidates$Hugo_Symbol, RefSeq$external_gene_name)]
candidates$mouse.end=RefSeq$mmusculus_homolog_chrom_end[match(candidates$Hugo_Symbol, RefSeq$external_gene_name)]
candidates$mouse.strand= RefSeq$mmusculus_strand[match(candidates$Hugo_Symbol, RefSeq$external_gene_name)]
candidates$mouse.percent.ID=RefSeq$mmusculus_homolog_perc_id[match(candidates$Hugo_Symbol, RefSeq$external_gene_name)]
mouse.ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host = "www.ensembl.org") #for mouse

# define a function to fetch mouse sequence for biomaRT
find_gene_sequence <- function(query_gene){
require(biomaRt)
if(length(query_gene) != 1){stop("Query a single gene")}
genes=NA
#genes = getBM(attributes=c("gene_exon_intron", "ensembl_gene_id"),filters="ensembl_gene_id",values=query_gene, mart=mouse.ensembl)
try({genes = getBM(attributes=c("gene_exon_intron", "ensembl_gene_id"),filters="ensembl_gene_id",values=query_gene, mart=mouse.ensembl)},silent=TRUE)
if(nrow(genes) != 0){return(genes)}
warning(query_gene, " does not match any mouse MGI gene symbol or ensembl gene ID")
return(NULL)
}#end function to fetch mouse seq
	



#cycle through 'candidates' and include mouse guides that target mutations uneditable in human genome
#first retreive analogous flanking sequence from mouse genome
candidates$mouse_flank=NA
ai=which(is.na(candidates$human.strand)==F & is.na(candidates$mouse.symbol)==F & is.na(candidates$mouse.ensg)==F & is.na(candidates$mouse.strand)==F)
for (m in ai){
	print(paste("Analyzing mutation",m,"of",nrow(candidates)))
	if (m==1){
	gene=candidates$Hugo_Symbol[m]	
	mouseGeneSeq=DNAString(as.character(find_gene_sequence(candidates$mouse.ensg[m])[1]))
	humanTranscriptStrand=candidates$human.strand[m]
	mouseTranscriptStrand=candidates$mouse.strand[m]
print("first gene")
	}else 	{if (candidates$Hugo_Symbol[m]!=candidates$Hugo_Symbol[m-1]){
	gene=candidates$Hugo_Symbol[m]	
	mouseGeneSeq=DNAString(as.character(find_gene_sequence(candidates$mouse.ensg[m])[1]))
	humanTranscriptStrand=candidates$human.strand[m]
	mouseTranscriptStrand=candidates$mouse.strand[m]
	print("new gene")
} }
	#if the human gen is on -1 strand, we must RC to align with mouse reference transcript, which is always delivered in sense direction
				query=DNAString(candidates$flank[m]) #look for region of interest in mouse genome
				if (humanTranscriptStrand== -1){mouseref =reverseComplement(mouseGeneSeq)}else mouseref= mouseGeneSeq
				a=b=c()
				a=pairwiseAlignment(mouseref,query, type = "local-global", gapExtension = -8, gapOpening = -150)
				b=pairwiseAlignment(reverseComplement(mouseref),query, type = "local-global", gapExtension = -8, gapOpening = -150)
			candidates$mouse_flank[m]=as.character(a)
			#mark if target base is intact in mouse and alignment is good 
			candidates$mouse_alignmentScore[m]=score(a)
			candidates$mouse_alignmentScoreRC[m]=score(b)
			candidates$mouse_alignmentScorePass[m]= (score(a)>score(b) & score(a) > 0)
			candidates$mouse_targetBaseIntact[m]=substr(as.character(a),31,31)==substr(as.character(query),31,31)
}#end retreive mouse flank
candidates$mouse_targetOK= (candidates$mouse_alignmentScorePass & candidates$mouse_targetBaseIntact)



#saveRDS(candidates,file=paste0(substr(Sys.time(),0,10),"_from",unlist(strsplit(maffile,split="[.]"))[1],".mousealigned.AA.candidates.rds"))



#AAconsequence
#specify AA consequences in human - derived from input HGVSp_Short call
candidates$human_ref_AA=substr(candidates$HGVSp_Short,3,3)
candidates$human_ref_AA[which(candidates$human_ref_AA %in%   c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")==F)]=NA

candidates$human_var_AA=NA
mni=which(candidates$Variant_Classification  %in% c("Missense_Mutation", "Nonsense_Mutation"))
candidates$human_var_AA[mni]=substr(candidates$HGVSp_Short[mni],nchar(candidates$HGVSp_Short[mni]),nchar(candidates$HGVSp_Short[mni]))
candidates$human_var_AA[which(candidates$human_var_AA %in%   c("*","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")==F)]=NA
candidates$human_pos_AA=parse_number(gsub("\\.","",candidates$HGVSp_Short)) #requires readr package


#function to pull top priority CDS. pre-processed in another script based on Gencode https://www.gencodegenes.org/
find_cds_sequence <- function(query_gene){
require("Biostrings")
if(length(query_gene) != 1){stop("Query a single gene")}
select_transcript=select_mouse_transcript_annotation$transcript_id[which(substr(select_mouse_transcript_annotation$gene_id,1,18)== query_gene)]
start=mfanno$CDSstart[which(mfanno$transcript_id == select_transcript)]
stop=mfanno$CDSend[which(mfanno$transcript_id == select_transcript)]
cds=c()
cds = subseq(mf[which(mfanno$transcript_id == select_transcript)],start,stop)
return(cds)
}#end function to fetch mouse seq
	





#main loop to locate position in mouse cDNA and determine edited AA consequence in mouse
candidates$mouse_HGVSp_Short=candidates$mouse_var_AA=candidates$mouse_ref_AA=candidates$mouse_target_codon_var=candidates$mouse_target_codon_ref=candidates$mouse_target_in_CDS=candidates$mouse_target_in_codon=candidates$mouse_pos_AA=candidates$mouse_transcript=NA
#candidates$msFlank_from_cdna=candidates$b=candidates$a=NA
target_in_query=31
ti=which(candidates$mouse_targetOK ==T & is.na(candidates$mouse.ensg)==F & candidates$Variant_Classification  %in% c("Missense_Mutation","Nonsense_Mutation"))

for (t in ti){
	#if advancing to a new gene, find cds
	if(t>1){if (candidates$mouse.ensg[t]!=candidates$mouse.ensg[max(ti[which(ti<t)])]){mouse_cds= find_cds_sequence(candidates$mouse.ensg[t])}  } else mouse_cds= find_cds_sequence(candidates$mouse.ensg[1])
	if (length(mouse_cds)==1){
	sub=a=b=c()	
	query=DNAString(candidates$mouse_flank[t])
	candidates$mouse_transcript[t]=substring(names(mouse_cds),1,18)
	
	#align to mouse cds and its reverse complement - find position in cds
	a=pairwiseAlignment(mouse_cds,query, type = "local-global", gapExtension = -8, gapOpening = -50, substitutionMatrix =nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = TRUE))
	#candidates$a[t]=score(a)
	b=pairwiseAlignment(mouse_cds,reverseComplement(query), type = "local-global", gapExtension = -8, gapOpening = -50, substitutionMatrix =nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = TRUE))
	#candidates$b[t]=score(b)
	if (score(a)>0|score(b)>0){
		
		#use best alignment to define nt substitution (with sense of cds) and location in cds
	if (score(a)>score(b) & (substr(as.character(pattern(a)), target_in_query, target_in_query)==substr(as.character(subject(a)), target_in_query, target_in_query)) ){
		
		#candidates$msFlank_from_cdna[t]="top"
		if( substr(candidates$change[t],1,1) == substr(as.character(pattern(a)), target_in_query, target_in_query) ){sub=substr(candidates$change[t],2,2)}
		if( substr(candidates$change[t],1,1) == as.character(reverseComplement(DNAString(substr(pattern(a), target_in_query, target_in_query)))) ){sub=  as.character(reverseComplement(DNAString( substr(candidates$change[t],2,2) )))  }		
		candidates$mouse_target_in_CDS[t]=start(pattern(a))+ target_in_query-1
	}	

	if (score(b)>score(a) & (substr(as.character(pattern(b)), target_in_query, target_in_query)==substr(as.character(subject(b)), target_in_query, target_in_query))){
		
		#candidates$msFlank_from_cdna[t]="bot"
		if( substr(candidates$change[t],1,1) == substr(as.character(pattern(b)), target_in_query, target_in_query) ){sub=substr(candidates$change[t],2,2)}
		if( substr(candidates$change[t],1,1) == as.character(reverseComplement(DNAString(substr(pattern(b), target_in_query, target_in_query)))) ){sub=  as.character(reverseComplement(DNAString( substr(candidates$change[t],2,2) )))  }		
		candidates$mouse_target_in_CDS[t]=start(pattern(b))+ target_in_query-1
	}		
	
	#skip if there was a problem identifying target_in_CDS  and substitution
	if (length(sub)==1& is.na(candidates$mouse_target_in_CDS[t])==F) {
	candidates$mouse_pos_AA[t]=ceiling(candidates$mouse_target_in_CDS[t]/3)
	candidates$mouse_target_in_codon[t]=candidates$mouse_target_in_CDS[t]%%3   
	if (candidates$mouse_target_in_codon[t]==0){
		candidates$mouse_target_in_codon[t]=3
		candidates$mouse_target_codon_ref[t]= substr(mouse_cds,candidates$mouse_target_in_CDS[t]-2,candidates$mouse_target_in_CDS[t])
		candidates$mouse_target_codon_var[t]=candidates$mouse_target_codon_ref[t]
		substring(candidates$mouse_target_codon_var[t],first=3) <- sub
		}
	if (candidates$mouse_target_in_codon[t]==1){
		candidates$mouse_target_codon_ref[t]= substr(mouse_cds,candidates$mouse_target_in_CDS[t],candidates$mouse_target_in_CDS[t]+2)
		candidates$mouse_target_codon_var[t]=candidates$mouse_target_codon_ref[t]
		substring(candidates$mouse_target_codon_var[t],first=1) <- sub
		}
	if (candidates$mouse_target_in_codon[t]==2){
		candidates$mouse_target_codon_ref[t]= substr(mouse_cds,candidates$mouse_target_in_CDS[t]-1,candidates$mouse_target_in_CDS[t]+1)
		candidates$mouse_target_codon_var[t]=candidates$mouse_target_codon_ref[t]
		substring(candidates$mouse_target_codon_var[t],first=2) <- sub
		}	
	candidates$mouse_ref_AA[t]= GENETIC_CODE[[  candidates$mouse_target_codon_ref[t]  ]]
	candidates$mouse_var_AA[t]= GENETIC_CODE[[  candidates$mouse_target_codon_var[t]  ]]
	#candidates$mouse_HGVSp_Short[t]=paste0("p.",candidates$mouse_ref_AA[t],candidates$mouse_pos_AA[t],candidates$mouse_var_AA[t])	
	
	
	}#end target_in_CDS check 
	}#end alignment check
	}#end cds check
	}#end AA annotate loop	

candidates$mouse_HGVSp_Short=NA
calls=which(is.na(candidates$mouse_ref_AA)==F&is.na(candidates$mouse_pos_AA)==F&is.na(candidates$mouse_var_AA)==F)
candidates$mouse_HGVSp_Short[calls]=paste0("p.",candidates$mouse_ref_AA[calls],candidates$mouse_pos_AA[calls],candidates$mouse_var_AA[calls])

candidates$Mm_orthology=NA
candidates$Mm_orthology[which(candidates$mouse_alignmentScorePass==FALSE | is.na(candidates$mouse_HGVSp_Short) )]="Low homology in target region"
candidates$Mm_orthology[which(candidates$mouse_alignmentScorePass==TRUE & candidates$mouse_targetBaseIntact ==FALSE)]="Target base not conserved"
candidates$Mm_orthology[which(candidates$mouse_targetOK==TRUE & candidates$concordant==FALSE)]="Target conserved - discordant AA"
candidates$Mm_orthology[which(candidates$concordant)]="Concordant_AA"
candidates$Mm_orthology[which(candidates$Variant_Classification =="Splice_Site" & candidates$mouse_targetOK==TRUE)]="Splice, conserved nt"
candidates$Mm_orthology[which(candidates$Variant_Classification =="Splice_Site" & candidates$mouse_targetOK==FALSE)]="Low homology in target region"
candidates$Mm_orthology[which(candidates$Variant_Classification =="Splice_Site" & candidates$mouse_alignmentScorePass==TRUE & candidates$mouse_targetBaseIntact ==FALSE)]="Splice, Target base not conserved"
candidates$Mm_orthology[which(is.na(candidates$mouse.ensg))]="No orthologous gene"


	if ("maffile" %in% ls()){
	saveRDS(candidates,file=paste0(substr(Sys.time(),0,10),"_from",unlist(strsplit(maffile,split="[.]"))[1],".mousealigned.AA.candidates.rds"))
} else saveRDS(candidates,file=paste0(substr(Sys.time(),0,10),".mousealigned.AA.candidates.rds"))

		

###########################################################################################################
#mouse guide search - if prealigned candidates file available. requires cas9 and editor characteristics inputs
#candidates=readRDS(url("https://wcm.box.com/shared/static/utvxwzt26fstjs0f4zmq02xbroe2p6ho.rds"))



ptm <- proc.time()
run.idx <- which(candidates$mouse_targetOK==TRUE)

output <- foreach(m = run.idx) %dopar% {
#make a new table for mouse outputs (number of mouse guides per variant is independent of human guides)
mouse_targets=candidates[1,]
mouse_targets$sgbot=mouse_targets$sgtop=mouse_targets$sensorv4=mouse_targets$mouse_targetBase_position=mouse_targets$mouse_adjacentEdit=mouse_targets$mouse_PAMseq=mouse_targets$mouse_window=mouse_targets$mouse_sgRNA=mouse_targets$mouse_cas=mouse_targets$mouse_editor=mouse_targets$mouse_PAMIntact=NA 
mouse_targets= mouse_targets[-1,]	
#new loop for guide design
#having looked up flanking seq for each candidate, look for possible guides in mouse genome
ms=0
#for (m in which(candidates$mouse_targetOK==TRUE)){ #replaced by %dopar%
	#print(m)
		#cycle through editors, cas9s
		for (c in 1:nrow(cas9s)){
		#define rules for this specific cas9
		cas=cas9s$Cas_ortholog[c] #which cas9 ortholog?
		PAM=unlist(strsplit(cas9s$PAM[c],split=",")) #cas9-specific PAM required
		PL=nchar(PAM) #length of the PAM
		GL=cas9s$gRNA_Length[c] #length of the guide
			for (e in 1:nrow(edit)){
				#define rules for this specific editor
				editor=edit$Editor[e] #pick editor
				PDmin=GL-edit$window_end[e]+1 #editor-specific window	
				if(edit$window_end[e]<0){PDmin =GL-edit$window_end[e]}	
				PDmax=GL-edit$window_start[e]+1
				if(edit$window_start[e]<0){PDmax=GL-edit$window_start[e]} #fixed for -1 editing
				BE=unlist(strsplit(edit$BE[e],split="[|]")) #which bases can be edited by this enzyme?
				#does this mutation match this enzyme conversion?
				if (candidates$change[m] %in% BE){
				forward=BE[1]
				reverse=BE[2]
		
				tmpflank=NA
				if (candidates$change[m]==forward){tmpflank=candidates$mouse_flank[m]} 
				if (candidates$change[m]==reverse){tmpflank=as.character(reverseComplement(DNAString(candidates$mouse_flank[m])))} 
					#scan through the possible guide spacing and check for a PAM 				
					for (s in PDmin:PDmax){
						#scan through the possible guide spacing and check for a PAM
						for (p in 1:length(PAM)){ 
							PL=nchar(PAM[p]) #length of the PAM						
						if(countPattern(DNAString(PAM[p]),DNAString(substr(tmpflank,31+s,31+s+PL-1)),fixed=FALSE)>0){ #fixed=FALSE to include wildcards
						ms=ms+1  #if PAM is good, make a new entry in mouse_targets
						mouse_targets[ms,1:ncol(candidates)]=candidates[m,]
						mouse_targets$mouse_flank[ms]=tmpflank
						mouse_targets$mouse_PAMIntact[ms]=TRUE
						mouse_targets$mouse_cas[ms]=cas
						mouse_targets$mouse_editor[ms]=editor
						mouse_targets$mouse_sgRNA[ms]=substr(tmpflank,31+s-GL,30+s)
						mouse_targets$mouse_sgRNA_G[ms]=paste0("G",substr(mouse_targets$mouse_sgRNA[ms],2,nchar(mouse_targets$mouse_sgRNA[ms]))) #substitute first nt of guide and guide within sensor with 'G'
						mouse_targets$mouse_PAMseq[ms]=PAM[p] #exact: substr(tmpflank,31+s,31+s+PL[p]-1)
						mouse_targets$mouse_window[ms]=substr(mouse_targets$mouse_sgRNA[ms],edit$window_start[which(edit$Editor== editor)],edit$window_end[which(edit$Editor== editor)])
						
						mouse_targets $mouse_window[ms]=substr(mouse_targets $mouse_flank[ms],30+s-GL+edit$window_start[e],30+s-GL+edit$window_end[e]) #version for edits in/after guide
						if (edit$window_start[e]<0){mouse_targets $mouse_window[ms]=substr(mouse_targets $mouse_flank[ms],31+s-GL+edit$window_start[e],30+s-GL+edit$window_end[e])} #if edited outside of guide range
						
						mouse_targets$mouse_adjacentEdit[ms]=countPattern(substr(BE[1],0,1),DNAString(mouse_targets$mouse_window[ms]))-1
						
						mouse_targets$mouse_targetBase_position[ms]=GL-s+1
						if (mouse_targets $mouse_targetBase_position[ms]<=0){mouse_targets$mouse_targetBase_position[ms]=mouse_targets $mouse_targetBase_position[ms]-1} #if edited outside of guide range
			
						mouse_targets$sensorv4[ms]= paste0("CATAGCGTACACGTCTCACACC", mouse_targets$mouse_sgRNA_G[ms], "GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTT", substr(mouse_targets$mouse_flank[ms],31+s-GL-10,30+s+10), "GAATTCTAGATCCGGTCGTCAACGGCAC")
						mouse_targets$sgtop[ms]=paste0("CACC",as.character(mouse_targets $mouse_sgRNA_G[ms]))
						mouse_targets$sgbot[ms]=paste0("AAAC", as.character(reverseComplement(DNAString(mouse_targets $mouse_sgRNA_G[ms]))))
						
						} #end if statement for mouse guide definition if good PAM is found in range
						}#end for PAM loop
					}#end for loop scan through window looking for this cas-editor combo
				}# end if statement testing if target change is in mouse		
			}#end loop through candidates looking for this editor
		}#end loop for this cas
		mouse_targets
}#end loop for this human mutation in candidates object




### now cleanup the final list and only keeps the ones with entries
### new version
names(output) <- candidates$ID[run.idx]
final.res.ms <- output %>% 
  tibble(Name = names(.), content = .) %>%
  mutate(nrow = map_dbl(content, nrow)) %>%
  filter(nrow > 0) %>%
  dplyr::select(-nrow, -Name) %>%
  unnest(content) %>%
  as.data.frame


keepcols=c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Strand", "Variant_Classification", "Reference_Allele", "Variant_Allele", "HGVSp_Short", "ID", "ID_AA", "change", "human_ref_AA", "human_var_AA", "numberMutationsObserved", "flank", "human.strand", "mouse.symbol", "mouse.ensg", "mouse.chromosome", "mouse.start", "mouse.end", "mouse.strand", "mouse.percent.ID", "mouse_flank", "mouse_targetBaseIntact",  "human_pos_AA", "mouse_transcript", "mouse_pos_AA", "mouse_target_in_codon", "mouse_target_in_CDS", "mouse_target_codon_ref", "mouse_target_codon_var", "mouse_ref_AA", "mouse_var_AA", "mouse_HGVSp_Short",  "mouse_editor", "mouse_cas", "mouse_sgRNA", "mouse_window", "mouse_PAMseq", "mouse_adjacentEdit", "mouse_targetBase_position", "Mm_orthology", "sensorv4", "mouse_sgRNA_G","sgtop","sgbot","numberSubstitutionsObserved","mouse_alignmentScore")

final.res.ms = final.res.ms[,which(colnames(final.res.ms) %in% keepcols),]

proc.time() - ptm


						
		if ("maffile" %in% ls()){
saveRDS(final.res.ms,file=paste0(substr(Sys.time(),0,10),"_from",unlist(strsplit(maffile,split="[.]"))[1],".AMINEsearch_mouse.targets.rds"))
} else saveRDS(final.res.ms,file=paste0(substr(Sys.time(),0,10),".AMINEsearch_mouse.targets.rds"))

				




}#end mouse search


