#AMINEsearch (Annotated Mutation-Informed Nucleotide Editing sgRNA search)
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
install.packages('gProfileR')
}

#load packages
library(BSgenome.Hsapiens.UCSC.hg19) 
library(Biostrings)
library(biomaRt)
library(hom.Hs.inp.db) 
library(org.Mm.eg.db) 
library(AnnotationDbi)
library(gProfileR)

###############################################################################################################
#settings - to be specified by end user
folder="/Users/ted/Documents/GitHub/AMINEsearch"
setwd(folder)
maffile="data_mutations_extended_20180227_patched.txt"

#make sure that the following columns are present with these exact names: 
#Hugo_Symbol	HGVSp_Short	Reference_Allele		Tumor_Seq_Allele2		Variant_Classification	Start_Position		End_Position	Strand	
#check variant classifications to match line 48

#define Cas9 and editor space
cas9s=read.table(file="cas9_properties_Dow.txt",stringsAsFactors=FALSE,sep="\t",header=TRUE)
edit=read.table(file="editors_Dow.txt",stringsAsFactors=FALSE,sep="\t",header=TRUE)


humanBarcode="BC1"
mouseBarcode="BC3"
BC= humanBarcode
BClib=read.table("FSR_OligoBC_lookup.txt",stringsAsFactors=FALSE,header=T,sep="\t")
#worklist
mousesearch=TRUE
###############################################################################################################everything else should be universal

#load .maf file for mutations - can use any standard .maf, double check heading names are consistent for new inputs. CCLE has 122,566 mutations
read.delim(maffile,stringsAsFactors=FALSE)->maf

maf$ID=paste(maf$Hugo_Symbol,maf$HGVSp_Short,sep=" ")
maf$recurrent=duplicated(maf$ID)
maf$change=paste(maf$Reference_Allele,maf$Tumor_Seq_Allele2,sep="")
candidates0 =maf[which(maf$Variant_Classification %in% c("Nonsense_Mutation", "Missense_Mutation", "Splice_Site", "Nonstop_Mutation") 
	& maf$recurrent==TRUE
	   ),]

#condense to one entry per mutant (already selected for recurrent mutations only)
candidates0$duplicated=duplicated(candidates0$ID)   
candidates=candidates0[which(candidates0$duplicated==FALSE),]	
for (c in 1:nrow(candidates)){candidates$numberMutationsObserved[c]=length(which(maf$ID==candidates$ID[c]))}
keepcols=c("Hugo_Symbol","Chromosome","Start_Position","End_Position","Strand","Reference_Allele","Tumor_Seq_Allele2","Variant_Classification","HGVSp_Short","ID","change", "numberMutationsObserved","flank")
candidates=candidates[which(candidates$numberMutationsObserved>3),which(colnames(candidates) %in% keepcols)]



# Function to pull flanking sequence. Defaults to +/- 30 bp
getflank <- function(chr, position,allele, offset=30) {
leftflank <- getSeq(Hsapiens,paste("chr",chr,sep=""),position-offset,position-1)
rightflank <- getSeq(Hsapiens,paste("chr",chr,sep=""),position+1,position+offset)
paste(leftflank,allele,rightflank,sep="")
}


#initialize new fields, objects
candidates$flank=NA
#cycle through each mutation>cas9>editor>possible spacing range
for (m in 1:nrow(candidates)){
	#get flanking sequence for each candidate mutation
	candidates$flank[m]=getflank(candidates$Chromosome[m],candidates$Start_Position[m],candidates$Reference_Allele[m])
	}
save(candidates,file=paste0(substr(Sys.time(),0,10),unlist(strsplit(maffile,split="[.]"))[1],"_candidates.Rdata"))

#load(file="20190616-candidates-mousealigned.Rdata")

targets=candidates[0,]
t=0	
for (m in 1:nrow(candidates)){	
	for (c in 1:nrow(cas9s)){
		#define rules for this specific cas9
		cas=cas9s$Cas_ortholog[c] #which cas9 ortholog?
		PAM=unlist(strsplit(cas9s$PAM[c],split=",")) #cas9-specific PAM required
		GL=cas9s$gRNA_Length[c] #length of the guide
			for (e in 1:nrow(edit)){
				#define rules for this specific editor
				editor=edit$Editor[e] #which editor?
				PDmin=GL-edit$window_end[e]+1 #editor-specific window
				PDmax=GL-edit$window_start[e]+1
				BE=unlist(strsplit(edit$BE[e],split="[|]")) #which bases can be edited by this enzyme?
				
				#does this mutation match this enzyme conversion?
				if (candidates$change[m] %in% BE){
				forward=BE[1]
				reverse=BE[2]
				#is the desired editing on the reverse strand? If so look at the reverse complement <<<<<<<<<<<<<<<<<<<<<fix

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
		targets$sensorv4[t]=paste0("CATAGCGTACACGTCTCACACC",targets$sgRNA_G[t],"GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTT", substr(targets$flank[t],31+s-GL-10,30+s+10), "GAATTCTAGATCCGGTCGTCAACGGCAC")
						targets$sensorv5[t]=paste0("CATAGCGTACACGTCTCACACC",targets$sgRNA_G[t],"GTTTAAGTGCAGGTTGCGAAATGCAT", substr(targets$flank[t],31+s-GL-10,30+s+10), "GAATTCTAGATCCGGTCGTCAACGGCAC")
						targets$window[t]=substr(targets$sgRNA[t],edit$window_start[e],edit$window_end[e])
						targets$PAM[t]=as.character(PAM[p])
						targets$cas9[t]= cas
						targets$editor[t]= editor
						targets$adjacentEdit[t]=countPattern(substr(BE[1],0,1),DNAString(targets$window[t]))-1
						targets$targetBase_position[t]=GL-s+1
						}#end for if statement
						}#end for pam
					}#end of sliding window
				} #end of base change check 
			}#end for editor
	}#end for cas9 ortholog
}#end for candidate mutation site
row.names(targets)=seq(from=1,to=nrow(targets))





keepcols=c("HGVSp_Short","Variant_Classification", "Chromosome", "Start_Position", "End_Position","Hugo_Symbol", "ID", "change", "numberMutationsObserved", "flank", "sgRNA","sgRNA_G", "window", "PAM", "cas9","editor","adjacentEdit","targetBase_position","sensorv4","sensorv5")
targets=targets[,which(colnames(targets) %in% keepcols)]

#oligo's for 1-by-1 or library format
#can vary barcodes to make complex libraries

LBC=BClib$LBC[which(BClib$Barcode.ID==BC)]
LAD=BClib$LAD[which(BClib$Barcode.ID==BC)]
RAD=BClib$RAD[which(BClib$Barcode.ID==BC)]
RBC=BClib$RBC[which(BClib$Barcode.ID==BC)]


for (i in 1:nrow(targets)){
	targets$sgtop[i]=paste0("CACC",as.character(targets$sgRNA_G[i]))
	targets$sgbot[i]=paste0("AAAC", as.character(reverseComplement(DNAString(targets$sgRNA[i]))))
	targets$oligo.lib.BC[i]=BC
	targets$oligo.lib[i]=paste(LBC,LAD,targets$sg[i],RAD,RBC,sep="")
}


write.table(targets,file=paste0("IMPACT_",substr(Sys.time(),0,10),"AMINEsearch_human.txt"),sep="\t",row.names=FALSE)
#save(targets,file="IMPACT_HU_targets.Rdata")




############################end human
###########################################################################################################
if (mousesearch==TRUE){
library(biomaRt)
setwd(folder)
#load(file="IMPACT_HU_targets_patched.Rdata")
#load(file=paste0(unlist(strsplit(maffile,split="[.]"))[1],"_candidates.Rdata"))

candidates=candidates[order(candidates$Hugo_Symbol),]

#refer to mouse ensembl to get mouse orthologs
#load Refseq (mouse>human orthologs)

load("REFSEQ_IMPACT.Rdata") #replace with full refseq for general case/ when distributing
RefSeq=BMH

candidates$human.strand= RefSeq$strand[match(candidates$Hugo_Symbol, RefSeq$external_gene_name)]
candidates$mouse.ensg=candidates$mouse.symbol= FALSE
candidates$mouse.symbol= RefSeq$mmusculus_homolog_associated_gene_name[match(candidates$Hugo_Symbol, RefSeq$external_gene_name)]
candidates$mouse.ensg = RefSeq$mmusculus_homolog_ensembl_gene[match(candidates$Hugo_Symbol, RefSeq$external_gene_name)]
candidates$mouse.chromosome=RefSeq$mmusculus_homolog_chromosome[match(candidates$Hugo_Symbol, RefSeq$external_gene_name)]
candidates$mouse.start=RefSeq$mmusculus_homolog_chrom_start[match(candidates$Hugo_Symbol, RefSeq$external_gene_name)]
candidates$mouse.end=RefSeq$mmusculus_homolog_chrom_end[match(candidates$Hugo_Symbol, RefSeq$external_gene_name)]
candidates$mouse.strand= RefSeq$mmusculus_strand[match(candidates$Hugo_Symbol, RefSeq$external_gene_name)]
candidates$mouse.percent.ID=RefSeq$mmusculus_homolog_perc_id[match(candidates$Hugo_Symbol, RefSeq$external_gene_name)]
mouse.ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host = "www.ensembl.org") #for mouse

# define a function to fetch mouse sequence
find_gene_sequence <- function(query_gene){
require(biomaRt)
if(length(query_gene) != 1){stop("Query a single gene")}
genes = getBM(attributes=c("gene_exon_intron", "ensembl_gene_id","mgi_symbol"),filters="ensembl_gene_id",values=query_gene, mart=mouse.ensembl)
if(nrow(genes) != 0){return(genes)}
warning(query_gene, " does not match any mouse MGI gene symbol or ensembl gene ID")
return(NULL)
}#end function to fetch mouse seq
	



#cycle through 'candidates' and include mouse guides that target mutations uneditable in human genome
#first retreive analogous flanking sequence from mouse genome
for (m in 1:nrow(candidates)){
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
	
	
	#if the human gene is on -1 strand, we must RC to align with mouse reference transcript, which is always delivered in sense direction
	
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

##save(candidates,file="20190614-candidates_aligned.Rdata")



#make a new table for mouse outputs (there could be multiple mouse guides per human guide)
mouse_targets=candidates[1,]
mouse_targets$sensorv5=mouse_targets$sensorv4=mouse_targets$mouse_targetBase_position=mouse_targets$mouse_adjacentEdit=mouse_targets$mouse_PAMseq=mouse_targets$mouse_window=mouse_targets$mouse_sgRNA=mouse_targets$mouse_cas=mouse_targets$mouse_editor=mouse_targets$mouse_PAMIntact=NA 	
#new loop for guide design
#having looked up flanking seq for each candidate, look for possible guides in mouse genome
ms=0
for (m in which(candidates$mouse_targetOK==TRUE)){ 
	print(m)
		#cycle through editors, cas9s
		for (c in 1:nrow(cas9s)){
		#define rules for this specific cas9
		cas=cas9s$Cas_ortholog[c] #which cas9 ortholog?
		PAM=unlist(strsplit(cas9s$PAM[c],split=",")) #cas9-specific PAM required
		PL=nchar(PAM) #length of the PAM
		GL=cas9s$gRNA_Length[c] #length of the guide
			for (e in 1:nrow(edit)){
				#define rules for this specific editor
				editor=edit$Editor[e] #which editor?
				PDmin=GL-edit$window_end[e]+1 #editor-specific window
				PDmax=GL-edit$window_start[e]+1
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
						mouse_targets$mouse_adjacentEdit[ms]=countPattern(substr(BE[1],0,1),DNAString(mouse_targets$mouse_window[ms]))-1
						mouse_targets$mouse_targetBase_position[ms]=GL-s+1

						#mouse_targets$mouse_sensor[ms]=paste0("GGTACACGTCTCACACC",mouse_targets$mouse_sgRNA_G[ms],"GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTT", substr(mouse_targets$mouse_flank[ms],31+s-GL-10,30+s+10), "GCATGAATTCTAGATCCGG")
						
						mouse_targets$sensorv4[ms]= paste0("CATAGCGTACACGTCTCACACC", mouse_targets$mouse_sgRNA_G[ms], "GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTT", substr(mouse_targets$mouse_flank[ms],31+s-GL-10,30+s+10), "GAATTCTAGATCCGGTCGTCAACGGCAC")
						
						mouse_targets$sensorv5[ms]=paste0("CATAGCGTACACGTCTCACACC", mouse_targets$mouse_sgRNA_G[ms],"GTTTAAGTGCAGGTTGCGAAATGCAT", substr(mouse_targets$mouse_flank[ms],31+s-GL-10,30+s+10), "GAATTCTAGATCCGGTCGTCAACGGCAC")
							
						} #end if statement for mouse guide definition if good PAM is found in range
						}#end for PAM loop
					}#end for loop scan through window looking for this cas-editor combo
				}# end if statement testing if target change is in mouse		
			}#end loop through candidates looking for this editor
		}#end loop for this cas
}#end loop for this human mutation in candidates object
row.names(mouse_targets)=seq(from=1,to=nrow(mouse_targets))




write.table(mouse_targets,file=paste0("IMPACT","_",substr(Sys.time(),0,10),"AMINEsearch_mouse.txt"),sep="\t",row.names=FALSE)
#mouse_targets=read.table(file="",stringsAsFactors=FALSE,sep="\t",header=TRUE)


targets.abbrev$mouse.sg=targets$mouse.sg
targets.abbrev$mousePAMseq=targets$mousePAMseq

###make oligo library with barcodes for mouse
#can vary barcodes to make complex libraries
BC= mouseBarcode
LBC=BClib$LBC[which(BClib$Barcode.ID==BC)]
LAD=BClib$LAD[which(BClib$Barcode.ID==BC)]
RAD=BClib$RAD[which(BClib$Barcode.ID==BC)]
RBC=BClib$RBC[which(BClib$Barcode.ID==BC)]


targets.abbrev$mouse.sgtop=targets.abbrev$mouse.sgbot=targets.abbrev$mouse.oligo.lib=NA
for (i in which(targets$mouse.sg!=FALSE)){
	targets$mouse.sgtop[i]=paste("CACCG",as.character(targets$mouse.sg[i]),sep="")
	targets$mouse.sgbot[i]=paste("AAAC", as.character(reverseComplement(DNAString(targets$mouse.sg[i]))) ,"C",sep="")
	targets$mouse.oligo.lib[i]=paste(LBC,LAD,targets$mouse.sg[i],RAD,RBC,sep="")
}

targets.abbrev$mouse.sgtop=targets$mouse.sgtop
targets.abbrev$mouse.sgbot=targets$mouse.sgbot
targets.abbrev$mouse.oligo.lib=targets$mouse.oligo.lib
targets.abbrev$mousePAMseq=targets$mousePAMseq


}#end mousesearch
###########################################################################################################
#output unique sensor oligos Cas9 or xCas9 APOBEC or 2X
lib.hu=unique(targets$sgRNA_G[which( (targets$cas9 %in% c("Sp Cas9","xCas9","Sc Cas9")) & (targets$editor=="2X" | targets$editor=="APOBEC") )])
hi=match(lib.hu ,targets$sgRNA_G)
lib.hu.lab=data.frame(ID=targets$ID[hi],sgRNA=lib.hu,sensorv4=targets$sensorv4[hi],sensorv5=targets$sensorv5[hi])



lib.ms=unique(mouse_targets$mouse_sgRNA_G[which( (mouse_targets$mouse_cas  %in% c("Sp Cas9","xCas9","Sc Cas9"))& (mouse_targets$mouse_editor=="2X" | mouse_targets$mouse_editor=="APOBEC") )])
mi=match(lib.ms ,mouse_targets$mouse_sgRNA_G)
lib.ms.lab=data.frame(ID=mouse_targets$ID[mi],sgRNA= lib.ms,sensorv4=mouse_targets$sensorv4[mi],sensorv5=mouse_targets$sensorv5[mi])


write.table(lib.hu.lab,file="unique_sensor_human_v20190717.txt",sep="\t",row.names=FALSE)
write.table(lib.ms.lab,file="unique_sensor_mouse_v20190717.txt",sep="\t",row.names=FALSE)

