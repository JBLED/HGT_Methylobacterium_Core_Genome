#### Definition of Methylobacterium Core genome

#Package seqinr to play with nucleotide sequences
library(seqinr)

#Package to align sequences
library(msa)

#package for parallelization
library(parallel)

#Number of cores to use
NC=8

pD<-"/myfolder/2-Core-genome/"

###STEP A ### Genomes and annotation statistics
###STEP B ### Define CORE genome 
###STEP C ### Generate fasta files for core genes
###STEP D ### Align each core gene (do separately for protein-encoding genes and rna genes)
###STEP E ### calculation of Average nucleotide identity among genomes using fastANI
###STEP F ### Function (KEGG analysis with KASS)
###STEP G ### for each gene alignment, calculate dN/dS (dNdS) and number of polymorphic sites per gene alignment (nPOL) that will be used to correct RF distances

###STEP A ##
##Genomes and annotation statistics 

setwd(pD)

#List of genome assemblies (Groups and Species are according to Leducq et al 2022 GBE for st1 to st213; groups are left undetermined (X) for st214-319; species only indicated if newly described)
data<-read.table("data.txt",header=T)

#Type of annotations
TYPE<-c("peg","rna","repeat")

#This hidden code is to convert new RAST output files into readable format
#sapply(c(311:317),function(x){
#	NX<-data[x,6]	
#	ANNOx<-read.delim(paste(pD,NX,sep="/"))	
#	write.table(ANNOx ,paste(pD,"RAST_files",NX,sep="/"))
#})

STAT<-mclapply(c(1:length(data[,1])),function(x){
	print(x)
	STx<-data[x,2]
	NX<-data[x,6]
	
	#Genome characteristics
	seqx<-read.fasta(paste(pD,"Genome_files/",data[x,7],sep=""))
	
	#number of contigs
	NCx=length(seqx)
	#contig size
	SCx<-sort(as.numeric(summary(seqx)[,1]),decreasing=T)
	#Total size
	TSx=sum(SCx)
	#N50
	N50x=subset(SCx,(cumsum(SCx)/TSx)>0.5)[1]
	#GC
	GCx<-GC(unlist(seqx))
	
	#Annotations characteristics
	ANNOx<-read.table(paste(pD,"RAST_files",NX,sep="/"))
	
	#number and type of annotations
	ANx<-sapply(TYPE,function(x){
		return(length(subset(ANNOx[,3],ANNOx[,3]==x)))
	})
	ANx<-c(ANx,dim(ANNOx)[1]-sum(ANx))
	sort(summary(as.factor(ANNOx[,3])),decreasing=T)
	#number of unique annotations
	AUx<-length(unique(ANNOx[,8]))
	#number of hypothetical proteins
	HYPx<-length(subset(ANNOx[,8],ANNOx[,8]=="hypothetical protein"))
	
	return(c(STx, TSx, NCx, N50x, GCx, ANx, AUx, HYPx))
	
},mc.cores=NC)

STAT<-matrix(unlist(STAT),ncol=11,byrow=T)

colnames(STAT)<-c("isolate","size","contigs","N50","GC","peg","rna","repeat","other","unique","hypothetical")

setwd(pD)

write.table(STAT,"Stat-genomes.txt")

###STEP B ##
##Define CORE genome 
####From here, exclude Microvirga and Enterovirga (293 genomes)

setwd(pD)

STAT<-read.table("Stat-genomes.txt",header=T)
data<-read.table("data.txt",header=T)

iMETH<-sort(intersect(subset(c(1:length(data[,1])),data[,4]!="M"),
	subset(c(1:length(data[,1])),data[,4]!="E")))

#List of unique gene annotations (9261, including “hypothetical protein”, “repeat region” and “mobile element protein”)

ANNOf<-sort(unique(unlist(mclapply(iMETH,function(x){
	print(x)
	NX<-data[x,6]
	ANNOx<-read.table(paste(pD,"RAST_files",NX,sep="/"))
	return(unique(ANNOx[,8]))
},mc.cores=NC))))

#Frequence of each annotation

ANNOa<-sapply(iMETH,function(x){
	print(x)
	NX<-data[x,6]
	ANNOx<-summary(as.factor(as.vector(read.table(
	paste(pD,"RAST_files",NX,sep="/"),header=T)[,8])),
		maxsum=length(ANNOf))

	ANNOo<-setdiff(ANNOf,names(ANNOx))
	ABx<-seq(0,0,length.out=length(ANNOo))
	names(ABx)<-ANNOo
	
	ANNOx<-c(ANNOx, ABx)
	return(ANNOx[order(names(ANNOx))])
	
})

colnames(ANNOa)<-data[iMETH,2]

setwd(pD)

write.table(ANNOa,"Annotation-occurence-raw.txt")

###Definition of the core genome. Only keep genes present in at least PC% of assemblies  and present in average in less than COP copies among assemblies - Calculate for each group (A, B, D, C) separatly in order to avoid biases induced by gene loss/duplication in some groups

setwd(pD)

STAT<-read.table("Stat-genomes.txt",header=T)
data<-read.table("data.txt",header=T)
ANNOa<-read.table("Annotation-occurence-raw.txt",header=T)

PC=0.8
COP=1.5

#average gene occurence (old; only if we don't care about the group biase)
#OCC<-as.vector(apply(ifelse(ANNOa==0,0,1),1,sum))/length(iMETH) 

#average gene copy number (old; only if we don't care about the group biase)
#DUP<-as.vector(apply(ANNOa,1,mean))

#Average gene occurence per Group
OCCg<-sapply(c("A","B","D","C"),function(x){
	iX<-subset(c(1:length(iMETH)),data[iMETH,4]==x)
	return(as.vector(apply(ifelse(ANNOa[, iX]==0,0,1),1,sum))/length(iX))
})

#average gene copy number per group
DUPg<-sapply(c("A","B","D","C"),function(x){
	iX<-subset(c(1:length(iMETH)),data[iMETH,4]==x)
	return(as.vector(apply(ANNOa[,iX],1,mean)))
})


#index of candidate core genes (1492 for PC=0.8 and COP=1.5 old; only if we don't care about the group biase)
#z<-intersect(subset(c(1:length(ANNOa[,1])), OCC>=PC),
#	subset(c(1:length(ANNOa[,1])), DUP<= COP)) 

#index of candidate core genes (1419 for PC=0.8 and COP=1.5)
z<-intersect(subset(c(1:length(ANNOa[,1])), apply(OCCg,1,min)>=PC),
	subset(c(1:length(ANNOa[,1])), apply(DUPg,1,min)<= COP)) 


# size of each copy of each candidate core gene in each assemblie
Sizec<-sapply(iMETH,function(x){
	print(x)
	NX<-data[x,6]
	ANNOx<-read.table(paste(pD,"RAST_files",NX,sep="/"))
	sapply(z,function(x){
		return(paste(nchar(subset(
			as.vector(ANNOx[,12]),
			as.vector(ANNOx[,8])== rownames(ANNOa)[x])),
			collapse="_"))
	})
})

# index of assemblies with N50 > 0.5Mb
iC<-subset(c(1:length(iMETH)),STAT[iMETH,4]>=500000)

# size of each copy of each candidate core gene normalized by average size in assemblies with N50 > 0.5Mb
SizeM<-sapply(c(1:length(Sizec[,1])),function(x){
	sx<-Sizec[x,]
	#mean size in assemblies with N50 > 0.5Mb 
	mx<-mean(as.numeric(
		unlist(strsplit(sx[iC],split="_"))))
		#unlist(strsplit(sx[iC],split="_"))))
	#normalized copy size
	sapply(c(1:length(sx)),function(x){
		return(paste(round(
			as.numeric(unlist(strsplit(
				sx[x],split="_")))/mx,2),collapse="_"
		))
	})	
})

hist(as.numeric(unlist(strsplit(SizeM,split="_"))),breaks=1000)

#real copy number per candidate core gene
RC<-sapply(c(1:length(SizeM[1,])),function(x){
	sx<-SizeM[,x]
	return(sapply(c(1:length(sx)),function(x){
		return(sum(as.numeric(unlist(
			strsplit(sx[x],split="_")))))
	}))
})

hist(RC,breaks=1000)

#replace values outside the range 0.7-1.3 as missing data
RC <-ifelse(RC <0.7,0,ifelse(RC>1.3,0,1))

#count observed duplications as missing data
RC <-ifelse(t(ANNOa[z,])==1, RC,0)

#Number of candidate core genes in function of N50 (points colored according to group)

COLg<-ifelse(data[,4]=="A","red",
	ifelse(data[,4]=="B","purple",
	ifelse(data[,4]=="C","green",
	ifelse(data[,4]=="D","blue","yellow"))))[iMETH]

plot(as.numeric(STAT[iMETH,4])/1000,apply(RC,1,sum),
	log="x",pch=21,bg=COLg,
	las=1,cex.axis=0.8,xlab="N50 (kB,log scale)",
	ylab="Core genes")

#Proportion of available data for each candidate core gene and each group

RCg<-sapply(c("A","B","D","C"),function(x){
	iX<-subset(c(1:length(iMETH)),data[iMETH,4]==x)
	return(as.vector(apply(RC[iX,],2,mean)))
})

#keep genes with data in at least PC2% of genomes, estimate separately for each group (1122 genes with PC2=0.7)

PC2=0.7

z2<-subset(z,apply(RCg,1,min)>= PC2)
RC2<-t(subset(t(RC),apply(RCg,1,min)>= PC2))

plot(as.numeric(STAT[iMETH,4])/1000,apply(RC2,1,sum),
	log="x",pch=21,bg=COLg,
	las=1,cex.axis=0.8,xlab="N50 (kB,log scale)",
	ylab="Core genes")


setwd(pD)

write.table(RC2,"Methylobacterium-core-gene-occ.txt")

#Core gene type (rna or peg; will be important for the type of model to consider for phylogeny inference)
CGtype<-sapply(iMETH,function(x){
	print(x)
	NX<-data[x,6]
	ANNOx<-read.table(paste(pD,"RAST_files",NX,sep="/"))
	as.vector(sapply(colnames(RC2),function(x){
		return(c(subset(
			ANNOx[,3],ANNOx[,8]== x),NA)[1])
	}))
})

CGtype<-sapply(c(1:length(CGtype[,1])),function(x){
	return(setdiff(CGtype[x,],NA))
})

#Give a short name for each core gene
CGname<-paste("CG",c(1:length(CGtype)),sep="_")

CGstat<-cbind(CGname,colnames(RC2), CGtype,
	as.numeric(apply(RC2,2,sum)))

colnames(CGstat)<-c("Name","Annotation","Type","Occurence")

write.table(CGstat,
	"Core-gene-code.txt",row.names=F)
	
pdf("CoreGenome_N50.pdf",width=6,height=6)
plot(as.numeric(STAT[iMETH,4])/1000,apply(RC2,1,sum),
	log="x",pch=21,bg=COLg,
	las=1,cex.axis=0.8,xlab="N50 (kB,log scale)",
	ylab="Core genes")
dev.off()	
	
###STEP C ##
#Generate fasta files for core genes

setwd(pD)

STAT<-read.table("Stat-genomes.txt",header=T)
data<-read.table("data.txt",header=T)
ANNOa<-read.table("Annotation-occurence-raw.txt",header=T)
RC2<-read.table("Methylobacterium-core-gene-occ.txt",header=T)
CGstat<-read.table("Core-gene-code.txt",header=T)
CGname<-CGstat[,1]

iMETH<-sort(intersect(subset(c(1:length(data[,1])),data[,4]!="M"),
	subset(c(1:length(data[,1])),data[,4]!="E")))

#write a fasta file with core genes for each genome
setwd(pD)
dir.create("core-genes")


sapply(c(1:length(iMETH)),function(x){
	print(x)
	tx<-data[iMETH[x],2]
	rcx= RC2[x,]
	ix<-subset(c(1:length(rcx)), rcx!=0)
	print(tx)
	###Annotation file
	NX<-data[iMETH[x],6]
	ANNOx<-read.table(paste(pD,"RAST_files",NX,sep="/"))

	SEQx<-sapply(ix,function(x){
		return(list(unlist(strsplit(subset(
			as.vector(ANNOx[,12]),
			as.vector(ANNOx[,8])== CGstat[x,2]),
			split=""))))
	})
	setwd(paste(pD,"core-genes",sep="/"))
	write.fasta(SEQx,names= CGname[ix],
		file.out=paste(tx,".fas",sep=""))
	return("")	
})

#Write one fasta per core gene

setwd(paste(pD,"core-genes",sep="/"))

mclapply(c(1:length(CGname)),function(x){
	Gx= CGname[x]
	rcx= RC2[,x]
	ix<-subset(c(1:length(iMETH)), rcx!=0)
	print(x)
	SEQx<-sapply(ix,function(x){
		#print(x)
		tx<-data[iMETH[x],2]
		seqx<-read.fasta(paste(tx,".fas",sep="")	)
		return(list(unlist(
			seqx[[subset(c(1:length(seqx)),
		names(seqx)==Gx)]])))
	})
	write.fasta(SEQx,names= data[iMETH[ix],2],
		file.out=paste(Gx,".fas",sep=""))
	return("")
},mc.cores=8)

#STEP D
#Align each core gene (do separately for protein-encoding genes and rna genes)


setwd(pD)

STAT<-read.table("Stat-genomes.txt",header=T)
data<-read.table("data.txt",header=T)
ANNOa<-read.table("Annotation-occurence-raw.txt",header=T)
RC2<-read.table("Methylobacterium-core-gene-occ.txt",header=T)
CGstat<-read.table("Core-gene-code.txt",header=T)
CGname<-CGstat[,1]
CGtype<-CGstat[,3]

####Align non protein coding genes
Irna<-subset(c(1:length(CGname)), CGtype =="rna")

setwd(paste(pD,"core-genes",sep="/"))

mclapply(Irna,function(x){	
	Ix<-CGname[x] #short gene name
	print(Ix)
	SEQw<-msa(paste(Ix,
		".fas",sep=""),method="ClustalW",type="dna",order="input")
	#Retrieve aligned sequence in seqinr format
	SEQn<-msaConvert(SEQw,type="seqinr::alignment")$seq
	#Get sequence names from the unaligned nucleotide sequence
	NX=names(read.fasta(paste(Ix,
		".fas",sep="")))
	#Convert
	SEQn<-sapply(c(1:length(SEQn)),function(x){
		return(strsplit(SEQn[x],split=""))
	})
	write.fasta(SEQn, names=NX,
		file.out =paste(Ix,
		"_aligned.fas",sep=""))
		
	file.remove(paste(Ix,
		c(".dnd",".aln"),sep=""))
	
		
},mc.cores=8)




####Align protein encoding genes
Ipeg<-subset(c(1:length(CGname)), CGtype =="peg")

setwd(paste(pD,"core-genes",sep="/"))

mclapply(Ipeg,function(x){
	Ix<-CGname[x] #short gene name
	print(Ix)

	#Get the unaligned nucleotide sequence
	SEQ<-read.fasta(paste(Ix,
		".fas",sep=""))

	#Sequence names
	NX=names(SEQ)	
	
	#Remove gaps
	SEQ<-sapply(c(1:length(SEQ)),function(x){
		return(list(subset(unlist(SEQ[[x]]),
			unlist(SEQ[[x]])!="-")))
	})
	
	#Sequence sizes
	Sx<-as.numeric(summary(SEQ)[,1])
	Sx<-c(Sx,min(Sx)-3)
	
	#Add a vector of "Ns" for potential missing data at the end and begining of each sequence (length: maximum-minimum sequence size)
	nx<-sapply(c(1:(max(Sx)-min(Sx))),function(x){
		return("N")
	})
	
	SEQ<-sapply(c(1:length(SEQ)),function(x){
		return(list(c(nx,as.vector(SEQ[[x]]),nx)))
	})
	
	#Convert into amino-acid sequences
	SEQaa<-getTrans(SEQ)

	#Write a temporary file with amino acid sequences unaligned
	write.fasta(SEQaa,names=NX,paste(Ix,"_AA.fas",sep=""))


	#Align AA sequences together using ClustalW - IMPORTANT: keep same order as input, since sequence names are removed during analysis, and retrieve sequence names in the original unaligned fasta file


	SEQw<-msa(paste(Ix,"_AA.fas",sep=""),method="ClustalW",
		type="protein",order="input")

	#Retrieve aligned sequence in seqinr format
	SEQn<-msaConvert(SEQw,type="seqinr::alignment")$seq
	
	#Convert back in nucleotide sequence

	SEQn<-sapply(c(1:length(NX)),function(x){
		#Unaligned nucleotide sequence
		seqx<-as.vector(SEQ[[x]])
		#Unaligned protein sequence
		sequ<-unlist(strsplit(SEQaa[[x]],split=""))
		#Aligned protein sequence
		seqa<-unlist(strsplit(SEQn[[x]],split=""))
		#Unaligned nucleotide sequence coordinates
		coordx<-c(1:length(seqx))
		#Unaligned proteine sequence coordinates
		coordu<-matrix(coordx,ncol=3,byrow=T)
		#Removed stop codon (because it's automaticaly removed by msa alignement)
		coordu<-subset(coordu, sequ!="*")
		#Aligned protein sequence coordinates
		coorda<-matrix(c(1:(length(seqa)*3)),
			ncol=3,byrow=T)
		#Ungapped positions in Aligned protein sequence coordinates
		coordau<-sort(as.vector(subset(coorda, 
			seqa!="-")))
		#Gapped position in Aligned protein sequence coordinates
	coordag<-sort(as.vector(subset(coorda, 
		seqa=="-")))		
		#Nucleotide sequence (ungapped position)	
		seqx<-seqx[sort(as.vector(coordu))]
		#Nucleotide sequence (Gapped position)
		seqg<-ifelse(seq(0,0,length.out=
			length(coordag))==0,"-","")
		return(list(c(seqg, seqx)
			[order(c(coordag, coordau))]))
	})
	
	#Remove position that only contain missing data
	SEQn<-sapply(c(1:length(NX)),function(x){
		return(as.vector(unlist(	SEQn[[x]])))
	})

	SEQn<-subset(SEQn,sapply(c(1:length(SEQn[,1])),function(x){
		seqx<-SEQn[x,]
		return(length(subset(seqx, seqx=="-"))+
		length(subset(seqx, seqx=="N")))
	})!=length(NX))
	
	SEQn<-sapply(c(1:length(NX)),function(x){
		return(list(SEQn[,x]))
	})
		
	write.fasta(SEQn, names=NX,
		file.out =paste(Ix,
		"_aligned.fas",sep=""))
		
	file.remove(paste(Ix,
		c("_AA.fas","_AA.dnd","_AA.aln"),sep=""))		
		
},mc.cores=8)

#Review each protein-coding alignment and replace "NNN(---)NNN" motifs by "NNN(NNN)NNN"

setwd(paste(pD,"core-genes",sep="/"))

mclapply(c(Ipeg,Irna),function(x){
	
	Ix<-CGname[x] #short gene name
	print(Ix)

	#Get the aligned nucleotide sequence
	SEQ<-read.fasta(paste(Ix,
		"_aligned.fas",sep=""))
	
	SEQn<-sapply(c(1:length(SEQ)),function(x){
		#print(x)
		seqx<-c(as.vector(SEQ[[x]]),"x","n","-","n","x","n","-","n")
		lx=length(seqx)
		z<-c(0,sapply(c(2: lx),function(x){
			ifelse(seqx[x-1]==seqx[x],1,0)
		}))
		seqxx<-subset(seqx,z==0)
		z<-subset(c(1:length(seqx)),z==0)
		
		y<-t(sapply(c(3:length(seqxx)),function(x){
			return(c(
			paste(z[(x-2):x],collapse="_"),
			paste(seqxx[(x-2):x],collapse="_")
			))
		}))
		
		y<-subset(y,y[,2]=="n_-_n")
		
		#positions to replace by 'Ns'
		pn<-unique(unlist(sapply(c(1:length(y[,1])),function(x){
			yx<-as.numeric(unlist(strsplit(y[x,1],split="_")))
			return(c(min(yx):max(yx)))
		})))
		sn<-sapply(pn,function(x){return("n")})
		#positions to keep
		pk<-setdiff(c(1:lx),pn)
		sk<-seqx[pk]
		
		seqn<-c(sn, sk)[order(c(pn,pk))]
		return(list(toupper(seqn[1:(lx-8)])))
	})
	
	write.fasta(SEQn,paste(Ix,"_alignedN.fas",sep=""),names=names(SEQ))
	
},mc.cores=6)


#Check that gap replacement by Ns did not affected alignments

setwd(paste(pD,"core-genes",sep="/"))

#only "n" should be returned

unique(unlist(
	sapply(c(1:length(CGname)),function(x){
	
	Ix<-CGname[x] #short gene name
	print(Ix)
	#Get the aligned nucleotide sequence (before and after N correcting)
	SEQ<-read.fasta(paste(Ix,
		"_aligned.fas",sep=""))	
	SEQn<-read.fasta(paste(Ix,
		"_alignedN.fas",sep=""))
	sumx<-unlist(sapply(c(1:length(SEQ)),function(x){
		seq<-as.vector(SEQ[[x]])
		seqn<-as.vector(SEQn[[x]])
		return(subset(seqn,seq!=seqn))
	}))
	print(c(x,unique(sumx))	)
	return(sumx)
})))


#Proportion of missing data (only for protein coding)
setwd(paste(pD,"core-genes",sep="/"))

PPX<-unlist(sapply(Ipeg,function(x){
	
	Ix<-CGname[x] #short gene name
	print(Ix)
	#Get the aligned nucleotide sequence
	SEQn<-read.fasta(paste(Ix,
		"_alignedN.fas",sep=""))
	NX<-names(SEQn)
	
	#Convert into amino-acid sequences
	SEQaa<-getTrans(SEQn)
	
	#Convert in matrix
	SEQn<-sapply(c(1:length(SEQn)),function(x){
		return(as.vector(SEQaa[[x]]))
	})
	
	#Count the proportion of X per position
	pn<-sapply(c(1:length(SEQn[,1])),function(x){
		length(subset(SEQn[x,],SEQn[x,]=="X"))
	})/length(NX)
	return(pn)
}))

#Remove aa positions with more than ppn% of missing data (only for protein-coding)
ppn=0.9

setwd(paste(pD,"core-genes",sep="/"))


sumNEW<-mclapply(Ipeg,function(x){
	
	Ix<-CGname[x] #short gene name
	print(Ix)
	#Get the aligned nucleotide sequence
	SEQn<-read.fasta(paste(Ix,
		"_alignedN.fas",sep=""))
	NX<-names(SEQn)
	
	#Convert into amino-acid sequences
	SEQaa<-getTrans(SEQn)
		
	#Convert both in matrix
	SEQn<-sapply(c(1:length(SEQn)),function(x){
		return(as.vector(SEQn[[x]]))
	})
	SEQaa <-sapply(c(1:length(SEQaa)),
		function(x){
		return(as.vector(SEQaa[[x]]))
	})
	#Count the proportion of X per position
	pn<-sapply(c(1:length(SEQaa[,1])),
		function(x){
		length(subset(SEQaa[x,], SEQaa[x,]=="X"))
	})/length(NX)
	
	pn <-as.vector(matrix(t(cbind(pn,pn,pn)),
		ncol=1,byrow=T))
	
	SEQn<-subset(SEQn, pn<(1-ppn))
	
	SEQn<-sapply(c(1:length(NX)),function(x){
		return(list(SEQn[,x]))
	})
	
	write.fasta(SEQn,paste(Ix,
		"_alignedNF.fas",sep=""),names=NX)
	return(c(x,Ix,length(SEQn),length(pn),
		length(SEQn[[1]])))	
		
},mc.cores=4)

sumNEW <-matrix(unlist(sumNEW),ncol=5,byrow=T)

#Same with rna

sumNEWrna<-mclapply(Irna,function(x){
	
	Ix<-CGname[x] #short gene name
	print(Ix)
	#Get the aligned nucleotide sequence
	SEQn<-read.fasta(paste(Ix,
		"_alignedN.fas",sep=""))
	NX<-names(SEQn)
	
	write.fasta(SEQn,paste(Ix,
		"_alignedNF.fas",sep=""),names=NX)
	return(c(x,Ix,length(SEQn),length(SEQn[[1]]),
		length(SEQn[[1]])))	
		
},mc.cores=4)

sumNEWrna <-matrix(unlist(sumNEWrna),ncol=5,byrow=T)

sumNEWrna<-cbind(sumNEWrna,"rna")
sumNEW<-cbind(sumNEW,"peg")
sumNEW<-rbind(sumNEWrna, sumNEW)

#remove intermediate files

sapply(c(Ipeg,Irna),function(x){
	
	Ix<-CGname[x] #short gene name
	print(Ix)
	file.remove(paste(Ix,
		c("_aligned.fas","_alignedN.fas"),sep=""))
})	

setwd(pD)
colnames(sumNEW)<-c("Name1","Name2","Sequences","Size1(bp)","Size2(bp)","Type")
write.table(sumNEW,"Summary-per-core-gene-after-alignment.txt",row.names=F)

#STEP E: calculation of Average nucleotide identity among genomes using fastANI

#Conversion of sequence file names in simple names (stXXX.fas) to be liste in a "genome-paths.txt" file

setwd(pD)

data<-read.table("data.txt",header=T)

mclapply(c(1:length(data[,1])),function(x){
	print(x)
	STx<-data[x,2]
	NX<-data[x,6]
	
	#Genome characteristics
	seqx<-read.fasta(paste(pD,"Genome_files/",data[x,7],sep=""))
	
	setwd(paste(pD,"ANI",sep="/"))
	
	write.fasta(seqx,names=names(seqx),paste(STx,".fna",sep=""))
},mc.cores=8)

setwd(paste(pD,"ANI",sep="/"))
write.table(t(t(paste(data[,2],".fna",sep=""))),"genome-paths.txt",quote=F,col.names=F,row.names=F)

#FastANI ran on server with following:

###		#!/bin/bash
###		#SBATCH --mem=16GB
###		#SBATCH -N 1
###		#SBATCH --time=1-00:00:00
###		#SBATCH -p batch_72h

###		module load StdEnv/2020  gcc/9.3.0
###		module load StdEnv/2020  intel/2020.1.217
###		module load fastani/1.32

###		fastANI --ql genome-paths.txt --rl genome-paths.txt -t 4 -o fastANI-out-raw.tsv



###FUNCTIONS
library(parallel)

num.strsplit.u<-function(x){return(as.numeric(unlist(strsplit(x,split= "_"))))}

reduce.group<-function(groups,ncor) {
	reduce.gx<-function(i){
		gx<-num.strsplit.u(groups[i])
		z<-sapply(groups,function(x){
			return(length(intersect(
				gx,num.strsplit.u(x))))
		})
		return(paste(sort(unique(
			num.strsplit.u(names(subset(z,z!=0))))),
			collapse="_"))
	}
	return(unlist(mclapply(c(1:length(groups)), 
		reduce.gx,mc.cores=ncor)))
}


setwd(paste(pD,"ANI",sep="/"))

#FastANI output
ANI<-read.table("fastANI-out-raw.tsv")

ANI<-cbind(gsub(".fna","",ANI[,1]),gsub(".fna","",ANI[,2]),ANI[,3])


iMETH<-sort(intersect(subset(c(1:length(data[,1])),data[,4]!="M"),
	subset(c(1:length(data[,1])),data[,4]!="E")))
	
#convert into a matrix	
ANImat<-sapply(iMETH,function(x){
	stx=data[x,2]
	ANIx<-rbind(subset(ANI[,c(1,3)],ANI[,2]==stx),
		subset(ANI[,c(2,3)],ANI[,1]==stx))
	sapply(iMETH,function(x){
		return(mean(as.numeric(subset(ANIx[,2],ANIx[,1]==data[x,2]))))
	})
})	


COLg<-ifelse(data[,4]=="A","red",
	ifelse(data[,4]=="B","purple",
	ifelse(data[,4]=="C","green",
	ifelse(data[,4]=="D","blue","yellow"))))[iMETH]

hist(ANImat,breaks=1000)


heatmap(ANImat,scale='none',ColSideColors=COLg,
	RowSideColors=COLg,distfun = function(x) dist(x, method = "maximum"))

#Group genomes by species (use 95% for ANI threshold)

TH=95
NCORE=8

ANIspe<-subset(ANI[,c(1,2)],as.numeric(ANI[,3])>=TH)
ANIspe<-cbind(as.numeric(gsub("st","", ANIspe[,1])),
	as.numeric(gsub("st","", ANIspe[,2])))
ANIspe<-subset(ANIspe, ANIspe[,1]!=ANIspe[,2])
ANIspe<-sapply(c(1:length(ANIspe[,1])),function(x){
	paste(sort(ANIspe[x,]),collapse="_")
})

GROUPx <-unique(reduce.group(
	unique(reduce.group(
	unique(reduce.group(
	unique(reduce.group(
	unique(reduce.group(
	unique(reduce.group(
	unique(reduce.group(
	unique(reduce.group(
	unique(reduce.group(
	unique(reduce.group(
	unique(reduce.group(
	unique(reduce.group(
ANIspe,ncor=NCORE)),ncor= NCORE)),
	ncor= NCORE)),ncor= NCORE)),
	ncor= NCORE)),ncor= NCORE)),
	ncor= NCORE)),ncor= NCORE)),
	ncor= NCORE)),ncor= NCORE)),
	ncor= NCORE)),ncor= NCORE))


GROUPx <-c(GROUPx,
	setdiff(as.numeric(gsub("st","",data[iMETH,2])),
	num.strsplit.u(GROUPx)))
print(length(GROUPx))

SPE<-sapply(c(1:length(GROUPx)),function(x){
	gx<-paste("st",num.strsplit.u(GROUPx[x]),sep="")
	return(paste(setdiff(sapply(gx,function(x){
		return(subset(data[,5],data[,2]==x))	
	}),"X"),collapse="="))
})

summary(as.factor(SPE),maxsum=length(GROUPx)) #ok, no major conflict with previous classifications from Leducq et al 2022 GBE

#In case of several names, keep in priority the official name (e.g. Meth_alsaeris replaces Meth_sp_045); in case of conflict for undescribed species, keep the lower number (eg Meth_sp_042 for Meth_sp_042=Meth_sp_043); in case of no name, give a number, starting from the lower attributed in the last classification, +1

#Undescribed Species (with numbers according to Leducq et al 2022 GBE)
NUM<-sapply(c(1:length(SPE)),function(x){
	spex<-unlist(strsplit(SPE[x],split="="))
	spex<-spex[grep("Meth_sp_",spex)]
	spex<-min(c(100000,as.numeric(gsub("Meth_sp_","", spex))))	
})

#Described species
REF<-sapply(c(1:length(SPE)),function(x){
	spex<-unlist(strsplit(SPE[x],split="="))
	spex<-spex[setdiff(c(1:length(spex)),grep("Meth_sp_",spex))]
	spex<-c(spex,NA)[1]
})

#In case of conflict, choose the described name (or lower number if undescribed)
SPEnew<-ifelse(is.na(REF)==T,NUM, REF)

#Index of species without name or number
UNK<-subset(c(1:length(SPEnew)),as.numeric(SPEnew)==100000)

#New number for these species
z<-c(1:length(UNK))+max(setdiff(NUM,100000))

#Index of species with a name or number
KNO<-setdiff(c(1:length(SPEnew)), UNK)

#Format new species names
SPEnew<-c(SPEnew[KNO],z)[order(c(KNO, UNK))]

SPEnew<-paste(ifelse(nchar(SPEnew)==1,"Meth_sp_00",
	ifelse(nchar(SPEnew)==2,"Meth_sp_0","")),SPEnew,sep="")
	

#Assign species name to Methylobacterium Assemblies
SPEx<-sapply(iMETH,function(x){
	X<-paste("_",gsub("st","",data[x,2]),"_",sep="")
	SPEnew[grep(X,paste("_",GROUPx,"_",sep=""))]	
})

#Format a species map - works with ASTRAL3 but not with ASTRAL4 :/

#species

SPEn<-sort(unique(SPEx))

#tips names

MAPspe<-as.vector(sapply(SPEn,function(x){
	ix<-subset(iMETH, SPEx==x)
	#Tn<-paste(data[ix,2], data[ix,4], data[ix,5],sep="-")
	Tn<-data[ix,2]
	return(paste(c(x,paste("[",length(ix),"]",sep=""), Tn),collapse=" "))
}))

setwd(pD)

write.table(MAPspe,"Methylobacterium-ListSpecies-Sept-2025.txt",row.names=F,col.names=F,quote=F)

#######Groupe genomes assemblies for different ANI thresholds

THX=c(84,85,88,92,95)

setwd(pD)

pdf("DistANI.pdf")

hist(ANImat,breaks=1000,las=1,cex.axis=0.8,xlab="ANI")

NCORE=8

sapply(THX,function(x){
segments(x,0,x,100000,col="red",lty=2,lwd=1.5)
})

dev.off()


GROUPth<-sapply(THX,function(x){

print(c(date(),x))
TH=x


ANIspe<-subset(ANI[,c(1,2)],as.numeric(ANI[,3])>=TH)
ANIspe<-cbind(as.numeric(gsub("st","", ANIspe[,1])),
	as.numeric(gsub("st","", ANIspe[,2])))
ANIspe<-subset(ANIspe, ANIspe[,1]!=ANIspe[,2])
ANIspe<-sapply(c(1:length(ANIspe[,1])),function(x){
	paste(sort(ANIspe[x,]),collapse="_")
})

GROUPx <-unique(reduce.group(
	unique(reduce.group(
	unique(reduce.group(
	unique(reduce.group(
	unique(reduce.group(
	unique(reduce.group(
	unique(reduce.group(
	unique(reduce.group(
	unique(reduce.group(
	unique(reduce.group(
	unique(reduce.group(
	unique(reduce.group(
ANIspe,ncor=NCORE)),ncor= NCORE)),
	ncor= NCORE)),ncor= NCORE)),
	ncor= NCORE)),ncor= NCORE)),
	ncor= NCORE)),ncor= NCORE)),
	ncor= NCORE)),ncor= NCORE)),
	ncor= NCORE)),ncor= NCORE))


GROUPx <-c(GROUPx,
	setdiff(as.numeric(gsub("st","",data[iMETH,2])),
	num.strsplit.u(GROUPx)))
print(c(date(),length(GROUPx)))

return(GROUPx)

})

GROUPsum<-sapply(c(1:length(THX)),function(x){
	X=x
	grx<-GROUPth[[X]]
	z<-matrix(unlist(sapply(c(1:length(grx)),function(x){
		t(cbind(paste("G",x,sep="-"),paste("st",num.strsplit.u(grx[x]),sep="")))
	})),ncol=2,byrow=T)
	sapply(data[iMETH,2],function(x){
		return(subset(z[,1],z[,2]==x))	
	})	
})

colnames(GROUPsum)<-paste("ANI>",THX,sep="")

#Built a summary table with short genome names and species names

setwd(pD)

#For each tip label, retrieve the species names
MAPspe<-as.vector(unlist(read.delim("Methylobacterium-ListSpecies-Sept-2025.txt",header=F)))
MAPspe<-unlist(strsplit(as.character(MAPspe),split="[[]"))
MAPspe<-unlist(strsplit(as.character(MAPspe),split="[]]"))
MAPspe<-matrix(MAPspe,byrow=T,ncol=3)
MAPspe<-sapply(c(1:length(MAPspe[,1])),function(x){
	spx<-setdiff(unlist(strsplit(MAPspe[x,1],split=" ")),"")
	stx<-setdiff(unlist(strsplit(MAPspe[x,3],split=" ")),"")
	return(list(t(cbind(spx,stx))))
})
MAPspe<-matrix(unlist(MAPspe),ncol=2,byrow=T)

MAPspe<-t(sapply(data[iMETH,2],function(x){
		return(c(subset(MAPspe[,1], MAPspe[,2]==x),
			subset(data[,c(3,4)], data[,2]==x)))
		
	}))
	
colnames(MAPspe)<-c("Species","Strain","Group")	

GROUPsum<-cbind(MAPspe,GROUPsum)

setwd(pD)

write.table(GROUPsum,"Classif-Methylo-ANI-sept2025.txt")

############################################################
######STEP F - Function (KEGG analysis with KASS)

##### SECTION 2: GO analysis ######
###Retrieve protein sequence of each core gene from the reference genome (CBMB20 = st36) for GO prediction KASS (https://www.genome.jp/tools/kaas/), with default parameters + add Methylobacterium/Methylorubrum, Beijerinkia, Methylocella. If A core gene is absent from the reference genome, take the sequence from the first one in the list with this gene

library(seqinr)

setwd(pD)

STAT<-read.table("Stat-genomes.txt",header=T)
data<-read.table("data.txt",header=T)
ANNOa<-read.table("Annotation-occurence-raw.txt",header=T)
RC2<-read.table("Methylobacterium-core-gene-occ.txt",header=T)
CGstat<-read.table("Core-gene-code.txt",header=T)
CGname<-CGstat[,1]
CGtype<-CGstat[,3]
sumNEW<-read.table("Summary-per-core-gene-after-alignment.txt",header=T)

#Generate a single fasta file with protein sequences from all core genes retrieved in the reference genome REFx

REFx="st36" #CBMB20 code


	Lx<-c(1:length(CGstat[,1]))

	#Index of peg
	iP<-subset(Lx,CGstat[,3]=="peg")

	setwd(paste(pD,"core-genes",sep="/"))
	
	#Genome assembly from which each sequence can be retrieved (st36 in priority, otherwise the first one in the list)
	Nref<-sapply(CGstat[iP,1],function(x){
		#print(x)
		seqx<-read.fasta(paste(x,
			"_alignedNF.fas",sep=""))
		return(c(subset(names(seqx),names(seqx)== REFx),
			names(seqx)[1])[1])
	})
	
	SEQrefx<-sapply(c(1:length(Nref)),function(x){
		#print(x)
		seqx<-read.fasta(paste(names(Nref)[x],
			"_alignedNF.fas",sep=""))
		seqx<-unlist(seqx[[grep(Nref[x],names(seqx))[1]]])
		seqx<-getTrans(subset(seqx, seqx!="-"))
		seqx<-subset(seqx,seqx!="X")
		return(list(as.vector(unlist(seqx))))
	})

	setwd(paste(pD,"KEGG",sep="/"))

	write.fasta(SEQrefx,names= paste(names(Nref), Nref,sep="-"),
		paste("SEQaa-CoreMeth-CBMB20.fa",sep=""))

#Also, prepare a map file locating Core genes on the reference genome chromosome (return NA if not identified)

ANNOref=subset(data[,6],data[,2]== REFx)

ANNOx<-read.table(paste(pD,"RAST_files", ANNOref,sep="/"))

#All core gene are actually present in the reference genome, but some are duplicated; report all positions

MAP<-sapply(CGstat[iP,1],function(x){
	CGx=x
	NX=subset(CGstat[,2], CGstat[,1]==x)
	posx<-subset(ANNOx[,c(5,6)],ANNOx[,8]==NX)
	posx<-paste(apply(posx,1,paste,collapse="-"), collapse =";")
	return(posx)
})

setwd(pD)

write.table(cbind(Nref,MAP),"MAP-core-genes-CBMB20.txt")

#Then submit the file to KASS (https://www.genome.jp/tools/kaas/)

#The output is not user friendly - contains KO ids for each gene when available, and pathways when identified - each gene can be associated with several pathway. I retrieved raw outpup from KASS and needed to convert it manually in excell.

setwd(paste(pD,"KEGG",sep="/"))

KEGG<-read.table("KEGG.txt",header=T)


############################################################
######STEP G - for each gene alignment, calculate dN/dS (dNdS) and number of polymorphic sites per gene alignment (nPOL) that will be used to correct RF distances


setwd(pD)

STAT<-read.table("Stat-genomes.txt",header=T)
data<-read.table("data.txt",header=T)
ANNOa<-read.table("Annotation-occurence-raw.txt",header=T)
RC2<-read.table("Methylobacterium-core-gene-occ.txt",header=T)
CGstat<-read.table("Core-gene-code.txt",header=T)
CGname<-CGstat[,1]
CGtype<-CGstat[,3]
sumNEW<-read.table("Summary-per-core-gene-after-alignment.txt",header=T)
MAP<-read.table("MAP-core-genes-CBMB20.txt",header=T)

#Retrieve alignments for each PEG core gene and calcultate the number of polymorphic sites per alignments and dN/dS ratio

Lx<-c(1:length(CGstat[,1]))

#Index of peg
iP<-subset(Lx,CGstat[,3]=="peg")

setwd(paste(pD,"core-genes",sep="/"))
	

STATn<-t(sapply(CGstat[iP,1],function(x){
	print(x)
	seqx <- read.alignment(paste(x,
		"_alignedNF.fas",sep=""),
		format = "fasta")
	dNdS <- kaks(seqx)
	dNdS<-median(dNdS $ka/dNdS $ks,na.rm=T)	

	#calculate the number of polymorphic sites
	seqx<-read.fasta(paste(x,
		"_alignedNF.fas",sep=""))
	#seq names	
	stx<-names(seqx)
		
	#convert alignment in matrix
	seqx <-sapply(stx,function(x){
		return(unlist(seqx[[
			subset(c(1:length(seqx)),
			names(seqx)==x)]]))
	})
	#number of variants per sequences
	mx<-sapply(c(1:length(stx)),function(x){
		return(length(unique(seqx[,x])))
	})
	#Summary per position (type and number of variants)
	sumx<-sapply(c(1:length(seqx[,1])),function(x){
		return(list(summary(as.factor(seqx[x,]))))
	})
	sumx <-sapply(c(1:length(sumx)),function(x){
		return(paste(sort(names(sumx[[x]])),
			collapse=""))
	})	
	#remove positions with "n" (not informative)
	sumx <-sumx[setdiff(c(1:length(sumx)),
		grep("n", sumx))]
	# remove gaps 
	sumx <-subset(sumx, sumx!="-")
	# remove Ns
	sumx <-subset(sumx, sumx!="n")	
	
	#Count of the number of variable and invariable site considering gaps
	VARg<-length(subset(sumx, nchar(sumx)!=1))
	INVg<-length(subset(sumx, nchar(sumx)==1))
	#Count of the number of variable and invariable site excluding gaps
	sumx <-gsub("-","", sumx)
	VAR<-length(subset(sumx, nchar(sumx)!=1))
	INV<-length(subset(sumx, nchar(sumx)==1))	
	#alignement size
	Sx=VAR+INV
	#number of genomes
	Mx=dim(seqx)[2]
	#Number of "n" position removed
	Nx=dim(seqx)[1]-Sx
	
	return(c(x, Mx, Nx,Sx, VAR,VARg, dNdS))

}))

STATn<-as.data.frame(STATn)

colnames(STATn)<-c("Gene","Genomes","N-removed","Align-size","nVarWOgap","nVarWgap","MEDdNdS")

setwd(pD)

write.table(STATn,"CoreGenes_nPOL_dNdS.txt",row.names=F,quote=T)
