#Script to identify cluster of genes that are conserved among different genomes

###FUNCTIONS
library(parallel)
library(gplots)
library(seqinr)
library(msa)

#install.packages("colorRamps")
library(colorRamps)

ncore=8

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
	
	#ALTERNATIVE? TO TEST
	#Pu<-sort(unique(num.strsplit.u(Pclu)))
	#Pclu<-unique(sapply(Pu,function(x){
	#	paste(sort(unique(num.strsplit.u(
	#	Pclu[grep(paste("_",x,"_",sep=""),
	#	paste("_", Pclu,"_",sep=""))]))),
	#	collapse="_")
	#}))

#Path

pP<-"/myFolder/3-Phylogeny/"

pD<-"/myFolder/2-Core-genome/"

setwd(pP)

#List of core genes #This file was generated in script 3-Phylogeny.R
KEGGres<-read.table("Summary-KEGG-RF-Pos-per-PEG.txt",header=T)

Nu<-as.vector(KEGGres[,2])

#List of genomes

setwd(pD)

#taxonomy table of strains (Groups, Species and Lineages) #those files were generated in script 3-Phylogeny.R
ORD<-read.table("Taxonomy.txt",header=T)
data<-read.table("data.txt",header=T)
STAT<-read.table("Stat-genomes.txt",header=T)

#Index of Methylobacterium genomes
iMETH<-sort(intersect(subset(c(1:length(data[,1])),data[,4]!="M"),
	subset(c(1:length(data[,1])),data[,4]!="E")))

#Index of Microvirga genomes
iMICR<-subset(c(1:length(data[,1])),data[,4]=="M")

#Entovirga genomes will be took from the Hyphomicrobiales dataset

###List of Hyphomicrobiales genomes and metagenomes (remove methylobacterium, enterovirga and microvirga)

#work directory for genomes and annotation files
pG<-"/myFolder/Genomes"

setwd(pG)

data1<-read.table("data4.txt",header=T) #contain information and path for genomes and metagenome files (RAST annotation file and fasta files)

#Don't take into account Methylobacterium and Microvirga genomes (keep metagenomes)

iHYP<-sort(unique(c(intersect(setdiff(c(1:length(data1[,1])),
	grep("Microvirga",data1[,3])),
	setdiff(c(1:length(data1[,1])),
	grep("Methylobacterium",data1[,3]))),
	grep("Metagenomics",data1[,7]))))

###Combine genome lists

DATA<-rbind(
	cbind(data[iMETH,2],"Methylobacteriaceae","Methylobacterium",
		apply(ORD[sapply(data[iMETH,2],function(x){
		grep(paste("_",x,"_",sep=""),paste("_",ORD[,5],"_",sep=""))
		}),c(3,1)],1,paste,collapse="="),paste(pD,"RAST_files",data[iMETH,6],sep="/")),
	cbind(data[iMICR,2],"Methylobacteriaceae","Microvirga",
		data[iMICR,5],
		paste(pD,"RAST_files",data[iMICR,6],sep="/")),
	cbind(gsub("seq","hyp",data1[iHYP,1]),
		data1[iHYP,2],"",data1[iHYP,3],
		gsub("baptisteleducq","baptiste",data1[iHYP,6])))

DATA<-cbind(DATA[,c(1,2)],
	ifelse(DATA[,3]=="",sapply(DATA[,4],function(x){
	unlist(strsplit(x,split="_"))[1]}),DATA[,3]),DATA[,c(4,5)])
	
DATA<-cbind(DATA[,c(1,2)],
	ifelse(DATA[,2]=="Rhabdaerophilaceae","Rhabdaerophilum",DATA[,3]),DATA[,c(4,5)])

#Phylogenomic order of Hyphomicrobiales genus from SYAM paper (Leducq 2026)
setwd(pP)
HYPorder<-read.table("Hyphomicrobiales-Genus-Order.txt",header=T)



#STEP A - BCHL gene occurence in Hyphomicrobiales

###Retrieve Methylobacterium core gene abundance in each genome

COLfam<-rbind(cbind(grep("Methylobacterium",DATA[,3]),"red"),
	cbind(grep("Enterovirga",DATA[,3]),"pink"),
	cbind(grep("Lichenihabitantaceae",DATA[,2]),"green"),
	cbind(grep("RH_AL",DATA[,2]),"blue"),
	cbind(grep("JAJXWB01",DATA[,2]),"yellow3"))

COLfam<-rbind(COLfam,
	cbind(setdiff(c(1:length(DATA[,1])),as.numeric(COLfam[,1])),"grey"))
COLfam<-COLfam[order(as.numeric(COLfam[,1])),2]

COLrf<-KEGGres[,5]-min(KEGGres[,5])
COLrf<-COLrf/max(COLrf)
COLrf<-rgb(COLrf,0, 1-COLrf)

#Gene abundance per genome
ANNOcore<-sapply(c(1:length(DATA[,1])),function(x){
	print(x)
	ANNOx<-read.table(DATA[x,5])
	sapply(KEGGres[,2],function(x){
		return(paste(nchar(subset(
			as.vector(ANNOx[,12]),
			as.vector(ANNOx[,8])== x)),
			collapse="_"))
	})
})

#Gene occurence per genome
ANNOcoreB<-sapply(c(1:length(ANNOcore[,1])),function(x){
	as.numeric(summary(strsplit(ANNOcore[x,],split="_"))[,1])
})


#################################
####Cluster enriched in outlier genes (higher residual RF distance) 

setwd(pP)

####List of genes
#Bachteriochlorophyll  and carotenoid biosynthesis
BCHL<-read.table("List-BCHL-genes_extended.txt",header=T)

pdf("HGTlikeliness-distribution.pdf",height=3.5,width=2.5)

hist(KEGGres[,5],breaks=25,las=1,
	xlab="HGT likelihood",
	ylab="Core genes (n = 1,099)",
	cex.axis=0.8,main="",border="grey",
	xlim=c(min(KEGGres[,5]),max(KEGGres[,5])),
	ylim=c(0,180))

RESbchl<-unlist(sapply(BCHL[,3],function(x){
	return(subset(KEGGres[,5],KEGGres[,1]==x))
}))

par(new=T)

hist(RESbchl,breaks=15,las=1,xlab="",
	ylab="",xaxt='n',yaxt="n",main="",
	border="green4",col="green3",
	xlim=c(min(KEGGres[,5]),max(KEGGres[,5])),
	ylim=c(0,180))
	
	

dev.off()

#Average occurence of each gene in each genus


BCHLhypOcc<-sapply(HYPorder[,2],function(x){
	print(x)	
	genx<-subset(DATA[,5],paste(DATA[,2],DATA[,3],sep="_")==x)
	as.vector(apply(ifelse(sapply(c(1:length(genx)),function(x){
		ANNOx<-read.table(genx[x])
		sapply(BCHL[,2],function(x){
			return(length(subset(
				as.vector(ANNOx[,12]),
				as.vector(ANNOx[,8])== x)))
		})
	
	})==0,0,1),1,mean))
})

rownames(BCHLhypOcc)<-BCHL[,1]

setwd(pP)

pdf("BCHLhypOcc.pdf")
heatmap(t(BCHLhypOcc),scale="none",Rowv=NA,
	margin=c(6,6),cexRow=0.6,cexCol=0.6)
dev.off()	
	
write.table(cbind(BCHL ,round(BCHLhypOcc,2)),"BCHLhypOcc.txt")	
	
##STEP B

BCHL1<-c("CG_798","CG_800","CG_747","CG_223","CG_400","CG_1041","CG_1083","CG_599","CG_1084","CG_24","CG_166","CG_515","CG_744","CG_840","CG_514","CG_516","CG_600","CG_799","CG_804")

BCHL2<-c("CG_224","CG_225","CG_401","CG_512","CG_513","CG_591","CG_743","CG_745","CG_746","CG_16","CG_88","CG_226","CG_478")

BCHLout<-setdiff(intersect(KEGGres[,1],BCHL[,3]),c(BCHL1, BCHL2))

# Osmoregulation (Sodium-Potassium Pump) 
OSMOK<-c("CG_677","CG_758","CG_759","CG_760")
#Tripartit ATP-independant perplasmic transposter
TRAP<-c("CG_990","CG_991","CG_992")
#One carbon (clusters containing outliers)
ONEC1<-c("CG_576","CG_586")
ONEC2<-c("CG_172","CG_996","CG_62","CG_629","CG_754","CG_391","CG_393","CG_392")
ONEC3<-c("CG_248","CG_336")
ONEC4<-c("CG_830","CG_831","CG_209","CG_652","CG_332","CG_21")
ONEC5<-c("CG_60","CG_61","CG_980")
#Iron-sulfur cluster assembly
SUF<-c("CG_158","CG_499","CG_257","CG_475","CG_498","CG_500","CG_501","CG_684")
#alkanesulfonate
ALKASU<-c("CG_387","CG_111","CG_112")
#Nicotinamidase
NICO<-c("CG_117","CG_688","CG_860","CG_653")
#Carbon starvation
STARV<-c("CG_197","CG_1091")
#NADH-ubiquinone oxidoreductase 
NADH<-c("CG_777","CG_523","CG_524","CG_851","CG_593","CG_643","CG_644","CG_645","CG_646","CG_647","CG_648","CG_649")
#Other outliers (not in clusters)
OUT<-c("CG_335","CG_446","CG_674","CG_807","CG_1099","CG_983","CG_468","CG_189","CG_264","CG_436")

LIST<-rbind(cbind("BCHL1", BCHL1),
	cbind("BCHL2", BCHL2),
	cbind("BCHLout", BCHLout),
	cbind("OSMOK", OSMOK),
	cbind("TRAP", TRAP),
	cbind("ONEC1", ONEC1),
	cbind("ONEC2", ONEC2),
	cbind("ONEC3", ONEC3),
	cbind("ONEC4", ONEC4),
	cbind("ONEC5", ONEC5),
	cbind("SUF", SUF),
	cbind("NICO", NICO),
	cbind("STARV", STARV),
	cbind("NADH", NADH),
	cbind(paste("OUT",c(1:length(OUT)),sep=""), OUT),
	cbind("ALKASU", ALKASU))


LIST <-cbind(LIST ,sapply(LIST[,2],function(x){
	return(subset(KEGGres[,2], KEGGres[,1]==x))
}))

colnames(LIST)<-c("CLUSTER","CG","RAST")

#iCRT<-grep("crt",BCHL[,1])
#CRT<-cbind("CRT", BCHL[iCRT,3], BCHL[iCRT,2])

#iBCHL<-setdiff(c(1:length(BCHL[,1])), iCRT)
#BCHL<-cbind("BCHL", BCHL[iBCHL,3], BCHL[iBCHL,2])

#BCHL<-cbind("BCHL", BCHL[,3], BCHL[,2])

#colnames(BCHL)<-colnames(LIST)
#colnames(CRT)<-colnames(LIST)

#LIST<-rbind(BCHL, CRT, LIST)
#LIST<-rbind(BCHL, LIST)

LIST<-cbind(LIST,sapply(LIST[,2],function(x){
	return(c(subset(KEGGres[,5],KEGGres[,1]==x),NA)[1])
}))

colnames(LIST)<-c(colnames(LIST)[1:3],"RF.residual")

setwd(pP)

write.table(LIST,"Outliers-for-Phylo-Hyphomicrobiales.txt",quote=T)



###################################################################################################
### STEP C
###bchG genomic context (to compare with its phylogeny) "Chlorophyll a synthase ChlG (EC 2.5.1.62)"

LIST<-read.table("Outliers-for-Phylo-Hyphomicrobiales.txt",header=T)

##For BCHL cluster, pool together genes that are located within 25KB
b1<-grep("BCHL1", LIST[,1])
b2<-grep("BCHL2", LIST[,1])


#occurence of BCHL genes per genome

oBCHL<-sapply(c(1:length(DATA[,1])),function(x){
	print(x)
	ANNOx<-read.table(DATA[x,5])
		
	p1<-sapply(b1,function(x){
		paste(subset(ANNOx[,4],ANNOx[,8]==LIST[x,3]),collapse="=")
	})
	p2<-sapply(b2,function(x){
		paste(subset(ANNOx[,4],ANNOx[,8]==LIST[x,3]),collapse="=")
	})
	
	return(sapply(c(p1,p2),function(x){
		ifelse(x=="",0,length(grep("=",x))+1)
		
	}))
})
rownames(oBCHL)<-LIST[c(b1,b2),6]

setwd(pP)

write.table(cbind(DATA[,c(1:4)],t(oBCHL)),"BCHL-gene-occurence-per-genome.txt")

#Average copy number per gene per genus

aBCHL<-sapply(unique(paste(DATA[,2],DATA[,3])),function(x){
	z<-subset(t(oBCHL),paste(DATA[,2],DATA[,3])==x)
	return(apply(z,2,mean))
	
})
rownames(aBCHL)<-LIST[c(b1,b2),6]

#number of BCHL genes per genome

nBCHL<-apply(ifelse(oBCHL==0,0,1),2,sum)


plot(as.factor(paste(DATA[,2],DATA[,3])), nBCHL,las=2,cex.axis=0.4)

#Estimate of the number of clusters per genus 

CLUgen<-t(sapply(unique(paste(DATA[,2],DATA[,3])),function(x){
	z<-subset(nBCHL,paste(DATA[,2],DATA[,3])==x)
	return(c(length(subset(z,z<10)),length(subset(z,z>=10))))
	
}))

iBCHL<-subset(c(1:length(nBCHL)),nBCHL>1)

#BCHL genes cluster architecture in genomes with at least 2 BCHL genes
posBCHL<-t(sapply(iBCHL,function(x){
	print(x)
	ANNOx<-read.table(DATA[x,5],header=T)
	
	
	p1<-sapply(b1,function(x){
		paste(subset(ANNOx[,4],ANNOx[,8]==LIST[x,3]),collapse="=")
	})
	p2<-sapply(b2,function(x){
		paste(subset(ANNOx[,4],ANNOx[,8]==LIST[x,3]),collapse="=")
	})
	
	#Pairwise minimum distance between each gene (return NA if on different scaffolds)
	bx<-subset(c(b1,b2),c(p1,p2)!="")
	
	distx<-unlist(sapply(c(1:(length(bx)-1)),function(x){
		#print(x)
		X=x
		X1=bx[x]
		S1=subset(ANNOx[,1],ANNOx[,8]==LIST[X1,3])
		s1=subset(ANNOx[,5],ANNOx[,8]==LIST[X1,3])
		e1=subset(ANNOx[,6],ANNOx[,8]==LIST[X1,3])
		sapply(c((X+1):(length(bx))),function(x){
			X2=bx[x]
			S2=subset(ANNOx[,1],ANNOx[,8]==LIST[X2,3])
			s2=subset(ANNOx[,5],ANNOx[,8]==LIST[X2,3])
			e2=subset(ANNOx[,6],ANNOx[,8]==LIST[X2,3])
			#print(x)
			min(sapply(c(1:length(S1)),function(x){
				Sx=S1[x]
				sx=s1[x]
				ex=e1[x]
				sapply(c(1:length(S2)),function(x){
					Sy=S2[x]
					sy=s2[x]
					ey=e2[x]
					z<-min(c(abs(sx-sy),abs(sx-ey),abs(ex-sy),abs(ex-ey)))
					paste(X1,X2,ifelse(Sx==Sy,z,NA),sep="_")
				})
			}),na.rm=T)	
		})
	}))
	distx<-matrix(as.numeric(unlist(strsplit(distx,split="_"))),ncol=3,byrow=T)
	distx<-subset(distx,is.na(distx[,3])==F)
	distx<-subset(distx,distx[,3]<=25000)
	distx<-paste(distx[,1], distx[,2],sep="_")
	distx <-unique(reduce.group(	unique(reduce.group(	unique(reduce.group(unique(reduce.group(
		unique(reduce.group(unique(reduce.group(unique(reduce.group(unique(reduce.group(
		unique(reduce.group(unique(reduce.group(unique(reduce.group(unique(reduce.group(
		distx,ncor=NCORE)),ncor= NCORE)),ncor= NCORE)),ncor= NCORE)),
		ncor= NCORE)),ncor= NCORE)),ncor= NCORE)),ncor= NCORE)),
		ncor= NCORE)),ncor= NCORE)),ncor= NCORE)),ncor= NCORE))	
	distx <-c(distx,setdiff(bx ,as.numeric(num.strsplit.u(distx))))
	distx <-distx[order(as.numeric(summary(strsplit(distx,split="_"))[,1]),decreasing=T)]
	
	sapply(c(b1,b2),function(x){
		c(grep(paste("_",x,"_",sep=""),paste("_", distx,"_",sep="")),NA)[1]
	})
	
	
}))


colnames(posBCHL)<-paste(LIST[c(b1,b2),1],LIST[c(b1,b2),2],sep="=")

posBCHL<-cbind(DATA[iBCHL,c(1:4)], posBCHL)

setwd(pP)

write.table(posBCHL,"BCHL-arch-Hyphomicr.txt")

########


###Retrieve outlier gene size in each genomes

ANNOa<-sapply(c(1:length(DATA[,1])),function(x){
	print(x)
	ANNOx<-read.table(DATA[x,5])
	sapply(LIST[,3],function(x){
		return(paste(nchar(subset(
			as.vector(ANNOx[,12]),
			as.vector(ANNOx[,8])== x)),
			collapse="_"))
	})
})

# STEP D

###Retrieve sequences in all genomes


#size of each gene copy normalized by median size accross genomes 
SizeM<-sapply(c(1:length(ANNOa[,1])),function(x){
	sx<-ANNOa[x,]
	#median size
	mx<-median(as.numeric(
		unlist(strsplit(sx,split="_"))))
		#unlist(strsplit(sx[iC],split="_"))))
	#normalized copy size
	sapply(c(1:length(sx)),function(x){
		return(paste(round(
			as.numeric(unlist(strsplit(
				sx[x],split="_")))/mx,2),collapse="_"
		))
	})			
})



#real copy number per gene
RC<-sapply(c(1:length(SizeM[1,])),function(x){
	sx<-SizeM[,x]
	return(sapply(c(1:length(sx)),function(x){
		return(sum(as.numeric(unlist(
			strsplit(sx[x],split="_")))))
	}))
})

setwd(pP)

pdf("Figure-Outliers-Genes-Hypomicrobiales.pdf")

hist(as.numeric(unlist(strsplit(SizeM,split="_"))),breaks=1000, main="relavite size")

hist(RC,breaks=1000, main="real copy number")

plot(as.factor(LIST[,1]),as.numeric(LIST[,4]),las=2,cex.axis=0.8,ylab="Residual RF",xlab="")


heatmap(t(RC),scale='none',labRow=paste(LIST[,2],LIST[,1]),cexRow=0.4,labCol=paste(DATA[,2], DATA[,3]),margin=c(10,3),cexCol=0.4,ColSideColors=COLfam)

heatmap(t(log10(RC+0.1)),scale='none',labRow=paste(LIST[,2],LIST[,1]),cexRow=0.4,labCol=paste(DATA[,2], DATA[,3]),margin=c(10,3),cexCol=0.4,ColSideColors=COLfam)

#Average copy number per genera per gene

RCgen<-t(sapply(HYPorder[order(HYPorder[,1],decreasing=T),2],function(x){
	
	apply(subset(RC,paste(DATA[,2],
		DATA[,3],sep="_")==x),2,mean)
}))

heatmap((RCgen),scale='none',labRow=HYPorder[order(HYPorder[,1],decreasing=T),2],cexRow=0.6,
	Rowv=NA,labCol=paste(LIST[,2],LIST[,1]),margin=c(10,10),cexCol=0.3,Colv=NA)

#Average copy number per genera per cluster

RCclu<-sapply(unique(LIST[,1]),function(x){
	
	apply(t(subset(t(RCgen),LIST[,1]==x)),1,mean)
})

heatmap(RCclu,scale='none',cexRow=0.6,
	Rowv=NA,margin=c(10,10),cexCol=1,Colv=NA)

dev.off()

#####For each genome and each gene, only keep copies that have relative size within the range 0.7-1.3 and write a fasta file for each genome

SUMall<-sapply(c(1:length(DATA[,1])),function(x){
	print(x)
	X=x
	SUMx<-sapply(c(1:length(LIST[,1])),function(x){
		sizex<-c(as.numeric(unlist(strsplit(SizeM[X,x],split="_"))),NA) #relative size
		annoax<-c(as.numeric(unlist(strsplit(ANNOa[x,X],split="_"))),NA) #real size
		ix<-order(annoax,decreasing=T) #index for copy 
		return(list(t(cbind(x,LIST[x,1],LIST[x,2], ix,annoax, sizex))))
	})
	SUMx<-matrix(unlist(SUMx),ncol=6,byrow=T)
	SUMx<-subset(SUMx,is.na(SUMx[,5])==F)

	ANNOx<-read.table(DATA[X,5])
	seqx<-unlist(sapply(unique(LIST[as.numeric(SUMx[,1]),3]),function(x){
		list(return(subset(
			as.vector(ANNOx[,12]),
			as.vector(ANNOx[,8])== x)))
	}))	
	print(c("are all sequences the expected length?", unique(as.numeric(SUMx[,5])==nchar(seqx))))
	#Remove gene copies outside the range 0.7-1.3
	seqx<-subset(seqx,as.numeric(SUMx[,6])>=0.7)
	SUMx<-subset(SUMx,as.numeric(SUMx[,6])>=0.7)
	seqx <-subset(seqx,as.numeric(SUMx[,6])<=1.3)
	SUMx<-subset(SUMx,as.numeric(SUMx[,6])<=1.3)	
	
	#sequence names (cluster, core gene code, gene copy number)
	nx<-paste(SUMx[,2],SUMx[,3],SUMx[,4],sep="-")
	
	#file name
	fx=paste(paste(DATA[X,c(1,2,4)],collapse="-"),".fas",sep="")
	
	seqx<-sapply(c(1:length(seqx)),function(x){
		return(list(unlist(strsplit(seqx[x],split=""))))
	})
	
	setwd(paste(pP,"Fasta-Outliers-Hyphomicrobiales",sep="/"))
	write.fasta(seqx,names= nx,file.out=fx)
	
	return(list(t(cbind(fx, SUMx))))
	
})

SUMall<-matrix(unlist(SUMall),ncol=7,byrow=T)

###Produce one sequence file per gene

setwd(paste(pP,"Fasta-Outliers-Hyphomicrobiales",sep="/"))


sapply(unique(LIST[,2]),function(x){
	print(x)
	X=x
	#List of genomes where sequences for this gene are available
	gx<-subset(SUMall[,1],SUMall[,4]==X)
	
	#names of sequences from this genes in those files
	nx<-subset(paste(SUMall[,3], SUMall[,4], SUMall[,5],sep="-"),SUMall[,4]==X)
	
	
	seqx<-sapply(c(1:length(gx)),function(x){
		#print(x)
		sx<-read.fasta(gx[x])
		return(list(unlist(
			sx[[subset(c(1:length(sx)),
		names(sx)==nx[x])]])))	
	})
	
	#sequence names
	NX<-subset(paste(gsub(".fas","",SUMall[,1]), SUMall[,5],sep="-"),SUMall[,4]==X)
	#file names
	FX<-paste(unique(subset(paste(SUMall[,3], SUMall[,4],sep="-"),SUMall[,4]==X)),".fas",sep="")
	
	print(c("are all sequences the expected length?",unique(as.numeric(summary(seqx)[,1])==as.numeric(subset(SUMall[,6],SUMall[,4]==X)))))
	
	
	write.fasta(seqx,names= NX,file.out=FX)
	
})

#Align each core gene based on amino acid sequence (and convert back to nucleotide)

mclapply(c(1:length(LIST[,1])),function(x){
	Ix<-paste(LIST[x,1],LIST[x,2],sep="-") #sequence file name
	print(Ix)
	
	setwd(paste(pP,"Fasta-Outliers-Hyphomicrobiales",sep="/"))
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

mclapply(c(1:length(LIST[,1])),function(x){
	Ix<-paste(LIST[x,1],LIST[x,2],sep="-") #sequence file name
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
	
},mc.cores=4)



#Check that gap replacement by Ns did not affected alignments


#only "n" should be returned

unique(unlist(
	sapply(c(1:length(LIST[,1])),function(x){
	Ix<-paste(LIST[x,1],LIST[x,2],sep="-") #sequence file name
	print(Ix)
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



#Proportion of missing data 

PPX<-unlist(sapply(c(1:length(LIST[,1])),function(x){
	
	Ix<-paste(LIST[x,1],LIST[x,2],sep="-") #sequence file name
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
#alkx
#Ony "n" should be returned


sumNEW<-mclapply(c(1:length(LIST[,1])),function(x){
	
	Ix<-paste(LIST[x,1],LIST[x,2],sep="-") #sequence file name
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

#remove intermediate files

sapply(c(1:length(LIST[,1])),function(x){
	
	Ix<-paste(LIST[x,1],LIST[x,2],sep="-") #sequence file name
	print(Ix)
	file.remove(paste(Ix,
		c("_aligned.fas","_alignedN.fas"),sep=""))
})	


sumNEW <-matrix(unlist(sumNEW),ncol=5,byrow=T)

colnames(sumNEW)<-c("Name1","Name2","Sequences","Size1(bp)","Size2(bp)")
sumNEW<-cbind(LIST, sumNEW)
write.table(sumNEW,"Summary-per-outlier-core-gene-Hyphomycrobiales.txt",row.names=F)



#Rename each alignment to be easily submitted for analysis with an array

sapply(c(1:length(LIST[,1])),function(x){
	
	Ix<-paste(LIST[x,1],LIST[x,2],sep="-") #sequence file name
	print(Ix)
	seqx<-read.fasta(paste(Ix,"_alignedNF.fas",sep=""))
	write.fasta(seqx,names=names(seqx),file.out=paste("CG-Hypho_",x,".fas",sep=""))
	
})	


###Run each core gene phylogeny separately with RAxML


#!/bin/bash
#SBATCH --mem=16GB
#SBATCH -N 1
#SBATCH --time=2-00:00:00
#SBATCH -p batch_72h
#SBATCH --ntasks-per-node=4
#SBATCH --array=1-109%8

###		module load StdEnv/2020  gcc/9.3.0 
###		module load raxml/8.2.12

###		raxmlHPC-MPI -T 4 -f a -s CG-Hypho_${SLURM_ARRAY_TASK_ID}.fas -x 12345 -p 12345 -# 1000 -m GTRGAMMA -n CG-Hypho_${SLURM_ARRAY_TASK_ID}-GTRGAMMA-1000

## STEP E
###Convert and draw RAxML output trees

pP<-"/myFolder/3-Phylogeny/"

setwd(paste(pP,"Fasta-Outliers-Hyphomicrobiales",sep="/"))
LIST<-read.table("Summary-per-outlier-core-gene-Hyphomycrobiales.txt",header=T)


library(ape)
library(ips)

#Color codes for Lineage
COLlin<-t(cbind(
	c("A1",rgb(1,0.6,0.6)),
	c("A2",rgb(1,0,0)),
	c("A3",rgb(0.7,0,0)),
	c("A4",rgb(0.35,0,0)),
	c("A5",rgb(0.5,0.25,0)),
	c("A6",rgb(0.75,0.375,0)),
	c("A7",rgb(1,0.5,0)),
	
	c("B1",rgb(1,0.6,1)),
	c("B2",rgb(1,0,1)),
	c("B3",rgb(0.7,0,0.7)),
	c("B4",rgb(0.4,0,0.4)),
	
	c("D1",rgb(0.7,0.7,1)),
	c("D2",rgb(0.35,0.35,1)),
	c("D3",rgb(0,0,1)),
	c("D4",rgb(0,0,0.5)),
	c("D5",rgb(0,0.4,0.75)),
	c("D6",rgb(0,0.8,1)),
	
	c("C1",rgb(0.7,1,0.7)),
	c("C2",rgb(0.35,1,0.35)),
	c("C3",rgb(0,0.7,0)),
	c("C4",rgb(0,0.35,0)),
	c("C5",rgb(0.4,0.5,0)),
	c("C6",rgb(0.6,0.75,0)),
	c("C7",rgb(0.8,1,0))
))


###For some genes, long branched cause by missalignments due to too short sequences in some genomes: problematic sequences were manually removed and phylogeny redone for the following genes: 


###### "CG_799" (bchH)
grep("CG_799",LIST[,2])

#RAxML_bipartitionsBranchLabels.BCHL1-CG_799_alignedNF-manual-GTRGAMMA-1000 replaces RAxML_bipartitionsBranchLabels.CG-Hypho_18-GTRGAMMA-1000 


pdf("Plot-Phylo-Outliers-Hyphomicrobiales.pdf",width=12,height=12)

sapply(c(1: length(LIST[,1])),function(x){
	print(x)
	#read tree as character
	treex<-as.vector(read.table(paste(
		"RAxML_bipartitionsBranchLabels.CG-Hypho_",x,"-GTRGAMMA-1000",sep="")))
	#split nodes
	treex<-unlist(strsplit(
		as.character(treex),split=":"))
	
	#identify nodes with brackets (node supports)
	bi<-grep("[[]",treex)
	bo<-setdiff(c(1:length(treex)), bi)
	#extract node support and convert in readable tree format: #BranchLength[NodeSupport])# converted in #NodeSupport:BranchLength)# 
	ns<-sapply(bi,function(x){
		nx<-unlist(strsplit(treex[x],split="[[]"))
		nx<-unlist(strsplit(nx,split="[]]"))
		return(paste(nx[2],":",nx[1],nx[3],sep=""))
	})
	treex<-c(ns, treex[bo])[order(c(bi,bo))]
	
	#if last sign is not ")" or ";", add ":"
	bi<-grep("[)]", treex)
	bo<-setdiff(c(1:length(treex)), bi)
	
	ns<-paste(treex[bo],":",sep="")
	
	treex<-c(ns, treex[bi])[order(c(bo,bi))]
	
	treex<-paste(treex,collapse="")
	write.table(treex,
		paste("RAxML_tree-",LIST[x,1],"-",LIST[x,2],"_reformated",sep=""),
		row.names=F,col.names=F,quote=F)
		
	treex<-read.tree(paste("RAxML_tree-",LIST[x,1],"-",LIST[x,2],"_reformated",sep=""))	

	#Color for tips (according to Methylo lineage; grey if other Hyphomicrobiales)
	COLx<-matrix(unlist(sapply(COLlin[,1],function(x){
		t(cbind(x,c(NA,grep(x,treex $tip))))
	})),ncol=2,byrow=T)
	COLx<-rbind(COLx,t(sapply(setdiff(c(1:length(treex $tip)),as.numeric(COLx[,2])),function(x){
		return(c(paste(unlist(strsplit(treex $tip[x],split="-"))[2:3],collapse="-"),x))
	})))
	
	COLx<-subset(COLx,is.na(COLx[,2])==F)
	COLx<-COLx[order(as.numeric(COLx[,2])),1]
	
	COLx<-cbind(COLx,sapply(COLx,function(x){
		return(c(subset(COLlin[,2], COLlin [,1]==x),"yellow2")[1])
	}))
	COLx<-cbind(COLx[,1],ifelse(COLx[,1]=="Methylobacteriaceae-Methylobacterium_sp","grey",COLx[,2]))
	
	#Root on other hyphomicrobiales
	OUT<-treex$tip[grep("yellow2",COLx[,2])]

	#plot.phylo(root(treex,OUT),type="phylogram",cex=0.3,main=paste(LIST[x,1],"\n",LIST[x,2],"\n",
	#	LIST[x,3],"\n resRF=",round(as.numeric(LIST[x,4]),3)),cex.main=0.8,
	#	tip.color= COLx[,2],lab4ut="axial",font=2, 	
	#	open.angle=20,
	#	node.color=rgb(0,0,0),
	#	edge.color=rgb(0,0,0))	
	
	plot.phylo(unroot(treex),type="unroot",cex=0.3,main=paste(LIST[x,1],"\n",LIST[x,2],"\n",
		LIST[x,3],"\n resRF=",round(as.numeric(LIST[x,4]),3)),cex.main=0.8,
		tip.color= COLx[,2],lab4ut="axial",font=2, 	
		open.angle=20,
		node.color=rgb(0,0,0),
		edge.color=rgb(0,0,0))	
	
})

dev.off()



######################################################################################################################################################################################################

##################STEP F- BCHL cluster architecture in other Hyphomicrobiales

pP<-"/myFolder/3-Phylogeny/"

setwd(paste(pP,"Fasta-Outliers-Hyphomicrobiales",sep="/"))
LIST<-read.table("Summary-per-outlier-core-gene-Hyphomycrobiales.txt",header=T)

setwd(pP)

##For BCHL cluster, pool together genes that are located within 25KB
b1<-grep("BCHL1", LIST[,1])
b2<-grep("BCHL2", LIST[,1])

#number of BCHL genes per genome

nBCHL<-sapply(c(1:length(DATA[,1])),function(x){
	print(x)
	ANNOx<-read.table(DATA[x,5])
		
	p1<-sapply(b1,function(x){
		paste(subset(ANNOx[,4],ANNOx[,8]==LIST[x,3]),collapse="=")
	})
	p2<-sapply(b2,function(x){
		paste(subset(ANNOx[,4],ANNOx[,8]==LIST[x,3]),collapse="=")
	})
	
	return(length(subset(c(p1,p2),c(p1,p2)!="")))
})


sapply(unique(paste(DATA[,2],DATA[,3])),function(x){
	return(subset(nBCHL,paste(DATA[,2],DATA[,3])==x ))
})

#Extract clusters for genomes with at least 25 BCHL genes, excluding Methylobacterium
iBCHL<-setdiff(subset(c(1:length(nBCHL)),nBCHL>=25),
 grep("Methylobacterium",DATA[,3]))

summary(as.factor(subset(paste(DATA[,2],DATA[,3]), nBCHL>=25)))

#As a reference, ad four Methylobacterium genomes (one from each group) : A1=Meth_radiotolerans=st68-BCHL.gbk, B1=Meth_extorquens=st49-BCHL.gbk, D1=Meth_sp_012=st1-BCHL.gbk, C1=Meth_alsaeris=st314-BCHL.gbk

iBCHL2 <-c(sapply(c("st68","st49","st1","st314"),function(x){
	return(subset(c(1:length(DATA[,1])),DATA[,1]==x))
}), iBCHL)


#Genome path for each genome

GP<-c(sapply(c("st68","st49","st1","st314"),function(x){
	return(paste(pD,"Genome_files/",subset(data[,7],data[,2]==x),sep=""))
}), sapply(iBCHL,function(x){
	stx<-gsub("hyp","seq",DATA[x,1])
	ix<-subset(c(1:length(data1[,1])),data1[,1]== stx)

	paste(gsub("leducq","",data1[ix,9]),data1[ix,5],sep="")
}))

NCORE=8

#BCHL genes cluster architecture 
POSbchl<-t(sapply(iBCHL2,function(x){
	print(x)
	ANNOx<-read.table(DATA[x,5],header=T)
	
	
	p1<-sapply(b1,function(x){
		c(apply(subset(ANNOx[,c(1,5,6)],ANNOx[,8]==LIST[x,3]),1,paste,collapse="="),"")[1]
	})
	p2<-sapply(b2,function(x){
		c(apply(subset(ANNOx[,c(1,5,6)],ANNOx[,8]==LIST[x,3]),1,paste,collapse="="),"")[1]
	})
	
	#Pairwise minimum distance between each gene (return NA if on different scaffolds)
	bx<-subset(c(b1,b2),c(p1,p2)!="")
	
	distx<-unlist(sapply(c(1:(length(bx)-1)),function(x){
		#print(x)
		X=x
		X1=bx[x]
		S1=subset(ANNOx[,1],ANNOx[,8]==LIST[X1,3])
		s1=subset(ANNOx[,5],ANNOx[,8]==LIST[X1,3])
		e1=subset(ANNOx[,6],ANNOx[,8]==LIST[X1,3])
		sapply(c((X+1):(length(bx))),function(x){
			X2=bx[x]
			S2=subset(ANNOx[,1],ANNOx[,8]==LIST[X2,3])
			s2=subset(ANNOx[,5],ANNOx[,8]==LIST[X2,3])
			e2=subset(ANNOx[,6],ANNOx[,8]==LIST[X2,3])
			#print(x)
			min(sapply(c(1:length(S1)),function(x){
				Sx=S1[x]
				sx=s1[x]
				ex=e1[x]
				sapply(c(1:length(S2)),function(x){
					Sy=S2[x]
					sy=s2[x]
					ey=e2[x]
					z<-min(c(abs(sx-sy),abs(sx-ey),abs(ex-sy),abs(ex-ey)))
					paste(X1,X2,ifelse(Sx==Sy,z,NA),sep="_")
				})
			}),na.rm=T)	
		})
	}))
	distx<-matrix(as.numeric(unlist(strsplit(distx,split="_"))),ncol=3,byrow=T)
	distx<-subset(distx,is.na(distx[,3])==F)
	distx<-subset(distx,distx[,3]<=25000)
	distx<-paste(distx[,1], distx[,2],sep="_")
	distx <-unique(reduce.group(	unique(reduce.group(	unique(reduce.group(unique(reduce.group(
		unique(reduce.group(unique(reduce.group(unique(reduce.group(unique(reduce.group(
		unique(reduce.group(unique(reduce.group(unique(reduce.group(unique(reduce.group(
		distx,ncor=NCORE)),ncor= NCORE)),ncor= NCORE)),ncor= NCORE)),
		ncor= NCORE)),ncor= NCORE)),ncor= NCORE)),ncor= NCORE)),
		ncor= NCORE)),ncor= NCORE)),ncor= NCORE)),ncor= NCORE))	
	distx <-c(distx,setdiff(bx ,as.numeric(num.strsplit.u(distx))))
	
	scax<-c(p1,p2)
	scax<-subset(scax, scax!="")
	scax<-matrix(unlist(strsplit(scax,split="=")),ncol=3,byrow=T)
	
	return(c(DATA[x,c(1:4)],
		length(unique(scax[,1])), #number of scaffolds (could be higher than the number of clusters if duplication)
		length(bx),    #number of genes
		length(distx),   #number of clusters
		paste(distx,collapse="===")))
	
	
}))

colnames(POSbchl)<-c("Strain","Family","Genus","Species","Scaffolds","Genes","Cluster","Structure")


#For the following genomes, scafold names in genome files do not match scaffold name in RAST (probably renammed during RAST analysis) ; the last number in RAST scaffold name (SPECIES_NAME_1) however corresponds to the index of the scaffold in genome file: need to extract this number for those particular genomes

iBCHL3<-c(5,6,11,13,14,15,16,18,19,20)
#iBCHL3, first Scaffold name in RAST file
#5, "Methylocapsa_palsarum_1"
#6, "Methylocella_silvestris_1"
#11, "Methylocystis_echinoides_1"
#13, "RH_AL1_1"
#14, "RH_AL8_1"
#15, "Rhodoblastus_acidophilus_1"
#16, "Rhodoblastus_sphagnicola_1"
#18, "Salinarimonas_ramus_1"
#19, "Salinarimonas_rosea_1"
#20, "Salinarimonas_soli_1"

####For genomes with scaffold names matching between genome and RAST files, retrieve complete genomic sequence from the original genomes file

SEQhyp2<-sapply(setdiff(c(1:length(iBCHL2)), iBCHL3),function(x){
	print(x)
	X=iBCHL2[x]
	POSx<-POSbchl[x,]
	#annotation file
	ANNOx<-read.table(DATA[X,5],header=T)

	#genome file
	SEQ<-read.fasta(GP[x])
	
	#BCHL gene position
	px<-unlist(sapply(c(b1,b2),function(x){
		list(apply(subset(ANNOx[,c(1,5,6)],ANNOx[,8]==LIST[x,3]),1,paste,collapse="="))
	}))
	px<-matrix(unlist(strsplit(px,split="=")),ncol=3,byrow=T)
	
	#Pairwise minimum distance between each gene (return NA if on different scaffolds)
	
	distx<-unlist(sapply(c(1:(length(px[,1])-1)),function(x){
		#print(x)
		X1=x
		S1= px[X1,1]
		s1=as.numeric(px[X1,2])
		e1=as.numeric(px[X1,3])
		sapply(c((X1+1):(length(px[,1]))),function(x){
			X2=x
			S2=px[X2,1]
			s2=as.numeric(px[X2,2])
			e2=as.numeric(px[X2,3])
			#print(x)
			z<-min(c(abs(s1-s2),abs(e1-e2),abs(s1-e2),abs(e1-s2)))
			paste(X1,X2,ifelse(S1==S2,z,NA),sep="_")	
		})
	}))
	distx<-matrix(as.numeric(unlist(strsplit(distx,split="_"))),ncol=3,byrow=T)
	distx<-subset(distx,is.na(distx[,3])==F)
	distx<-subset(distx,distx[,3]<=25000)
	distx<-paste(distx[,1], distx[,2],sep="_")
	distx <-unique(reduce.group(	unique(reduce.group(	unique(reduce.group(unique(reduce.group(
		unique(reduce.group(unique(reduce.group(unique(reduce.group(unique(reduce.group(
		unique(reduce.group(unique(reduce.group(unique(reduce.group(unique(reduce.group(
		distx,ncor=NCORE)),ncor= NCORE)),ncor= NCORE)),ncor= NCORE)),
		ncor= NCORE)),ncor= NCORE)),ncor= NCORE)),ncor= NCORE)),
		ncor= NCORE)),ncor= NCORE)),ncor= NCORE)),ncor= NCORE))	
	distx <-c(distx,setdiff(c(1:length(px[,1])) ,as.numeric(num.strsplit.u(distx))))

	seqc<-sapply(c(1:length(distx)),function(x){	
		pxx<-num.strsplit.u(distx[x])
		scax<-unique(px[pxx,1])
		posmin=min(as.numeric(px[pxx,c(2:3)]))
		posmax=max(as.numeric(px[pxx,c(2:3)]))
		seqx<-paste(SEQ[[subset(c(1:length(SEQ)),names(SEQ)== scax)]][posmin: posmax],collapse="")
	})
	
	scac<-sapply(c(1:length(distx)),function(x){	
		pxx<-num.strsplit.u(distx[x])
		scax<-unique(px[pxx,1])
	})

	t(cbind(paste(POSx[1:4],collapse=" "),c(1:length(seqc)), scac, seqc))

})

SEQhyp2<-matrix(unlist(SEQhyp2),ncol=4,byrow=T)
colnames(SEQhyp2)<-c("genome","cluster","scaffold","sequence")

####For genomes with scaffold names different between genome and RAST files, retrieve complete genomic sequence from the original genomes file

SEQhyp3<-sapply(iBCHL3,function(x){
	print(x)
	X=iBCHL2[x]
	POSx<-POSbchl[x,]
	#annotation file
	ANNOx<-read.table(DATA[X,5],header=T)

	#genome file
	SEQ<-read.fasta(GP[x])
	
	#BCHL gene position
	px<-unlist(sapply(c(b1,b2),function(x){
		list(apply(subset(ANNOx[,c(1,5,6)],ANNOx[,8]==LIST[x,3]),1,paste,collapse="="))
	}))
	px<-matrix(unlist(strsplit(px,split="=")),ncol=3,byrow=T)
	
	#scaffold index
	scai<-as.numeric(matrix(unlist(strsplit(px[,1],split="_")),ncol=3,byrow=T)[,3])
	
	#Pairwise minimum distance between each gene (return NA if on different scaffolds)
	
	distx<-unlist(sapply(c(1:(length(px[,1])-1)),function(x){
		#print(x)
		X1=x
		S1= px[X1,1]
		s1=as.numeric(px[X1,2])
		e1=as.numeric(px[X1,3])
		sapply(c((X1+1):(length(px[,1]))),function(x){
			X2=x
			S2=px[X2,1]
			s2=as.numeric(px[X2,2])
			e2=as.numeric(px[X2,3])
			#print(x)
			z<-min(c(abs(s1-s2),abs(e1-e2),abs(s1-e2),abs(e1-s2)))
			paste(X1,X2,ifelse(S1==S2,z,NA),sep="_")	
		})
	}))
	distx<-matrix(as.numeric(unlist(strsplit(distx,split="_"))),ncol=3,byrow=T)
	distx<-subset(distx,is.na(distx[,3])==F)
	distx<-subset(distx,distx[,3]<=25000)
	distx<-paste(distx[,1], distx[,2],sep="_")
	distx <-unique(reduce.group(	unique(reduce.group(	unique(reduce.group(unique(reduce.group(
		unique(reduce.group(unique(reduce.group(unique(reduce.group(unique(reduce.group(
		unique(reduce.group(unique(reduce.group(unique(reduce.group(unique(reduce.group(
		distx,ncor=NCORE)),ncor= NCORE)),ncor= NCORE)),ncor= NCORE)),
		ncor= NCORE)),ncor= NCORE)),ncor= NCORE)),ncor= NCORE)),
		ncor= NCORE)),ncor= NCORE)),ncor= NCORE)),ncor= NCORE))	
	distx <-c(distx,setdiff(c(1:length(px[,1])) ,as.numeric(num.strsplit.u(distx))))

	seqc<-sapply(c(1:length(distx)),function(x){	
		pxx<-num.strsplit.u(distx[x])
		scax<-unique(scai[pxx])
		posmin=min(as.numeric(px[pxx,c(2:3)]))
		posmax=max(as.numeric(px[pxx,c(2:3)]))
		seqx<-paste(SEQ[[scax]][posmin: posmax],collapse="")
	})
	
	scac<-sapply(c(1:length(distx)),function(x){	
		pxx<-num.strsplit.u(distx[x])
		scax<-unique(scai[pxx])
	})

	t(cbind(paste(POSx[1:4],collapse=" "),c(1:length(seqc)), scac, seqc))

})

SEQhyp3<-matrix(unlist(SEQhyp3),ncol=4,byrow=T)
colnames(SEQhyp3)<-c("genome","cluster","scaffold","sequence")

###Write a single fasta file to annotate with RAST

SEQhyp<-rbind(SEQhyp2, SEQhyp3)


setqd(pP)
write.table(SEQhyp,"BCHL-cluster-sequences-Hyphomicrobiales.txt",quote=T)

#Names: seq ID + higher taxonomy level + cluster index
SEQhyNames<-matrix(unlist(strsplit(SEQhyp[,1],split=" ")),ncol=4,byrow=T)
SEQhyNames<-paste(SEQhyNames [,1],ifelse(SEQhyNames[,3]=="NA", SEQhyNames[,2],ifelse(SEQhyNames[,4]=="NA", SEQhyNames[,3],SEQhyNames[,4])),SEQhyp[,2],sep="-")


write.fasta(sapply(c(1:length(SEQhyp[,1])),function(x){
	
	return(list(unlist(strsplit(SEQhyp[x,4],split=""))))
}),names= SEQhyNames,"BCHL-cluster-sequences-Hyphomicrobiales.fas")

###Annotate in RAST

#Sparse gbk file and for each species pool clusters from the same species in a new gbk file

setwd(pP)

SEQallgbk<-c("//",as.vector(read.delim("6666666.1487196.gbk",header=F)[,1]))

sp<-grep("//", SEQallgbk)

NAMESbgk<-t(sapply(c(2:length(sp)),function(x){
	seqx<-SEQallgbk[sp[x-1]:sp[x]]
	nx<-setdiff(unlist(strsplit(seqx[2],
		split=" ")),"")[2]
	return(c(unlist(strsplit(nx,split="-"))[1:2],sp[x-1],sp[x]))
}))

NAMESbgk<-cbind(apply(NAMESbgk[,c(1:2)],1,paste,collapse="-"), 
	NAMESbgk[,c(3:4)])

NAMESu<-sort(unique(NAMESbgk[,1]))

sapply(NAMESu,function(x){
	print(x)
	X=x
	nx<-subset(NAMESbgk, NAMESbgk[,1]==X)
	SPEgbk<-unlist(sapply(
		c(1:length(nx[,1])),function(x){
		return(list(SEQallgbk[as.numeric(
			nx[x,2]):as.numeric(nx[x,3])]))
	
	}))
	setwd(paste(pP,"For_Clinker",sep="/"))
	write.table(SPEgbk[2:length(SPEgbk)],
		paste(X,"-BCHL.gbk",sep=""),
		row.names=F,col.names=F,quote=F)	
})

#Align clusters with CLINKER https://cagecat.bioinformatics.nl/tools/clinker
	
###Do the same by removing genomes that are too fragmented
###Remove genome with n50 (calculated on BCHL cluster only) <20000 (too fragmented )

SUM<-t(sapply(unique(SEQhyp[,1]),function(x){
	
	seqx<-subset(SEQhyp, SEQhyp[,1]==x)
	sizex<-sort(nchar(seqx[,4]),decreasing=T)
	n50<-subset(sizex ,(cumsum(sizex)/sum(sizex))>=0.5)[1]
	return(c(sum(sizex), n50,length(seqx[,3]),length(unique(seqx[,3]))))
}))
plot(SUM[,2], SUM[,4])
SUM<-subset(SUM,SUM[,2]>=20000)

NAMESuF<-as.vector(sapply(rownames(SUM),function(x){
	stx<-unlist(strsplit(x,split=" "))[1]
	NAMESu[grep(paste(stx,"-",sep=""), NAMESu)]
}))

sapply(NAMESuF,function(x){
	print(x)
	X=x
	nx<-subset(NAMESbgk, NAMESbgk[,1]==X)
	SPEgbk<-unlist(sapply(
		c(1:length(nx[,1])),function(x){
		return(list(SEQallgbk[as.numeric(
			nx[x,2]):as.numeric(nx[x,3])]))
	
	}))
	setwd(paste(pP,"For_Clinker",sep="/"))
	write.table(SPEgbk[2:length(SPEgbk)],
		paste(X,"-BCHL.gbk",sep=""),
		row.names=F,col.names=F,quote=F)	
})

######MARCH 25th 2026: After a comparison with CLINKER (three different sessions), refining of genomes to choose for the sumary figure. 
#CLINKER SESSIONS: 
# W678X869S241X23: Methylobacterium dec 2025; 
# G035U634M234V83: Hyphomicrobiales extended (with 4 representative Methylobacterium); # N467B133V970G77 : Hyphomicrobiales restricted to good quality genomes and MAGs (with 4 representative Methylobacterium)

#METHYLOBACTERIUM: best genome for each lineage (excluding C3 (M. nodulans): no BCHL cluster: same genomes than  W678X869S241X23)

# D6=Meth_iners=st135# D5=Meth_sp_017=st107# D4=Meth_sp_001=st33# D3=Meth_sp_011=st120# D2=Meth_flocculans=st226# D1=Meth_sp_012=st1# C7=Meth_carpenticola=st315# C6=Meth_oryzihabitans=st185# C5=Meth_oryzisoli=st313# C4=Meth_sp_044=st16# C3 : no# C2=Meth_crusticola=st151# C1=Meth_alsaeris=st314# B4=Meth_organophilum=st216# B3=Meth_brachythecii=st55# B2=Meth_suomiense=st155# B1=Meth_extorquens=st49# A7=Meth_jeotgali=st137# A6=Meth_trifolii=st147# A5=Meth_sp_018=st62# A4=Meth_radiodurans=st11# A3=Meth_durans=st7# A2=Meth_sp_046=st236# A1=Meth_radiotolerans=st68

####IN G035U634M234V83 only (too fragmented, excluded for N467B133V970G77) - After this round, genera Enterovirga, JAJXWB01 and Rhabdaerophilum are thus excluded

# hyp52-Enterovirga_sp-BCHL.gbk
# hyp38-Enterovirga_sp-BCHL.gbk
# hyp165-Salinarimonas_rosea-BCHL.gbk
# hyp30-JAJXWB01-BCHL.gbk
# hyp79-RH_AL1-BCHL.gbk# hyp76-RH_AL1-BCHL.gbk# hyp74-RH_AL1-BCHL.gbk# hyp72-RH_AL1-BCHL.gbk# hyp70-RH_AL1-BCHL.gbk# hyp68-RH_AL1-BCHL.gbk# hyp67-RH_AL1-BCHL.gbk# hyp66-RH_AL1-BCHL.gbk# hyp64-RH_AL1-BCHL.gbk# hyp47-RH_AL1-BCHL.gbk# hyp45-RH_AL1-BCHL.gbk# hyp44-RH_AL1-BCHL.gbk# hyp35-RH_AL1-BCHL.gbk
# hyp75-Rhabdaerophilum-BCHL.gbk

####IN G035U634M234V83 and N467B133V970G77

#Methylobacterium: same genomes as for W678X869S241X23
# st68-A1=Meth_radiotolerans-BCHL.gbk
# st49-B1=Meth_extorquens-BCHL.gbk
# st1-D1=Meth_sp_012-BCHL.gbk
# st314-C1=Meth_alsaeris-BCHL.gbk

#Salinarimonas: the two clusters have same architecture but one is more complete# hyp166-Salinarimonas_soli-BCHL.gbk : THIS ONE TO CONSERVE (all genes in a single cluster)# hyp164-Salinarimonas_ramus-BCHL.gbk : to discard (too fragmented)

#Bosea: dubious (fragmented): BCHL cluster is only present in one of the 14 genomes
# hyp117-Bosea_massiliensis-BCHL.gbk : THIS ONE TO CONSERVE but fragmented (dubious structure)

#Roseiarcus: only one genome, with a BCHL cluster
# hyp1-Roseiarcus_fermantans-BCHL.gbk  : THIS ONE TO CONSERVE but fragmented (dubious structure)

#Rhodoplastus: the two BCHL clusters have nearly identical architecture# hyp163-Rhodoblastus_sphagnicola-BCHL.gbk : to discard (lower N50 for Rhodoblastus)# hyp162-Rhodoblastus_acidophilus-BCHL.gbk : THIS ONE TO CONSERVE (all genes in a single cluster and higher N50 for Rhodoblastus))

#Methylocystis: the two BCHL clusters have nearly identical architecture (all genes in a single cluster and high quality genomes)
# hyp83-Methylocystis_rosea-BCHL.gbk : THIS ONE TO CONSERVE (genome with 2 scaffolds)
# hyp167-Methylocystis_echinoides-BCHL.gbk: to discard (genome with 11 scaffolds)

#Methylocella
# hyp161-Methylocella_silvestris-BCHL.gbk : THIS ONE TO CONSERVE (all genes in a single cluster)

#Methylocapsa
# hyp160-Methylocapsa_palsarum-BCHL.gbk : THIS ONE TO CONSERVE (all genes in a single cluster)

#RH_AL1 : nearly identical and high quality genomes but cluster in four parts: fragmentation confirmed. Keep the genome with the most complete sequence (hyp168 is lacking complete bchD and rpoB sequence)# hyp169-RH_AL1-BCHL.gbk: THIS ONE TO CONSERVE (bchD sequence is complete)# hyp168-RH_AL1-BCHL.gbk: to discard (bchD sequence fragmented)

#Lichenifustis# hyp8-Lichenifustis_flavocetrariae-BCHL.gbk: THIS ONE TO CONSERVE (all genes in a single cluster, BUT bchH is duplicated: dubious because the extra copy is in a different scaffold)

#Lichenihabitans# hyp22-Lichenihabitans_sp-BCHL.gbk: THIS ONE TO CONSERVE but fragmented (dubious structure)

#Lichenibacterium: keep hyp21 (all genes in a signe cluster), hyp3 (two clusters) + hyp3 (two clusters, alternative bchG positions)# hyp21-Lichenibacterium_sp-BCHL.gbk: THIS ONE TO CONSERVE (all genes in a single cluster)# hyp19-Lichenibacterium_sp-BCHL.gbk: THIS ONE TO CONSERVE (genes in two clusters)# hyp3-Lichenibacterium_minor-BCHL.gbk: THIS ONE TO CONSERVE (genes in two clusters, alternative bchG positions)# hyp6-Lichenibacterium_dinghuense-BCHL.gbk: to discard (almost identical to hyp3 but more fragmented)
#### SUMMARY Hyphomicrobiales genomes to keep (in addition to Methylobacterium representatives for eac lineage):

#Salinarimonas: hyp166-Salinarimonas_soli-BCHL.gbk (all genes in a single cluster)
#Bosea: hyp117-Bosea_massiliensis-BCHL.gbk (dubious cluster fragmentation but can be fused with overlaping sequence)
#Roseiarcus: hyp1-Roseiarcus_fermantans-BCHL.gbk (several clusters, to confirm)
#Rhodoplastus: hyp162-Rhodoblastus_acidophilus-BCHL.gbk (all genes in a single cluster)
#Methylocystis: hyp83-Methylocystis_rosea-BCHL.gbk (all genes in a single cluster)
#Methylocella: hyp161-Methylocella_silvestris-BCHL.gbk (all genes in a single cluster)
#Methylocapsa: hyp160-Methylocapsa_palsarum-BCHL.gbk (all genes in a single cluster)
#RH_AL1 : hyp169-RH_AL1-BCHL.gbk (cluster fragmentation confirmed)
#Lichenifustis: hyp8-Lichenifustis_flavocetrariae-BCHL.gbk (all genes in a single cluster, BUT bchH is duplicated: dubious because the extra copy is in a different scaffold)
#Lichenihabitans: hyp22-Lichenihabitans_sp-BCHL.gbk (dubious cluster fragmentation but can be fused with overlaping sequence)
#Lichenibacterium (3 genomes): hyp21-Lichenibacterium_sp-BCHL.gbk (all genes in a single cluster), hyp19-Lichenibacterium_sp-BCHL.gbk (genes in two clusters) and hyp3-Lichenibacterium_minor-BCHL.gbk (genes in two clusters, alternative bchG positions) 

#hyphomicrobiale genomes with unusual bchG localisation in the cluster: hyp19, hyp22, hyp160, hyp161, hyp1, hyp117, hyp168,169 (these, together with Methylobacterium lineages D3, D5 and B3, form a monophyletic group in bchG phylogeny)

############################################################################################################################################################################################################
######STEP G - Format the tree for figure displaying CBHL cluster architecture in Methylobacterium and related Hyphomycrobiales

#Combine "Methylobacterium-AsterSpeTree-PEG-SPE-sept-2025_mBS=10.tre"

#and Lichenihabitantaceae-AsterLinTree-nov-2025_mBS=50.tre



library(ape)
library(ips)
setwd(pP)

TREElin<-read.tree("Methylobacterium-AsterSpeTree-PEG-SPE-sept-2025_mBS=10.tre")

linx<-sapply(TREElin$tip,function(x){
	return(unique(unlist(strsplit(DATA[grep(x,DATA[,4]),4],split="="))[1]))
})

plot(root(TREElin,outgroup=c(TREElin$tip[grep("A", linx)],TREElin$tip[grep("B", linx)],TREElin$tip[grep("D", linx)])),cex=0.5)

plot(root(TREElin,outgroup=TREElin$tip[grep("C", linx)]),cex=0.5)


# Lichenihabitantaceae-AsterLinTree-nov-2025_mBS=50.tre (Clades manually simplified at the genus level)

### ==> HYPHOMICROBIALES-Manually-Simplified.tre
(((((((((((((((((Methylocapsa)1:0.0296423,(Methylocella)1:0.100233)1:0.0928316,(Beijerinckia)1:0.216918)1:0.0170314,Methyloferula)1:0.0857297,Methylovirgula)1:0.0968992,(RH_AL1)1:0.129063)1:0.0913055,((((Methylocystis)1:0.0469977,(Methylosinus)1:0.121849)1:0.164874,(JAJXWB01)1:0.103225)0.999995:0.0312205,(Rhodoblastus)1:0.213136)1:0.0869375)1:0.103233,((Lichenibacterium)1:0.131331,((Lichenihabitans)1:0.0841063,(Lichenifustis)1:0.120374)0.999689:0.0379801)1:0.116523)0.0595489:0.00588331,Roseiarcus)1:0.126702,(Rhabdaerophilum)1:0.0918673)1:0.0515366,(Alsobacter)1:0.102946)1:0.0705957,(Bosea)1:0.187638)0.999923:0.0286404,((METHYLOBACTERIUM)1:0.0781796,(((Microvirga)1:0.106257,Hohaiivirga)1:0.0551153,(Salinarimonas)1:0.0629015)0.98569:0.0877014)1:0.093443)1:0.0684973,((Pseudochelatococcus)0.999992:0.0396635,(Chelatococcus)1:0.0348014)1:0.0471207)1:0.142027,seq97_Chelatococcaceae_Chelatococcus_reniformis:0.18207)1:0.0657717,(seq119_Chelatococcaceae_Camelimonas_lactis:0.020459,seq98_Chelatococcaceae_Camelimonas_fluminis:0.0286223)1:0.0882976):0.0791238,seq116_Chelatococcaceae_Camelimonas_abortus:0.0791238);

#then then open in figtree, root on rhabdaerophilum and export as newick (save as currently displayed); manually add branch labels (lost in the process)

### HYPHOMICROBIALES-Manually-Simplified-Rooted.tre
((((((((((Methylocapsa:1.0296423,Methylocella:1.100233)1:0.0928316,Beijerinckia:1.216918)1:0.0170314,Methyloferula:1.0)1:0.0857297,Methylovirgula:1.0)1:0.0968992,RH_AL1:1.129063)1:0.0913055,(((Methylocystis:1.0469977,Methylosinus:1.121849)1:0.164874,JAJXWB01:1.103225)1:0.0312205,Rhodoblastus:1.213136)1:0.0869375)1:0.103233,(Lichenibacterium:1.131331,(Lichenihabitans:1.0841063,Lichenifustis:1.120374)1:0.0379801)1:0.116523)0.06:0.00588331,Roseiarcus:1.0)1:0.126702,(Alsobacter:1.102946,(Bosea:1.187638,((METHYLOBACTERIUM:1.0781796,((Microvirga:1.106257,Hohaiivirga:1.0)1:0.0551153,Salinarimonas:1.0629015)0.99:0.0877014)1:0.093443,((Pseudochelatococcus:1.0396635,Chelatococcus:1.0348014)1:0.0471207,(Camelimonas):0.142027)1:0.0684973)1:0.0286404)1:0.0705957)1:0.0515366)1:0.04593365,Rhabdaerophilum)1;

#Open Methylobacterium-AsterSpeTree-PEG-SPE-sept-2025_mBS=10-FIGTREE-Lineages in figtree, root on group C and export as newick tree (save as currently displayed); manually add branch labels (lost in the process)

### Methylobacterium-Species-Manually-Simplified-Rooted.tre

(((((((((((((((Meth_sp_020:0.034346,Meth_guangdongense:0.029036)1:0.006974,Meth_sp_019:0.030816)1:0.022691,Meth_phyllostachyos:0.034216)1:0.004341,Meth_sp_023:0.049483)1:0.002127,(Meth_sp_021:0.015257,Meth_sp_027:0.015591)1:0.040805)1:0.033535,Meth_sp_022:0.056765)1:0.025714,(((((Meth_sp_028:0.019266,Meth_radiotolerans:0.016251)1:0.021522,Meth_sp_031:0.041033)0.98:0.001334,Meth_fujisawaense:0.047852)1:0.016033,Meth_tardum:0.050151)1:0.00778,(Meth_longum:0.068892,Meth_sp_025:0.056442)1:0.0177937)1:0.0085919)1:0.014797,((Meth_sp_032:0.016546,Meth_mesophilicum:0.0090173)1:0.035720,(((((Meth_pseudosasicola:0.031846,Meth_sp_029:0.027506)1:0.004127,Meth_sp_024:0.034571)0.53:4.949E-4,Meth_sp_026:0.0055572)1:0.014782,Meth_sp_030:0.049144)1:0.0078178,Meth_brachiatum:0.049588)1:5.894E-4)1:0.0411986)1:0.037596,(((Meth_sp_051:0.053655,Meth_sp_049:0.05282)1:0.03547,Meth_aerolatum:0.086317)1:0.022944,((Meth_sp_034:0.050064,Meth_persicinum:0.082557)1:0.0256813,(((Meth_sp_046:0.020293,Meth_komagatae:0.0184766)1:0.067270,Meth_litchii:0.057231)1:0.0154090,Meth_sp_052:0.088450)1:0.017940)1:0.0211540)1:0.0159713)1:0.075894,((Meth_symbioticum:0.027374,Meth_dankookense:0.0130460)1:0.059927,((Meth_gregans:0.039628,Meth_hispanicum:0.03321)0.96:0.0078676,Meth_radiodurans:0.031067)1:0.0539)1:0.0411964)1:0.012954,(Meth_trifolii:0.083523,(Meth_sp_018:0.020214,Meth_cerastii:0.0291480)1:0.109211)1:0.025529)1:0.042161,(((Meth_nigriterrae:0.084043,(Meth_planeticum:0.04296,Meth_soli:0.051043)1:0.040192)1:0.020058,((Meth_durans:0.046953,Meth_segetis:0.039424)1:0.0153539,Meth_oxalidis:0.047271)1:0.052568)1:0.027285,Meth_jeotgali:0.169275)1:0.0070275)1:0.026054,((((((((Meth_sp_037:0.027697,(Meth_extorquens:0.0242025,Meth_sp_040:0.0225413)1:0.005103)1:0.00319504,Meth_sp_036:0.0267687)1:0.033297,((((Meth_aminovorans:0.051095,(Meth_sp_041:0.0205577,Meth_rhodesianum:0.0260558)1:0.021027)0.58:5.9102E-4,Meth_sp_039:0.052955)1:0.004265,Meth_podarium:0.0278956)1:0.010245,Meth_thiocyanatum:0.043158)0.78:0.006466)1:0.025243,Meth_sp_035:0.064162)1:0.008904,Meth_zatmanii:0.072411)1:0.033378,(Meth_suomiense:0.079577,(Meth_rhodinum:0.042613,((Meth_sp_038:0.0300254,(Meth_sp_047:0.0160013,Meth_subtropicum:0.0153348)1:0.0164275)1:0.0120273,(Meth_sp_048:0.0259317,Meth_salsuginis:0.0197660)1:0.015896)1:0.011616)1:0.0145685)1:0.042928)1:0.062535,Meth_organophilum:0.102709)1:0.0295714,((Meth_brachythecii:0.054746,Meth_gnaphalii:0.068363)1:0.051653,(Meth_thuringiense:0.049700,Meth_haplocladii:0.053974)1:0.030575)1:0.033424)1:0.045433)1:0.0085270,(Meth_iners:0.141913,(Meth_sp_017:0.167726,(((Meth_sp_013:0.062254,Meth_marchantiae:0.073353)1:0.0064789,(Meth_sp_016:0.046505,((Meth_sp_014:0.0266253,Meth_sp_012:0.028291)1:0.006357,Meth_bullatum:0.024475)1:0.0206185)1:0.036810)1:0.088916,(((((Meth_sp_011:0.0283368,Meth_sp_004:0.0273990)1:0.018683,Meth_sp_007:0.057703)1:0.015135,((Meth_sp_053:0.020892,Meth_gossipiicola:0.024168)1:0.034148,Meth_adhaesivum:0.050958)1:0.023390)1:0.064278,Meth_sp_001:0.085727)1:0.0514154,((((((Meth_sp_010:0.039986,Meth_sp_008:0.032729)1:0.00210085,Meth_sp_006:0.0290032)1:0.0221035,Meth_sp_005:0.0391104)1:0.0241355,Meth_sp_009:0.090856)1:0.0196357,(Meth_sp_003:0.0057616,Meth_goesingense:0.0120836)1:0.05343095)1:0.0155035,(Meth_sp_002:0.069606,Meth_flocculans:0.078791)1:0.0108300)1:0.06952)1:0.034373)1:0.025525)1:0.0091620)1:0.058248)1:0.0917205,(((Meth_sp_044:0.078861,(Meth_carpenticola:0.051005,Meth_nodulans:0.054660)1:0.0206114)1:0.0255344,(Meth_oryzisoli:0.0195427,Meth_isbiliense:0.0185376)1:0.0582016)1:0.0549086,(Meth_oryzihabitans:0.084343,(Meth_crusticola:0.081242,(Meth_variabile:0.034259,(Meth_nonmethylotrophicum:0.053180,((((Meth_sp_042:0.035237,Meth_currus:0.0296652)1:0.0064229,Meth_frigidaeris:0.035165)1:0.0093214,(Meth_aquaticum:0.050692,Meth_tarhaniae:0.0339556)1:0.0068242)1:0.0143216,(((Meth_indicum:0.0197964,Meth_ajmalii:0.020337)1:0.031827,Meth_sp_050:0.042023)1:0.0053915,((Meth_terricola:0.0328626,Meth_alsaeris:0.030549)1:0.0136117,(Meth_terrae:0.0500074,Meth_platani:0.033241)0.59:5.6642E-4)1:0.0185617)1:0.00148133)1:0.0102805)0.97:0.0018655)1:0.069791)1:0.037849)1:0.0116930)1:0.0917204)1;

#Insert Methylobacterium-Species-Manually-Simplified-Rooted.tre in HYPHOMICROBIALES-Manually-Simplified-Rooted.tre (in clade "METHYLOBACTERIUM")

### HYPHOMICROBIALES+METHYLOBACERIUM-Manually-Simplified-Rooted.tre
((((((((((Methylocapsa:1.0296423,Methylocella:1.100233)1:0.0928316,Beijerinckia:1.216918)1:0.0170314,Methyloferula:1.0)1:0.0857297,Methylovirgula:1.0)1:0.0968992,RH_AL1:1.129063)1:0.0913055,(((Methylocystis:1.0469977,Methylosinus:1.121849)1:0.164874,JAJXWB01:1.103225)1:0.0312205,Rhodoblastus:1.213136)1:0.0869375)1:0.103233,(Lichenibacterium:1.131331,(Lichenihabitans:1.0841063,Lichenifustis:1.120374)1:0.0379801)1:0.116523)0.06:0.00588331,Roseiarcus:1.0)1:0.126702,(Alsobacter:1.102946,(Bosea:1.187638,(((((((((((((((((Meth_sp_020:0.034346,Meth_guangdongense:0.029036)1:0.006974,Meth_sp_019:0.030816)1:0.022691,Meth_phyllostachyos:0.034216)1:0.004341,Meth_sp_023:0.049483)1:0.002127,(Meth_sp_021:0.015257,Meth_sp_027:0.015591)1:0.040805)1:0.033535,Meth_sp_022:0.056765)1:0.025714,(((((Meth_sp_028:0.019266,Meth_radiotolerans:0.016251)1:0.021522,Meth_sp_031:0.041033)0.98:0.001334,Meth_fujisawaense:0.047852)1:0.016033,Meth_tardum:0.050151)1:0.00778,(Meth_longum:0.068892,Meth_sp_025:0.056442)1:0.0177937)1:0.0085919)1:0.014797,((Meth_sp_032:0.016546,Meth_mesophilicum:0.0090173)1:0.035720,(((((Meth_pseudosasicola:0.031846,Meth_sp_029:0.027506)1:0.004127,Meth_sp_024:0.034571)0.53:4.949E-4,Meth_sp_026:0.0055572)1:0.014782,Meth_sp_030:0.049144)1:0.0078178,Meth_brachiatum:0.049588)1:5.894E-4)1:0.0411986)1:0.037596,(((Meth_sp_051:0.053655,Meth_sp_049:0.05282)1:0.03547,Meth_aerolatum:0.086317)1:0.022944,((Meth_sp_034:0.050064,Meth_persicinum:0.082557)1:0.0256813,(((Meth_sp_046:0.020293,Meth_komagatae:0.0184766)1:0.067270,Meth_litchii:0.057231)1:0.0154090,Meth_sp_052:0.088450)1:0.017940)1:0.0211540)1:0.0159713)1:0.075894,((Meth_symbioticum:0.027374,Meth_dankookense:0.0130460)1:0.059927,((Meth_gregans:0.039628,Meth_hispanicum:0.03321)0.96:0.0078676,Meth_radiodurans:0.031067)1:0.0539)1:0.0411964)1:0.012954,(Meth_trifolii:0.083523,(Meth_sp_018:0.020214,Meth_cerastii:0.0291480)1:0.109211)1:0.025529)1:0.042161,(((Meth_nigriterrae:0.084043,(Meth_planeticum:0.04296,Meth_soli:0.051043)1:0.040192)1:0.020058,((Meth_durans:0.046953,Meth_segetis:0.039424)1:0.0153539,Meth_oxalidis:0.047271)1:0.052568)1:0.027285,Meth_jeotgali:0.169275)1:0.0070275)1:0.026054,((((((((Meth_sp_037:0.027697,(Meth_extorquens:0.0242025,Meth_sp_040:0.0225413)1:0.005103)1:0.00319504,Meth_sp_036:0.0267687)1:0.033297,((((Meth_aminovorans:0.051095,(Meth_sp_041:0.0205577,Meth_rhodesianum:0.0260558)1:0.021027)0.58:5.9102E-4,Meth_sp_039:0.052955)1:0.004265,Meth_podarium:0.0278956)1:0.010245,Meth_thiocyanatum:0.043158)0.78:0.006466)1:0.025243,Meth_sp_035:0.064162)1:0.008904,Meth_zatmanii:0.072411)1:0.033378,(Meth_suomiense:0.079577,(Meth_rhodinum:0.042613,((Meth_sp_038:0.0300254,(Meth_sp_047:0.0160013,Meth_subtropicum:0.0153348)1:0.0164275)1:0.0120273,(Meth_sp_048:0.0259317,Meth_salsuginis:0.0197660)1:0.015896)1:0.011616)1:0.0145685)1:0.042928)1:0.062535,Meth_organophilum:0.102709)1:0.0295714,((Meth_brachythecii:0.054746,Meth_gnaphalii:0.068363)1:0.051653,(Meth_thuringiense:0.049700,Meth_haplocladii:0.053974)1:0.030575)1:0.033424)1:0.045433)1:0.0085270,(Meth_iners:0.141913,(Meth_sp_017:0.167726,(((Meth_sp_013:0.062254,Meth_marchantiae:0.073353)1:0.0064789,(Meth_sp_016:0.046505,((Meth_sp_014:0.0266253,Meth_sp_012:0.028291)1:0.006357,Meth_bullatum:0.024475)1:0.0206185)1:0.036810)1:0.088916,(((((Meth_sp_011:0.0283368,Meth_sp_004:0.0273990)1:0.018683,Meth_sp_007:0.057703)1:0.015135,((Meth_sp_053:0.020892,Meth_gossipiicola:0.024168)1:0.034148,Meth_adhaesivum:0.050958)1:0.023390)1:0.064278,Meth_sp_001:0.085727)1:0.0514154,((((((Meth_sp_010:0.039986,Meth_sp_008:0.032729)1:0.00210085,Meth_sp_006:0.0290032)1:0.0221035,Meth_sp_005:0.0391104)1:0.0241355,Meth_sp_009:0.090856)1:0.0196357,(Meth_sp_003:0.0057616,Meth_goesingense:0.0120836)1:0.05343095)1:0.0155035,(Meth_sp_002:0.069606,Meth_flocculans:0.078791)1:0.0108300)1:0.06952)1:0.034373)1:0.025525)1:0.0091620)1:0.058248)1:0.0917205,(((Meth_sp_044:0.078861,(Meth_carpenticola:0.051005,Meth_nodulans:0.054660)1:0.0206114)1:0.0255344,(Meth_oryzisoli:0.0195427,Meth_isbiliense:0.0185376)1:0.0582016)1:0.0549086,(Meth_oryzihabitans:0.084343,(Meth_crusticola:0.081242,(Meth_variabile:0.034259,(Meth_nonmethylotrophicum:0.053180,((((Meth_sp_042:0.035237,Meth_currus:0.0296652)1:0.0064229,Meth_frigidaeris:0.035165)1:0.0093214,(Meth_aquaticum:0.050692,Meth_tarhaniae:0.0339556)1:0.0068242)1:0.0143216,(((Meth_indicum:0.0197964,Meth_ajmalii:0.020337)1:0.031827,Meth_sp_050:0.042023)1:0.0053915,((Meth_terricola:0.0328626,Meth_alsaeris:0.030549)1:0.0136117,(Meth_terrae:0.0500074,Meth_platani:0.033241)0.59:5.6642E-4)1:0.0185617)1:0.00148133)1:0.0102805)0.97:0.0018655)1:0.069791)1:0.037849)1:0.0116930)1:0.0917204)1:1.0781796,((Microvirga:1.106257,Hohaiivirga:1.0)1:0.0551153,Salinarimonas:1.0629015)0.99:0.0877014)1:0.093443,((Pseudochelatococcus:1.0396635,Chelatococcus:1.0348014)1:0.0471207,(Camelimonas):0.142027)1:0.0684973)1:0.0286404)1:0.0705957)1:0.0515366)1:0.04593365,Rhabdaerophilum)1;

### HYPHOMICROBIALES+METHYLOBACERIUM-Manually-Simplified-Rooted_adj.tre #terminal branches removed 
((((((((((Methylocapsa:0.0001,Methylocella:0.0001)1:0.0928316,Beijerinckia:0.0001)1:0.0170314,Methyloferula:0.0001)1:0.0857297,Methylovirgula:0.0001)1:0.0968992,RH_AL1:0.0001)1:0.0913055,(((Methylocystis:0.0001,Methylosinus:0.0001)1:0.164874,JAJXWB01:0.0001)1:0.0312205,Rhodoblastus:0.0001)1:0.0869375)1:0.103233,(Lichenibacterium:0.0001,(Lichenihabitans:0.0001,Lichenifustis:0.0001)1:0.0379801)1:0.116523)0.06:0.00588331,Roseiarcus:0.0001)1:0.126702,(Alsobacter:0.0001,(Bosea:0.0001,(((((((((((((((((Meth_sp_020:0.034346,Meth_guangdongense:0.029036)1:0.006974,Meth_sp_019:0.030816)1:0.022691,Meth_phyllostachyos:0.034216)1:0.004341,Meth_sp_023:0.049483)1:0.002127,(Meth_sp_021:0.015257,Meth_sp_027:0.015591)1:0.040805)1:0.033535,Meth_sp_022:0.056765)1:0.025714,(((((Meth_sp_028:0.019266,Meth_radiotolerans:0.016251)1:0.021522,Meth_sp_031:0.041033)0.98:0.001334,Meth_fujisawaense:0.047852)1:0.016033,Meth_tardum:0.050151)1:0.00778,(Meth_longum:0.068892,Meth_sp_025:0.056442)1:0.0177937)1:0.0085919)1:0.014797,((Meth_sp_032:0.016546,Meth_mesophilicum:0.0090173)1:0.035720,(((((Meth_pseudosasicola:0.031846,Meth_sp_029:0.027506)1:0.004127,Meth_sp_024:0.034571)0.53:4.949E-4,Meth_sp_026:0.0055572)1:0.014782,Meth_sp_030:0.049144)1:0.0078178,Meth_brachiatum:0.049588)1:5.894E-4)1:0.0411986)1:0.037596,(((Meth_sp_051:0.053655,Meth_sp_049:0.05282)1:0.03547,Meth_aerolatum:0.086317)1:0.022944,((Meth_sp_034:0.050064,Meth_persicinum:0.082557)1:0.0256813,(((Meth_sp_046:0.020293,Meth_komagatae:0.0184766)1:0.067270,Meth_litchii:0.057231)1:0.0154090,Meth_sp_052:0.088450)1:0.017940)1:0.0211540)1:0.0159713)1:0.075894,((Meth_symbioticum:0.027374,Meth_dankookense:0.0130460)1:0.059927,((Meth_gregans:0.039628,Meth_hispanicum:0.03321)0.96:0.0078676,Meth_radiodurans:0.031067)1:0.0539)1:0.0411964)1:0.012954,(Meth_trifolii:0.083523,(Meth_sp_018:0.020214,Meth_cerastii:0.0291480)1:0.109211)1:0.025529)1:0.042161,(((Meth_nigriterrae:0.084043,(Meth_planeticum:0.04296,Meth_soli:0.051043)1:0.040192)1:0.020058,((Meth_durans:0.046953,Meth_segetis:0.039424)1:0.0153539,Meth_oxalidis:0.047271)1:0.052568)1:0.027285,Meth_jeotgali:0.169275)1:0.0070275)1:0.026054,((((((((Meth_sp_037:0.027697,(Meth_extorquens:0.0242025,Meth_sp_040:0.0225413)1:0.005103)1:0.00319504,Meth_sp_036:0.0267687)1:0.033297,((((Meth_aminovorans:0.051095,(Meth_sp_041:0.0205577,Meth_rhodesianum:0.0260558)1:0.021027)0.58:5.9102E-4,Meth_sp_039:0.052955)1:0.004265,Meth_podarium:0.0278956)1:0.010245,Meth_thiocyanatum:0.043158)0.78:0.006466)1:0.025243,Meth_sp_035:0.064162)1:0.008904,Meth_zatmanii:0.072411)1:0.033378,(Meth_suomiense:0.079577,(Meth_rhodinum:0.042613,((Meth_sp_038:0.0300254,(Meth_sp_047:0.0160013,Meth_subtropicum:0.0153348)1:0.0164275)1:0.0120273,(Meth_sp_048:0.0259317,Meth_salsuginis:0.0197660)1:0.015896)1:0.011616)1:0.0145685)1:0.042928)1:0.062535,Meth_organophilum:0.102709)1:0.0295714,((Meth_brachythecii:0.054746,Meth_gnaphalii:0.068363)1:0.051653,(Meth_thuringiense:0.049700,Meth_haplocladii:0.053974)1:0.030575)1:0.033424)1:0.045433)1:0.0085270,(Meth_iners:0.141913,(Meth_sp_017:0.167726,(((Meth_sp_013:0.062254,Meth_marchantiae:0.073353)1:0.0064789,(Meth_sp_016:0.046505,((Meth_sp_014:0.0266253,Meth_sp_012:0.028291)1:0.006357,Meth_bullatum:0.024475)1:0.0206185)1:0.036810)1:0.088916,(((((Meth_sp_011:0.0283368,Meth_sp_004:0.0273990)1:0.018683,Meth_sp_007:0.057703)1:0.015135,((Meth_sp_053:0.020892,Meth_gossipiicola:0.024168)1:0.034148,Meth_adhaesivum:0.050958)1:0.023390)1:0.064278,Meth_sp_001:0.085727)1:0.0514154,((((((Meth_sp_010:0.039986,Meth_sp_008:0.032729)1:0.00210085,Meth_sp_006:0.0290032)1:0.0221035,Meth_sp_005:0.0391104)1:0.0241355,Meth_sp_009:0.090856)1:0.0196357,(Meth_sp_003:0.0057616,Meth_goesingense:0.0120836)1:0.05343095)1:0.0155035,(Meth_sp_002:0.069606,Meth_flocculans:0.078791)1:0.0108300)1:0.06952)1:0.034373)1:0.025525)1:0.0091620)1:0.058248)1:0.0917205,(((Meth_sp_044:0.078861,(Meth_carpenticola:0.051005,Meth_nodulans:0.054660)1:0.0206114)1:0.0255344,(Meth_oryzisoli:0.0195427,Meth_isbiliense:0.0185376)1:0.0582016)1:0.0549086,(Meth_oryzihabitans:0.084343,(Meth_crusticola:0.081242,(Meth_variabile:0.034259,(Meth_nonmethylotrophicum:0.053180,((((Meth_sp_042:0.035237,Meth_currus:0.0296652)1:0.0064229,Meth_frigidaeris:0.035165)1:0.0093214,(Meth_aquaticum:0.050692,Meth_tarhaniae:0.0339556)1:0.0068242)1:0.0143216,(((Meth_indicum:0.0197964,Meth_ajmalii:0.020337)1:0.031827,Meth_sp_050:0.042023)1:0.0053915,((Meth_terricola:0.0328626,Meth_alsaeris:0.030549)1:0.0136117,(Meth_terrae:0.0500074,Meth_platani:0.033241)0.59:5.6642E-4)1:0.0185617)1:0.00148133)1:0.0102805)0.97:0.0018655)1:0.069791)1:0.037849)1:0.0116930)1:0.0917204)1,((Microvirga:0.0001,Hohaiivirga:0.0001)1:0.0551153,Salinarimonas:0.0001)0.99:0.0877014)1:0.093443,((Pseudochelatococcus:0.0001,Chelatococcus:0.0001)1:0.0471207,Camelimonas:0.0001)1:0.0684973)1:0.0286404)1:0.0705957)1:0.0515366)1:0.04593365,Rhabdaerophilum:0.0001)1;


TREElin<-read.tree("HYPHOMICROBIALES-Manually-Simplified.tre")

plot(root(TREElin,outgroup= "Rhabdaerophilum"))


#Hyphomicrobiales+Methylob-rooted.tre (combine Lichenihabitantaceae-AsterLinTree-nov-2025_mBS=50.tre rooted on Rhabdaerophylum and Methylobacterium-Species-Manually-Simplified-Rooted.tre)

(((((((((((seq160_Beijerinckiaceae_Methylocapsa_palsarum:0.1775469999999999,(seq128_Beijerinckiaceae_Methylocapsa_aurea:0.11263899999999993,(seq131_Beijerinckiaceae_Methylocapsa_polymorpha:0.08400019999999997,seq126_Beijerinckiaceae_Methylocapsa_acidiphila:0.12403600000000004):0.014256999999999964):0.02340830000000005):0.0296422999999999,(seq130_Beijerinckiaceae_Methylocella_tundrae:0.09174729999999998,seq161_Beijerinckiaceae_Methylocella_silvestris:0.11345399999999994):0.10023300000000002):0.09283160000000001,(seq124_Beijerinckiaceae_Beijerinckia_indica:0.09837849999999992,seq127_Beijerinckiaceae_Beijerinckia_mobilis:0.11612500000000003):0.21691799999999994):0.01703140000000003,seq125_Beijerinckiaceae_Methyloferula_stellata:0.16832500000000006):0.08572970000000002,seq129_Beijerinckiaceae_Methylovirgula_ligni:0.1980320000000001):0.09689919999999996,(((((((((((seq44_Meta_Phyllosphere:0.0012810899999999847,(seq45_Meta_Phyllosphere:2.930080000000945E-4,seq35_Meta_Phyllosphere:4.0226900000006616E-4):6.239820000000229E-4):0.002425349999999993,seq41_Meta_Phyllosphere:0.0014475199999999244):0.002651339999999891,seq47_Meta_Phyllosphere:0.008601029999999898):0.06945420000000002,seq79_Meta_Moss:0.045791899999999996):0.015206300000000006,seq65_Meta_Lichen:0.07499440000000002):0.0028376400000000856,seq76_Meta_Lichen:0.06483419999999995):0.030755199999999983,seq73_Meta_Lichen:0.06949870000000002):0.12234999999999996,((seq74_Meta_Lichen:0.10763100000000003,seq72_Meta_Lichen:0.08873449999999994):0.03623920000000003,(seq68_Meta_Lichen:0.11383199999999993,(seq78_Meta_Lichen:0.08055670000000004,seq66_Meta_Lichen:0.1417060000000001):0.052978699999999934):0.0032834399999999153):0.05220960000000008):0.05142880000000005,((seq70_Meta_Lichen:0.055880699999999894,seq77_Meta_Lichen:0.05491489999999999):9.214340000001098E-4,seq64_Meta_Lichen:0.04918160000000005):0.09216399999999991):0.02327100000000004,(((seq81_Meta_Moss:0.04489910000000008,(seq168_RH_AL1_NA:0.006559020000000082,seq169_RH_AL1_NA:0.006637039999999983):0.07222649999999997):0.05135690000000004,seq67_Meta_Lichen:0.10139000000000009):0.02685139999999997,seq69_Meta_Lichen:0.04900149999999992):0.00844661999999996):0.09576059999999997,seq63_Meta_Moss:0.1768590000000001):0.12906300000000004):0.09130550000000004,((((((seq113_Methylocystaceae_Methylocystis_borbori:0.07904679999999997,(((seq90_Methylocystaceae_Methylocystis_hirsuta:0.015599799999999941,seq93_Methylocystaceae_Methylocystis_silviterrae:0.01825220000000005):4.495550000001014E-4,seq104_Methylocystaceae_Methylocystis_suflitae:0.021580900000000014):0.0026438300000000137,seq83_Methylocystaceae_Methylocystis_rosea:0.016557799999999956):0.08662820000000004):0.02435389999999993,((seq106_Methylocystaceae_Methylocystis_parvus:0.07986730000000009,seq107_Methylocystaceae_Methylocystis_iwaonis:0.0761658999999999):0.024549699999999897,seq167_Methylocystaceae_Methylocystis_echinoides:0.07928249999999992):0.05664350000000007):0.07316459999999991,(seq108_Methylocystaceae_Methylocystis_bryophila:0.17290399999999995,seq91_Methylocystaceae_Methylocystis_heyeri:0.13837499999999991):0.0680419000000001):0.04699769999999992,(seq82_Methylocystaceae_Methylosinus_trichosporium:0.051187400000000105,seq88_Methylocystaceae_Methylosinus_sporium:0.061960900000000096):0.1218490000000001):0.16487400000000008,(seq132_JAJXWB01_NA:0.167937,((seq30_Meta_Phyllosphere:0.024251800000000046,seq48_Meta_Phyllosphere:0.023000200000000026):0.2851189999999999,seq71_Meta_Lichen:0.13440300000000005):0.149362):0.10322500000000001):0.031220499999999984,(seq102_Rhodoblastaceae_Rhodoblastus_alkanivorans:0.10277800000000004,(seq163_Rhodoblastaceae_Rhodoblastus_sphagnicola:0.05506820000000001,seq162_Rhodoblastaceae_Rhodoblastus_acidophilus:0.04691880000000004):0.10085600000000006):0.213136):0.0869375):0.10323300000000002,((((((((seq27_Meta_Lichen:0.010346200000000083,seq29_Meta_Phyllosphere:0.012936800000000082):0.042856500000000075,seq21_Meta_Lichen:0.039110200000000095):0.03232919999999995,seq42_Meta_Phyllosphere:0.036199400000000104):0.027776500000000093,seq31_Meta_Phyllosphere:0.054398500000000016):0.1215520000000001,seq25_Meta_Lichen:0.10189099999999995):0.01802640000000011,((((seq3_Lichenihabitantaceae_Lichenibacterium_minor:0.046435900000000085,seq13_Meta_Lichen:0.0625325000000001):0.007176800000000094,seq6_Lichenihabitantaceae_Lichenibacterium_dinghuense:0.04246700000000003):0.05736249999999998,((seq11_Meta_Sediment:0.07369899999999996,seq50_Meta_Moss:0.049107800000000035):0.05427129999999991,seq2_Lichenihabitantaceae_Lichenibacterium_ramalinae:0.057916999999999996):0.0139876000000001):0.05839430000000001,(seq14_Meta_Lichen:0.011103499999999933,seq17_Meta_Lichen:0.017493699999999945):0.07641830000000005):0.06092539999999991):0.10202200000000006,(((seq10_Meta_Soil:0.11012399999999989,seq12_Meta_Lichen:0.167945):0.019405400000000017,seq23_Meta_Lichen:0.12170599999999998):0.06623049999999997,seq19_Meta_Lichen:0.1459140000000001):0.028294599999999948):0.13133099999999998,((((seq16_Meta_Lichen:0.09220410000000001,(seq18_Meta_Lichen:0.06493440000000006,seq15_Meta_Lichen:0.06155850000000007):0.10049700000000006):0.06813800000000003,((((seq4_Lichenihabitantaceae_Lichenihabitans_psoromatis:0.0029111600000000237,seq5_Lichenihabitantaceae_Lichenihabitans_psoromatis:0.002671330000000083):0.0089748300000001,seq7_Lichenihabitantaceae_Lichenihabitans_sp:0.01569799999999999):2.3668900000006765E-4,seq20_Meta_Lichen:0.12050399999999994):0.060103099999999965,seq22_Meta_Lichen:0.10438500000000006):0.05466720000000003):0.05486150000000001,(seq24_Meta_Lichen:0.12492999999999999,(seq33_Meta_Phyllosphere:0.09487360000000011,seq26_Meta_Lichen:0.1028420000000001):0.08721500000000004):0.05669139999999995):0.08410629999999997,((((seq49_Meta_Phyllosphere:0.042442700000000055,(seq34_Meta_Phyllosphere:6.913710000000517E-4,seq43_Meta_Phyllosphere:0.0024448800000000936):2.3342199999998314E-4):0.17281499999999994,seq32_Meta_Phyllosphere:0.0940319999999999):0.06894820000000002,(seq36_Meta_Phyllosphere:0.020159899999999897,seq40_Meta_Phyllosphere:0.020161899999999955):0.2563700000000001):0.05605020000000005,(seq8_Lichenihabitantaceae_Lichenifustis_flavocetrariae:0.14171899999999993,seq9_Lichenihabitantaceae_Lichenifustis_sp:0.11122600000000005):0.08244580000000001):0.12037399999999998):0.03798009999999996):0.11652300000000004):0.005883309999999975,seq1_Roseiarcaceae_Roseiarcus_fermantans:0.29206200000000004):0.12670199999999998,(((seq96_Alsobacteraceae_Alsobacter_metallidurans:0.0747198,seq87_Alsobacteraceae_Alsobacter_soli:0.06884190000000001):0.09084219999999998,seq103_Alsobacteraceae_Alsobacter_ponti:0.10612299999999997):0.10294600000000009,((((((seq120_Boseaceae_Bosea_robiniae:0.026509000000000005,seq111_Boseaceae_Bosea_beijingensis:0.027131400000000028):0.03877520000000001,seq123_Boseaceae_Bosea_thiooxidans:0.04258519999999999):0.027865600000000046,(seq105_Boseaceae_Bosea_minatitlanensis:0.04465779999999997,seq99_Boseaceae_Bosea_spartocytisi:0.05256709999999998):0.009755969999999947):0.0657451,((((seq86_Boseaceae_Bosea_psychrotolerans:0.0014730499999999758,seq85_Boseaceae_Bosea_vaviloviae:0.0014887399999999884):0.08192429999999995,seq121_Boseaceae_Bosea_lathyri:0.05852670000000004):0.04416690000000001,seq117_Boseaceae_Bosea_massiliensis:0.07089990000000002):0.03266190000000002,seq89_Boseaceae_Bosea_caraganae:0.044343699999999986):0.05031540000000001):0.03640549999999998,(((seq118_Boseaceae_Bosea_eneae:0.028064199999999984,seq122_Boseaceae_Bosea_lupini:0.027298299999999998):0.011655700000000047,seq110_Boseaceae_Bosea_vestrisii:0.041422700000000034):0.04975839999999998,seq112_Boseaceae_Bosea_rubneri:0.053489200000000015):0.07179400000000002):0.18763799999999997,((((((((((((((((((Meth_sp_020:0.034346,Meth_guangdongense:0.029036)1:0.006974,Meth_sp_019:0.030816)1:0.022691,Meth_phyllostachyos:0.034216)1:0.004341,Meth_sp_023:0.049483)1:0.002127,(Meth_sp_021:0.015257,Meth_sp_027:0.015591)1:0.040805)1:0.033535,Meth_sp_022:0.056765)1:0.025714,(((((Meth_sp_028:0.019266,Meth_radiotolerans:0.016251)1:0.021522,Meth_sp_031:0.041033)0.98:0.001334,Meth_fujisawaense:0.047852)1:0.016033,Meth_tardum:0.050151)1:0.00778,(Meth_longum:0.068892,Meth_sp_025:0.056442)1:0.0177937)1:0.0085919)1:0.014797,((Meth_sp_032:0.016546,Meth_mesophilicum:0.0090173)1:0.035720,(((((Meth_pseudosasicola:0.031846,Meth_sp_029:0.027506)1:0.004127,Meth_sp_024:0.034571)0.53:4.949E-4,Meth_sp_026:0.0055572)1:0.014782,Meth_sp_030:0.049144)1:0.0078178,Meth_brachiatum:0.049588)1:5.894E-4)1:0.0411986)1:0.037596,(((Meth_sp_051:0.053655,Meth_sp_049:0.05282)1:0.03547,Meth_aerolatum:0.086317)1:0.022944,((Meth_sp_034:0.050064,Meth_persicinum:0.082557)1:0.0256813,(((Meth_sp_046:0.020293,Meth_komagatae:0.0184766)1:0.067270,Meth_litchii:0.057231)1:0.0154090,Meth_sp_052:0.088450)1:0.017940)1:0.0211540)1:0.0159713)1:0.075894,((Meth_symbioticum:0.027374,Meth_dankookense:0.0130460)1:0.059927,((Meth_gregans:0.039628,Meth_hispanicum:0.03321)0.96:0.0078676,Meth_radiodurans:0.031067)1:0.0539)1:0.0411964)1:0.012954,(Meth_trifolii:0.083523,(Meth_sp_018:0.020214,Meth_cerastii:0.0291480)1:0.109211)1:0.025529)1:0.042161,(((Meth_nigriterrae:0.084043,(Meth_planeticum:0.04296,Meth_soli:0.051043)1:0.040192)1:0.020058,((Meth_durans:0.046953,Meth_segetis:0.039424)1:0.0153539,Meth_oxalidis:0.047271)1:0.052568)1:0.027285,Meth_jeotgali:0.169275)1:0.0070275)1:0.026054,((((((((Meth_sp_037:0.027697,(Meth_extorquens:0.0242025,Meth_sp_040:0.0225413)1:0.005103)1:0.00319504,Meth_sp_036:0.0267687)1:0.033297,((((Meth_aminovorans:0.051095,(Meth_sp_041:0.0205577,Meth_rhodesianum:0.0260558)1:0.021027)0.58:5.9102E-4,Meth_sp_039:0.052955)1:0.004265,Meth_podarium:0.0278956)1:0.010245,Meth_thiocyanatum:0.043158)0.78:0.006466)1:0.025243,Meth_sp_035:0.064162)1:0.008904,Meth_zatmanii:0.072411)1:0.033378,(Meth_suomiense:0.079577,(Meth_rhodinum:0.042613,((Meth_sp_038:0.0300254,(Meth_sp_047:0.0160013,Meth_subtropicum:0.0153348)1:0.0164275)1:0.0120273,(Meth_sp_048:0.0259317,Meth_salsuginis:0.0197660)1:0.015896)1:0.011616)1:0.0145685)1:0.042928)1:0.062535,Meth_organophilum:0.102709)1:0.0295714,((Meth_brachythecii:0.054746,Meth_gnaphalii:0.068363)1:0.051653,(Meth_thuringiense:0.049700,Meth_haplocladii:0.053974)1:0.030575)1:0.033424)1:0.045433)1:0.0085270,(Meth_iners:0.141913,(Meth_sp_017:0.167726,(((Meth_sp_013:0.062254,Meth_marchantiae:0.073353)1:0.0064789,(Meth_sp_016:0.046505,((Meth_sp_014:0.0266253,Meth_sp_012:0.028291)1:0.006357,Meth_bullatum:0.024475)1:0.0206185)1:0.036810)1:0.088916,(((((Meth_sp_011:0.0283368,Meth_sp_004:0.0273990)1:0.018683,Meth_sp_007:0.057703)1:0.015135,((Meth_sp_053:0.020892,Meth_gossipiicola:0.024168)1:0.034148,Meth_adhaesivum:0.050958)1:0.023390)1:0.064278,Meth_sp_001:0.085727)1:0.0514154,((((((Meth_sp_010:0.039986,Meth_sp_008:0.032729)1:0.00210085,Meth_sp_006:0.0290032)1:0.0221035,Meth_sp_005:0.0391104)1:0.0241355,Meth_sp_009:0.090856)1:0.0196357,(Meth_sp_003:0.0057616,Meth_goesingense:0.0120836)1:0.05343095)1:0.0155035,(Meth_sp_002:0.069606,Meth_flocculans:0.078791)1:0.0108300)1:0.06952)1:0.034373)1:0.025525)1:0.0091620)1:0.058248)1:0.0917205,(((Meth_sp_044:0.078861,(Meth_carpenticola:0.051005,Meth_nodulans:0.054660)1:0.0206114)1:0.0255344,(Meth_oryzisoli:0.0195427,Meth_isbiliense:0.0185376)1:0.0582016)1:0.0549086,(Meth_oryzihabitans:0.084343,(Meth_crusticola:0.081242,(Meth_variabile:0.034259,(Meth_nonmethylotrophicum:0.053180,((((Meth_sp_042:0.035237,Meth_currus:0.0296652)1:0.0064229,Meth_frigidaeris:0.035165)1:0.0093214,(Meth_aquaticum:0.050692,Meth_tarhaniae:0.0339556)1:0.0068242)1:0.0143216,(((Meth_indicum:0.0197964,Meth_ajmalii:0.020337)1:0.031827,Meth_sp_050:0.042023)1:0.0053915,((Meth_terricola:0.0328626,Meth_alsaeris:0.030549)1:0.0136117,(Meth_terrae:0.0500074,Meth_platani:0.033241)0.59:5.6642E-4)1:0.0185617)1:0.00148133)1:0.0102805)0.97:0.0018655)1:0.069791)1:0.037849)1:0.0116930)1:0.0917204)1:0.14144299999999999,((((seq52_Meta_Lichen:0.10626900000000006,seq53_Meta_Lichen:0.10271799999999998):0.021070699999999998,seq147_Methylobacteriaceae_Enterovirga_rhinocerotis:0.12271299999999996):0.08969369999999999,((seq54_Meta_Lichen:0.13803100000000001,((seq92_Methylobacteriaceae_Enterovirga_aerilata:0.010499000000000036,seq151_Methylobacteriaceae_Enterovirga_sp_001:0.010499000000000036):0.11260199999999998,seq56_Meta_Lichen:0.10873):0.04759879999999994):0.051192899999999986,(seq55_Meta_Lichen:0.013611599999999946,seq51_Meta_Lichen:0.019416600000000006):0.18135500000000004):0.03235080000000001):0.020415300000000025,seq38_Meta_Phyllosphere:0.15263099999999996):0.12849699999999997):0.07817960000000002,(((seq144_Methylobacteriaceae_Microvirga_ossetica:0.06292609999999998,seq149_Methylobacteriaceae_Microvirga_sp_001:0.06972080000000003):0.10625700000000005,seq114_Methylobacteriaceae_Hohaiivirga_grylli:0.5101929999999999):0.05511529999999998,((seq165_Methylobacteriaceae_Salinarimonas_rosea:0.047891499999999976,seq164_Methylobacteriaceae_Salinarimonas_ramus:0.053181900000000004):0.09586039999999996,seq166_Methylobacteriaceae_Salinarimonas_soli:0.11382700000000001):0.06290149999999994):0.08770140000000004):0.09344299999999994,((((seq115_Chelatococcaceae_Pseudochelatococcus_lubricantis:0.06515300000000002,seq94_Chelatococcaceae_Pseudochelatococcus_contaminans:0.16230199999999995):0.19251299999999993,seq100_Chelatococcaceae_Chelatococcus_asaccharovorans:0.16662899999999992):0.039663500000000074,(seq109_Chelatococcaceae_Chelatococcus_albus:0.0974238999999999,(seq101_Chelatococcaceae_Chelatococcus_composti:0.03325210000000001,(seq84_Chelatococcaceae_Chelatococcus_sambhunathii:0.020946800000000043,seq95_Chelatococcaceae_Chelatococcus_caeni:0.0049481300000000505):0.03162860000000001):0.1385979999999999):0.03480140000000009):0.047120700000000015,(seq97_Chelatococcaceae_Chelatococcus_reniformis:0.18206999999999995,((seq119_Chelatococcaceae_Camelimonas_lactis:0.020459000000000005,seq98_Chelatococcaceae_Camelimonas_fluminis:0.02862229999999988):0.08829759999999998,seq116_Chelatococcaceae_Camelimonas_abortus:0.15824760000000015):0.06577169999999999):0.1420269999999999):0.06849729999999998):0.02864040000000001):0.07059569999999993):0.051536600000000043):0.04593365,(seq75_Meta_Lichen:0.279103,seq170_Rhabdaerophilaceae_Rhabdaerophilum_calidifontis:0.207858):0.045933649999999986);

