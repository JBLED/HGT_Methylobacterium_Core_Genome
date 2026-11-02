#Script to identify cluster of genes that are conserved among different genomes

## STEP A - Calcultation of Average cluster conservation (ACC) and average ANI per species

###FUNCTIONS
library(parallel)
library(gplots)

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

#Only look in methylobacterium genomes
iMETH<-sort(intersect(subset(c(1:length(data[,1])),data[,4]!="M"),
	subset(c(1:length(data[,1])),data[,4]!="E")))


#For each genome, identify pairs of core genes that are located within 25Kb 
CluF<-mclapply(iMETH,function(x){
	print(date())

	print(x)
	stx=data[x,2]
	linx<-ORD[grep(paste("_", stx,"_",sep=""),paste("_",ORD[,5],"_",sep="")),3]
	NX<-subset(data[,6],data[,2]==stx)
	ANNOx<-read.table(paste(pD,"RAST_files",NX,sep="/"))	
	#print(date())
	#Position of each gene if present in this genome
	z<-as.vector(unlist(sapply(c(1:length(Nu)),function(x){
		
		#subset(c(1:length(ANNOx[,1])),ANNOx[,8]==x)
		paste(x,
		apply(subset(ANNOx[,c(1,5,6)],ANNOx[,8]==Nu[x]),1,paste,collapse="="),sep="=")
	})))
	#print(date())
	
	#Remove genes that are absent in this genoms	
	z<-subset(z,summary(strsplit(z,split="="))[,1]=="4")
	
	#pool together genes that are located in the same scaffold and within 25Kb
	
	Pclu<-subset(unlist(sapply(c(1:(length(z)-1)),function(x){
		X=x
		sapply(c((X+1):length(z)),function(x){
			paste(X,x,sep="_")
		})
	})),
	unlist(sapply(c(1:(length(z)-1)),function(x){
		X=x
		z1=unlist(strsplit(z[X],split="="))
		sapply(c((X+1):length(z)),function(x){
			z2=unlist(strsplit(z[x],split="="))
			ifelse(z1[2]!=z2[2],0,
			ifelse(min(c(abs(as.numeric(z2[c(3:4)])-as.numeric(z1[c(3)])),abs(as.numeric(z2[c(3:4)])-as.numeric(z1[c(4)]))))<=25000,1,0))
		})
	}))==1)
	#print(date())
	

	z<-	unique(sapply(c(1:length(Pclu)),function(x){
		cx<-num.strsplit.u(Pclu[x])
		paste(sapply(cx,function(x){
			return(as.numeric(unlist(
			strsplit(z[x],split="="))[1]))
		}),collapse="_")
		
	})	)
	setwd(paste(pP,"pairs-per-genome",sep="/"))
	write.table(z,paste("pairs_",stx,".txt",sep=""),row.names=F,col.names=F,quote=F)
	
	list(z)		
	#print(date())
		
},mc.cores= ncore)

###Safe point
CluF<-sapply(iMETH,function(x){
	print(x)
	stx=data[x,2]
	setwd(paste(pP,"pairs-per-genome",sep="/"))
	z<-as.vector(read.table(paste("pairs_",stx,".txt",sep=""))[,1])
	return(list(z))
})
	
#COnvert in Matrix of pairwise similarity in cluster conservation among genomes (Average Cluster Conservation = ACC)

MATgen<-sapply(c(1:length(CluF)),function(x){
	print(x)
	z1<-unlist(CluF[[x]])
	l1=length(z1)
	sapply(c(1:length(CluF)),function(x){
		z2<-unlist(CluF[[x]])
		return(length(intersect(z1,z2))/min(c(length(z2),l1)))
	})	
})

#average per species (ordered according to the species tree)

MATspe<-sapply(c(1:length(ORD[,1])),function(x){
	print(x)
	spex=ORD[x,1]
	stx<-unlist(strsplit(ORD[x,5],split="_"))
	ix<-sapply(stx,function(x){
		z<-subset(c(1:length(data[,1])),data[,2]==x)
		subset(c(1:length(iMETH)), iMETH==z)
	})
	sapply(c(1:length(ORD[,1])),function(x){
		print(x)
		spex=ORD[x,1]
		stx<-unlist(strsplit(ORD[x,5],split="_"))
		iy<-sapply(stx,function(x){
			z<-subset(c(1:length(data[,1])),data[,2]==x)
			subset(c(1:length(iMETH)), iMETH==z)
		})
		mean(MATgen[ix,iy])
	})
})

setwd(pP)
colnames(MATspe)<-ORD[,1]
rownames(MATspe)<-ORD[,1]
write.table(MATspe,"Average-cluster-cons-per-species.txt")

setwd(pP)
MATspe <-read.table("Average-cluster-cons-per-species.txt",header=T)

COLclu1<-seq(1,0,length.out=6)
COLclu2<-seq(1,0.5,length.out=6)
COLclu <-rgb(COLclu2, COLclu2, COLclu1)

pdf("Average-cluster-per-species_yellow.pdf",height=10,width=10)
heatmap.2(as.matrix(MATspe),scale="none",
	col= COLclu,
	Rowv=NA,Colv=NA,cexRow=0.4,cexCol=0.4,trace="none",
	key.title="Cluster conservation",breaks=length(COLclu)+1)	
dev.off()


pdf("Average-cluster-cons-per-species.pdf",height=10,width=10)
heatmap.2(as.matrix(MATspe),scale="none",
	col=heat.colors(n = 50),
	Rowv=NA,Colv=NA,cexRow=0.4,cexCol=0.4,trace="none",
	key.title="Cluster conservation")
dev.off()

#Do the same with ANI (see script 2-Core-genome.R)
setwd(paste(pD,"ANI",sep="/"))

ANI<-read.table("fastANI-out-raw.tsv")
ANI<-as.matrix(cbind(gsub(".fna","",ANI[,1]),gsub(".fna","",ANI[,2]),as.numeric(ANI[,3])))


#convert into a matrix	
ANImat<-sapply(iMETH,function(x){
	print(x)
	stx=data[x,2]
	ANIx<-rbind(subset(ANI[,c(1,3)],ANI[,2]==stx),
		subset(ANI[,c(2,3)],ANI[,1]==stx))
	sapply(iMETH,function(x){
		return(mean(as.numeric(subset(ANIx[,2],ANIx[,1]==data[x,2]))))
	})
})	

#Mean per species
ANIspe<-sapply(c(1:length(ORD[,1])),function(x){
	print(x)
	spex=ORD[x,1]
	stx<-unlist(strsplit(ORD[x,5],split="_"))
	ix<-sapply(stx,function(x){
		return(subset(c(1:length(iMETH)),data[iMETH,2]==x))
	})	
	sapply(c(1:length(ORD[,1])),function(x){
		#print(x)
		spey=ORD[x,1]
		sty<-unlist(strsplit(ORD[x,5],split="_"))
		iy<-sapply(sty,function(x){
			return(subset(c(1:length(iMETH)),data[iMETH,2]==x))
		})		
		mean(as.numeric(ANImat[ix,iy]))
	})
})


setwd(pP)
colnames(ANIspe)<-ORD[,1]
rownames(ANIspe)<-ORD[,1]
write.table(ANIspe,"Average-ANI-per-species.txt")

setwd(pP)
ANIspe<-read.table("Average-ANI-per-species.txt",header=T)

pdf("Average-ANI-per-species.pdf",height=10,width=10)
heatmap.2(as.matrix(ANIspe),scale="none",
	col=matlab.like2(n = 50),
	Rowv=NA,Colv=NA,cexRow=0.4,cexCol=0.4,trace="none",
	key.title="ANI")		
dev.off()

COLani<-seq(1,0,length.out=6)
COLani<-rgb(COLani, COLani, COLani)

pdf("Average-ANI-per-species_GREY.pdf",height=10,width=10)
heatmap.2(as.matrix(ANIspe),scale="none",
	col= COLani,
	Rowv=NA,Colv=NA,cexRow=0.4,cexCol=0.4,trace="none",
	key.title="ANI",breaks=length(COLani)+1)	
dev.off()


### STEP B : Agregate pairs of core genes in clusters , according to different levels of conservation accross lineages. For each level (range 8-24 lineages), do permutation tests to identifiy clusters that have higher residual RF values than expected


Slim=0.4 #for each species, keep pairs with relative frequency above this threshold

Llim=0.2 #for each lineage, keep pairs with relative frequency above this threshold

#Agregate core genes by species
CluSPE<-sapply(c(1:length(ORD[,1])),function(x){
	print(x)
	spex=ORD[x,1]
	stx<-unlist(strsplit(ORD[x,5],split="_"))
	clux<-unlist(sapply(stx,function(x){
		ix<-subset(c(1:length(data[,1])),data[,2]==x)
		CluF[[subset(c(1:length(iMETH)), iMETH==ix)]]
	}))
	clux<-summary(as.factor(clux),maxsum=length(unique(clux)))/length(stx)
	return(list(subset(names(clux),clux>= Slim)))
})

#Agregate core genes by lineage (weigted by the number of species per lineage)
CluLIN<-sapply(unique(ORD[,3]),function(x){
	print(x)
	ix<-subset(c(1:length(ORD[,1])),ORD[,3]==x)
	
	clux<-unlist(sapply(ix,function(x){
		return(unlist(CluSPE[[x]]))
	}))
	clux<-summary(as.factor(clux),maxsum=length(unique(clux)))/length(ix)
	return(list(subset(names(clux),clux>= Llim)))
})


CluLIN <-summary(as.factor(unlist(CluLIN)),maxsum=length(unlist(CluLIN)))

Zlim<-c(length(unique(ORD[,3])):8) #levels to test (conserved cluster of genes in the range 8-24 lineage; 24 means that cluster ar conserved accross all lineages)

NCORE=8

CLUall<-sapply(Zlim,function(x){
	Mlim=x #after filtering per lineage, keep pairs with absolute frequency above this threshold
	CluLINx <-names(CluLIN[subset(c(1:length(CluLIN)), CluLIN>= Mlim)])
	CluLINx <-unique(reduce.group(
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
	CluLINx,ncor=NCORE)),ncor= NCORE)),
		ncor= NCORE)),ncor= NCORE)),
		ncor= NCORE)),ncor= NCORE)),
		ncor= NCORE)),ncor= NCORE)),
		ncor= NCORE)),ncor= NCORE)),
		ncor= NCORE)),ncor= NCORE))	
	print(c(x,length(CluLINx),length(num.strsplit.u(CluLINx)),
		length(unique(num.strsplit.u(CluLINx)))))
	return(list(CluLINx))
})

CLUall2 <-sapply(c(1:length(Zlim)),function(x){
	clux<-c(CLUall[[x]],setdiff(c(1:length(KEGGres[,1])),num.strsplit.u(CLUall[[x]])))
	sapply(c(1:length(KEGGres[,1])),function(x){
		return(c(grep(paste("_",x,"_",sep=""),paste("_", clux,"_",sep="")),NA)[1])	
	})	
})

rownames(CLUall2)<-KEGGres[,1]
colnames(CLUall2)<-paste("L", Zlim,sep="")

setwd(pP)

write.table(CLUall2,"clustering-genes-lineages.txt")

CLUall2<-read.table("clustering-genes-lineages.txt",header=T) #Each line correspon d to a core gene; each column correspon to a level; for a given column, genes with the same index belong to the same cluster

###For each hierarchical level and each cluster, calculate average residual RF value

RFclu<-sapply(c(1:length(Zlim)),function(x){
	clux<-CLUall2[,x]
	list(sapply(sort(setdiff(clux,NA)),function(x){
		mean(subset(KEGGres[,5], clux==x))	
	}))
})

###Recalculate mean RF residual after K permutations of genes

k=1000
nrep=100

setwd(paste(pP,"perm-cluster",sep="/"))

mclapply(c(1:nrep),function(x){

X=x
Krf<-sapply(c(1:k),function(x){
	#print(x)
	rfx<-sample(KEGGres[,5])
	RFcluK<-sapply(c(1:length(Zlim)),function(x){
		clux<-CLUall2[,x]
		list(sapply(sort(setdiff(clux,NA)),function(x){
			mean(subset(rfx, clux==x))	
		}))	
	})
	ifelse(unlist(RFclu)>unlist(RFcluK),1,0)
})

Krf<-apply(Krf,1,sum)

write.table(Krf,paste("Perm_k=",k,"_X=",X,".txt",sep=""),row.names=F,col.names=F)

},mc.cores=8)

Krf<-apply(sapply(c(1: nrep),function(x){
	return(as.vector(read.table(paste("Perm_k=",k,"_X=",x,".txt",sep=""))[,1]))
}),1,sum)

#Pvalue for each level

LEVEL<-unlist(sapply(c(1:length(Zlim)),function(x){
	seq(x,x,length.out=length(unlist(RFclu[[x]])))
}))

K=nrep*k

pVAL<-sapply(c(1:length(Zlim)),function(x){
	X=x
	Kx<-subset(Krf, LEVEL==X)
	Mx= RFclu[[x]]
	
	Px<-(K-Kx+1)/(K+1)
	Padjx<-p.adjust(Px,method="BH")
	
	paste(CLUall2[,X],
	ifelse(Padjx[CLUall2[,X]]<=0.001,"***",
		ifelse(Padjx[CLUall2[,X]]<=0.01,"**",
		ifelse(Padjx[CLUall2[,X]]<=0.05,"*",
		ifelse(Px[CLUall2[,X]]<=0.001,"°°°",
		ifelse(Px[CLUall2[,X]]<=0.01,"°°",
		ifelse(Px[CLUall2[,X]]<=0.05,"°","")))))),sep="")	
	
	#plot(Mx, Padjx,log="y")
	
	#SIG<-subset(c(1:length(Mx)), Padjx<=0.05)
	#CLUsig<-sapply(SIG ,function(x){
	#	list(KEGGres[subset(c(1:length(CLUall2[,X])),CLUall2[,X]==x),c(1,2,5)])
	#})
	#names(CLUsig)<-	paste("TH=", Zlim[X],"_clu=", SIG,
	#	"_meanRF=",round(RFclu[[X]][SIG],4),
	#	"_pvalBH=",round(subset(Padjx ,Padjx<=0.05),4),sep="")
	
})

rownames(pVAL)<-KEGGres[,1]
colnames(pVAL)<-paste("L", Zlim,sep="")

setwd(pP)

write.table(cbind(KEGGres, pVAL),"clustering-genes-lineages-pVal.txt")


##STEP C - For each possible pair of gene, calculate the CCL value (in the range 1-24, this will be used as input for cytoscape, each node being a core gene; each edge being the CCL value) #This measure should be weighted by the number of genomes within species and the number of species per lineages, and should take into account that some genes are missing in some genomes. CCL = 24: this pair of gene is conserved accross all lienages

####Compared pairwise RF (pRF) residual between genes (NOT BETWEEN GENE AND THE CONSENSUS PHYLOGENY) and proportion of species/lineage for wich each pair of gene is located winthin 25Kb




#All possible pairs between genes

ALLp<-c("0_0",unlist(sapply(c(1:(length(Nu)-1)),function(x){
	return(paste(x,c((x+1):length(Nu)),sep="_"))
})))


#Average pair occurence By lineage, weighted by species and by genomes within species - also report pRF (corrected RF distance) between genes

CLUav<-sapply(unique(ORD[,3]),function(x){
	LINx=x
	print(x)
	lx<-subset(c(1:length(ORD[,1])),ORD[,3]== LINx)
	
	CluSPEav<-sapply(lx,function(x){
	
		spex=ORD[x,1] #species name
		print(spex)
		stx<-unlist(strsplit(ORD[x,5],split="_")) #short genome names
		#absolute abundance of each pair in this species
		clux<-unlist(sapply(stx,function(x){
			#index of this genome in the data file
			ix<-subset(c(1:length(data[,1])),data[,2]==x)
			#pairs of genes located within 25Kb window in this genomes
			CluF[[subset(c(1:length(iMETH)), iMETH==ix)]]
		}))
		#remove pairs of the same gene
		clux<-matrix(unlist(strsplit(clux,split="_")),ncol=2,byrow=T)
		clux<-subset(clux ,clux[,1]!=clux[,2])
		clux<-paste(clux[,1],clux[,2],sep="_")
		clux<-summary(as.factor(clux),maxsum=length(unique(clux)))

	
	
		#number of times each gene is missing in genomes of this species
		misx<-as.vector(unlist(sapply(stx,function(x){
			NX<-subset(data[,6],data[,2]==x)
			ANNOx<-read.table(paste(pD,"RAST_files",NX,sep="/"))	
			#Occurence each gene in this genome
			z<-sapply(c(1:length(Nu)),function(x){
				length(subset(c(1:length(ANNOx[,1])),ANNOx[,8]==Nu[x]))
			})
			subset(c(1:length(Nu)),z==0)
		})))
		misx <-c(summary(as.factor(misx),maxsum=length(unique(misx))))
	#Calculate the number of time each pair of gene could be in theory observable, takin into account missing gene
	
		TOTx<-setdiff(c("0_0=1",unlist(sapply(c(1:length(misx)),function(x){
			paste(ALLp[grep(paste("_",names(misx)[x],"_",sep=""),
			paste("_", ALLp,"_",sep=""))], misx[x],sep="=")
		}))),"=NA")
		TOTx<-matrix(unlist(strsplit(TOTx,split="=")),ncol=2,byrow=T) 
		TOTx<-sapply(unique(TOTx[,1]),function(x){
			return(mean(subset(as.numeric(TOTx[,2]), TOTx[,1]==x)))
		})
		TOTo<-seq(0,0,length.out=length(setdiff(ALLp,names(TOTx))))
		names(TOTo)<-setdiff(ALLp,names(TOTx))
		TOTx<-c(TOTo,TOTx)
		rm(TOTo)
		TOTx<-length(stx)-TOTx
	
		#pairs that are not observed in this species
		cluo<-seq(0,0,length.out=length(setdiff(ALLp,names(clux))))
		names(cluo)<-setdiff(ALLp,names(clux))
	
		clux<-c(clux,cluo)
		rm(cluo)
	
		#relative abudance of each pair in this species, taking into accout impossible pairs
		clur<-clux[order(names(clux))]/TOTx[order(names(TOTx))]
	
	
	})

	return(apply(CluSPEav,1,mean,na.rm=T))

})


setwd(pP)

write.table(CLUav,"AverageCLU-per-LIN.txt")

#############


pD<-"/myFolder/2-Core-genome/"

setwd(pD)

ORD<-read.table("Taxonomy.txt",header=T)
data<-read.table("data.txt",header=T)
STAT<-read.table("Stat-genomes.txt",header=T)


pP<-"/myFolder/3-Phylogeny/"

setwd(pP)

KEGGres<-read.table("Summary-KEGG-RF-Pos-per-PEG.txt",header=T)


#Retrieve average colocalization file per core gene per
CLUav<-apply(read.table("AverageCLU-per-LIN.txt",header=T),1,mean,na.rm=T)
CLUav<-subset(CLUav,names(CLUav)!="0_0")

#pRFgene file was generated in script 3-Phylogeny.R
pRFgene<-read.table("pDistanceRF.txt",header=T)

#Number of KEGGs 
pKEGG<-cbind(paste(gsub("CG_","", pRFgene[,1]),
	gsub("CG_","", pRFgene[,2]),sep="_"), pRFgene[,6])

#RF distance
pRFgene<-cbind(paste(gsub("CG_","", pRFgene[,1]),
	gsub("CG_","", pRFgene[,2]),sep="_"), pRFgene[,5])
	
np<-pRFgene[,1]
pRFgene<-as.numeric(pRFgene[,2])
names(pRFgene)<-np
rm(np)

np<-pKEGG[,1]
pKEGG <-as.numeric(pKEGG[,2])
names(pKEGG)<-np
rm(np)

#Convert relative abundance of CLUav a the number of lineages
#CLUav<-round(length(unique(ORD[,3]))* CLUav)

CLUav<-cbind(CLUav[order(names(CLUav))] ,pRFgene[order(names(pRFgene))], pKEGG[order(names(pKEGG))])
colnames(CLUav)<-c("CluLin","RFres","pKEGG")

####For Cytoscape: network according to the average conservation. of gene pairs (remove rare pairs)

minLin=1 #minimum average number of Lineage in which the pair is conserved (filter rare colocalizations)


CLUavF<-subset(CLUav, round(length(unique(ORD[,3]))* CLUav[,1])> minLin)

CLUavF<-cbind(
	KEGGres[as.numeric(matrix(unlist(strsplit(rownames(CLUavF),split="_")),
	ncol=2,byrow=T)[,1]),1],
	KEGGres[as.numeric(matrix(unlist(strsplit(rownames(CLUavF),split="_")),
	ncol=2,byrow=T)[,2]),1], round(length(unique(ORD[,3]))* CLUavF[,1]))
	
colnames(CLUavF)<-c("Gene1","Gene2","Lin")

write.table(CLUavF,paste("Pairs-CG-minLin=", minLin,".txt",sep=""),
	row.names=F,quote=F)

##Display core genes with cytoscape (input file: Pairs-CG-minLin=XXX.txt), linking genes according to the number of times they are colocolizing in genomes and color them according to RF dist (value of each gene, not pRF!).

##STEP D - Compare CCL with pRF (the more frequently two genes are found in the same cluster, the more they are expected to coevolve) - test for regression

###############Previsualizing data

CLUavo<-as.data.frame(CLUav)
CLUavo<-as.data.frame(subset(CLUav, CLUav[,1]!=0))

pK<-ifelse(CLUavo[,3]==0,0,1)

plot(as.factor(CLUavo[,3]), CLUavo[,1])
plot(as.factor(CLUavo[,3]), CLUavo[,2])


plot(CLUavo[,2], CLUavo[,1], pch=19,cex=0.2,col=rgb(0,0,0,0.05),log="")
plot(as.factor(CLUavo[,3]), CLUavo[,2], pch=19,cex=0.2,col=rgb(0,0,0,0.05),log="")
plot(as.factor(CLUavo[,3]), CLUavo[,1], pch=19,cex=0.2,col=rgb(0,0,0,0.05),log="")


#############

library(ggplot2)

df <- data.frame(
  x = round(CLUavo[,1],3),
  y = round(CLUavo[,2],3)
)

#Non-linear model
ggplot(df, aes(x, y)) +
  geom_point(alpha = 0.1, color =rgb(0,0,0),cex=0.2) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cr"), color = "red", size = 1.2) +
  labs(title = "Non-linear regression (Splines)",
       subtitle = "Modèle asymétrique robuste sans erreur de convergence") +
  theme_minimal()


# Model 1 : simple linear model (H0 : slope is constant)
mod_lin <- lm(y ~ x, data = df)

# Model 2 : GAM (H1 : slight curve)
#mod_gam <- gam(y ~ s(x), data = df)

# Comparison by AIC (the smallest, the best)
#AIC(mod_lin, mod_gam)


# Slope test
res <- summary(mod_lin)
print(res$coefficients)

library(dplyr)

# 1. Binning (to avoid overepresentation of rare pairs)
# Divide data in 20 segments (0.05 width)
df_bins <- df %>%
  mutate(bin = cut(x, breaks = seq(0, 1, 0.02), include.lowest = TRUE)) %>%
  group_by(bin) %>%
  summarise(
    x_mean = mean(x, na.rm = TRUE),
    y_mean = mean(y, na.rm = TRUE),
    n = n()
  ) %>%
  filter(!is.na(bin))

stats_brutes <- summary(mod_lin)

r2_val <- round(stats_brutes$r.squared, 4)        # R² (ex: 0.003)
p_val  <- stats_brutes$coefficients[2, 4]       # slope p-value 
p_text <- ifelse(p_val < 0.001, "< 0.001", round(p_val, 3))

print(paste("R2 =", r2_val, "| P-value =", p_text))

label_stats <- paste0("R² = ", r2_val, "\np-value ", p_text)


setwd(pP)

# 2. disoplay raw data

pdf("Pairwise_RF-vs_colocalization-raw.pdf",width=5,height=10)


ggplot() +
  # raw data
  geom_point(data = df, aes(x = x, y = y), alpha = 0.05, color = "grey", size = 0.5) +

#mean per bin
  geom_point(data = df_bins, aes(x = x_mean, y = y_mean, size = n), 
             color = "black", alpha = 0.8,shape=1) +
 geom_smooth(data = df, aes(x = x, y = y, color = "Linear regression (global)"), 
              method = "lm", se = T, size = 1,fill="black") +
            
  geom_smooth(data = df_bins, aes(x = x_mean, y = y_mean, color = "Linear regression (bins)"), 
              method = "lm", se = T, linetype = "dashed", size = 1,fill="grey") +
  # Configuration esthétique
  scale_color_manual(name = "", values = c("Linear regression (global)" = "red", 
                                                 "Linear regression (bins)" = "blue")) +
  scale_size_continuous(name = "Pairs per bin",breaks = c(200,2000,20000)) +
  theme_minimal() +
  labs(title = "Gene colocalization and pairwise phylogenetic dissimilarity",
       subtitle = "113,169 pairs; 1,099 core genes",
       x = "Gene colocalization within 25Kb (prop. of lineages)",
       y = "Pairwise phylogenetic dissimilarity (pRFr)")
      
dev.off()


#Bins only
pdf("Pairwise_RF-vs_colocalization-bins.pdf",width=4.5,height=5.5)

ggplot() +
  geom_point(data = df_bins, aes(x = x_mean, y = y_mean, size = n), 
             color = "black", alpha = 0.8,shape=1) +
             
  geom_smooth(data = df, aes(x = x, y = y, color = "Global"), 
              method = "lm", se = T, size = 1,fill="black") +
            
  geom_smooth(data = df_bins, aes(x = x_mean, y = y_mean, color = "Bins"), 
              method = "lm", se = T, linetype = "dashed", size = 1,fill="grey") +
  

  scale_color_manual(name = "Linear regression", values = c("Global" = "red", 
                                                 "Bins" = "blue")) +
  scale_size_continuous(name = "Pairs per bin",breaks = c(200,20000)) +
  theme_minimal() +
  labs(title = "Gene colocalization and pairwise phylogenetic dissimilarity",
       subtitle = "113,169 pairs; 1,099 core genes",
       x = "Conservation of gene colocalization\n among lineages (within 25Kb)",
       y = "Pairwise phylogenetic dissimilarity between genes (pRFr)")
      


dev.off()


#Statistics and pvalues for slopes

summary(lm(y_mean ~ x_mean, data = df_bins))
summary(lm(y ~ x, data = df))

###STEP E : calculate average + SD ANI and ACC for each group, lineage and species, compare ANI and ACC (model and regression)

setwd(pP)

ANI<-read.table("Average-ANI-per-species.txt",header=T)
ACC<-read.table("Average-cluster-cons-per-species.txt",header=T)

#number of core genes with data for phylogeny per genome
library(seqinr)

setwd(paste(pD,"core-genes",sep="/"))
CGx<-sapply(KEGGres[,1],function(x){
	names(read.fasta(paste(x,"_alignedNF.fas",sep="")))
})
CGx<-summary(as.factor(unlist(CGx)),maxsum=length(unique(unlist(CGx))))

#Statistics per species (mean among genomes within species)
SPEstat<-t(sapply(c(1:length(ORD[,1])),function(x){
	#Species
	sx=ORD[x,1]
	#Group
	Gx=ORD[x,2]
	#Lineage
	lx=ORD[x,3]
	#strains
	stx<-unlist(strsplit(ORD[x,5],split="_"))
	#number of strains
	nx=length(stx)
	#Average ANI
	anix<-ANI[x,x]
	#Average ACC
	accx<-ACC[x,x]	
	#average genome size, number of contigs, N50 and GC content
	gx<-apply(matrix(as.numeric(sapply(c(stx,stx),function(x){
		subset(STAT[,c(2:5)],STAT[,1]==x)
	})),ncol=4,byrow=T),2,mean)
	#average number of core genes
	cx<-mean(sapply(c(stx),function(x){
		subset(CGx,names(CGx)==x)
	}))
	c(sx,Gx,lx,nx,round(anix,3),round(accx,3),
		round(gx[1:3],0),round(gx[4],3),round(cx,3))
}))

colnames(SPEstat)<-c("Species","Group","Lineage","Strains","ANI","ACC","Size","contigs","N50","GC","Core-genes")

setwd(pP)
write.table(SPEstat,"Stat-per-species.txt",row.names=F)

#Statistics per lineage (mean among species within lineages)
LINstatMEAN <-t(sapply(unique(ORD[,3]),function(x){
	#Lineage
	lx=x	
	#Group
	Gx=unique(subset(ORD[,2],ORD[,3]==x))
	#species
	SPEx<-subset(ORD[,1],ORD[,3]==x)
	#number of species
	nx=length(SPEx)	
	#average and sd genome size, number of contigs, N50, GC content and number of core genes
	gx<-matrix(as.numeric(sapply(c(SPEx, SPEx),function(x){
		subset(SPEstat[,c(7:11)], SPEstat[,1]==x)
	})),ncol=5,byrow=T)
	gx<-round(apply(gx,2,mean),3)
	#Index of species within this lineage
	ix<-subset(c(1:length(ORD[,1])),ORD[,3]==x)
	#Average and sd ANI (among species)
	anix<-as.numeric(sapply(ix,function(x){
		ANI[x, setdiff(ix,x)]
	}))
	anix <-round(mean(anix),3)
	#Average and sd ACC (among species)
	accx<-as.numeric(sapply(ix,function(x){
		ACC[x, setdiff(ix,x)]
	}))
	accx <-round(mean(accx),3)
	
	c(lx,Gx,nx,anix,accx,gx)	
}))	
	
colnames(LINstatMEAN)<-c("Lineage","Group","Species","ANI","ACC","Size","contigs","N50","GC","Core-genes")

setwd(pP)
write.table(LINstatMEAN,"Stat-per-lineage.txt",row.names=F)

#### Average ANI and ACC for each pairwise comparison among lineages 

ANI.ACClin<-matrix(unlist(sapply(c(1:(length(unique(ORD[,3]))-1)),function(x){
	X=x
	#Lineage
	lx=unique(ORD[,3])[x]
	#Group
	gx<-unique(subset(ORD[,2],ORD[,3]==lx))
	#Index of species within this lineage
	ix<-subset(c(1:length(ORD[,1])),ORD[,3]==lx)
	sapply(c((X+1):(length(unique(ORD[,3])))),function(x){
		Y=x
		#Lineage
		ly=unique(ORD[,3])[x]
		#Group
		gy<-unique(subset(ORD[,2],ORD[,3]==ly))
		#Index of species within this lineage
		iy<-subset(c(1:length(ORD[,1])),ORD[,3]==ly)
		t(c(gx,gy,lx,ly,round(mean(as.numeric(unlist(ANI[ix, iy]))),3),
			round(mean(as.numeric(unlist(ACC[ix, iy]))),3)))
	})
})),ncol=6,byrow=T)	

colnames(ANI.ACClin)<-c("Group","Group","Lineage","Lineage","ANI","ACC")
ANI.ACClin<-rbind(LINstatMEAN[,c(2,2,1,1,4,5)], ANI.ACClin)

setwd(pP)

write.table(ANI.ACClin,"Mean_ANI_ACC_among_lineage.txt",row.names=F)

####

pdf("ANI_vs_ACC.pdf",width=5,height=5)

x1=as.numeric(unlist(ANI))
x2=as.numeric(unlist(ACC))

plot(x1, x2,
	las=1,cex.axis=0.8,xlab="ANI",ylab="ACC",
	pch=19,cex=0.4,col=rgb(0.9,0.9,0.9,0.5))
points(as.numeric(SPEstat[,5]), as.numeric(SPEstat[,6]),
	pch=19,cex=0.4,col="black")
	
points(as.numeric(ANI.ACClin[,5]), as.numeric(ANI.ACClin[,6]),
	pch=ifelse(ANI.ACClin[,3]==ANI.ACClin[,4],19,
	ifelse(ANI.ACClin[,1]==ANI.ACClin[,2],4,5)),
	col=ifelse(ANI.ACClin[,3]==ANI.ACClin[,4],NA,
		ifelse(ANI.ACClin[,1]!=ANI.ACClin[,2],rgb(0.6,0.6,0.6,0.7),NA)),
	cex=0.6,lwd=1.5
)

MOD<-lm(x2 ~ x1 + I(x1^2))

summary(MOD)

# Call:
# lm(formula = x2 ~ x1 + I(x1^2))
#
# Residuals:
#       Min        1Q    Median        3Q       Max 
# -0.310032 -0.059944 -0.009262  0.055070  0.289885 
#
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -2.051e+01  3.247e-01  -63.17   <2e-16 ***
# x1           4.419e-01  7.511e-03   58.84   <2e-16 ***
# I(x1^2)     -2.271e-03  4.331e-05  -52.44   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.08944 on 13453 degrees of freedom
# Multiple R-squared:  0.7139,	Adjusted R-squared:  0.7139 
# F-statistic: 1.679e+04 on 2 and 13453 DF,  p-value: < 2.2e-16

# ACC = -20.509175 + 0.441930.ANI - 0.002271.ANI2

grille <- data.frame(x1 = seq(min(x1), max(x1), length.out = 200))

grille$x2_pred <- predict(MOD, newdata = grille, type = "response")

lines(grille$x1, grille$x2_pred, col = "black", lwd = 1.5,lty=2)	
		
points(as.numeric(ANI.ACClin[,5]), as.numeric(ANI.ACClin[,6]),
	pch=ifelse(ANI.ACClin[,3]==ANI.ACClin[,4],19,
	ifelse(ANI.ACClin[,1]==ANI.ACClin[,2],4,5)),
	col=ifelse(ANI.ACClin[,3]==ANI.ACClin[,4],"orange",
		ifelse(ANI.ACClin[,1]!=ANI.ACClin[,2],NA,
		ifelse(ANI.ACClin[,1]=="A","red",
		ifelse(ANI.ACClin[,1]=="B","purple2",
		ifelse(ANI.ACClin[,1]=="D","blue",
		ifelse(ANI.ACClin[,1]=="C","green4","white")))))),
	cex=0.8,lwd=1.5
)

legend("bottomright",cex=0.7, horiz=F,
	legend= c("Among species",
		"Within species",
		"lm(ACC ~ ANI + I(ANI^2))",
		"Within lineages",
		"Within group A",
		"Within group B",
		"Within group D",
		"Within group C",
		"Among groups"),
	pt.cex=c(0.4,0.4,0,0.8,0.8,0.8,0.8,0.8,		0.6),
	pt.lwd=1.5,
	text.col= "black",box.col=NA,
	pch=c(19,19,0,19,4,4,4,4,5),
	lty=c(0,0,2,0,0,0,0,0,0),	
	col=c(rgb(0.9,0.9,0.9,0.5),
		"black","black","orange","red","purple2","blue","green4",
		rgb(0.6,0.6,0.6,0.7)),
	border=NA,
	title="",
	title.col="black",ncol=1)	
	
dev.off()




