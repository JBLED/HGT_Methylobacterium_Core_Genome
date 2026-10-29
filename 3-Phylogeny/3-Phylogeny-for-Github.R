###FUNCTIONS

#package for parallelization
library(parallel)

#Homemade functions
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


#Package seqinr to play with nucleotide sequences
library(seqinr)

#Package to align sequences
library(msa)

#Number of cores to use
NC=8

pD<-"myFolder/2-Core-genome/"
pP<-"myFolder/3-Phylogeny/"

###STEP A # Phylogenies : first, run each core gene separatly in RAxML; second: perform a consensus tree (lineage tree) with ASTRAL (combine core gene trees); third, same with a species tree

###Run each core gene phylogeny separately with RAxML (do it on a server)


###		#!/bin/bash
###		#SBATCH --mem=128GB
###		#SBATCH -N 1
###		#SBATCH --time=2-00:00:00
###		#SBATCH -p batch_72h
###		#SBATCH --ntasks-per-node=4
###		#SBATCH --array=1-1122%8

###		module load StdEnv/2020  gcc/9.3.0 
###		module load raxml/8.2.12

###		raxmlHPC-MPI -T 4 -f a -s CG_${SLURM_ARRAY_TASK_ID}_alignedNF.fas -x 12345 -p 12345 -# 1000 -m GTRGAMMA -n CG_${SLURM_ARRAY_TASK_ID}-GTRGAMMA-1000


#Reopen each tree in R and format them for ASTRAL analysis

library(ape)
library(ips)

setwd(pD)

STAT<-read.table("Stat-genomes.txt",header=T)
data<-read.table("data.txt",header=T)
ANNOa<-read.table("Annotation-occurence-raw.txt",header=T)
RC2<-read.table("Methylobacterium-core-gene-occ.txt",header=T)
CGstat<-read.table("Core-gene-code.txt",header=T)
CGname<-CGstat[,1]
CGtype<-CGstat[,3]
sumNEW<-read.table("Summary-per-core-gene-after-alignment.txt",header=T)

#New tree labels
LABnew<-paste(data[,2], data[,4], data[,5],sep="-")

mBS=10 #minimum node support to keep node

LX=length(CGstat[,1])

setwd(paste(pP,"RAxML-trees",sep="/"))

#######Multi gene tree (input for ASTRAL) to estimate the lineage tree with all core genes

TREE<-sapply(c(1:LX),function(x){
	print(x)
	#read tree as character
	treex<-as.vector(read.table(paste(
		"RAxML_bipartitionsBranchLabels.CG_",x,"-GTRGAMMA-1000",sep="")))
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
		paste("RAxML_tree.",x,"_reformated",sep=""),
		row.names=F,col.names=F,quote=F)
	
	treex<-read.tree(paste("RAxML_tree.",
		x,"_reformated",sep=""))
		
	#replace tree tip labels	
	#TIP<-treex$tip.label
	#nTIP<-as.vector(sapply(TIP,function(x){
	#	return(subset(LABnew ,data[,2]==x))
	#}))		
	#treex$tip.label <-nTIP
	
	treex <-collapseUnsupportedEdges(treex, 
		value = "node.label", cutoff= mBS)	
	write.tree(treex,paste("RAxMLtree-CG_",
		x,"-reformated_mBS=",mBS,sep=""))
	return(list(unroot(treex)))
})

class(TREE) <- "multiPhylo"

setwd(pP)

write.tree(TREE,paste("Methylobacterium-Multi-coregene-sept-2025_mBS=",mBS,".tre",sep=""))

#######Multi gene tree (input for ASTRAL) to estimate the lineage tree only with protein-encoding core genes (PEG)
#Version without rna genes (no same model of evolution, too short and low number of polymorphic sites)

setwd(paste(pP,"RAxML-trees",sep="/"))

TREEpeg<-sapply(setdiff(c(1:LX),subset(sumNEW[,1], sumNEW[,6]=="rna")),function(x){
	print(x)
	#read tree as character
	treex<-as.vector(read.table(paste(
		"RAxML_bipartitionsBranchLabels.CG_",x,"-GTRGAMMA-1000",sep="")))
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
		paste("RAxML_tree.",x,"_reformated",sep=""),
		row.names=F,col.names=F,quote=F)
	
	treex<-read.tree(paste("RAxML_tree.",
		x,"_reformated",sep=""))
		
	#replace tree tip labels	
	#replace tree tip labels	
	#TIP<-treex$tip.label
	#nTIP<-as.vector(sapply(TIP,function(x){
	#	return(subset(LABnew ,data[,2]==x))
	#}))		
	#treex$tip.label <-nTIP
	
	treex <-collapseUnsupportedEdges(treex, 
		value = "node.label", cutoff= mBS)	
	write.tree(treex,paste("RAxMLtree-CG_",
		x,"-reformated_mBS=",mBS,sep=""))
	return(list(unroot(treex)))
})

class(TREEpeg) <- "multiPhylo"

setwd(pP)

write.tree(TREEpeg,paste("Methylobacterium-Multi-coregene-PEG-sept-2025_mBS=",mBS,".tre",sep=""))


#######Multi gene tree (input for ASTRAL) to estimate the species tree only with protein-encoding core genes (PEG)
#Version without rna genes AND BY REPLACING DIRECTLY STRAIN NAME BY SPECIES NAME (this is an alternative way to force ASTRAL to group several sequences from the same species; because the species map file does work with ASTRAL 3 but not ASTRAL 4)

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

setwd(paste(pP,"RAxML-trees",sep="/"))

TREEpegSPE<-sapply(setdiff(c(1:LX),subset(sumNEW[,1], sumNEW[,6]=="rna")),function(x){
	print(x)
	#read tree as character
	treex<-as.vector(read.table(paste(
		"RAxML_bipartitionsBranchLabels.CG_",x,"-GTRGAMMA-1000",sep="")))
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
		paste("RAxML_tree.",x,"_reformated",sep=""),
		row.names=F,col.names=F,quote=F)
	
	treex<-read.tree(paste("RAxML_tree.",
		x,"_reformated",sep=""))
		
	#replace tree tip labels	
	#TIP<-treex$tip.label
	#nTIP<-as.vector(sapply(TIP,function(x){
	#	return(subset(LABnew ,data[,2]==x))
	#}))		
	#treex$tip.label <-nTIP

	#replace tree tip labels	
	TIP<-treex$tip.label
	nTIP<-as.vector(sapply(TIP,function(x){
		return(subset(MAPspe[,1] , MAPspe[,2]==x))
	}))		
	treex$tip.label <-nTIP	
	
	treex <-collapseUnsupportedEdges(treex, 
		value = "node.label", cutoff= mBS)	
	write.tree(treex,paste("RAxMLtree-CG_",
		x,"-reformated_mBS=",mBS,sep=""))
	return(list(unroot(treex)))
})

class(TREEpegSPE) <- "multiPhylo"

setwd(pP)

write.tree(TREEpegSPE,paste("Methylobacterium-Multi-coregene-PEG-SPE-sept-2025_mBS=",mBS,".tre",sep=""))


#STEP B - COMPUTE LINEAGE AND SPECIES TREE WITH ASTRAL 4 IMPLEMENTED IN ASTER (https://doi.org/10.1093/molbev/msaf172; different version than used in Leducq et al 2022 GBE paper; because this new version estimates terminal branch lengths)

# 1 lineage tree with all genes
# bin/astral4 -i ../Methylobacterium-Multi-coregene-sept-2025_mBS=10.tre -o ../Methylobacterium-AsterLinTree-sept-2025_mBS=10.tre

# 2 lineage tree only with protein-encoding genes (PEG)
# bin/astral4 -i ../Methylobacterium-Multi-coregene-PEG-sept-2025_mBS=10.tre -o ../Methylobacterium-AsterLinTree-PEG-sept-2025_mBS=10.tre

# 3 Species tree only with protein-encoding genes (PEG)
# bin/astral4 -i ../Methylobacterium-Multi-coregene-PEG-SPE-sept-2025_mBS=10.tre -o ../Methylobacterium-AsterSpeTree-PEG-SPE-sept-2025_mBS=10.tre

#STEP C - Calculate RF distance between each PEG core gene phylogeny and the lineage consensus tree built from PEGs. For each gene, report statistics that will we used to correct RF distances:   number of polymorphic sites per gene alignment (nPOL), dN/dS (dNdS), number of tips in the gene tree (nTIP), and median nodal support (Mns),

#Calculate RF distance with RF.dist function in phangorn package: do it separately for each gene tree by pruning the lineage tree (remove any tip absent in each gene tree) - Also return the median branch support for each gene tree

library(phangorn)

setwd(pD)

STAT<-read.table("Stat-genomes.txt",header=T)
data<-read.table("data.txt",header=T)
ANNOa<-read.table("Annotation-occurence-raw.txt",header=T)
RC2<-read.table("Methylobacterium-core-gene-occ.txt",header=T)
CGstat<-read.table("Core-gene-code.txt",header=T)
CGname<-CGstat[,1]
CGtype<-CGstat[,3]
sumNEW<-read.table("Summary-per-core-gene-after-alignment.txt",header=T)

#Location of core genes onto the reference genome
MAP<-read.table("MAP-core-genes-CBMB20.txt",header=T)

#Number of poymorphic sites and dNdS per gene
STATn<-read.table("CoreGenes_nPOL_dNdS.txt",header=T)

#Consensus Lineage tree (ASTRAL) only with PEG

setwd(pP)
LINtree<-read.tree("Methylobacterium-AsterLinTree-PEG-sept-2025_mBS=10.tre")


#Retrieve RAxML tree for each PEG core gene and calcultate RF distance and median nodal support

Lx<-c(1:length(CGstat[,1]))

#Index of peg
iP<-subset(Lx,CGstat[,3]=="peg")

setwd(paste(pP,"RAxML-trees",sep="/"))
	
STATtree<-t(sapply(CGstat[iP,1],function(x){
	print(x)

	treex<-read.tree(paste("RAxML_tree.",gsub("CG_","",x),"_reformated",sep=""))
	
	#treex<-unroot(read.tree(paste(
	#	"RAxML_bipartitionsBranchLabels.",x,"-GTRGAMMA-1000",sep="")))
	
	stx<-treex$tip #tips in this tree
	
	#in the lineage consensus tree, remove the missing tips
	treel10<-unroot(keep.tip(LINtree,tip= stx))
	
	#RF distance between the gene tree and the consensus lineage tree
	RF10x<-RF.dist(treex, treel10)
	
	#Median nodal support
	Mnsx<-median(as.numeric(treex$node),na.rm=T)
	
	#Retrieve the number of polymorphic sites 
	nPOLx<-as.numeric(subset(STATn$nVarWgap, STATn[,1]==x))
	#Retrieve the median dNdS for this gene
	dNdSx<-as.numeric(subset(STATn$MEDdNdS, STATn[,1]==x))
	
	return(c(length(stx),RF10x, Mnsx, nPOLx, dNdSx))

}))

STATtree<-as.data.frame(STATtree)

colnames(STATtree)<-c("TIPs","RF10","Mns","nPOL","dNdS")	

setwd(pP)

pdf("RFdist-Correction.pdf",height=8,width=4)

par(mar=c(4,4,1,1),mfcol=c(4,2))

#RF distance (log transformation to meet normality)
RF<-STATtree$RF10
#hist(RF,breaks=100)
RFlog<-log10(RF)
#hist(RFlog,breaks=100)
LogRF<-log10(STATtree$RF10)

###Model on RF distances for each scenario according to 
#dNdS per gene (normal distribution)
dNdS <-STATtree$dNdS
#hist(dNdS,breaks=100)
plot(RF, dNdS,pch=19,cex=0.5,col=rgb(0,0,0,0.2),las=1,cex.axis=0.8,log="")

#Number of polymorphic sites per gene #(log transformation to meet normality)
nPOL<-STATtree$nPOL 
#hist(nPOL,breaks=100)
nPOLlog<-log10(nPOL)
#hist(nPOLlog,breaks=100)
plot(RF, nPOL,pch=19,cex=0.5,col=rgb(0,0,0,0.2),las=1,cex.axis=0.8,log="")

#median nodal support per gene
Mns<-STATtree$Mns
#hist(Mns,breaks=100)
nMns<-Mns^2
nMns <-nMns/max(nMns)
#hist(nMns,breaks=100)
plot(RF, Mns,pch=19,cex=0.5,col=rgb(0,0,0,0.2),las=1,cex.axis=0.8,log="")

#Number of genomes per gene
TIP<-STATtree $TIPs #number of tips in the tree
#hist(TIP,breaks=100)
nTIP<-log10(max(TIP)-TIP +1) #1 - log transformation to be closer to normality
#hist(nTIP,breaks=100)
plot(RF, TIP,pch=19,cex=0.5,col=rgb(0,0,0,0.2),las=1,cex.axis=0.8,log="")


#Start with the most complex model: all factors, all interactions
MOD10<-lm(RFlog ~ Mns* dNdS* nPOL* TIP)

anova(MOD10) # remove not significant factors (p>0.1)

MOD11<-lm(RFlog ~ Mns + dNdS + nPOL + TIP + Mns:dNdS + Mns:nPOL + dNdS:nPOL + nPOL:TIP + Mns:dNdS:nPOL + Mns:nPOL:TIP +dNdS:nPOL:TIP)

anova(MOD11)

MOD12<-lm(RFlog ~ Mns + dNdS + nPOL + TIP + Mns:dNdS + Mns:nPOL + dNdS:nPOL + Mns:dNdS:nPOL + Mns:nPOL:TIP +dNdS:nPOL:TIP)

anova(MOD12)

MOD13<-lm(RFlog ~ Mns + dNdS + nPOL + TIP + Mns:dNdS + Mns:nPOL + dNdS:nPOL + Mns:dNdS:nPOL+dNdS:nPOL:TIP)

anova(MOD13)

anova(MOD10,MOD11,MOD12, MOD13)

anova(MOD10, MOD13)

c(AIC(MOD10),AIC(MOD11),AIC(MOD12),AIC(MOD13))

#Keep model MOD13: extract residuals ("modified" RF values)

RES10<-MOD13 $residuals

#hist(RES10,breaks=100)

plot(RES10, dNdS,pch=19,cex=0.5,col=rgb(0,0,0,0.2),las=1,cex.axis=0.8,log="")
plot(RES10, nPOL,pch=19,cex=0.5,col=rgb(0,0,0,0.2),las=1,cex.axis=0.8,log="")
plot(RES10, Mns,pch=19,cex=0.5,col=rgb(0,0,0,0.2),las=1,cex.axis=0.8,log="")
plot(RES10, TIP,pch=19,cex=0.5,col=rgb(0,0,0,0.2),las=1,cex.axis=0.8,log="")

dev.off()

write.table(anova(MOD13),"ANOVA-RFdist-Correction.txt")

STATtree<-cbind(STATtree, RES10)

colnames(STATtree)<-c("TIPs","RF10","Mns","nPOL","dNdS","Residuals")	

setwd(pP)

write.table(STATtree,"GeneTreeStat.txt")

