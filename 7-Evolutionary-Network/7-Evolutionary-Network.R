####################################
###STEP A: For each protein-encoding core gene, report (1) If groups A, B, C and D are monophyletic; (2) If lineages are monophyletic; (3) The proportion of species that are monophyletic. At the same time, draw each gene tree (color tips according to Groups/Lineages), as well as the species and the lineage tree 

pD<-"/myFolder/2-Core-genome/"
pP<-"/myFolder/3-Phylogeny/"

library(ape)
library(ips)

setwd(pD)

#taxonomy table of strains (Groups, Species and Lineages)
TAXA<-read.table("Taxonomy.txt",header=T)
data<-read.table("data.txt",header=T)

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

#Number of species and strains per lineage

COLlin <-cbind(COLlin ,t(sapply(COLlin[,1],function(x){
	c(length(grep(x,TAXA[,3])),
		length(unlist(strsplit(subset(TAXA[,5],TAXA[,3]==x),split="_")))
		)
})))

#Check colors

Z1<-sample(c(1:length(COLlin[,1])))
Z2<-sample(c(1:length(COLlin[,1])))

plot(Z1,Z2,
	col=COLlin[,2],pch=19,cex=4)
text(Z1,Z2,
	COLlin[,1],cex=0.8,font=2,col="white")	

#Species tree

setwd(pP)

SPEtree<-read.tree("Methylobacterium-AsterSpeTree-PEG-SPE-sept-2025_mBS=10.tre")

#Color for tips (according to lineage)
COLx<-sapply(SPEtree $tip,function(x){
	LINx<-TAXA [grep(x, TAXA[,1]),3]
	return(subset(COLlin[,2], 
		COLlin[,1]== LINx))
})

#Outgroup (group C)
OUT<-setdiff(TAXA[,1],TAXA[grep("C",TAXA[,2]),1])
#OUT<-TAXA[grep("C",TAXA[,2]),1]

setwd(paste(pP,"PlotTrees",sep="/"))

pdf("Species_tree.pdf")

plot.phylo(root(SPEtree,outgroup= OUT),cex=0.5,font=2,
		tip.color= COLx,type="radial",
		show.node.label=F,cex.main=1,
		main= "Species tree", 
		lab4ut="axial", open.angle=1,
		node.color=rgb(0.7,0.7,0.7),
		edge.color=rgb(0.7,0.7,0.7))

legend("bottomleft",cex=0.4, horiz=F,
	legend= paste(COLlin[,1]," (",COLlin[,3],")",sep=""),
	pt.cex=1,
	text.col= "black",box.col=NA,
	pch=19,	col=COLlin[,2],
	border=NA,
	title="Lineages",
	title.col="black",ncol=4)

plot.phylo(root(SPEtree,outgroup= OUT),cex=0.5,font=2,
		tip.color= COLx,type="fan",
		show.node.label=F,cex.main=1,
		main= "Species tree", 
		lab4ut="axial", open.angle=10,
		node.color=rgb(0.7,0.7,0.7),
		edge.color=rgb(0.7,0.7,0.7))

axisPhylo(side = 1,cex=0.5)

legend("bottomleft",cex=0.4, horiz=F,
	legend= paste(COLlin[,1]," (",COLlin[,3],")",sep=""),
	pt.cex=1,	
	text.col= "black",box.col=NA,
	pch=19,	col=COLlin[,2],
	border=NA,
	title="Lineages",
	title.col="black",ncol=4)


dev.off()

#Lineage tree

setwd(pP)

LINtree<-read.tree("Methylobacterium-AsterLinTree-PEG-sept-2025_mBS=10.tre")

#Color for tips (according to lineage)
COLx<-sapply(LINtree $tip,function(x){
	LINx<-TAXA [grep(paste("_",x,"_",sep=""),
		paste("_", TAXA[,5],"_",sep="")),3]
	return(subset(COLlin[,2], 
		COLlin[,1]== LINx))
})

#Outgroup (group C; actually ABD)
OUT<-setdiff(unlist(strsplit(TAXA[,5],split="_")),
	unlist(strsplit(
	TAXA[grep("C",TAXA[,2]),5],split="_")))

#Root the tree on C (actually ABD) and collapse unsupported nodes
LINtreeR<-root(collapseUnsupportedEdges(
	LINtree,value = "node.label", 
	cutoff= 0.95),outgroup= OUT)

#real strain name
STn<-sapply(LINtreeR $tip,function(x){
	subset(data[,3],data[,2]==x)
})

#real strain name + species
STs<-sapply(LINtreeR $tip,function(x){
	paste(subset(data[,3],data[,2]==x),
	TAXA [grep(paste("_",x,"_",sep=""),
		paste("_", TAXA[,5],"_",sep="")),1],sep=" - ")
})

LINtreeR$tip.label<-STs

setwd(paste(pP,"PlotTrees",sep="/"))

pdf("Lineage_tree.pdf")

plot.phylo(LINtreeR,cex=0.2,
		tip.color= COLx,type="radial",
		show.node.label=F,cex.main=1,
		main= "Lineage tree", 
		lab4ut="axial",font=2, 	
		node.color=rgb(0.7,0.7,0.7),
		edge.color=rgb(0.7,0.7,0.7))

legend("bottomleft",cex=0.4, horiz=F,
	legend= paste(COLlin[,1]," (",COLlin[,4],")",sep=""),
	pt.cex=1,
	text.col= "black",box.col=NA,
	pch=19,	col=COLlin[,2],
	border=NA,
	title="Lineages",
	title.col="black",ncol=4)

plot.phylo(LINtreeR,cex=0.2,
		tip.color= COLx,type="fan",
		show.node.label=F,cex.main=1,
		main= "Lineage tree", 
		lab4ut="axial",font=2, 	
		open.angle=20,
		node.color=rgb(0.7,0.7,0.7),
		edge.color=rgb(0.7,0.7,0.7))

axisPhylo(side = 1,cex=0.5)

legend("bottomleft",cex=0.4, horiz=F,
	legend= paste(COLlin[,1]," (",COLlin[,4],")",sep=""),
	pt.cex=1,
	text.col= "black",box.col=NA,
	pch=19,	col=COLlin[,2],
	border=NA,
	title="Lineages",
	title.col="black",ncol=4)


dev.off()
		
#Each gene tree (also report (1) If groups A, B, C and D are monophyletic; (2) If lineages are monophyletic; (3) The proportion of species that are monophyletic.)

setwd(pD)

CGstat<-read.table("Core-gene-code.txt",header=T)
CGname<-CGstat[,1]
CGtype<-CGstat[,3]
sumNEW<-read.table("Summary-per-core-gene-after-alignment.txt",header=T)

setwd(pP)
KEGGres<-read.table("Summary-KEGG-RF-Pos-per-PEG.txt",header=T)		
STATtree<-read.table("GeneTreeStat.txt")

#######Multi gene tree (input for ASTRAL) to estimate the lineage tree with all core genes


#Cutoff on bootstrapp support 
mBS=10

#Outgroup (group C; actually ABD)
OUT<-setdiff(unlist(strsplit(TAXA[,5],split="_")),
	unlist(strsplit(
TAXA[grep("C",TAXA[,2]),5],split="_")))

#Code for PEG core genes
CG<-subset(CGname, CGtype=="peg")

#Optional: only trees with the 50 top RF residual values
#t(sapply(rownames(STATtree[order(STATtree[,6],decreasing=T),])[1:50],function(x){

#Optional: only trees with the 100 top ratio between max branch length et mean branch length
#t(sapply(rownames(STATtree[order(as.numeric(MONO[,33])/as.numeric(MONO[,31]),decreasing=T),])[1:100],function(x){

MONO<-t(sapply(CG,function(x){
	print(x)
	
	Fx=subset(CGstat[,2], CGstat[,1]==x)
	
	setwd(paste(pP,"RAxML-trees",sep="/"))
	treex<-read.tree(paste("RAxML_tree.",
		gsub("CG_","",x),"_reformated",sep=""))

	#treex<-collapseUnsupportedEdges(treex, 
	#	value = "node.label", cutoff= mBS)

	#Color for tips (according to lineage)
	COLx<-sapply(treex$tip,function(x){
		LINx<-TAXA [grep(paste("_",x,"_",sep=""),
			paste("_", TAXA[,5],"_",sep="")),3]
		return(subset(COLlin[,2], 
			COLlin[,1]== LINx))
	})
		
	#Root the tree on C (actually ABD) 
	#treeR<-root(treex ,
	#	outgroup= intersect(OUT,treex$tip))
	
	#tree tips
	STn<-treex $tip

	#real strain name + species
	#STs<-sapply(STn,function(x){
	#	paste(subset(data[,3],data[,2]==x),
	#	TAXA [grep(paste("_",x,"_",sep=""),
	#		paste("_", TAXA[,5],"_",
	#		sep="")),1],sep=" - ")
	#})

	#treex $tip.label<-STs
	
	#Stats for this gene
	Sx<-paste(paste(colnames(STATtree),
		subset(round(STATtree,3),row.names(STATtree)==x),sep=" = "),collapse=" ; ")
	#KEGGs for this gene
	Kx<-paste("KEGG =",subset(KEGGres[,6], KEGGres [,1]==x))

	setwd(paste(pP,"PlotTrees",sep="/"))		
	
	pdf(paste(x,"-tree.pdf",sep=""))
	
	par(mar=c(4,4,4,4),bty="n")
	
	plot.phylo(treex,cex=0.5,
		tip.color= COLx,type="unroot",
		show.node.label=F,cex.main=1,
		main= paste(x,"\n",Fx,"\n",Sx,"\n",Kx,sep=""),
		cex.main=0.6, 
		lab4ut="axial",font=2, 	
		open.angle=20,
		node.color=rgb(0.7,0.7,0.7),
		edge.color=rgb(0.7,0.7,0.7))
	
	#Number of tips and bootstrapp support for each group and lineage (if monophyletic)
	
	MONOg<-t(sapply(c("A","B","D","C"),function(x){
		stx<-intersect(STn,unlist(strsplit(subset(
			TAXA[,5],TAXA[,2]==x),split="_")))
		c(length(stx),is.monophyletic(treex, stx))
	}))
	
	MONOl<-t(sapply(COLlin[,1],function(x){
		stx<-intersect(STn,unlist(strsplit(subset(
			TAXA[,5],TAXA[,3]==x),split="_")))
		z<-length(stx)
		c(z,ifelse(z!=0,is.monophyletic(treex, stx),NA))
	}))	
	#axisPhylo(side = 1,cex=0.5)
	
	par(mar=c(1,1,1,1),bty="n",new=T)
	plot(-10,-10,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1))
	
	legend("bottom",cex=0.4, horiz=F,
		legend= paste(c("A","B","D","C",COLlin[,1])," (", c(MONOg[,1],MONOl[,1]),")",sep=""),
		pt.cex=1,
		text.col= "black",box.col=NA,
		pch=ifelse(c(MONOg[,2],MONOl[,2])==0,1,19),	
		col=c("red","purple","blue","green",COLlin[,2]),
		border=NA,
		title="Groups/Lineages",
		title.col="black",ncol=14)	
		
	dev.off()
	
	
	#median, average, SD and max branch length
	BRANCH<-round(c(median(treex[[2]]),mean(treex[[2]]),
		sd(treex[[2]]),max(treex[[2]])),3)
	names(BRANCH)<-c("median","mean","sd","max")
	return(c(MONOg[,2], MONOl[,2],BRANCH))
}))	

setwd(pP)

write.table(MONO,"Monophyly+BranchLength-per-PEG.txt")		

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### STEP B ### Methylobacterium evolutionary network: reconstruction on the network by combining gene tree topologies simplified at the lineage level (gene swaps prior lineage radiation are not considered). In the network, each node will correspond to the consensus tree internal node, and each edge will correspond to pairwise combinations among those nodes (partitions), which frequency accross core gene trees will be calculated (edge width proportionnal to partition frequency) - only gene trees that have at have sequences for at least 23 out of the 24 Methylobacterium lineages are considered.

library(ape)
library(ips)

#Retrieve groups of lineages forming monophyletic groups for each gene 
pD<-"/myFolder/2-Core-genome/"
pP<-"/myFolder/3-Phylogeny/"

setwd(pD)
#taxonomy table of strains (Groups, Species and Lineages)
TAXA<-read.table("Taxonomy.txt",header=T)
data<-read.table("data.txt",header=T)


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

setwd(pP)
MONO<-read.table("Monophyly+BranchLength-per-PEG.txt",header=T)
LINtree<-read.tree("Methylobacterium-AsterLinTree-PEG-sept-2025_mBS=10.tre")
KEGGres<-read.table("Summary-KEGG-RF-Pos-per-PEG.txt",header=T)


#MONOc<-subset(rownames(MONO),MONO[,4]==1)
MONOc<-rownames(MONO)

#number of missing lineages per gene 
z<-sapply(c(1:length(MONO[,1])),function(x){
	z<-unlist(MONO[x,c(5:28)])
	return(length(subset(z,is.na(z)==T)))
})
#Number of missing genes per lineage
w<-sapply(c(5:28),function(x){
	z<-unlist(MONO[,x])
	return(length(subset(z,is.na(z)==T)))
})


#Only consider genes that have one missing lineage
MONO1<-subset(rownames(MONO),z==1)

#Partitions of gene trees (genes that have one missing lineage)
PARTLINs1<-sapply(MONO1,function(x){
	print(x)
		
	setwd(paste(pP,"RAxML-trees",sep="/"))
	treex<-read.tree(paste("RAxML_tree.",
		gsub("CG_","",x),"_reformated",sep=""))	
	#tree tips
	STn<-treex $tip
	
	#Lineages that are monophyletic for this tree
	MONOx<-subset(MONO,rownames(MONO)==x)
	LINm<-intersect(COLlin[,1],
		subset(colnames(MONO),as.vector(unlist(MONOx))==1))
	
	#Root tree on the first monophyletic lineage and return partitions
	treex<-root(treex,intersect(STn,unlist(strsplit(
		subset(TAXA[,5],TAXA[,3]== LINm[1]),split="_"))))
	#Retrieve monophyletic groups with complete lineages
	partitions <- prop.part(treex)
	lins<-sapply(STn,function(x){
		return(TAXA[grep(paste("_",x,"_",sep=""),
			paste("_",TAXA[,5],"_",sep="")),3])
	})
	partlins<-sapply(c(1:length(partitions)),function(x){
		return(list(lins[partitions[[x]]]))
	})
	lins<-summary(as.factor(lins),maxsum=length(unique(lins)))
	partlins<-sapply(c(1:length(partlins)),function(x){
		linx= partlins[[x]]
		paste(setdiff(sapply(c(1:length(lins)),function(x){
			ifelse(length(grep(names(lins[x]), linx))== lins[x],names(lins[x]),'')
		}),''),collapse="_")
	})
	partlins<-setdiff(partlins,'')
	#For each group, generate the miror clade and only return the smallest 
	partlins2<-as.vector(sapply(partlins,function(x){
		return(paste(sort(setdiff(LINm,
			unlist(strsplit(x,split="_")))),collapse="_"))
	}))
	
	partlins<-unique(ifelse(nchar(partlins)<=nchar(partlins2), partlins, partlins2))
	
	#Missing lineage for this gene
	MISS<-setdiff(COLlin [,1],
		unlist(strsplit(partlins,split="_")))
	
	#Artificially ad those lineages to partitions containing the closest lineage according to the consensus tree 
	setwd(pP)
	CONS<-read.table("consensus.txt",header=T)
	#Grep the frequently associated lineages in the consensus tree 
	COMP<-summary(as.factor(unlist(strsplit(
		CONS[grep(MISS, CONS[,1]),1],split="_"))),
		maxsum=length(COLlin[,1]))
	COMP<-subset(COMP,names(COMP)!= MISS)
	COMP<-sort(names(subset(COMP, COMP==max(COMP))))
	
	partlins<-unlist(c(gsub(
		paste(COMP,collapse="_"),
		paste(sort(c(MISS,COMP)),collapse="_"), 
		partlins), paste(COMP,collapse="_")))
	
	partlins <-unique(as.vector(sapply(c(1:length(partlins)),function(x){
		return(paste(sort(			unlist(strsplit(partlins[x],split="_"))),collapse="_"))
	})))
		
	return(partlins)
		
})	

#Only consider genes that have no missing lineage
MONOc<-subset(rownames(MONO),z==0)

#Partitions of gene trees (genes that have no missing lineage)
PARTLINs<-sapply(MONOc,function(x){
	print(x)
		
	setwd(paste(pP,"RAxML-trees",sep="/"))
	treex<-read.tree(paste("RAxML_tree.",
		gsub("CG_","",x),"_reformated",sep=""))	
	#tree tips
	STn<-treex $tip
	
	#Lineages that are monophyletic for this tree
	MONOx<-subset(MONO,rownames(MONO)==x)
	LINm<-intersect(COLlin[,1],
		subset(colnames(MONO),as.vector(unlist(MONOx))==1))
	
	#Root on the first monophyletic lineage and return partitions

	treex<-root(treex,intersect(STn,unlist(strsplit(
		subset(TAXA[,5],TAXA[,3]== LINm[1]),split="_"))))
	#Retrieve monophyletic groups with complete lineages
	partitions <- prop.part(treex)
	lins<-sapply(STn,function(x){
		return(TAXA[grep(paste("_",x,"_",sep=""),
			paste("_",TAXA[,5],"_",sep="")),3])
	})
	partlins<-sapply(c(1:length(partitions)),function(x){
		return(list(lins[partitions[[x]]]))
	})
	lins<-summary(as.factor(lins),maxsum=length(unique(lins)))
	partlins<-sapply(c(1:length(partlins)),function(x){
		linx= partlins[[x]]
		paste(setdiff(sapply(c(1:length(lins)),function(x){
			ifelse(length(grep(names(lins[x]), linx))== lins[x],names(lins[x]),'')
		}),''),collapse="_")
	})
	partlins<-setdiff(partlins,'')
	#For each group, generate the miror clade and only return the smallest 
	partlins2<-as.vector(sapply(partlins,function(x){
		return(paste(sort(setdiff(LINm,
			unlist(strsplit(x,split="_")))),collapse="_"))
	}))
	
	partlins<-unique(ifelse(nchar(partlins)<=nchar(partlins2), partlins, partlins2))
		
	partlins <-unique(as.vector(sapply(c(1:length(partlins)),function(x){
		return(paste(sort(			unlist(strsplit(partlins[x],split="_"))),collapse="_"))
	})))
		
	return(partlins)
		
})	

###Trees with more than one missing lineages will not be considered for this analysis (too complicated scenarios to consider)


####To retrieve gene names for a particular partition
t(t(subset(MONOc,sapply(c(1:length(PARTLINs)),function(x){
	return(length(subset(PARTLINs[[x]], PARTLINs[[x]]=="A7_D1_D2_D3_D4_D5_D6")))
})==1)))

t(t(subset(MONO1,sapply(c(1:length(PARTLINs1)),function(x){
	return(length(subset(PARTLINs1[[x]], PARTLINs1[[x]]=="A7_D1_D2_D3_D4_D5_D6")))
})==1)))

#### Pool partition

MONOc<-c(MONOc,MONO1)
PARTLINs<-c(PARTLINs, PARTLINs1)


PARTLINs2<-summary(as.factor(unlist(PARTLINs)),maxsum=length(unique(unlist(PARTLINs))))

####### Manually generate partition of the consensus tree (simplified at the lineage level) and corresponding node names (for the network): file "nodes.txt"

#First, import the list of nodes in the network (node names and corresponding partitions)

setwd(pP)

NODE<-read.table("nodes.txt",header=T) #first column: node name; secund column: corresponding partition in the consensus tree

#Second, import the structure of the consensus tree manually simplified at the lineage level

PARTcons<-read.table("consensus.txt",header=T) #first column: name of the consensus partition; secund column: corresponding daughter node; third column: corresponding mother node #For tips of the consensus tree (lineages), counts will correspond to the number of PEG tree in which the lineage appear as monophyletic (return 0 for lienages with only one representative) 

##for each posible pairwise combination of nodes in the consensus tree (combination of nodes from the NODE file), return the resulting partition and determine whether it corresponds to a partition of the consensus tree (in the PARTcons file) or an alternative partition

PART<-matrix(unlist(strsplit(unlist(sapply(c(1:(length(NODE[,1])-1)),function(x){
	X=x
	n1= NODE[x,1] #name of this node
	p1= NODE[x,2] #partition of this node
	sapply(c((X+1):(length(NODE[,1]))),function(x){
		Y=x
		n2= NODE[x,1] #name of this node
		p2= NODE[x,2] #partition of this node
		
		return(paste(c(n1, n2,
			paste(sort(unique(unlist(strsplit(
			c(p1,p2),split="_")))),
			collapse="_")),collapse="="))
	})
})),split="=")),ncol=3,byrow=T)

colnames(PART)<-c("NODE1","NODE2","PART.EDGE")

###For each partition (group of lineages that are monophyletic), generate the miror partition and only return the smallest group
PART2<-as.vector(sapply(PART[,3],function(x){
	return(paste(sort(setdiff(COLlin[,1],
		unlist(strsplit(x,split="_")))),collapse="_"))
}))

PART2 <-ifelse(nchar(PART[,3])<=nchar(PART2), PART[,3], PART2)

#Remove partitions that correspond to the consensus tree

z<-sapply(c(1:length(PART[,1])),function(x){
	return(length(intersect(PART2[x],PARTcons[,1])))
})

PART <-subset(PART,z==0)
PART2<-subset(PART2,z==0)
PART <-subset(PART, PART2 !="")
PART2 <-subset(PART2, PART2 !="")

#Some alternative partitions correspond to different scenarios (swaps of deep nodes); mostly because the network in unrooted. In those cases, just consider them twice 

DUP<-sapply(c(1:length(PART2)),function(x){
	return(length(subset(PART2,PART2== PART2[x])))
}) #return 2 if the partition has an alternative scenario

#Index of genes that were conserved for this analysis
iMONOc<-sapply(MONOc,function(x){
	return(subset(c(1:length(MONO[,1])),rownames(MONO)==x))
})

#Count the occurence of consensus partition in gene trees - #For lineages, return total of observed monophyly (NA excluded) in MONO 
CONScount<-sapply(PARTcons[,1],function(x){
	z<-subset(PARTLINs2, names(PARTLINs2)==x)
	z2=sum(MONO[iMONOc,subset(c(1:length(MONO[1,])),colnames(MONO)==x)],na.rm=T)
	return(max(c(z,z2)))
})

names(CONScount)<-PARTcons[,1]

#Count the occurence of alternative partitions in gene trees
ALTcount<-sapply(PART2,function(x){
	z<-c(subset(PARTLINs2, names(PARTLINs2)==x),0)[1]
})
names(ALTcount)<-PART2

#Other partitions (that won't be displayed in the network because more complex history)

OTHER<-setdiff(names(PARTLINs2),c("",PART2,PARTcons[,1]))

OTHERcount<-sapply(OTHER,function(x){
	z<-c(subset(PARTLINs2, names(PARTLINs2)==x),0)[1]
})
names(OTHERcount)<-OTHER

#Summary file for partition that can be displayed in the network

PARTsum1<-cbind(PARTcons[,c(2,3)],CONScount,"Cons")
PARTsum2<-cbind(PART[,c(1,2)], ALTcount,"Alt")
	
colnames(PARTsum1)<-c("Node1","Node2","Occurrence","Type")
colnames(PARTsum2)<-c("Node1","Node2","Occurrence","Type")

PARTsum<-rbind(PARTsum1, PARTsum2)
PARTsumO<-subset(PARTsum, PARTsum[,3]=="0") #partition that were not observed
PARTsum<-subset(PARTsum, PARTsum[,3]!="0")


###Write a network file

setwd(pP)

write.table(PARTsum,"Partitions-for-Cyto_Sept2026-2.txt",quote=F)

write.table(PARTsum,"Partitions-for-Cyto_Sept2026.txt",row.names=F,quote=F) #input for CYTOSCAPE

write.table(subset(PARTsum,PARTsum[,4]=="Cons"),"Partitions-for-Cyto_CONS_Sept2026.txt",row.names=F,quote=F) #input for consensus partition only CYTOSCAPE

write.table(subset(PARTsum,PARTsum[,4]=="Alt"),"Partitions-for-Cyto_ALT_Sept2026.txt",row.names=F,quote=F) #input for alternative partition only CYTOSCAPE

write.table(t(t(OTHERcount)),"Other-Partitions_ALT_Sept2026.txt")
setwd(pP)

###Step C: statistics on partition categories
####Statistics per partition category: number of partitions, average, sd, min and max number of genes with this partition; 

fx=as.numeric(PARTsum[,3])

###Categories

#0 Unsolved partitions
N0=length(OTHERcount)
A0=mean(OTHERcount)
S0=sd(OTHERcount)
I0=min(OTHERcount)
M0=max(OTHERcount)

#1L Consensus partitions (lineages only)
i1L<-sapply(COLlin[,1],function(x){
	return(subset(c(1:length(PARTsum[,1])),
		rownames(PARTsum)==x))
})

N1L=length(i1L)
A1L=mean(fx[i1L])
S1L=sd(fx[i1L])
I1L=min(fx[i1L])
M1L=max(fx[i1L])

#1I Consensus partitions (internal nodes only)
i1I<-setdiff(grep("Cons", PARTsum[,4]), i1L)

N1I=length(i1I)
A1I=mean(fx[i1I])
S1I=sd(fx[i1I])
I1I=min(fx[i1I])
M1I=max(fx[i1I])

#2 Alternative partitions (all, remove duplicated partition)

i2<-setdiff(grep("Alt", PARTsum[,4]), 
	grep(".1",rownames(PARTsum),fixed=T)) #".1" for duplicated partition names 
	
N2=length(i2)
A2=mean(fx[i2])
S2=sd(fx[i2])
I2=min(fx[i2])
M2=max(fx[i2])

STAT<-cbind(c(0, 0, 0,length(PARTsumO[,1])),
	c(N0, N1L, N1I,N2),
	c(A0, A1L, A1I,A2),
	c(S0, S1L, S1I,S2),
	c(I0, I1L, I1I,I2),
	c(M0, M1L, M1I,M2))

colnames(STAT)<-c("non-observed","observed","Mean","SD","min","max")

rownames(STAT)<-c("other","lineages","consensus","alternative")
#In alternative partitions, identify partion among groups and within groups = remove alternative partitions involving nodes (A, B, C, D, AB and ABD) 

#observed alternative partitions not involving deep nodes 
i3<-setdiff(grep("Alt", PARTsum[,4]),unlist(
	sapply(c("A","B","C","D","AB","ABD"),function(x){
	return(c(subset(c(1:length(PARTsum[,1])), PARTsum[,1]==x),
		subset(c(1:length(PARTsum[,1])), PARTsum[,2]==x)))
})))

#non-observed alternative partitions not involving deep nodes
i3o<-setdiff(c(1:length(PARTsumO[,1])),unlist(
	sapply(c("A","B","C","D","AB","ABD"),function(x){
	return(c(subset(c(1:length(PARTsumO[,1])), PARTsumO[,1]==x),
		subset(c(1:length(PARTsumO[,1])), PARTsumO[,2]==x)))
})))



#stats of alternative partitions within and among groups

setwd(pP)

STATw<-t(sapply(c("A-A","B-B","D-D","C-C","A-B","A-D","B-D","A-C","B-C","C-D"),function(x){
	#groups
	gx<-unlist(strsplit(x,split="-"))
	#observed partitions among thoses groups
	ix<-intersect(i3, c(intersect(grep(gx[1],PARTsum[,1]),
		grep(gx[2],PARTsum[,2])),intersect(grep(gx[2],PARTsum[,1]),
		grep(gx[1],PARTsum[,2]))))
	#non-observed partitions among those groups
	iy<-intersect(i3o, c(intersect(grep(gx[1],PARTsumO[,1]),
		grep(gx[2],PARTsumO[,2])),intersect(grep(gx[2],PARTsumO[,1]),
		grep(gx[1],PARTsumO[,2]))))	
	Nx=length(ix)
	Ax=mean(fx[ix])
	Sx=sd(fx[ix])
	Ix=min(fx[ix])
	Mx=max(fx[ix])
	Ny=length(iy)
	return(c(Ny, Nx,Ax,Sx,Ix,Mx))
}))


TYPEw<-matrix(unlist(strsplit(unlist(sapply(c("A-A","B-B","D-D","C-C","A-B","A-D","B-D","A-C","B-C","C-D"),function(x){
	#groups
	gx<-unlist(strsplit(x,split="-"))
	#observed partitions among thoses groups
	ix<-intersect(i3, c(intersect(grep(gx[1],PARTsum[,1]),
		grep(gx[2],PARTsum[,2])),intersect(grep(gx[2],PARTsum[,1]),
		grep(gx[1],PARTsum[,2]))))
	return(paste(fx[ix],x,sep="="))
})),split="=")),ncol=2,byrow=T)

colnames(STATw)<-c("non-observed","observed","Mean","SD","min","max")

write.table(rbind(STAT, STATw),"Stat-network.txt")

pdf("Statistics-Networks.pdf")

barplot(STATw[,2]/apply(STATw[,c(1:2)],1,sum),
	ylim=c(0,1),las=1,cex.axis=0.8,
	col=c("red","purple","blue","green3",
		"cyan","cyan","cyan",
		"yellow3","yellow3","yellow3"),border=NA,
	ylab="proportion of expected partitions",
	names=paste(rownames(STATw),apply(STATw[,c(1:2)],1,sum),
		sep="\n"),
	cex.names=0.8,xlab="Type of partition (expected observations)")

plot(as.factor(TYPEw[,2]),as.numeric(TYPEw[,1]),log="y",
	ylab="# PEG tree (log scale)",las=1,cex.axis=0.8)


########################
####Tukey test on partition categories
 
 # install.packages(c("ggplot2", "dplyr", "forcats", "multcompView"))
 library(ggplot2)
 library(dplyr)
 #install.packages("forcats")
 library(forcats)       # Pour réordonner facilement les facteurs
 #install.packages("multcompView")
 library(multcompView)  # Pour générer automatiquement les lettres du test de Tukey
 
 TYPEt <-data.frame(CAT= factor(gsub("-","<=>",TYPEw[,2])),OCC=(as.numeric(TYPEw[,1])))
 
 TYPEt <- TYPEt %>%
   mutate(CAT = fct_reorder(CAT, OCC, .fun = median))
 
 modele_anova <- aov(OCC ~ CAT, data = TYPEt)
 test_tukey <- TukeyHSD(modele_anova)
 
 lettres_tukey <- multcompLetters4(modele_anova, test_tukey)
 df_lettres <- data.frame(lettres = lettres_tukey$CAT$Letters) %>%
   tibble::rownames_to_column("CAT")
 
 df_complet <- TYPEt %>%
   group_by(CAT) %>%
   summarise(ymax = max(OCC) + (max(OCC) * 0.05)) %>% # Décale la lettre de 5% au-dessus du max
   left_join(df_lettres, by = "CAT")
 
 ggplot(TYPEt, aes(x = CAT, y = OCC)) +
   geom_boxplot(fill = "skyblue", alpha = 0.7, outlier.alpha = 0.5) +
   # Ajout des lettres de Tukey (elles s'aligneront parfaitement sur l'axe) 
   scale_y_log10()+
   geom_text(data = df_complet, aes(x = CAT, y = ymax, label = lettres), 
             hjust = -0.2, color = "darkred", fontface = "bold", size = 3.5) +
   # Inversion des axes
   coord_flip() + 
   theme_minimal() +
   labs(
     title = "PEG tree occurence per partition category (Tukey test)",
     x = "Partition categories",
     y = "PEG trees (log scale)"
   )

#Quick test with 3 polled categories: 1) within groups, 2) between groups A, B D and 3) between group C and groups A, B and D

POOLt<-matrix(unlist(strsplit(as.vector(TYPEt[,1]),
	split="<=>")),ncol=2,byrow=T)
POOLt <-ifelse(POOLt[,1]== POOLt[,2],"Within",
	ifelse(POOLt[,1]=="C","C<=>X",
	ifelse(POOLt[,2]=="C","C<=>X","A<=>B<=>D")))
POOLtSUM <-t(sapply(unique(POOLt),function(x){
	z<-subset(as.numeric(TYPEt[,2]),POOLt==x)
	return(c(length(z),mean(z),sd(z),min(z),max(z)))
}))

modele_anovaPOOL <- aov(TYPEt$OCC ~ POOLt)
 test_tukeyPOOL <- TukeyHSD(modele_anovaPOOL)

################
####Proportion test on observed vs non-observed partitions 

library(ggplot2)
library(dplyr)
library(forcats)
library(rstatix) # Pour les tests par paires robustes sur matrice
library(cld)     # Pour générer les lettres de significativité automatiquement

data<-as.data.frame(STATw[,c(1,2)])
rownames(data)<-gsub("-","<=>",rownames(data))
colnames(data)<-c("nOBS","OBS")

test_global <- chisq.test(data)
print(test_global)

test_paires <- pairwise_prop_test(data, p.adjust.method = "holm")

lettres_prop <- make_cld(test_paires)

df_graph <- as.data.frame(data) %>%
  tibble::rownames_to_column("Classe") %>%
  mutate(
    Total = nOBS + OBS,
    Proportion_Observed = OBS / Total
  ) %>%
  # Liaison avec les lettres
  left_join(lettres_prop, by = c("Classe" = "group")) %>%
  # Tri des classes de façon définitive par leur proportion constatée
  mutate(Classe = fct_reorder(Classe, Proportion_Observed))


ggplot(df_graph, aes(x = Classe, y = Proportion_Observed)) +
  geom_col(fill = "steelblue", alpha = 0.8, width = 0.7) +
  # Placement de la lettre juste au-dessus du sommet de la barre
  geom_text(aes(label = cld, y = Proportion_Observed + 0.01), 
            hjust = 0, color = "darkred", fontface = "bold", size = 3) +
  coord_flip() + # Basculement à l'horizontale pour une lisibilité parfaite des axes
  theme_minimal() +
  labs(
    title = "Proportion of expected partition (proportion test)",
    x = "Type of partition",
    y = "Proportion of expected partitions"
  )

dev.off()


#Quick test with 3 polled categories: 1) within groups, 2) between groups A, B D and 3) between group C and groups A, B and D

POOL<-matrix(unlist(strsplit(rownames(data),
	split="<=>")),ncol=2,byrow=T)
POOL<-ifelse(POOL[,1]==POOL[,2],"Within",
	ifelse(POOL[,1]=="C","C<=>X",
	ifelse(POOL[,2]=="C","C<=>X","A<=>B<=>D")))
POOL<-sapply(unique(POOL),function(x){
	return(apply(subset(data,POOL==x),2,sum))
})

pairwise_prop_test(POOL, p.adjust.method = "holm")

###########STEP D - Permutation test of KEGGs among genes to identify partitions enriched in some pathways/functions (to display in the network)

KEGGres<-read.table("Summary-KEGG-RF-Pos-per-PEG.txt",header=T)		


PartKEGG<-unlist(sapply(c(1:length(MONOc)),function(x){
	#print(x)
	#KEGG for this gene
	kx<-subset(KEGGres[,6], KEGGres[,1]== MONOc[x])
	kx<-ifelse(kx=="","NA",kx)
	kx<-unlist(strsplit(kx,split="_"))
	#Partitions for this gene
	px<-PARTLINs[[x]]
	px<-setdiff(px,"")
	
	return(sapply(kx,function(x){
		paste(x,px,sep="=")
	}))
}))

PartKEGG<-summary(as.factor(PartKEGG),
	maxsum=length(unique(PartKEGG)))
	
PartKEGG<-cbind(matrix(unlist(
	strsplit(names(PartKEGG),split="=")),
	ncol=2,byrow=T),as.vector(PartKEGG))	

#COnvert into matrix; only keep partition also present in the summary file (PARTsum) 
PartKEGG<-sapply(sort(unique(PartKEGG[,1])),function(x){
	px<-subset(PartKEGG[,c(2,3)], PartKEGG[,1]==x)
	return(sapply(rownames(PARTsum),function(x){
		as.numeric(c(subset(px[,2],px[,1]==x),0)[1])
	}))
})

####Do K*k permutations of Keggs among genes and recalculate expected frequencies of each KEGG for each part

library(parallel)
NCORE=8

K= 100
k= 1000

setwd(paste(pP,"/Perm_Part_KEGG",sep=""))

#Relative abundance of PART

#subset(names(SUMpart),names(SUMpart)!=rownames(PartKEGG))###check names correspondance

#PartKEGGr<-PartKEGG/SUMpart

sapply(c(1:K),function(x){
	print(c(x,date()))
	Kx=x
	mclapply(c(1:k),function(x){
		MONOk<-sample(MONOc)
		PartKEGGk<-unlist(sapply(c(1:length(MONOk)),function(x){
			#print(x)
			#KEGG for this gene
			kx<-subset(KEGGres[,6], KEGGres[,1]== MONOk[x])
			kx<-ifelse(kx=="","NA",kx)
			kx<-unlist(strsplit(kx,split="_"))
			#Partitions for this gene
			px<-PARTLINs[[x]]
			px<-setdiff(px,"")
			return(sapply(kx,function(x){
				paste(x,px,sep="=")
			}))
		}))
		PartKEGGk <-summary(as.factor(PartKEGGk),
			maxsum=length(unique(PartKEGGk)))
		PartKEGGk <-cbind(matrix(unlist(
			strsplit(names(PartKEGGk),split="=")),
			ncol=2,byrow=T),as.vector(PartKEGGk))	
		PartKEGGk <-sapply(colnames(PartKEGG),function(x){
			px<-subset(PartKEGGk[,c(2,3)], PartKEGGk[,1]==x)
			return(sapply(rownames(PARTsum),function(x){
				as.numeric(c(subset(px[,2],px[,1]==x),0)[1])
			}))
		})
		PartKEGGk<-ifelse(PartKEGG> PartKEGGk,1,0)
		write.table(PartKEGGk,paste("perm_k=",x,sep=""),
			quote=F,row.names=F,col.names=F)
		print(max(PartKEGGk))
	},mc.cores= NCORE)
	setwd(paste(pP,"/Perm_Part_KEGG",sep=""))
	PERM<-PartKEGG*0
	write.table(PERM,paste("PartKEGG_",Kx,".txt",sep=""))
	sapply(c(1:k),function(x){
		Zx<-read.table(paste("PartKEGG_",Kx,".txt",sep=""),header=T)
		zx<-read.table(paste("perm_k=",x,sep=""))
		write.table(Zx+ zx,paste("PartKEGG_",Kx,".txt",sep=""))
		file.remove(paste("perm_k=",x,sep=""))
	})
	return("")
})	

setwd(paste(pP,"/Perm_Part_KEGG",sep=""))
PERM<-PartKEGG*0
write.table(PERM,paste("PartKEGG_K=",K*k,".txt",sep=""))

sapply(c(1:K),function(x){
	Zx<-read.table(paste("PartKEGG_",x,".txt",sep=""),header=T)
	zx<-read.table(paste("PartKEGG_K=",K*k,".txt",sep=""))
	write.table(Zx+ zx,paste("PartKEGG_K=",K*k,".txt",sep=""))
})


#For each PART, calculate P-value associated to each KEGG (number of times KEGG is observed in higher relative frequencies than expected)

PERM<-read.table(paste("PartKEGG_K=",K*k,".txt",sep=""),header=T)

#partition
#Pf<-subset(c(1:length(rownames(PartKEGG))),apply(PartKEGG,1,sum)>=5)
Pf<-c(1:length(PARTsum[,1]))

SUMx<-sapply(Pf,function(x){
	#print(x)
	X=x
	PARTx=rownames(PARTsum)[X]
	OBS<-PartKEGG[X,]
	Kx<-PERM[X,]
	
	#Occurence 
	Fpart<-subset(PARTsum[,3],rownames(PARTsum)== PARTx)
	
	Px<-(K*k)-Kx
	Px <-(Px +1)/((K*k)+1)
	Px <- p.adjust(Px,method="BH")
	Sx<-ifelse(Px <=0.001,"***",ifelse(Px <=0.01,"**",
		ifelse(Px <=0.05,"*",ifelse(Px <=0.1,"°",""))))
		
	t(subset(cbind(PARTx ,paste("K",colnames(PartKEGG),sep=""), Fpart, OBS, round(Px,4), Sx),Sx!="")	)
})

SUMx<-matrix(unlist(SUMx),ncol=6,byrow=T)
colnames(SUMx)<-c("Partition","KEGG","PartitionOcc","KEGGocc","AdjPval","Sig")

setwd(paste(pP,"/Perm_Part_KEGG",sep=""))

write.table(SUMx,paste("Sum_PartKEGG_K=",K*k,".txt",sep=""),row.names=F)

###return the list of genes for a particular partition

#Glutamine synthetase types I (typical of procaryotes) and II (typical of plants) both found in some hyphomicrobiales; case for Methylobacterium 

#Il all partitions: 
#CG_414 "CG_414" "Glutamine synthetase type I (EC 6.3.1.2)"                                     
#CG_415 "CG_415" "Glutamine synthetase type II, eukaryotic (EC 6.3.1.2)" 

#CG_412 ("Glutaminase (EC 3.5.1.2)"), CG_737 ("Phosphoribosylformylglycinamidine synthase, glutamine amidotransferase subunit (EC 6.3.5.3)"), CG_100 ("Adenylosuccinate lyase (EC 4.3.2.2) @ SAICAR lyase (EC 4.3.2.2)") only in partition "A7_B1_B2_B4" 

#Why do these specific bacteria keep both types?
#Bacteria that possess GSII—like Rhizobium, Frankia, and Streptomyces—almost always #have Type I (GSI) as well. Retaining both enzymes provides a distinct physiological #advantage: 
	•	#	•	Adapting to Extreme Nitrogen Shifts: These bacteria are heavily involved in complex nitrogen metabolisms (like nitrogen fixation or specialized soil decomposition). 
	•	#	•	Dual Control Systems: GSI has a highly sensitive, energy-expensive "on/off" switch regulated by adenylation. GSII lacks this specific covalent control. By utilizing GSI under certain metabolic conditions and GSII under others, the bacteria can fine-tune their ammonia assimilation depending on whether they are living freely in the soil or symbiotically inside a plant host.
# REF: 
# https://pubmed.ncbi.nlm.nih.gov/6107288/
# https://pmc.ncbi.nlm.nih.gov/articles/PMC7105534/
# https://pubmed.ncbi.nlm.nih.gov/2570059/

###"Glutamine synthetase"                   
partx="A5_D1_D2_D3_D4" 
partx="A6_C2"
partx="A7_B1_B2_B4" 

#BCHL
partx="A7_D2"
partx="A3_A5"
partx="A7_D2_D3_D4"
partx="A7_D1_D2_D3_D4_D5_D6"
partx="C1_C6"

#Propanoate
partx="A3_A7_B1_B2_B3_B4"

#Photosynthesis
partx="A4_B1"

CGx<-subset(names(PARTLINs),
	sapply(c(1:length(PARTLINs)),function(x){
	length(intersect(partx, PARTLINs[[x]]))	
})==1)

paste(unlist(t(sapply(CGx,function(x){
	return(subset(KEGGres, KEGGres[,1]==x))
}))[,1]),collapse=";")

#Generate a network file with partitions enriched in KEGGs

setwd(pP)
PARTsum<-read.table("Partitions-for-Cyto_Sept2026-2.txt",
	header=T)


##Remove partitions corresponding to lineages (false positives)

SUMxFILT<-subset(SUMx,nchar(SUMx[,1])>2)
SUMxFILT<-subset(SUMxFILT,as.numeric(SUMxFILT[,5])<=0.05)
#SUMxFILT<-subset(SUMxFILT, SUMxFILT[,2]!="K03010") ##Remove "K03010_Ribosome" (false positive?)

#Use threshold on pvalue: 0.05
PARTsig<-summary(as.factor(SUMxFILT[,1]))
KEGGsig<-summary(as.factor(SUMxFILT[,2]))

CYTOsign<-t(sapply(names(PARTsig),function(x){
	subset(PARTsum,rownames(PARTsum)==x)
}))

#Add corresponding KEGGs

TYPE<-sapply(rownames(CYTOsign),function(x){
	paste(sort(subset(SUMxFILT[,2], SUMxFILT[,1]==x)),collapse="_")
})

CYTOsign<-cbind(CYTOsign[,c(1:3)], TYPE)
colnames(CYTOsign)<-c("Node1","Node2","Occurrence","Type")

#Give a new name to partitions corresponding to lineages with only one species (not relevant for monophyly)

#LINu<-summary(as.factor(TAXA[,3]))
#LINu<-subset(names(LINu),LINu==1)

#CYTOlinu<-t(sapply(LINu,function(x){
#	subset(PARTsum,rownames(PARTsum)==x)
#}))

#CYTOlinu <-cbind(CYTOlinu[,c(1:3)],"Unique")
#colnames(CYTOlinu)<-c("Node1","Node2","Occurrence","Type")

#Remove CYTOlinu and CYTOsign from the original network file

#PARTsumN<-rbind(t(sapply(setdiff(rownames(PARTsum),
#	c(rownames(CYTOlinu),rownames(CYTOsign))),function(x){
#	subset(PARTsum,rownames(PARTsum)==x)	
#})), CYTOlinu, CYTOsign)

PARTsumN<-rbind(t(sapply(setdiff(rownames(PARTsum),
	rownames(CYTOsign)),function(x){
	subset(PARTsum,rownames(PARTsum)==x)	
})), CYTOsign)

###Write a new network file

setwd(pP)

write.table(PARTsumN,"Partitions-for-Cyto_Sept2026-New2.txt",quote=F)

write.table(PARTsumN,"Partitions-for-Cyto_Sept2026-New.txt",row.names=F,quote=F)

write.table(subset(PARTsumN,PARTsumN[,4]!="Cons"),"Partitions-for-Cyto_Sept2026-ALT_New.txt",row.names=F,quote=F)

write.table(subset(PARTsumN,PARTsumN[,4]=="Cons"),"Partitions-for-Cyto_Sept2026-CONS_New.txt",row.names=F,quote=F)

##Count the number of ALT partition between and among groups

###OBSERVED
PARTsumN12<-subset(PARTsumN, nchar(PARTsumN[,1])!=1)
PARTsumN12<-subset(PARTsumN12, nchar(PARTsumN12[,2])!=1)
PARTsumN12<-subset(PARTsumN12, PARTsumN12[,4]!="Cons")

PARTsumN1<-gsub("1","", PARTsumN12[,1])
PARTsumN1<-gsub("2","", PARTsumN1)
PARTsumN1<-gsub("3","", PARTsumN1)
PARTsumN1<-gsub("4","", PARTsumN1)
PARTsumN1<-gsub("5","", PARTsumN1)
PARTsumN1<-gsub("6","", PARTsumN1)
PARTsumN1<-gsub("7","", PARTsumN1)
PARTsumN2<-gsub("1","", PARTsumN12[,2])
PARTsumN2<-gsub("2","", PARTsumN2)
PARTsumN2<-gsub("3","", PARTsumN2)
PARTsumN2<-gsub("4","", PARTsumN2)
PARTsumN2<-gsub("5","", PARTsumN2)
PARTsumN2<-gsub("6","", PARTsumN2)
PARTsumN2<-gsub("7","", PARTsumN2)

PARTsumN12<-cbind(paste(PARTsumN1, PARTsumN2,sep="_"),PARTsumN12[,c(3:4)])

PARTsumN12<-subset(PARTsumN12, nchar(PARTsumN12[,1])==3)

PARTsumN12<-cbind(as.vector(sapply(as.vector(unlist(PARTsumN12[,1])),function(x){
	return(paste(sort(unlist(strsplit(x,split="_"))),collapse="_"))
})), PARTsumN12[,2])

###EXPECTED

PART12<-subset(PART, nchar(PART[,1])!=1)
PART12 <-subset(PART12, nchar(PART12[,2])!=1)

PART1 <-gsub("1","", PART12[,1])
PART1 <-gsub("2","", PART1)
PART1 <-gsub("3","", PART1)
PART1 <-gsub("4","", PART1)
PART1 <-gsub("5","", PART1)
PART1 <-gsub("6","", PART1)
PART1 <-gsub("7","", PART1)
PART2 <-gsub("1","", PART12[,2])
PART2 <-gsub("2","", PART2)
PART2 <-gsub("3","", PART2)
PART2 <-gsub("4","", PART2)
PART2 <-gsub("5","", PART2)
PART2 <-gsub("6","", PART2)
PART2 <-gsub("7","", PART2)

PART12 <-paste(PART1, PART2,sep="_")

PART12 <-subset(PART12, nchar(PART12)==3)

PART12 <-as.vector(sapply(as.vector(PART12),function(x){
	return(paste(sort(unlist(strsplit(x,split="_"))),collapse="_"))
}))


PARTsumN12<-t(sapply(sort(unique(as.vector(unlist(PARTsumN12[,1])))),function(x){
	obs<-as.vector(unlist((subset(PARTsumN12[,2], PARTsumN12[,1]==x))))
	exp<-length(subset(PART12, PART12==x))
	
	return(c(exp ,length(obs),sum(obs),sum(obs)/exp))
	
}))

colnames(PARTsumN12)<-c("MaxExpected","nObserved","Total","Total/MaxExpected")

barplot(t(PARTsumN12[c(1,5,10,8,2,4,7,3,6,9),c(1:2)]),beside=T,
	legend=T,las=1,cex.axis=0.8)

barplot(PARTsumN12[c(1,5,10,8,2,4,7,3,6,9),4],
	las=1,cex.axis=0.8)


