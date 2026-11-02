#install.packages("ggplot2", repos = "https://cloud.r-project.org")
#install.packages(c("patchwork", "ggpubr", "ape"), repos = "https://cloud.r-project.org")
#install.packages("rhierbaps")
library(ggplot2)
library(rhierbaps)
library(parallel)

#Retrieve the list of core genes
pP<-"/myFolder/3-Phylogeny/"
setwd(pP)
KEGGres<-read.table("Summary-KEGG-RF-Pos-per-PEG.txt",header=T)

#Output folder
pB<-"/myFolder/4-BAPS"

#Input folder
pD<-"/myFolder/2-Core-genome"

setwd(pD)

ORD<-read.table("Taxonomy.txt",header=T) #taxonomy and species order according to the species tree
ORD<-ORD[(order(ORD[,4])),]

ncore=8


#STEP A
#Classify strains for each gene in nLin lineages (when the number of lineages is defined): vary the number of lineages between 2 and nVAR
nVAR=60
nLEV=1

mclapply(KEGGres[setdiff(setdiff(c(1:190),Temp)),1],function(x){
  	X=x
	# Load the multiple sequence alignment
	snp_matrix <- load_fasta(paste(pD,"/core-genes/",
		paste(x,"_alignedNF.fas",sep=""),sep=""))		
	sapply(c(2: nVAR),function(x){
 	 	hb_results <- hierBAPS(snp_matrix, max.depth = nLEV, 
 	 		n.pops = x, quiet = TRUE,
  			assignment.probs=FALSE)    
  		# Extract the sequence classification table
	 	classification_table <- hb_results$partition.df	
	 	# Extract the posterior probabilities table
 		#ASS1<-hb_results$cluster.assignment.prob[[1]]
		setwd(paste(pB,"/LEVEL", nLEV,"_nLIN=", paste("2-",nVAR,sep=""),sep=""))
		write.table(classification_table,paste(X,"-rhierbaps_nPOP=",x,".txt",
			sep=""),row.names=F)
		#write.table(ASS1,paste(X,"-rhierbaps_nPOP=",x,"_prob.txt",
		#	sep=""),row.names=F)			
	})	
  		
  	CLASSall<-sapply(c(2: nVAR),function(x){
  		CLASS<-read.table(paste(X,"-rhierbaps_nPOP=",x,".txt",sep=""),header=T)
  		return(CLASS[order(CLASS[,1]),2])
  	})
  	CLASSo<-read.table(paste(X,"-rhierbaps_nPOP=",2,".txt",sep=""),header=T)
  	linx<-sapply(sort(CLASSo[,1]),function(x){
		return(ORD[grep(paste("_",x,"_",sep=""),
			paste("_",ORD[,5],"_",sep="")),3])
	})
	colnames(CLASSall)<-paste("B",c(2: nVAR),sep="")
	rownames(CLASSall)<-paste(sort(CLASSo[,1]),linx,sep="_")
	
	setwd(paste(pB,"/LEVEL", nLEV,"_nLIN=", paste("2-",nVAR,sep=""),sep=""))
	write.table(CLASSall,paste(X,"-rhierbaps_ALL.txt",sep=""))
	
  	
  	#pdf(paste(X,"-rhierbaps.pdf",sep=""),height=4,width=4)
  	#plot(c(2: nLIN),as.numeric(summary(apply(CLASSall,2,unique))[,1]),
  	#	xlab="Pred Npop",ylab="Post Npop",
  	#	cex.axis=0.8,las=1,type="l",col="blue")
  	#dev.off()
  	
  	sapply(c(2: nVAR),function(x){
  		file.remove(paste(X,"-rhierbaps_nPOP=",x,".txt",sep=""))
  	})
  	
 },mc.cores= ncore)
  
 ####Sumarize result: number of populations in function of number tested
setwd(paste(pB,"/LEVEL", nLEV,"_nLIN=", paste("2-",nVAR,sep=""),sep=""))
  	
  Temp<-c(1:190)  	
  	
pdf("Npop_Combined.pdf",height=4,width=4)	
 
par(mar=c(4,4,1,1),bty="n") 
  	
plot(-10,10,xlim=c(0, nVAR),ylim=c(0, nLIN),
 	xlab="n.pops",ylab="BAPS lineages",
  		cex.axis=0.8,las=1) 	
 
 Z<-sapply(Temp,function(x){
   	X=paste("CG",x,sep="_")
   	CLASSall<-read.table(paste(X,"-rhierbaps_ALL.txt",sep=""),header=T)
  	z<-	as.numeric(summary(apply(CLASSall,2,unique))[,1])
 	points(c(2: nVAR),z,
  		type="l",col=rgb(0,0,1,0.02),lwd=11,lend=1)
  	list(t(cbind(c(2: nVAR),z)))
 })
 
 Z <-matrix(unlist(Z),ncol=2,byrow=T)
 
sumZ<-t(sapply(sort(unique(Z[,1])),function(x){
	return(c(x,mean(subset(Z[,2],Z[,1]==x)),
	sd(subset(Z[,2],Z[,1]==x))))
})) 

points(sumZ[,1], sumZ[,2],type="l",col="black",lwd=3,lend=1) 

plot(as.factor(Z[,1]),Z[,2],ylim=c(0, nLIN),
	xlab="n.pops",ylab="BAPS lineages",
  		cex.axis=0.8,las=1,pch=19,cex=0.5) 	
 

  	
 dev.off() 	
  

#Classify strains for each gene in nLin lineages (correspond to the number of Methylobacterium lineages)


#Number of lineages
nLIN=length(unique(ORD[,3]))
#Number of levels
nLEV=2


ncore=8


mclapply(KEGGres[,1],function(x){

  	X=x
	# Load the multiple sequence alignment
	snp_matrix <- load_fasta(paste(pD,"/core-genes/",
		paste(x,"_alignedNF.fas",sep=""),sep=""))
	# Run the BAPS clustering model
  	# n.pops: number of known lineages
  	# max.depth: Number of hierarchical nesting levels (usually 2 is standard)
  	hb_results <- hierBAPS(snp_matrix, max.depth = nLEV, 
  		n.pops = nLIN, quiet = TRUE,
  		assignment.probs=FALSE)
 	# Extract the sequence classification table
 	classification_table <- hb_results$partition.df	
 	# Extract the posterior probabilities table (NO)
 	#ASS1<-hb_results$cluster.assignment.prob[[1]]
	setwd(paste(pB,"/LEVEL", nLEV,"_nLIN=", nLIN,sep=""))
	write.table(classification_table,paste(x,"-rhierbaps.txt",sep=""),
		row.names=F)

},mc.cores=ncore)

###############################################################
### Retrieve rhierbaps output for each gene and classify according to BAPS results


pD<-"/myFolder/2-Core-genome/"

pP<-"/myFolder/3-Phylogeny/"

# rhierbaps Output folder
pB<-"/myFolder/4-BAPS"

num.strsplit.u<-function(x){return(as.numeric(unlist(strsplit(x,split= "_"))))}

library(ape)
library(ips)

#Species tree order
setwd(pD)

ORD<-read.table("Taxonomy.txt",header=T) #taxonomy and species order according to the species tree
ORD<-ORD[(order(ORD[,4])),]

###Color codes for graphique
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

plot(c(1:length(COLlin[,1])),pch=19,cex=5,col=COLlin[,2])
text(c(1:length(COLlin[,1])),c(1:length(COLlin[,1])),COLlin[,1])

COLgroup<-t(cbind(
	c("A","red"),
	c("B","purple2"),
	c("C","green4"),
	c("D","blue")
))

plot(c(1:length(COLgroup[,1])),pch=19,cex=5,col= COLgroup[,2])
text(c(1:length(COLgroup[,1])),c(1:length(COLgroup[,1])), COLgroup[,1])


#color codes for hybrid lineages (orange and yellow scale)

COLhyb<-c("yellow2","orange2","yellow4","orange4")


#Consensus Lineage tree (ASTRAL) only with PEG

setwd(pP)
LINtree<-read.tree("Methylobacterium-AsterLinTree-PEG-sept-2025_mBS=10.tre")

#Core gene statistics

KEGGres<-read.table("Summary-KEGG-RF-Pos-per-PEG.txt",header=T)

#Statistics for a subset of genes (nPOP vary)

nVAR=60
nLEV=1
Temp<-c(1:190)


VATstat<-sapply(KEGGres[Temp,1],function(x){
  	X=x
  	print(X)
   	setwd(paste(pB,"/LEVEL", nLEV,"_nLIN=", paste("2-",nVAR,sep=""),sep=""))
   	linx<-read.table(paste(X,"-rhierbaps_ALL.txt",sep=""),header=T)
   	linSP<-rownames(linx)
   	linSP<-matrix(unlist(strsplit(linSP,split="_")),ncol=2,byrow=T)
	#Gene code
	CGx= X
	#phylogenetic tree for this gene
	setwd(paste(pP,"RAxML-trees",sep="/"))
	treex<-read.tree(paste("RAxML_tree.",
		gsub("CG_","", CGx),"_reformated",sep="")) 
		
	#in the lineage consensus tree, only keep tips from this tree
	treel10<-unroot(keep.tip(LINtree,tip= treex$tip))
	
	Z<-sapply(c(2: nVAR),function(x){
		NX=x
		liny<-linx[,(NX-1)]

	#Convert it so it fit the species tree (for each species, report the BAP lineage(s) )	
		LNx<-sapply(c(1:length(ORD[,4])),function(x){
			stx<-unlist(strsplit(ORD[x,5],split="_"))	
			paste(c(unique(sapply(stx,function(x){
				lx<-c(subset(liny, linSP[,1]==x),NA)[1]
				return(lx)
			})),NA),collapse="_")
		})
		
		LNx<-gsub("_NA","",LNx)
		LNx<-gsub("NA_","",LNx)

		#For each BAP lineage, return: 1) corresponding Methylo groups; 2) the corresponding methylo group with the larger number of species in this BAP lineage; 3) corresponding Methylo lineages; 3) the corresponding methylo lineage with the larger number of species in this BAP lineage; 4) Whether this BAP is monophyletic or not in the corresponding tree; 5) is it monophyletic in the consensus lineage tree
		LINbap<-t(sapply(sort(unique(liny)),function(x){
			stx<-subset(linSP[,1],liny==x) #strains assigned to this LINbap 
			z<-unique(sapply(stx,function(x){
				grep(paste("_",x,"_",sep=""),
					paste("_",ORD[,5],"_",sep=""))
			}))
			#corresponding Methylo groups (number of species)
			gx<-sort(summary(as.factor(ORD[z,2])),decreasing=T)
			#corresponding Methylo lineages (number of species)
			lx<-sort(summary(as.factor(ORD[z,3])),decreasing=T)
			gy<-paste(names(gx),collapse="_")
			ly<-paste(names(lx),collapse="_")
		
			return(c(x,length(stx),sum(lx),names(gx)[1],
				names(lx)[1],gy,ly,
				is.monophyletic(unroot(treex), stx), #is this BAP monophyletic in the gene tree?
				is.monophyletic(unroot(treel10), stx) #is this BAP monophyletic in the lineage consensus tree?
			))
		}))
		
		colnames(LINbap)<-c("BAPS","strains","species",
			"Main_group","Main_lineage","Groups",
			"Lineages","MonoGeneTree?","MonoLinTree?")

		STATx<-ifelse(LINbap[,8]==F,"ignore",   #this BAP is not monophyletic in the gene tree: ignore
			ifelse(LINbap[,9]==T,"conserved", #this BAP is monophyletic both in the gene tree and the consensus lineage tree: not affected by recombination
			"HGT?" #this BAP is monophyletic in the gene tree but not in the lineage consensus tree: potentialy affected by HGT
		))
		
		#return the number of BAPS lineage that are conserved, putative HGT or ignored and the number of embeded species 
		
		c(sapply(c("conserved","HGT?","ignore"),function(x){
			return(length(subset(STATx, STATx==x)))
		}),
		sapply(c("conserved","HGT?","ignore"),function(x){
			return(sum(subset(as.numeric(LINbap[,3]), STATx==x)))
		}))		
		  
	})

	#proportion of BAPS lineage that are conserved for analysis
	LINcons<-(Z[1,]+ Z[2,])/(Z[1,]+Z[2,]+Z[3,])
	#proportion of species that are conserved for analysis
	SPEcons<-(Z[4,]+ Z[5,])/(Z[4,]+Z[5,]+Z[6,])

	return(round(c(LINcons, SPEcons),3))
})




setwd(paste(pB,"/LEVEL", nLEV,"_nLIN=", paste("2-",nVAR,sep=""),sep=""))
  	
  Temp<-c(1:190)  	
  	
pdf("Npop_Combined_VATstat.pdf",height=4,width=4)	
 
par(mar=c(4,4,1,1),bty="n") 
  	
plot(-10,10,xlim=c(0, nVAR),ylim=c(0, 1),
 	xlab="n.pops",ylab="prop. of conserved BAPS lineages",
  		cex.axis=0.8,las=1) 	
 
 Z<-sapply(Temp,function(x){
   	X=x
  	z<-	VATstat[1:(nVAR-1),x]
 	points(c(2: nVAR),z,
  		type="l",col=rgb(0,0,1,0.02),lwd=3,lend=1)
  	list(t(cbind(c(2: nVAR),z)))
 })
 
 Z <-matrix(unlist(Z),ncol=2,byrow=T)
 
sumZ<-t(sapply(sort(unique(Z[,1])),function(x){
	return(c(x,mean(subset(Z[,2],Z[,1]==x)),
	sd(subset(Z[,2],Z[,1]==x))))
})) 

points(sumZ[,1], sumZ[,2],type="l",col="black",lwd=3,lend=1) 

plot(as.factor(Z[,1]),Z[,2],ylim=c(0, 1),
	xlab="n.pops",ylab="prop. of conserved BAPS lineages",
  		cex.axis=0.8,las=1,pch=19,cex=0.5) 	
 
  	
plot(-10,10,xlim=c(0, nVAR),ylim=c(0, 1),
 	xlab="n.pops",ylab="prop. of conserved species",
  		cex.axis=0.8,las=1) 	
 
 Z<-sapply(Temp,function(x){
   	X=x
  	z<-	VATstat[c(nVAR :(2*nVAR-2)),x]
 	points(c(2: nVAR),z,
  		type="l",col=rgb(0,0,1,0.02),lwd=3,lend=1)
  	list(t(cbind(c(2: nVAR),z)))
 })
 
 Z <-matrix(unlist(Z),ncol=2,byrow=T)
 
sumZ<-t(sapply(sort(unique(Z[,1])),function(x){
	return(c(x,mean(subset(Z[,2],Z[,1]==x)),
	sd(subset(Z[,2],Z[,1]==x))))
})) 

points(sumZ[,1], sumZ[,2],type="l",col="black",lwd=3,lend=1) 

plot(as.factor(Z[,1]),Z[,2],ylim=c(0, 1),
	xlab="n.pops",ylab="prop. of conserved BAPS species",
  		cex.axis=0.8,las=1,pch=19,cex=0.5) 
  	
 dev.off() 	
  


############################################################################################################################################
####################
#STEP B
#rhierbaps results per core gene (nLIN fixed) + generate color code for graphic


num.strsplit.u<-function(x){return(as.numeric(unlist(strsplit(x,split= "_"))))}

library(ape)
library(ips)

#Retrieve the list of core genes
pP<-"/myFolder/3-Phylogeny/"
setwd(pP)
KEGGres<-read.table("Summary-KEGG-RF-Pos-per-PEG.txt",header=T)

#Output folder
pB<-"/myFolder/4-BAPS"

#Input folder
pD<-"/myFolder/2-Core-genome"

setwd(pD)

ORD<-read.table("Taxonomy.txt",header=T) #taxonomy and species order according to the species tree
ORD<-ORD[(order(ORD[,4])),]
#Sort species by keeping order within groups but flip groups in hte order:
ORD<-matrix(unlist(sapply(c("D","B","A","C"),function(x){
	list(t(subset(ORD, ORD[,2]==x)))	
})),ncol=5,byrow=T)


###Color codes for graphique
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

plot(c(1:length(COLlin[,1])),pch=19,cex=5,col=COLlin[,2])
text(c(1:length(COLlin[,1])),c(1:length(COLlin[,1])),COLlin[,1])

COLgroup<-t(cbind(
	c("A","red"),
	c("B","purple2"),
	c("C","green4"),
	c("D","blue")
))

plot(c(1:length(COLgroup[,1])),pch=19,cex=5,col= COLgroup[,2])
text(c(1:length(COLgroup[,1])),c(1:length(COLgroup[,1])), COLgroup[,1])


#color codes for hybrid lineages (orange and yellow scale)

COLhyb<-c("yellow2","orange2","yellow4","orange4")


#Consensus Lineage tree (ASTRAL) only with PEG

setwd(pP)
LINtree<-read.tree("Methylobacterium-AsterLinTree-PEG-sept-2025_mBS=10.tre")

nLEV=2
nLIN=length(unique(ORD[,3]))

CG_rhierbaps<-sapply(c(1:length(KEGGres[,1])),function(x){
	print(x)
	X=x
	
	#Gene code
	CGx= KEGGres[X,1]
			
	#retrieve BAPS assignation
	setwd(paste(pB,"/LEVEL", nLEV,"_nLIN=", nLIN,sep=""))

	linx<-read.table(paste(CGx,"-rhierbaps.txt",sep=""),header=T) #Lineage file with how rhierbaps clustered strains for this gene
		
	#phylogenetic tree for this gene
	setwd(paste(pP,"RAxML-trees",sep="/"))
	treex<-read.tree(paste("RAxML_tree.",
		gsub("CG_","", CGx),"_reformated",sep="")) 
		
	#in the lineage consensus tree, only keep tips from this tree
	treel10<-unroot(keep.tip(LINtree,tip= treex$tip))
		
	#Convert it so it fit the species tree (for each species, report the BAP lineage(s) )	
	LNx<-sapply(c(1:length(ORD[,4])),function(x){
		stx<-unlist(strsplit(ORD[x,5],split="_"))	
		paste(c(unique(sapply(stx,function(x){
			lx<-c(subset(linx[,2],linx[,1]==x),NA)[1]
			return(lx)
		})),NA),collapse="_")
	})
		
	LNx<-gsub("_NA","",LNx)
	LNx<-gsub("NA_","",LNx)
		
	#For each BAP lineage, return: 1) corresponding Methylo groups; 2) the corresponding methylo group with the larger number of species in this BAP lineage; 3) corresponding Methylo lineages; 3) the corresponding methylo lineage with the larger number of species in this BAP lineage; 4) Whether this BAP is monophyletic or not in the corresponding tree; 5) is it monophyletic in the consensus lineage tree
	LINbap<-t(sapply(sort(unique(linx[,2])),function(x){
		stx<-subset(linx[,1],linx[,2]==x) #strains assigned to this LINbap 
		z<-unique(sapply(stx,function(x){
			grep(paste("_",x,"_",sep=""),
				paste("_",ORD[,5],"_",sep=""))
		}))
		#corresponding Methylo groups (number of species)
		gx<-sort(summary(as.factor(ORD[z,2])),decreasing=T)
		#corresponding Methylo lineages (number of species)
		lx<-sort(summary(as.factor(ORD[z,3])),decreasing=T)
		gy<-paste(names(gx),collapse="_")
		ly<-paste(names(lx),collapse="_")
		
		return(c(x,length(stx),sum(lx),names(gx)[1],
			names(lx)[1],gy,ly,
			is.monophyletic(unroot(treex), stx), #is this BAP monophyletic in the gene tree?
			is.monophyletic(unroot(treel10), stx) #is this BAP monophyletic in the lineage consensus tree?
			))
	}))
		
	colnames(LINbap)<-c("BAPS","strains","species",
		"Main_group","Main_lineage","Groups",
		"Lineages","MonoGeneTree?","MonoLinTree?")


	#Add a upper level (CLADES): test all pairs of LINbap and check if these pairs form a monophyletic group in the tree	
	cladx<-unlist(sapply(c(1:(length(LINbap[,1])-1)),function(x){
		X=x
		sapply(c((X+1):(length(LINbap[,1]))),function(x){
			 return(paste(X,x,sep="_"))
		})
	}))
		
	CLADbap<-subset(cladx,sapply(c(1:length(cladx)),function(x){
		stx1<-subset(linx[,1],linx[,2]==num.strsplit.u(cladx[x])[1]) #strains assigned to the first LINbap of this pair
		stx2<-subset(linx[,1],linx[,2]==num.strsplit.u(cladx[x])[2]) #strains assigned to the second LINbap of this pair			
		#if the two LINbap were already monophyletic, return 0; if these two LINbap form a monophyletic clade when merged together, return 1, otherwise return 0
		ifelse(
		(ifelse(LINbap[num.strsplit.u(cladx[x])[1],8]==T,1,0)+
		ifelse(LINbap[num.strsplit.u(cladx[x])[2],8]==T,1,0))==2,0,			ifelse(is.monophyletic(unroot(treex), c(stx1,stx2))==T,1,0))
			#	
	})==1)
		
		
	CLADbap<-as.character(c(setdiff(as.numeric(LINbap[,1]),
		num.strsplit.u(CLADbap)), CLADbap))

	CLADbap<-t(sapply(CLADbap,function(x){
		stx<-c(subset(linx[,1],linx[,2]==num.strsplit.u(x)[1]),
			subset(linx[,1],linx[,2]==num.strsplit.u(x)[2]))		
		z<-unique(sapply(stx,function(x){
			grep(paste("_",x,"_",sep=""),
			paste("_",ORD[,5],"_",sep=""))
		}))
		#corresponding Methylo groups (number of species)
		gx<-sort(summary(as.factor(ORD[z,2])),decreasing=T)
		#corresponding Methylo lineages (number of species)
		lx<-sort(summary(as.factor(ORD[z,3])),decreasing=T)		
		gy<-paste(names(gx),collapse="_")
		ly<-paste(names(lx),collapse="_")		
		return(c(x,length(stx),sum(lx),
		names(gx)[1],names(lx)[1],gy,ly,			
		is.monophyletic(unroot(treex), stx), #is this BAP monophyletic in the gene tree?
		is.monophyletic(unroot(treel10), stx) #is this BAP monophyletic in the lineage consensus tree?
			))
	}))
		
	colnames(CLADbap)<-c("BAPS","strains","species",
		"Main_group","Main_lineage","Groups",
		"Lineages","MonoGeneTree?","MonoLinTree?")

		
	#Different cases for each cluster/BAP_lineage
	# 1) if a BAP is composed by all strains from a same and unique group (monophyletic in the lineage consensus tree): color according to this group (in COLgroup but lighter)
	# 2) if a BAP is composed by all strains from a same and unique lineage  (monophyletic in the lineage consensus tree): color according to this group (in COLgroup but lighter)
	# 3) if a BAP is composed by different but not all strains/lineages from a same and unique group:
	# 3a: this BAP is not monophyletic in the gene tree: ignore
	# 3b: this BAP is monophyletic in the gene tree and the lineage consensus tree: ignore
	# 3c: this BAP is monophyletic in the gene tree BUT NOT in the lineage consensus tree: color according to the group (in COLgroup)
	# 4) if a BAP is composed by different lineage(s) from different groups: 
	# 4a: this BAP is not monophyletic in the gene tree: ignore
	# 4b: this BAP is monophyletic in the gene tree: use a special color (in COLhyb)	
		

	STATx<-ifelse(CLADbap[,8]==F,"ignore",   #this BAP is not monophyletic in the gene tree: ignore
		ifelse(CLADbap[,9]==T,"conserved", #this BAP is monophyletic both in the gene tree and the consensus lineage tree: not affected by recombination
		"HGT?" #this BAP is monophyletic in the gene tree but not in the lineage consensus tree: potentialy affected by HGT
	))
	
	#Color code for graphic representation according to PAB congruence with phylogenies			
	COLx<-t(sapply(c(1:length(CLADbap[,1])),function(x){	
		LX= CLADbap[x,7]
			
		#1 is this BAP composed by all strains from a same and unique group?if yes, 1
		C1=ifelse(CLADbap[x,4]!= CLADbap[x,6],0,
			ifelse(length(grep(CLADbap[x,4], CLADbap[,6]))==1,1,0))
		#2 is this BAP composed by all strains from a same and unique lineage? if yes, 1
		C2=ifelse(C1==1,0,ifelse(CLADbap[x,5]!= CLADbap[x,7],0,
			ifelse(length(grep(CLADbap[x,5], CLADbap[,7]))==1,1,0)))
		#3 is this BAP composed by different but not all strains/lineages from a same and unique group? if yes, 1
		C3=ifelse(C1==1,0,ifelse(C2==1,0,
			ifelse(
				sum(sapply(COLgroup[,1],function(x){
					length(grep(x, LX))				
				}))==1,1,0
			)))
		#4 is this BAP composed  by different lineage(s) from different groups? (0 at the last question)	
		colx<-ifelse(C1==1,subset(COLgroup[,2],
			COLgroup[,1]== CLADbap[x,4]),
			ifelse(C2==1,subset(COLgroup[,2],
			COLgroup[,1]== CLADbap[x,4]),
			ifelse(C3==1,subset(COLgroup[,2],
			COLgroup[,1]== CLADbap[x,4]),NA)))			
			
		colx<-ifelse(STATx[x]=="ignore",rgb(0.85,0.85,0.85),colx)
				return(c(x,colx))
	}))
	
	COLhybx<-c(subset(COLx[,1],is.na(COLx[,2])==T),1000,2000)
	COLhybx<-cbind(COLhybx,c(COLhyb, COLhyb)[c(1:length(COLhybx))])
	colnames(COLhybx)<-colnames(COLx)
	COLx<-rbind(COLx, COLhybx)
	COLx<-subset(COLx,is.na(COLx[,2])==F)
	COLx<-subset(COLx,as.numeric(COLx[,1])<1000)
	COLx<-COLx[order(as.numeric(COLx[,1])),]
	COLx<-cbind(t(col2rgb(COLx[,2])/255),
		ifelse(STATx=="HGT?",1,0.15))
	COLx<-rgb(COLx[,1], COLx[,2],COLx[,3],COLx[,4])
	LNx<-ifelse(LNx=="NA","0",LNx)
	CLAx<-as.vector(sapply(LNx,function(x){
		c(CLADbap[,1][grep(paste("_",x,"_",sep=""),
		paste("_",CLADbap[,1],"_",sep=""))],0)[1]
	}))
	#Color for each species according to assignation
	COLsp<-sapply(CLAx,function(x){	
		z<-subset(COLx, CLADbap[,1]==x)
		return(c(z,"black")[1])
	})	
	
	setwd(paste(pB,"/LEVEL", nLEV,"_nLIN=", nLIN,sep=""))

	write.table(cbind(CLADbap ,STATx, COLx),
		paste("SUM-", CGx,".txt",sep=""))
	
	
	return(COLsp)
})
	
rownames(CG_rhierbaps)<-ORD[,1]
colnames(CG_rhierbaps)<-KEGGres[,1]

setwd(paste(pB,"/LEVEL", nLEV,"_nLIN=", nLIN,sep=""))

write.table(CG_rhierbaps,"Colors-rhierbaps-species.txt")

########################################################################################################################
## STEP C
### Summary Graphic


#Retrieve the list of core genes
pP<-"/myFolder/3-Phylogeny/"
setwd(pP)
KEGGres<-read.table("Summary-KEGG-RF-Pos-per-PEG.txt",header=T)

#Output folder
pB<-"/myFolder/4-BAPS"

#Input folder
pD<-"/myFolder/2-Core-genome"

setwd(pD)

ORD<-read.table("Taxonomy.txt",header=T) #taxonomy and species order according to the species tree
ORD<-ORD[(order(ORD[,4])),]

#Sort species by keeping order within groups but flip groups in hte order:
ORD<-matrix(unlist(sapply(c("D","B","A","C"),function(x){
	list(t(subset(ORD, ORD[,2]==x)))	
})),ncol=5,byrow=T)


#Number of lineages
nLIN=length(unique(ORD[,3]))
#Number of levels
nLEV=2

setwd(paste(pB,"/LEVEL", nLEV,"_nLIN=", nLIN,sep=""))

#Color code per gene and per species

CG_rhierbaps <-read.table("Colors-rhierbaps-species.txt",header=T)

CG_rhierbaps<-CG_rhierbaps[,order(as.numeric(KEGGres[,3]))]
KEGGres<-KEGGres[order(as.numeric(KEGGres[,3])),]

###Structure of the graphic (circular):
# Keep the top left quarter to display the species tree
# Each species in row in the circular panel
#Each species in column (circle)


####Graphique parameters

##COORDINATES OF GENES (sorted according to position in the reference genome)

COORD<-cbind(c(0:(length(KEGGres[,1])-1)),c(1:(length(KEGGres[,1]))), KEGGres[,5],(KEGGres[,4]+KEGGres[,3])/2000000)
rownames(COORD)<-KEGGres[,1]
colnames(COORD)<-c("start","end","RES","POS")

#Openness of the circle
OPEN=0.75

##Angular coordinates of genes (order along hte reference genome)

angS<-2*pi*(COORD[,1]/length(COORD[,1]))*OPEN
angE<-2*pi*(COORD[,2]/length(COORD[,1]))*OPEN

##Angular coordinates of genes (real position along hte reference genome)

angR<-2*pi*(COORD[,4]/max(COORD[,4]))*OPEN


##Radial coordinates (Y axis)
#Limit of the graphic
rLIM=1.65
#Coordinates of the plot
re=1.42 #external radius of gene frames
ri=0.15 #internal radius of gene frames 
#Coordinates of RES values
rc1= re+0.02
rc2= rLIM-0.02

#Convert RES values

rc<-((COORD[,3]-min(COORD[,3]))/(max(COORD[,3])-min(COORD[,3]))*(rc2-rc1))+ rc1

#Y Axis (Y, position, graduations at -0.08,-0.04,0,0.04,0.08,0.12)
Yaxis<-c(-0.05,0,0.05,0.1)

rY<-((Yaxis-min(COORD[,3]))/(max(COORD[,3])-min(COORD[,3]))*(rc2-rc1))+ rc1

#Radius for species (reverse order if you want group C is inside)

rs<-seq(re,ri,length.out=(length(ORD[,1])+1))
#rs<-seq(ri,re,length.out=(length(ORD[,1])+1))

#Clusters of interest

setwd(pP)

CGf<-read.table("CG_FastGear.txt",header=T) #Positions of clusters of interest in the CBMB20 genome
CLU<-c("KDP","BCHL1","BCHL2")

#return names and positions (max and min) in the heatmap and in the reference genome

POSclu<-t(sapply(CLU,function(x){
	nx=x
	gcx<-subset(CGf[,3],CGf[,7]==x)
	rank<-sapply(gcx,function(x){
		return(subset(c(1:length(KEGGres[,1])), 
		KEGGres[,1]==x))
	})
	angx<-sapply(gcx,function(x){
		return(c(subset(angS,names(angS)==x),
			subset(angE,names(angE)==x),
			subset(angR,names(angR)==x)))
	})
	return(c(nx, 
	angx[,c(subset(c(1:length(gcx)),rank==min(rank)))],
	angx[,c(subset(c(1:length(gcx)),rank==max(rank)))],
	paste(gcx,collapse="-")))
}))
colnames(POSclu)<-c("Cluster","angS1","angE1","angR1","angS2","angE2","angR2","CG")

setwd(paste(pB,"/LEVEL", nLEV,"_nLIN=", nLIN,sep=""))

pdf("Circular_BAPS_nLIN=24.pdf")

###Frame of the plot
zones<-matrix(c(1),ncol=1)
layout(zones)
#layout.show(max((zones)))
par(bty="n",mar=c(1,1,1,1))
plot(-10,-10,xlim=c(-rLIM, rLIM),ylim=c(-rLIM, rLIM),xaxt="n",xlab="",yaxt="n",ylab="")

###Frames for clusters of interest (POSclu)

sapply(c(1:length(POSclu[,1])),function(x){
	nx= POSclu[x,1]
	#frame in the heatmap
	xx1<-min(rs)*sin(as.numeric(POSclu[x,2]))
	yy1<-min(rs)*cos(as.numeric(POSclu[x,2]))
	xx2<-min(rs)*sin(as.numeric(POSclu[x,6]))
	yy2<-min(rs)*cos(as.numeric(POSclu[x,6]))
	
	xx3<-max(rs)*sin(as.numeric(POSclu[x,2]))
	yy3<-max(rs)*cos(as.numeric(POSclu[x,2]))
	xx4<-max(rs)*sin(as.numeric(POSclu[x,6]))
	yy4<-max(rs)*cos(as.numeric(POSclu[x,6]))
	
	polygon(c(xx1,xx2,xx4,xx3),c(yy1,yy2,yy4,yy3),
		col=NA,border="yellow4",lwd=0.5)
	#frame in the genome position
	xx1<-min(rc)*sin(as.numeric(POSclu[x,4]))
	yy1<-min(rc)*cos(as.numeric(POSclu[x,4]))
	xx2<-min(rc)*sin(as.numeric(POSclu[x,7]))
	yy2<-min(rc)*cos(as.numeric(POSclu[x,7]))
	
	xx3<-max(rc)*sin(as.numeric(POSclu[x,4]))
	yy3<-max(rc)*cos(as.numeric(POSclu[x,4]))
	xx4<-max(rc)*sin(as.numeric(POSclu[x,7]))
	yy4<-max(rc)*cos(as.numeric(POSclu[x,7]))
	
	polygon(c(xx1,xx2,xx4,xx3),c(yy1,yy2,yy4,yy3),
		col="yellow",border=NA)	
		
	
})


#Gene frames
RESpt<-t(sapply(c(1:length(COORD[,1])),function(x){
	#print(x)
	ex<-angE[x]
	sx<-angS[x]
	rx<-angR[x]
	nx=rownames(COORD)[x]
	cx<-CG_rhierbaps[,x]
	rsx<-COORD[x,3]
	
	segx<-sapply(c(1:length(cx)),function(x){
		
		XXeS<-rs[x]*sin(sx)
		XXeE<-rs[x]*sin(ex)
		XXiE<-rs[x+1]*sin(ex)
		XXiS<-rs[x+1]*sin(sx)
	
		YYeS<-rs[x]* cos(sx)
		YYeE<-rs[x]*cos(ex)
		YYiE<-rs[x+1]* cos(ex)
		YYiS<-rs[x+1]* cos(sx)
		
		polygon(c(XXeS, XXeE, XXiE, XXiS),
		c(YYeS, YYeE, YYiE, YYiS),col= cx[x],border=NA)	
		return(c(mean(c(c(XXeS, XXeE, XXiE, XXiS))),
			mean(c(YYeS, YYeE, YYiE, YYiS))))
	})[,(ifelse(rs[1]>rs[2],1,length(ORD[,1])))]
	
	#to link each gene between its localization in the heatmap and in the genome
	XXrc<-rc[x]*sin(rx)
	YYrc<-rc[x]*cos(rx)		
	#XXrc<-rc[x]*sin((sx+ex)/2)
	#YYrc<-rc[x]*cos((sx+ex)/2)	
	
	#segments(segx[1], segx[2], XXrc, YYrc,col=rgb(0.8,0.8,0.8,0.5),lwd=0.5)
	
	#points(XXrc, YYrc,pch=19,cex=0.2,col="black")
	
	colx<-length(intersect(rownames(COORD)[x],
		unlist(strsplit(POSclu[,8],split="-"))))
	
	return(c(segx[1], segx[2], XXrc, YYrc,colx))
}))



###RES values

#Y axis

absy<-sapply(Yaxis,function(x){
	X=subset(rY, Yaxis==x)
	xx<-X*sin((angS+angE)/2)
	yy<-X*cos((angS+angE)/2)
	lines(xx,yy,lwd=0.5,col=ifelse(x==0,"black",rgb(0.9,0.9,0.9)))
	text(xx[1],yy[1],x,pos=2,cex=0.3)
	return(yy[1])
})

segments(-0.02,min(absy),-0.02,max(absy),lwd=0.5)
segments(-0.03, absy,-0.02, absy,lwd=0.5)

#Genome positions (every 5 Gb)
Xaxis<-seq(0,round(max(KEGGres[,c(3,4)])/10^6,0),length.out=13)

angPOS<-Xaxis/(max(KEGGres[,c(3:4)])/10^6)*2*pi*OPEN

segments(min(absy)*sin(angPOS),min(absy)*cos(angPOS),
	max(absy)*sin(angPOS),max(absy)*cos(angPOS),
	col="black",lwd=0.5)

text(rLIM*sin(angPOS),
 	rLIM*cos(angPOS),
 	Xaxis, cex=0.3)

#RES Values per gene
segments(RESpt[,1], RESpt[,2], RESpt[,3], RESpt[,4],
	lwd=0.5,
	col=ifelse(RESpt[,5]==1,"yellow3",rgb(0.7,0.7,0.7,0.2)))

points(RESpt[,3], RESpt[,4],pch=19,cex=0.2,
	col=ifelse(RESpt[,5]==1,"yellow4",rgb(0.2,0.2,0.2,0.6)))	



#Add a reference bar of color for groups and lineages

coordG<-c(-0.075, -0.06)
coordL<-c(-0.05, -0.02)
	
sapply(COLgroup[,1],function(x){
	colx<-col2rgb(subset(COLgroup[,2], 
		COLgroup[,1]==x))[,1]/255
	colx<-rgb(colx[1], colx[2], colx[3],1)
	z<-subset(cbind(rs[2:length(rs)],rs[1:(length(rs)-1)]),
		ORD[,2]==x)
	polygon(c(coordG[1], coordG[2], coordG[2], coordG[1]),
		c(min(z),min(z),max(z),max(z)),
		col= colx,border="white",lwd=0.2)
	text(max(coordG)+0.02,mean(z),x,col= colx,
	font=2,cex=0.5,pos=2)	
		
})

sapply(COLlin[,1],function(x){
	colx<-subset(COLlin[,2], COLlin[,1]==x)
	z<-subset(cbind(rs[2:length(rs)],rs[1:(length(rs)-1)]),
		ORD[,3]==x)
	polygon(c(coordL[1], coordL[2], coordL[2], coordL[1]),
	c(min(z),min(z),max(z),max(z)),
		col= colx,border="white",lwd=0.2)
	text(mean(coordL),mean(z),x,col="white",font=2,cex=0.15)		
})



dev.off()
	
### STEP D

##Proportion of each type per species
colz<-sort(unique(
		as.vector(unlist(CG_rhierbaps))))

PROPcg<-t(sapply(c(1:length(CG_rhierbaps[,1])),function(x){
	
	cx<-as.vector(unlist(CG_rhierbaps[x,]))
	
	return(sapply(colz,function(x){
		length(subset(cx, cx==x))
		
	}))
	
}))


colz<-cbind( c("No\nHGT","Within\ngroup","No\nHGT",
"Within\ngroup" ,"Among\ngroups","No\nHGT",
"Within\ngroup","Low\nresolution","Among\ngroups" ,
"Among\ngroups","No\nHGT","Within\ngroup","Missing"),
c("#0000FF26","#0000FFFF","#008B0026","#008B00FF" ,"#8B8B00FF","#912CEE26","#912CEEFF","#D9D9D926","#EE9A00FF", "#EEEE00FF","#FF000026","#FF0000FF","black"))

#barplot(t(PROPcg),col= colz,border=NA,legend=T)

PROPcg2<-sapply(sort(unique(colz[,1])),function(x){
	apply(subset(t(PROPcg), colz[,1]==x),2,sum)
})

rownames(PROPcg2)<-ORD[,1]

pAmong<-PROPcg2[,1]/(PROPcg2[,1]+PROPcg2[,4]+PROPcg2[,5])
pWithin<-PROPcg2[,5]/(PROPcg2[,1]+PROPcg2[,4]+PROPcg2[,5])

colPROP<-sapply(ORD[,3],function(x){
	subset(COLlin[,2],COLlin[,1]==x)
})




###Color codes for graphique
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


setwd(paste(pB,"/LEVEL", nLEV,"_nLIN=", nLIN,sep=""))

pdf("proportion_species_BAPS_nLIN=24.pdf",height=2,width=6)

par(mar=c(0.5,6,0.5,0.5),mfrow=c(5,1))
sapply(c(4,5,1,2,3),function(x){
	barplot(PROPcg2[,x],las=2,cex.names=0.01,
	ylab=colnames(PROPcg2)[x],
	col= colPROP,border= colPROP,cex.axis=0.7,
	cex.lab=0.8)
})

dev.off()



par(mar=c(0.5,6,0.5,0.5),mfrow=c(1,1))

barplot(t(cbind(PROPcg2[,4]/(PROPcg2[,1]+PROPcg2[,4]+PROPcg2[,5]),
	PROPcg2[,5]/(PROPcg2[,1]+PROPcg2[,4]+PROPcg2[,5]),
	PROPcg2[,1]/(PROPcg2[,1]+PROPcg2[,4]+PROPcg2[,5]))))

### STEP E
###average proportion of each type of BAPS lineage per Methylo lineages

PROPlinM<-t(sapply(unique(ORD[,3]),function(x){
	
	prox<-subset(PROPcg2,ORD[,3]==x)
	round(apply(prox,2,mean),0)
}))

PROPlinSD<-t(sapply(unique(ORD[,3]),function(x){
	
	prox<-subset(PROPcg2,ORD[,3]==x)
	round(apply(prox,2,sd),0)
}))


setwd(paste(pB,"/LEVEL", nLEV,"_nLIN=", nLIN,sep=""))

pdf("BAPSlin_mean-per-METHlin.pdf",height=6,width=10)

par(mar=c(4,4,1,4))

barplot(t(PROPlinM[order(PROPlinM[,5]+PROPlinM[,1],decreasing=T),c(5,1,4,2,3)]),
	col=c("yellow3","grey","black","white","red")[c(5,1,4,2,3)],
	las=1,cex.axis=0.8,cex.names=0.8,ylab="Core PEGs",
	xlab="Methylobacterium lineage")

barplot(t(PROPlinM[order(PROPlinM[,5]+PROPlinM[,1],decreasing=T),c(5,1,4)]/apply(PROPlinM[order(PROPlinM[,5]+PROPlinM[,1],decreasing=T),c(5,1,4)],1,sum)),
	las=1,cex.axis=0.8,cex.names=0.8,ylab="prop. Core PEGs (w/o missing and low res.)",
	col=c("red","yellow3","white"),
	xlab="Methylobacterium lineage")

dev.off()

PROPlin<-matrix(paste(PROPlinM[order(PROPlinM[,5]+PROPlinM[,1],decreasing=T),c(5,1,4,2,3)],
	PROPlinSD[order(PROPlinM[,5]+PROPlinM[,1],decreasing=T),c(5,1,4,2,3)],sep=" ± "),
	ncol=5)
colnames(PROPlin)<-colnames(PROPlinM)[c(5,1,4,2,3)]
rownames(PROPlin)<-rownames(PROPlinM)[order(PROPlinM[,5]+PROPlinM[,1],decreasing=T)]
PROPlin<-gsub(" ± NA","",PROPlin)

write.table(PROPlin,"BAPSlin_mean-per-METHlin.txt")

#proportion of each type of BAPS lineage per gene (and calculate per BCHL, KDP and others)


PROPcg3<-t(sapply(c(1:length(CG_rhierbaps[1,])),function(x){
	
	cx<-as.vector(unlist(CG_rhierbaps[,x]))
	
	return(sapply(colz,function(x){
		length(subset(cx, cx==x))
		
	}))
	
}))


PROPcg4<-sapply(sort(unique(colz[,1])),function(x){
	apply(subset(t(PROPcg3), colz[,1]==x),2,sum)
})


PROPcg4<-ifelse(PROPcg4==0,0,1)

summary(as.factor(apply(PROPcg4[,c(5,1,4)],1,paste,collapse="")))

# Genes with too much missing information / low resolution: 0
# Genes with no evidence of HGT: 514 ("001")
# Genes with evidence of HGT within groups only: 495 ("101")
# Genes with evidence of HGT among groups only: 33 ("011")
# Genes with evidence of HGT among and within groups: 57 ("111")


			#No HGT  #Within only #Among only #Among and within
			#001 #101 #011 #111
#In BCHL1.  #0   #6   #4   #9
#In BCHL2	#2   #1   #5   #5
#in KDP.    #0   #2   #0   #2

#### STEP F - Compare BAPS incongruence with RF residuals

		
###RES value in function of the proportion of BAPs lineages with putative HGT signal	(IF STATx is returned)

CGfsum<-sapply(c(1:length(KEGGres[,1])),function(x){
	zx<-CG_rhierbaps[[x]]
	return(sapply(sort(unique(unlist(CG_rhierbaps))),
	function(x){
		length(subset(zx, zx==x))
	}))
})

mis<-CGfsum[3,]  #no resolution
hgt<-CGfsum[2,] #potential HGT
no<-CGfsum[1,] #no HGT detected 

Phgt<-((hgt)/(hgt+no))
RES<-as.numeric(KEGGres[,5])

MOD<-lm(RES ~ Phgt)

summary(MOD)

#Call:
#lm(formula = RES ~ Phgt)

#Residuals:
#      Min        1Q    Median        3Q       Max 
#-0.093080 -0.016124  0.000424  0.015159  0.106380 

#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept) -0.0063249  0.0009894  -6.393 2.41e-10 ***
#Phgt         0.0657103  0.0064613  10.170  < 2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.02551 on 1097 degrees of freedom
#Multiple R-squared:  0.08616,	Adjusted R-squared:  0.08532 
#F-statistic: 103.4 on 1 and 1097 DF,  p-value: < 2.2e-16

anova(MOD)

cor.test(RES, Phgt, method = "pearson")
cor.test(RES, Phgt, method = "spearman")

library(ggplot2)

dat<-as.data.frame(cbind(RES, Phgt))

r2 <- round(summary(MOD)$r.squared, 3)
pval <- signif(summary(MOD)$coefficients[2,4], 3)
 
setwd(pP)

pdf("HGT_rhierbaps_vs_RF_linear.pdf",height=5.5,width=6.5)

ggplot(dat, aes(x = Phgt, y = RES)) +
	geom_point(size = 2, alpha = 0.2) +
	geom_smooth(method = "lm",
	color = "steelblue",
	fill = "lightblue") +
annotate(
"text",
x = max(dat$Phgt)*0.7,
y = max(dat$RES),
label = paste0("R² = ", r2,
"\np = ", pval),
hjust = 0
) +
theme_classic(base_size = 14) +
labs(
x = "Proportion of phylogenetically \n incongruent BAPS lineages",
y = "HGT likelihood (RES)"
)

dev.off()


#Log model
mod.log <- lm(RES ~ log(Phgt + 0.01), data = dat)
summary(mod.log)

#Call:
#lm(formula = RES ~ log(Phgt + 0.01), data = dat)

#Residuals:
#      Min        1Q    Median        3Q       Max 
#-0.093820 -0.015955  0.000408  0.016011  0.106437 

#Coefficients:
#                  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)      0.0136150  0.0018429   7.388 2.95e-13 ***
#log(Phgt + 0.01) 0.0043924  0.0005384   8.158 9.28e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.02591 on 1097 degrees of freedom
#Multiple R-squared:  0.05719,	Adjusted R-squared:  0.05633 
#F-statistic: 66.55 on 1 and 1097 DF,  p-value: 9.281e-16

AIC(MOD, mod.log)

# Extraction des statistiques
r2 <- summary(mod.log)$r.squared
p <- summary(mod.log)$coefficients[2,4]

# Grille de prédiction
newdat <- data.frame(
  Phgt = seq(min(dat$Phgt),
             max(dat$Phgt),
             length.out = 200)
)

pred <- predict(mod.log,
                newdata = newdat,
                interval = "confidence")

newdat <- cbind(newdat, pred)

# Texte des statistiques
lab <- paste0(
  "R² = ", round(r2, 3),
  "\np = ", signif(p, 3)
)

# Figure

setwd(pP)

pdf("HGT_rhierbaps_vs_RF_log.pdf",height=5,width=6)

ggplot(dat, aes(x = Phgt, y = RES)) +
  geom_point(size = 2, alpha = 0.2) +

  geom_ribbon(
    data = newdat,
    aes(x = Phgt, ymin = lwr, ymax = upr),
    fill = "steelblue",
    alpha = 0.2,
    inherit.aes = FALSE
  ) +

  geom_line(
    data = newdat,
    aes(x = Phgt, y = fit),
    colour = "steelblue",
    linewidth = 1.2,
    inherit.aes = FALSE
  ) +

  annotate(
    "text",
    x = max(dat$Phgt) * 0.7,
    y = max(dat$RES),
    label = lab,
    hjust = 0
  ) +

  theme_classic() +
  labs(x = "Proportion of phylogenetically \n incongruent BAPS lineages",
y = "HGT likelihood (RES)")
  
dev.off()
  
