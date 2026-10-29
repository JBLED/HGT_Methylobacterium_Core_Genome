##############################################################################
#################################################################
####Step A - Identify gene clusters that are significantly enriched for HGT (Higher RES10 values). We first need to determine a cluster size

#Cluster size will be defined according to the reference genome CBM20, using:
# 1 - the distance between two genes involved in the same function (KEGG). Within each window in the range 1-100 kB, we calcultate the proportion of gene pairs that share at least one KEGG.
# 2 - the similarity of evolutionary history between two geness (RF distance). Within each window in the range 1-100 kB, we calculate average RF distance between all pairs of gene trees (corrected as for the RF distance between the gene trees and consensus tree)

# each variable is calculated for different window sizes along the genome. A segment test is used to compare linear regressions between pairwise distance among  genes (log value) and average RF value between genes (1) or proportion of genes pairs shared at least one KEGG (2).

library(ape)
library(ips)
library(phangorn)
library(segmented) 

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

Lx<-c(1:length(CGstat[,1]))

#Index of peg
iP<-subset(Lx,CGstat[,3]=="peg")

setwd(pP)

STATtree<-read.table("GeneTreeStat.txt")

#####For each pair of genes, calculate distance in the reference genome, RF distance and the number of shared KEGG

#KEGG pathways
setwd(paste(pD,"KEGG",sep="/"))

KEGG<-read.table("KEGG.txt",header=T)

#Remove "Human_Diseases" and "Organismal_Systems"
KEGG<-subset(KEGG, KEGG[,4]!="Human_Diseases")
KEGG<-subset(KEGG, KEGG[,4]!="Organismal_Systems")

#KEGG pathway short names
FEGGcode<-sapply(unique(KEGG[,6]),function(x){
	return(unlist(strsplit(x,split="_"))[1])
})

#Simplified KESSres
KEGGres<-t(sapply(CGstat[iP,1],function(x){
	#print(x)
	#KEGG pathways(s) for this gene (only code)
	Kx<-KEGG[grep(paste(x,"-",sep=""),KEGG[,3]),6]
	Kx<-c(paste(sapply(Kx,function(x){
		return(unlist(strsplit(x,split="_"))[1])
	}),collapse="_"),NA)[1]
	#Residual value of RF distance for this gene
	RESx<-subset(STATtree[,6],rownames(STATtree)==x)
	#Position of this gene in the reference genome
	px<-unlist(strsplit(subset(MAP[,2],rownames(MAP)==x),split=";"))
	px<-as.numeric(unlist(strsplit(px,split="-")))
	#RAST annotation
	Rx<-subset(CGstat[,2], CGstat[,1]==x)
	return(c(x, Rx, px[1],px[2],RESx, Kx))
}))

colnames(KEGGres)<-c("Gene","RAST","Start","End","ResidualRF","KEGGs")

setwd(pP)

write.table(KEGGres,"Summary-KEGG-RF-Pos-per-PEG.txt")

#####For each pair of genes, calculate distance in the reference genome, RF and the number of shared KEGGS   (for the pairwise distance, take into account that the chromosome is circular)

#Reference genome size
Sref<-subset(STAT[,4], STAT[,1]=="st36")

ncor=6

setwd(paste(pP,"RAxML-trees",sep="/"))
	
pDIST<-sapply(c(1:(length(iP)-1)),function(x){
	print(x)
	X=x
	CGx=CGstat[iP,1][x]
	treex<-read.tree(paste("RAxML_tree.",gsub("CG_","", CGx),"_reformated",sep=""))
	stx<-treex$tip #tips in this tree
	posx<-unlist(strsplit(subset(MAP[,2],rownames(MAP)== CGx),split=";"))[1]
	posx<-mean(as.numeric(unlist(strsplit(posx,split="-")))) #average position of this gene in the reference genome
	kx=unlist(strsplit(KEGGres[x,6],split="_"))	#KEGGs for this gene


	mclapply(c((X+1):(length(iP))),function(x){
		#print(x)
		CGy=CGstat[iP,1][x]
		treey<-read.tree(paste("RAxML_tree.",gsub("CG_","", CGy),"_reformated",sep=""))
		sty<-treey$tip #tips in this tree
		stxy<-intersect(stx,sty) #shared tips
		#in the lineage consensus tree, remove the missing tips
		treelx2<-unroot(keep.tip(treex,tip= stxy))
		treely<-unroot(keep.tip(treey,tip= stxy))
		posy<-unlist(strsplit(subset(MAP[,2],rownames(MAP)== CGy),split=";"))[1]
		posy<-mean(as.numeric(unlist(strsplit(posy,split="-")))) #average position of this gene in the reference genome	
		ky=unlist(strsplit(KEGGres[x,6],split="_"))	#KEGGs for this gene
		#Distance (in pb) between both gene in the genome and RF distance between both trees (for the pairwise distance, take into account that the chromosome is circular)
		return(c(iP[X],iP[x], #gene number
			ifelse(abs(posx-posy)>(Sref/2),(abs(posx-posy)-(Sref/2)),abs(posx-posy)), #physical distance
			RF.dist(treelx2, treely), #RF distance
			length(intersect(kx,ky))))	#number of shared KEGGs	
	},mc.cores=ncor)
})
	
pDIST <-matrix(unlist(pDIST),ncol=5,byrow=T)	

colnames(pDIST)<-c("Gene1","Gene2","pDIST","RF","KEGGS")

	
#Corrected RF distance 

iPAIRS<-mclapply(c(1:(length(iP)-1)),function(x){
	print(x)
	X=x
	sapply(c((X+1):(length(iP))),function(x){
		return(c(X,x))	
	})
},mc.cores=ncor)

iPAIRS <-matrix(unlist(iPAIRS),ncol=2,byrow=T)	

MODdist<-lm(pDIST[,4] ~ STATtree$TIPs[iPAIRS[,1]]* STATtree$TIPs[iPAIRS[,2]]
	* STATtree$Mns[iPAIRS[,1]]* STATtree$Mns[iPAIRS[,2]]
	*STATtree$nPOL[iPAIRS[,1]]*STATtree$nPOL[iPAIRS[,2]]
	*STATtree$dNdS[iPAIRS[,1]]*STATtree$dNdS[iPAIRS[,2]])

RESdist<-MODdist $residuals	

SUMdist<-cbind(CGstat[iP,1][iPAIRS[,1]],CGstat[iP,1][iPAIRS[,2]] , pDIST[,3], pDIST[,4],round(RESdist,3),pDIST[,5])

###Replace the number of shared KEGGs by NA, when KEGG is unknow for at least one gene
KEGGrep<-ifelse(KEGGres[,6][iPAIRS[,1]]=="","NA",
	ifelse(KEGGres[,6][iPAIRS[,2]]=="","NA",
	SUMdist[,6]))

SUMdist<-cbind(SUMdist[,c(1:5)], KEGGrep)

colnames(SUMdist)<-c("Gene1","Gene2","pDist","RF","RFResiduals","KEGGs")

setwd(pP)

write.table(SUMdist,"pDistanceRF.txt",row.names=F,quote=F)


####For each Kb window, calculate the average and SD residual RF distance and the proportion of pairs of genes that share at least one KEGG

setwd(pP)

SUMdist<-read.table("pDistanceRF.txt",header=T)

#Define window bounadries to ensure about the same number of observations per window

Nobs=600 #number of observation per window

nW<-round(length(SUMdist[,1])/Nobs,0) #number of windows

z<-sort(as.numeric(SUMdist[,3]))

Wz<-ceiling(seq(0.000001, nW,length.out=length(z)))

Wz <-t(sapply(c(1: nW),function(x){
	return(c(x,min(subset(z, Wz==x)),max(subset(z, Wz==x))))
}))

#calculate values of each window size (discrete: as size increase, values below the lower bundary are not considered)
pW<-t(sapply(Wz[,1],function(x){
	print(x)
	px=subset(SUMdist,as.numeric(SUMdist[,3])>= Wz[x,2]) #keep only pairwise comparisons among genes distant by more than this value
	px=subset(px,as.numeric(px[,3])< Wz[x,3]) #keep only pairwise comparisons among genes distant by less than this value
	
	px2<-subset(px,is.na(px[,6])==F) #remove pairs for which at least one gene has no know KEGG
	return(c(
	x,
	round(mean(as.numeric(px[,3])),3),
	round(mean(as.numeric(px[,5])),3),
	round(sd(as.numeric(px[,5])),3),
	dim(px)[1],
	dim(px2)[1],
	round(length(subset(px2[,1],px2[,6]!="0"))/dim(px2)[1],3)
	))
}))

colnames(pW)<-c("W","pDIST","meanRF","sdRF","N","Nkegg","pKEGG")

#Same but cumulated (windows larger and larger)
pWc<-t(sapply(Wz[,1],function(x){
	print(x)
	px=subset(SUMdist,as.numeric(SUMdist[,3])< Wz[x,3]) #this is the part that change compared to pW: all pairwise comparisons below this distance are considered - so values are cumulating as window size is increasing
	px2<-subset(px,is.na(px[,6])==F) #remove pairs for which at least one gene has no know KEGG
	return(c(
	x,
	round(mean(as.numeric(px[,3])),3),
	round(mean(as.numeric(px[,5])),3),
	round(sd(as.numeric(px[,5])),3),
	dim(px)[1],
	dim(px2)[1],
	round(length(subset(px2[,1],px2[,6]!="0"))/dim(px2)[1],3)
	))
}))

colnames(pWc)<-c("W","pDIST","meanRF","sdRF","N","Nkegg","pKEGG")



setwd(pP)

write.table(pW,"PairwiseDistRF-KEGG-per-window.txt")
write.table(pWc,"PairwiseDistRF-KEGG-per-window-cumulated.txt")

setwd(pP)
pW<-read.table("PairwiseDistRF-KEGG-per-window.txt",header=T)
pWc <-read.table("PairwiseDistRF-KEGG-per-window-cumulated.txt",header=T)

pdf("PairwiseDistRF-KEGG.pdf",height=6,width=10)
	
par(mar=c(4,4,1,1),bty="n",mfrow=c(2,3))	
	
plot(pW[,2]/1000, pW[,3],
	type="l",col="red",lwd=1.5,las=1,cex.axis=0.8,log="x",
	xlab="Mean pairwise distance between genes (Kb, log scale)",
	ylab="Mean residual RF distance between genes",
	ylim=c(min(pW[,3]-pW[,4]),max(pW[,3]+pW[,4])))	
segments(pW[,2]/1000,pW[,3]+pW[,4],pW[,2]/1000,pW[,3]-pW[,4],
	col=rgb(0.8,0.8,0.8,0.8),lwd=2)
points(pW[,2]/1000, pW[,3],
	type="l",col="red",lwd=1.5)			
points(pW[,2]/1000, pW[,3],
	pch=19,cex=0.4)

plot(pW[,2]/1000, pW[,7],
	type="l",col="red",lwd=1.5,las=1,cex.axis=0.8,log="x",
	xlab="Mean pairwise distance between genes (Kb, log scale)",
	ylab="Prop. of genes with same KEGG")
			
points(pW[,2]/1000, pW[,7],
	pch=19,cex=0.4)	

plot(pW[,3], pW[,7],
	type="l",col="red",lwd=1.5,las=1,cex.axis=0.8,log="",
	xlab="Mean residual RF distance between genes",
	ylab="Prop. of genes with same KEGG")
	
points(pW[,3], pW[,7],
	pch=19,cex=0.4)	
	
plot(pWc[,2]/1000, pWc[,3],
	type="l",col="red",lwd=1.5,las=1,cex.axis=0.8,log="x",
	xlab="Mean pairwise distance between genes (Kb, log scale)",
	ylab="Mean residual RF distance between genes",
	ylim=c(min(pWc[,3]-pWc[,4]),max(pWc[,3]+pWc[,4])))	
segments(pWc[,2]/1000,pWc[,3]+pWc[,4],pWc[,2]/1000,pWc[,3]-pWc[,4],
	col=rgb(0.8,0.8,0.8,0.8),lwd=2)
points(pWc[,2]/1000, pWc[,3],
	type="l",col="red",lwd=1.5)			
points(pWc[,2]/1000, pWc[,3],
	pch=19,cex=0.4)

plot(pWc[,2]/1000, pWc[,7],
	type="l",col="red",lwd=1.5,las=1,cex.axis=0.8,log="x",
	xlab="Mean pairwise distance between genes (Kb, log scale)",
	ylab="Prop. of genes with same KEGG")
			
points(pWc[,2]/1000, pWc[,7],
	pch=19,cex=0.4)	

plot(pWc[,3], pWc[,7],
	type="l",col="red",lwd=1.5,las=1,cex.axis=0.8,log="",
	xlab="Mean residual RF distance between genes",
	ylab="Prop. of genes with same KEGG")
	
points(pWc[,3], pWc[,7],
	pch=19,cex=0.4)			

dev.off()


#####Detecting the breakpoint Using function segmented 


#Raw data

setwd(pP)

pdf("PairwiseDistRF-KEGG+Segmented-raw.pdf",height=5,width=5)

DISTscale=c(0.1,0.3,1,3,10,30,100,300,1000,3000)

plot(log(SUMdist[,3]),SUMdist[,5],pch=19,cex=0.4,
	log="",col=rgb(0,0,0,0.05),
	xlab="Pairwise distance between genes (Kb, log scale)",
	ylab="Residual RF distance between genes",
	xaxt="n")

axis(1,log(DISTscale*1000),DISTscale,cex.axis=0.8)

dataRF<-data.frame(x=log(SUMdist[,3]),y=SUMdist[,5])
#dataRF<-subset(dataRF, dataRF[,1]<=12.5)
mRF<-lm(y~x,data= dataRF)
zRF<-segmented(mRF, npsi=1)

summary(zRF)

exp(zRF $psi[2])


plot(zRF,add=T,col="red",lwd=2.5,lty=2)

sRF<-exp(segmented(mRF)$psi[2])

segments(log(sRF),-1000,log(sRF),1000,col="blue",lty=2,lwd=1.5)
text((log(sRF)),10,paste("size = ", round((sRF))," bp",sep=""),pos=4,cex=0.8,col="white")

#Per Kb window

dataRF<-data.frame(x=log(pW[,2]),y=pW[,3],z=pW[,7])
mRF<-lm(y~x,data= dataRF)
mKEGG<-lm(z~x,data= dataRF)

zRF<-segmented(mRF, npsi=1)
zKEGG<-segmented(mKEGG, npsi=1)

#Slopes
summary(zRF)
summary(zKEGG)

#Breakpoints
exp(zRF $psi[2])
exp(zKEGG $psi[2])


#Slopes
slope(zRF)
slope(zKEGG)

2 * pt(abs(slope(zKEGG)$x[, "t value"]), df = zKEGG$df.residual, lower.tail = FALSE)

2 * pt(abs(slope(zRF)$x[, "t value"]), df = zRF$df.residual, lower.tail = FALSE)

dev.off()

#davies.test to evaluate if the break in slope is significant (do not take the estimate of breakpoint; use the segemented test value instead: just use this test to see if there's a significant breakpoint in the slope)

dRF<-davies.test(mRF)
dKEGG<-davies.test(mKEGG)

#Breakpoints
exp(dRF$statistic)
exp(dKEGG $statistic)

#pvalues
dRF$p.value
dKEGG $p.value

#global regression

setwd(pP)

pdf("PairwiseDistRF-KEGG+Segmented.pdf",height=5,width=10)

DISTscale=c(1,3,10,30,100,300,1000,3000)

par(mar=c(4,4,1,1),bty="n",mfrow=c(1,2))	

#RF

plot(log(pW[,2]), pW[,3],
	type="l",col="white",lwd=1.5,las=1,cex.axis=0.8,log="",
	xlab="Mean pairwise distance between genes (Kb, log scale)",
	ylab="Mean residual RF distance between genes",
	ylim=c(min(pW[,3]-pW[,4]),max(pW[,3]+pW[,4])),
	xaxt="n")	

segments(log(pW[,2]),pW[,3]+pW[,4],log(pW[,2]),pW[,3]-pW[,4],
	col=rgb(0.9,0.9,0.9,0.8),lwd=2)		
points(log(pW[,2]), pW[,3],
	pch=19,cex=0.6,col=rgb(0.6,0.6,0.6))

axis(1,log(DISTscale*1000),DISTscale,cex.axis=0.8)

plot(zRF,add=T,col="red",lwd=2.5,lty=2)

sRF<-exp( zRF $psi[2])

segments(log(sRF),-1000,log(sRF),1000,col="blue",lty=2,lwd=1.5)
text((log(sRF)),10,paste("breakpoint = ", round((sRF))," bp",sep=""),pos=4,cex=1,col="blue")


###KEGG

plot(log(pW[,2]), pW[,7],
	type="l",col=rgb(0.9,0.9,0.9),lwd=1.5,las=1,cex.axis=0.8,log="x",
	xlab="Mean pairwise distance between genes (Kb, log scale)",
	ylab="Prop. of genes with same KEGG",
	xaxt="n")
			
points(log(pW[,2]), pW[,7],
	pch=19,cex=0.6,col=rgb(0.6,0.6,0.6)	)

axis(1,log(DISTscale*1000),DISTscale,cex.axis=0.8)


plot(zKEGG,add=T,col="red",lwd=2.5,lty=2)

sKEGG<-exp(zKEGG $psi[2])

segments(log(sKEGG),-1000,log(sKEGG),1000,col="blue",lty=2,lwd=1.5)
text((log(sKEGG)),0.65,paste("breakpoint = ", round((sKEGG))," bp",sep=""),pos=4,cex=1,col="blue")

dev.off()
