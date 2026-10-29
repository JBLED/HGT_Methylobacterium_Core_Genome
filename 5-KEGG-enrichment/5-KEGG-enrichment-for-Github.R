##############################################################################
#################################################################
####STEP A - Identify KEGGS  that are significantly enriched for HGT (Higher RES10 values) 

#Path

pP<-"/Users/jean-baptiste/Dropbox/Article-HGT/For-Github/3-Phylogeny/"

pD<-"/Users/jean-baptiste/Dropbox/Article-HGT/For-Github/2-Core-genome/"



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

#Retrieve KEGG information
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
setwd(pP)

KEGGres<-read.table("Summary-KEGG-RF-Pos-per-PEG.txt",header=T)

#For each KEGG pathway, gene indexes, the average and sd of residual values 
KEGGresSUM<-t(sapply(FEGGcode,function(x){
	ix<-grep(paste("_",x,"_",sep=""),
		paste("_", KEGGres[,6],"_",sep=""))
	resx<-as.numeric(KEGGres [ix,5])
	return(c(paste(ix,collapse="_"),mean(resx),sd(resx)))
}))

####Do K permutation of genes along the genome to estimate KEGGs that have significantly higher residual RF that expected

K=10000
Nrep=50 #number of replicates (to avoid loosing data) 

sapply(c(1:Nrep),function(x){
	print(x)
	RFperm<-mclapply(c(1:K),function(x){
		#print(x)
		RFk<-sample(KEGGres[,5]) #randomize RF values
		
		#Mean residual value recalculated after permutations
		#For KEGGs
		Kk<-as.numeric(sapply(FEGGcode,function(x){
			mean(RFk[grep(paste("_",x,"_",sep=""),
				paste("_", KEGGres[,6],"_",sep=""))])
		}))
		return(round(Kk,4))
	},mc.cores=8)

	RFperm <-t(sapply(c(1:K),function(x){
		return(as.vector(RFperm[[x]]))
	}))
	setwd(paste(pP,"PERM",sep="/"))
	write.table(RFperm,paste("PERM-RF-W+KEGG_rep=",x,".txt",sep=""),row.names=F,col.names=F)

})

#Observed values
OBSz<-as.numeric(KEGGresSUM[,2])

SUMperm<-sapply(c(1:Nrep),function(x){
	print(x)
	setwd(paste(pP,"PERM",sep="/"))
	RFperm<-read.table(paste("PERM-RF-W+KEGG_rep=",x,".txt",sep=""))

	#number of permutations with observed RF lower than expected
	RFexp1x<-sapply(c(1:length(OBSz)),function(x){
		z<-RFperm[,x]
		sum(ifelse(OBSz[x] <z,1,0))
	})
	#number of permutations with observed RF higher than expected
	RFexp2x<-sapply(c(1:length(OBSz)),function(x){
		z<-RFperm[,x]
		sum(ifelse(OBSz[x] >z,1,0))
	})
	return(list(t(cbind(RFexp1x, RFexp2x))))
})

RFexp1 <-apply(sapply(c(1:Nrep),function(x){
	return(SUMperm[[x]][1,])
}),1,sum)	

RFexp2 <-apply(sapply(c(1:Nrep),function(x){
	return(SUMperm[[x]][2,])
}),1,sum)	

	
#Calculate pvalues 
Prf2 =((K* Nrep)-RFexp2)
Prf1 =((K* Nrep)-RFexp1)

Prf<-(apply(cbind(Prf1, Prf2),1,min) +1)/((K* Nrep)+1)

#adjusted pvalues for KEGGs
PFadj<-p.adjust(Prf,method="bonferroni")

#Real gene indexes for KEGGs
CGk<-sapply(c(1:length(KEGGresSUM[,1])),function(x){
	Gx<-num.strsplit.u(KEGGresSUM[x,1])
	paste(c("CG",gsub("CG_","",KEGGres[Gx,1])),collapse="_")
})

KEGGresSUM2<-cbind(CGk , KEGGresSUM, Prf1 , 
	Prf2, Prf, PFadj)

colnames(KEGGresSUM2)<-c("Genes","Index","MeanResRF","SdResRF","obs>exp","obs<exp","p.value","adj.p.value")

setwd(pP)

write.table(KEGGresSUM2,"Perm-RF-Summary-Per-KEGG.txt")

############################################################################################################ STEP B : DISPLAY RESULTS FROM PERMUTATIONS TO VISUALIZE KEGGS WITH HIGHER RF THAN EXPECTED 


setwd(pP)

KEGGresSUM2 <-read.table("Perm-RF-Summary-Per-KEGG.txt",header=T)
KEGGres<-read.table("Summary-KEGG-RF-Pos-per-PEG.txt",header=T)


Nkegg<-gsub("One_carbon","One-carbon",rownames(KEGGresSUM2))
Nkegg<-gsub("_","\n", Nkegg)
Nkegg<-gsub("pool\nby\nfolate","pool by folate", Nkegg)

setwd(pP)

pT=0.05
pOUT<-pT/4
		
pdf("KEGG-enrichment.pdf",height=8,width=6)

par(mar=c(4,4,1,1),bty="n")
		
plot(as.numeric(KEGGresSUM2[,7]),as.numeric(KEGGresSUM2[,3]),log="x",
	xlim=c(0.5,min(as.numeric(KEGGresSUM2[,7]))/100),xlab="p.value",ylab="mean residual RF",
	cex.axis=0.8,las=1,
	cex= sqrt(as.numeric(summary(strsplit((KEGGresSUM2[,2]),split="_"))[,1])/pi),
	pch=19,
	col=ifelse(as.numeric(KEGGresSUM2[,3])<0,rgb(0.5,0.5,0.5,0.1),
		ifelse(as.numeric(KEGGresSUM2[,8]) <= pT,rgb(1,0,0,0.5),
		ifelse(as.numeric(KEGGresSUM2[,7]) <= (pOUT),rgb(0.5,0.5,0,0.5),
		rgb(0,0,0,0.1)))),
	ylim=c(min(as.numeric(KEGGresSUM2[,3])),max(as.numeric(KEGGresSUM2[,3]))))

plotNames<-gsub("_"," ", rownames(KEGGresSUM2))
plotNames<-sapply(c(1:length(plotNames)),function(x){
	pnx<-unlist(strsplit(plotNames[x],split=" "))
	return(paste(pnx[1],"\n",paste(pnx[2:length(pnx)],collapse=" "),sep=""))
})

text( as.numeric(KEGGresSUM2[,7]),as.numeric(KEGGresSUM2[,3]),
	plotNames,font=2,pos=4,
	cex=ifelse(as.numeric(KEGGresSUM2[,3])<= 0,0.0001,
		ifelse(as.numeric(KEGGresSUM2[,8])<= pT,0.9,
		ifelse(as.numeric(KEGGresSUM2[,7])<= pOUT,0.6,0.0001))),
	col=ifelse(as.numeric(KEGGresSUM2[,3])<0,rgb(0.6,0.6,0.6),
		ifelse( as.numeric(KEGGresSUM2[,8]) <= pT,
		"red3",rgb(0.5,0.5,0))),srt=0)

LEGx<-c(2,50,200)
LEGy<-c(25,25)

par(new=T,mar=c(1,4,18,4))

plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
	xaxt="n",yaxt="n",xlab="",ylab="")

legend("right",cex=0.8, horiz=T,
	legend= LEGx,pt.cex=sqrt(LEGx/pi),
	text.col= "black",box.col=NA,
	pch=19,	col=rgb(0.5,0.5,0.5,0.3),
	border=NA,
	title="Genes per KEGG",
	title.col="black",ncol=1)

par(new=T,mar=c(1,4,24,4))

plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
	xaxt="n",yaxt="n",xlab="",ylab="")

legend("right",cex=0.8, horiz=F,
	legend= c("Adjusted p-value \u2264 0.05",
	"Adjusted p-value > 0.05"),pt.cex=sqrt(LEGy/pi),
	text.col= "black",box.col=NA,
	pch=19,	col=c(rgb(1,0,0,0.5),rgb(0.5,0.5,0,0.5)),
	border=NA,y.intersp = 2,
	title="",
	title.col="black",ncol=1)	
		
dev.off()		
