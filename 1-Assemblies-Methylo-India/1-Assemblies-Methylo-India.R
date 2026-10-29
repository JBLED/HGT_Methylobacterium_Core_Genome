##STEP A

#Format and calculate statistic of assemblies

#install.packages("seqinr")
library(seqinr)

p0="my_folder"

setwd(paste(p0,"Final-assemblies",sep="/"))

#List of selected assemblies
ASS<-read.table("List-final-assemblies.txt",header=T)

#Remove uncorrect assemlies (too low quality)
ASS<-subset(ASS,ASS[,8]!="none")

STATS<-t(sapply(c(1:length(ASS[,1])),function(x){
	Nx<-ASS[x,8]
	print(Nx)
	seqx<-read.fasta(paste(p0,ASS[x,7],sep=""))
	
	#number of scaffolds
	nx<-length(seqx)
	#Scaffold size
	sx<-as.numeric(summary(seqx)[,1])
	#Coverage per scaffold (Megahit)
	covMeg<-as.numeric(sapply(c(1: nx),
		function(x){
		att<-strsplit(attributes(
			seqx[[x]])$Annot,
			split=" ")[[1]][3:4]
		cov<-as.numeric(strsplit(att[1],
			split="=")[[1]][2])
	}))	
	#Coverage per scaffold (Unicycler)	
	covUni<-as.numeric(sapply(c(1: nx),
		function(x){
		att<-strsplit(attributes(
			seqx[[x]])$Annot,
			split=" ")[[1]][3]
		cov<-as.numeric(strsplit(strsplit(att,
			split="=")[[1]][2],
			split="x")[[1]][1])
	}))
	#Number of reads per base (Racon)	
	covRac<-(sapply(c(1: nx),
		function(x){
		att<-strsplit(attributes(
			seqx[[x]])$Annot,
			split=" ")[[1]][3]
		cov<-as.numeric(strsplit(att,
			split=":")[[1]][3])/sx[x]
	}))
	
	#ANNOTATION
	ANNOT<-sapply(c(1: nx),
		function(x){
		return(attributes(
			seqx[[x]])$Annot)
	})
	
	
	COVx<-apply(cbind(covMeg, covUni, covRac),1,mean,na.rm=T)
		
	#Total size
	tx=sum(sx)
	#scaffold GC content
	GCx<-sapply(c(1: nx),function(x){
		return(GC(seqx[[x]]))
	})
	#average GC content
	gcx=sum(GCx* sx)/tx
	#average coverage
	covx=sum(COVx* sx)/tx
	#N50	
	n50x<-sort(sx,decreasing=T)
	n50x <-subset(n50x,
		cumsum(n50x)>= (tx/2))[1]
			
	sumx<-cbind(sx,COVx,GCx)	
	rownames(sumx)<-names(seqx)
	colnames(sumx)<-c("size","coverage","GC")
		
	setwd(paste(p0,"Final-assemblies",sep="/"))
		
	write.table(sumx,paste("sum_", Nx,
		".txt",sep=""))
	
	write.fasta(seqx,names= ANNOT,
		file.out=paste(Nx,".fa",sep=""))	
		
	pdf(paste(Nx,".pdf",sep=""))
	par(mar=c(4,4,1,5),bty="n",
		cex.axis=0.8,las=1)		
	
	plot(GCx[order(sx,decreasing=T)],
		log="x",xlab="",ylab="",
		col="cyan",type="l",lty=2,
		cex.axis=0.8,xaxt="n",yaxt="n")	
	points(GCx[order(sx,decreasing=T)],
		col="blue",pch=19,cex=0.2)

	axis(side=4,at=seq(min(GCx),
		max(GCx),length.out=5),
		labels=paste(100*round(seq(min(GCx),
		max(GCx),length.out=5),3),"%"),
		main="GC content")
	
	par(new=T)

	plot(cumsum(sort(sx,decreasing=T))/10^6,
		log="x",xlab="Scaffolds (log scale)",
		ylab="Cumulated size (Mb)",
		col="grey",type="l",
		cex.axis=0.8,las=1,main=Nx)
	points(cumsum(sort(sx,decreasing=T))/10^6,
		col=ifelse(sort(sx,decreasing=T)<1000,"red",
		ifelse(sort(sx,decreasing=T)<2000,"orange",
		"green3")),pch=19,cex=0.5)	
		
	par(mar=c(4,4,1,5),bty="n",
		cex.axis=0.8,las=1)	
		
	plot(COVx,GCx,log="x",pch=19,
		cex=5*sqrt((sx/max(sx)/pi)),
		col=ifelse(sx<1000,rgb(1,0,0,0.3),
		ifelse(sx<2000,rgb(0.5,0.5,0,0.3),
		rgb(0,0.5,0,0.3))),
		xlab="Coverage (log scale)",
		ylab="GC content")
	
	dev.off()	
			
	return(c(Nx,nx,tx/10^6,max(sx)/10^6,
		n50x/10^6,covx,gcx))
}))
	
	colnames(STATS)<-c("Isolate-Method-Library",
		"n.Scaffolds","Total.size(Mb)","Max.size(Mb)",
		"N50(Mb)","Coverage","GC")

setwd(paste(p0,"Final-assemblies",sep="/"))
	
write.table(STATS,paste("Sum-Final-assemblies.txt",sep=""),row.names=F)

###Scaffold size and coverage filtering to apply to all Megahit assemblies

cx=10
lx1=500
lx2=2000

Mi<-grep("Megahit",ASS[,8])

setwd(paste(p0,"Final-assemblies",sep="/"))

pdf("Megahit-Size-Cov-GC.pdf")

plot(1,1,col="white",log="x",
	ylab="Cumulated Size (Mb)",xlab="Scaffolds",
	cex.axis=0.8,las=1,
	xlim=c(1,10000),ylim=c(0.01,15))

SUMall<-matrix(unlist(sapply(Mi,function(x){
	Nx<-ASS[x,8]
	
	sumx<-read.table(paste("sum_", Nx,
		".txt",sep=""),header=T)
		
	ix<-intersect(subset(c(1:length(sumx[,1])),
		sumx[,1]<lx2),
		subset(c(1:length(sumx[,1])),
		sumx[,2]<cx))
	ix<-unique(c(ix,subset(c(1:length(sumx[,1])),
		sumx[,1]<lx1)))	
	
	sumf<-sumx[ix,]
	sumx<-sumx[setdiff(c(1:length(sumx[,1])),
		ix),]
	
	z<-order(sumx[,1],decreasing=T)

	points(cumsum(sumx[z,1])/10^6,
		type="l",col="grey",lwd=0.5)	
	points(cumsum(sumx[z,1])/10^6,
		pch=19,cex=0.2,
		col=ifelse(sumx[z,1]<lx1,"orange",
		ifelse(sumx[z,1]<lx2,"yellow3",
		"green3")))
	
	sumx<-cbind(1,sumx)
	sumf<-cbind(0,sumf)
	colnames(sumf)<-colnames(sumx)
	
	return(t(rbind(sumx,sumf)))
})),ncol=4,byrow=T)

colx=ifelse(SUMall [,1]==0,rgb(0,0,0,0.2),
		ifelse(SUMall [,2]<lx1,rgb(1,0.5,0,0.3),
		ifelse(SUMall [,2]<lx2,rgb(0.5,0.5,0,0.3),
		rgb(0,0.5,0,0.3))))


plot(SUMall[,2], SUMall[,3],log="xy",
	xlab="Size (bp)",ylab="Coverage",
	pch=19,cex=0.3,col= colx,
	cex.axis=0.8,las=1)
	
plot(SUMall[,2], SUMall[,4],log="x",
	xlab="Size (bp)",ylab="GC",
	pch=19,cex=0.3,col= colx,
	cex.axis=0.8,las=1)	
	
plot(SUMall[,3], SUMall[,4],log="x",
	xlab="Coverage",ylab="GC",
	pch=19,cex=0.3,col= colx,
	cex.axis=0.8,las=1)		

dev.off()


####Apply filters on Megahit assemblies

STATSmf<-t(sapply(Mi,function(x){
	Nx<-ASS[x,8]
	seqx<-read.fasta(paste(Nx,".fa",sep=""))
	nx=length(seqx)
	sumx<-read.table(paste("sum_", Nx,
		".txt",sep=""),header=T)
		
	ix<-intersect(subset(c(1:length(sumx[,1])),
		sumx[,1]<lx2),
		subset(c(1:length(sumx[,1])),
		sumx[,2]<cx))
	ix<-unique(c(ix,subset(c(1:length(sumx[,1])),
		sumx[,1]<lx1)))	
	
	iy<-setdiff(c(1:length(sumx[,1])),
		ix)
	
	sumx<-sumx[iy,]
		
	seqf<-sapply(iy,function(x){
		return(list(unlist(seqx[[x]])))
	})	
	
	ANNOT<-sapply(c(1: nx),
		function(x){
		return(attributes(
			seqx[[x]])$Annot)
	})
	
	write.fasta(seqf,names= ANNOT[iy],
		paste(Nx,"-filt.fa",sep=""))
	
	sx=sumx[,1]
	nx=length(sx)
	tx=sum(sx)
	
	#average GC content
	gcx=sum(sumx[,3]* sx)/tx
	#average coverage
	covx=sum(sumx[,2]* sx)/tx
	#N50	
	n50x<-sort(sx,decreasing=T)
	n50x <-subset(n50x,
		cumsum(n50x)>= (tx/2))[1]

	return(c(Nx,nx,tx/10^6,max(sx)/10^6,
		n50x/10^6,covx,gcx))
		
}))		

	colnames(STATSmf)<-c("Isolate-Method-Library",
		"n.Scaffolds","Total.size(Mb)","Max.size(Mb)",
		"N50(Mb)","Coverage","GC")

setwd(paste(p0,"Final-assemblies",sep="/"))
	
	write.table(STATSmf,paste("Sum-Final-assemblies-Megahit-filtered.txt",sep=""),row.names=F)

###STEP B ###
#"Manual" cleaning of contaminated Megahit assemblies


library(seqinr)
library(parallel)
library(ade4)

setwd(paste(p0,"Final-assemblies",sep="/"))

#List of selected assemblies
ASS<-read.table("List-final-assemblies.txt",header=T)

#Remove uncorrect assemlies (too low quality)
ASS<-subset(ASS,ASS[,8]!="none")

Mi<-grep("Megahit",ASS[,8])

STATSmf<-read.table("Sum-Final-assemblies-Megahit-filtered.txt",header=T)

Ci<-subset(Mi,as.numeric(STATSmf[,3])>8)

####Group all contaminated assemblies (filtered) together and perform scaffold filtering using
#1 coverage
#2 library id
#3 tetranucleotide frequency

#Kmer size (+1)
K=3

SEQall<-unlist(sapply(Ci,function(x){
	Nx<-ASS[x,8]
	print(Nx)
	seqx<-read.fasta(paste(Nx,"-filt.fa",sep=""))
	return(sapply(c(1:length(seqx)),function(x){
		return(paste(tolower(unlist(seqx[[x]])),
			collapse=""))
	}))
}))

#Assembly index per scaffold
SEQi<-unlist(sapply(c(1:length(Ci)),function(x){
	Nx<-ASS[Ci[x],8]
	print(Nx)
	seqx<-read.fasta(paste(Nx,"-filt.fa",sep=""))
	return(seq(x,x,length.out=length(seqx)))
}))

COLci<-c("black","grey","blue4","blue","cyan",
	"green4","green","yellow2","orange","red",
	"red4","purple","pink")
	
COLci2<-rgb(col2rgb(COLci)[1,]/255,
			col2rgb(COLci)[2,]/255,
			col2rgb(COLci)[3,]/255,0.2)	
	
#Coverage
COVi<-unlist(sapply(c(1:length(Ci)),function(x){
	Nx<-ASS[Ci[x],8]
	print(Nx)
	seqx<-read.fasta(paste(Nx,"-filt.fa",sep=""))
	nx=length(seqx)
	COVx<-as.numeric(sapply(c(1: nx),
		function(x){
		att<-strsplit(attributes(
			seqx[[x]])$Annot,
			split=" ")[[1]][3:4]
		cov<-as.numeric(strsplit(att[1],
			split="=")[[1]][2])
	}))	
	return(COVx)
}))	

#number of scaffolds
Ns<-length(SEQall)
	
#scaffold size
SC<-nchar(SEQall)

####scaffolds  Filtering according to coverage (Isolates  84 can't be saved: two genomes with too similar coverage and tetranuc. composition)

tCOV<-c(40,40,20,20,20,NA,40,50,40,"12;25",30,40,30)

Fi<-subset(c(1:length(Ci)),is.na(tCOV)==F)

STATfc<-t(sapply(Fi,function(x){
	
	NX=ASS[Ci[x],8]
	
	seqx<-read.fasta(paste(NX,"-filt.fa",sep=""))
	
	ANNOT<-sapply(c(1: length(seqx)),
		function(x){
		return(attributes(
			seqx[[x]])$Annot)
	})
	
	ix<-subset(c(1:length(SEQi)), SEQi==x)
	scx<-SC[ix]
	covx<-COVi[ix]
		
	colx<-COLci[x]
	colx2<-COLci2[x]
	
	tcov=as.numeric(unlist(strsplit(
		tCOV[x],split=";")))
		
	nz<-c("a","b","c")[c(1:(length(tcov)+1))]
		
	tcov<-c(tcov,max(covx))
	
	#assembly 1 (min coverage)
	A1<-subset(ix, covx<=tcov[1])
	#assembly 2 (mid coverage)
	A2<-intersect(subset(ix, covx>tcov[1]),
		subset(ix, covx<=tcov[2]))
	#assembly 3 (max coverage)
	A3<-subset(ix, covx>tcov[2])

	z<-ifelse( covx<=tcov[1],1,
		ifelse( covx>tcov[2],3,2))
			
	#scaffold GC content
	GCx<-sapply(c(1: length(seqx)),function(x){
		return(GC(seqx[[x]]))
	})

	pdf(paste("Filtering-Coverage-",ASS[Ci[x],8],
		".pdf",sep=""))
		plot(as.factor(z), covx,log="y",
			ylab="Coverage (log scale)",
			xlab="Putative genomes",
			col=colx2,border=colx,
			main= NX,names=nz,
			las=1,cex.axis=0.8)
			
		plot(as.factor(z), GCx,log="",
			ylab="GC content",
			xlab="Putative genomes",
			col=colx2,border=colx,
			main= NX,names=nz,
			las=1,cex.axis=0.8)			
	
	sapply(c(1: NF),function(x){
		
		Ax=x
		pcax<-PCA$li[ix,Ax]	
		plot(as.factor(z), pcax,log="",
			ylab=paste("Tetranuc. Comp. Axis ",
				Ax, " (", PART[Ax]," %)",sep=""),
			xlab="Putative genomes",
			col=colx2,border=colx,
			main= NX,names=nz,
			las=1,cex.axis=0.8)
	})
	
	dev.off()
	
	#Statistics per genomes
	
	STATs<-t(sapply(c(1:length(nz)),function(x){
		
		vx<-subset(c(1:length(seqx)),z==x)
		
		seqn<-sapply(vx,function(x){
			return(unlist(seqx[[x]]))
		})
		
		namex<-subset(ANNOT,z==x)
		
		nx=length(vx)
		gcn=subset(GCx,z==x)
		sx= subset(scx,z==x)
		tx=sum(sx)
		covn=subset(covx,z==x)
	
		#average GC content
		gcx=sum(gcn* sx)/tx
		#average coverage
		covx=sum(covn* sx)/tx
		#N50	
		n50x<-sort(sx,decreasing=T)
		n50x <-subset(n50x,
			cumsum(n50x)>= (tx/2))[1]

		NEW=paste(NX,"-filt=",nz[x],sep="")	
		
		write.fasta(seqn,names= namex,
			paste(NEW,".fa",sep=""))
		
		return(c(NEW,nx,tx/10^6,max(sx)/10^6,
			n50x/10^6,covx,gcx))
		
	}))
	
	#colnames(STATs)<-c("Putative-genome",
	#	"n.Scaffolds","Total.size(Mb)","Max.size(Mb)",
	#	"N50(Mb)","Coverage","GC")
	#write.table(STATs,paste("Sum-Filtering-cov-", 
	#	NX,".txt",sep=""),row.names=F)
	
	return(t(STATs))
	
}))	


STATfc <-matrix(unlist(STATfc),ncol=7,byrow=T)

colnames(STATfc)<-c("Putative-genome",
	"n.Scaffolds","Total.size(Mb)","Max.size(Mb)",
	"N50(Mb)","Coverage","GC")

write.table(STATfc,paste("Sum-Filtering-conta.txt",sep=""),row.names=F)

