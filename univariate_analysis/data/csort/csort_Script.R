library(stringr)
library(dplyr)
library(Rtsne)
library("RColorBrewer")
library(Hmisc)
library(viridis)
library(corrplot)

setwd("/Users/jdmills/Dropbox/Csort/")

####After correction: age corrected
expression<-read.delim(file="RNAseq_Results_AC.txt",header=TRUE)

####Before age correction
expression2<-read.delim(file="RNAseq_Results.txt",header=TRUE)
row.names(expression2)<-expression2$X
expression2<-expression2[,c(2:291)]

csort<-read_csv("csort/Cell_Proportions_RNASeq_James.csv")

###Filter the file out of the RNASeq sample_data

sample_data<-subset(sample_data,sample_data$GSID..RNA.seq.raw.file.!="NA")
sample_data$GSID..RNA.seq.raw.file.<-sub("^","X",sample_data$GSID..RNA.seq.raw.file.)
sample_data$GSID..RNA.seq.raw.file.<-gsub("-",".",sample_data$GSID..RNA.seq.raw.file.)

sample_data<-sample_data[order(sample_data$GSID..RNA.seq.raw.file.),]



sample_data<-subset(sample_data,sample_data$qualified_for_analysis!="no")
sample_data<-sample_data[-c(236:242),]
sample_data<-sample_data[-c(243:249),]
sample_data<-sample_data[-c(258:263),]
sample_data<-sample_data[-c(261:266),]

row.names(sample_data)<-sample_data$GSID..RNA.seq.raw.file.

temp<-colnames(expression2)
temp<-as.data.frame(temp)


sample_data<-merge(temp,sample_data,by.x=1,by.y=0)
row.names(sample_data)<-sample_data$temp

row.names(sample_data)==colnames(expression2)

temp<-row.names(sample_data)

csort<-select(csort,!!temp)

####Gene list for RNA-Seq

genes<-read.csv(file="names_of_selected_RNAs.csv",header=TRUE)
expression2<-merge(genes,expression2,by.x=1,by.y=0)
row.names(expression2)<-expression2$names
expression2<-expression2[,c(2:291)]

###Plot tSNE: colour by age

expression<-expression2
expression<-expression[complete.cases(expression),]

vst_remove<-t(expression)

vst_remove <- Rtsne(vst_remove,perplexity = 5)

dimRed<-vst_remove$Y

rownames(dimRed)<-rownames(sample_data)
x2<-merge(dimRed,sample_data,by=0)


#####First clustering by timepoint:
#plot rtsne
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

A<-x2[x2$Dev_Age_Grouping_new=="0-10 weeks",]
B<-x2[x2$Dev_Age_Grouping_new=="11-40 weeks",]
C<-x2[x2$Dev_Age_Grouping_new=="> 40 weeks",]


A<-A[,c(2:3)]
B<-B[,c(2:3)]
C<-C[,c(2:3)]

brewer.pal(n=8,"Dark2")


tiff("./TSNE_RNAexpression.tiff", res=300,width=3000,height=3000)
par(mai = c(0.8, 1, 0.8,2.1 )) #Bottom,left,top,right
plot(A,col="#1B9E77" ,pch=1,xlim=c(-50,50),ylim=c(-50,50),main="",cex=2.5,cex.lab=2,cex.axis=1.5,cex.main=2.5)
par(new=T)
plot(B,col="#7570B3",pch=2,xlim=c(-50,50),ylim=c(-50,50),cex=2.5,cex.lab=2,cex.axis=1.5)
par(new=T)
plot(C,col="#D95F02",pch=5,xlim=c(-50,50),ylim=c(-50,50),cex=2.5,cex.lab=2,cex.axis=1.5)
par(new=T)


par(fig=c(0,0.6,0,1),new=TRUE)
labels<-c("0-10 weeks","11-40 weeks","> 40 weeks")
add_legend("right",labels,col=c("#1B9E77", "#7570B3","#D95F02"),pch=c(1,2,5),box.col=NULL,bty="n",pt.cex=2,cex=1.5)


dev.off()
#######Lets do the same thing except with csort results


expression<-csort

vst_remove<-t(expression)

vst_remove <- Rtsne(vst_remove,perplexity = 5)

dimRed<-vst_remove$Y

rownames(dimRed)<-rownames(sample_data)
x2<-merge(dimRed,sample_data,by=0)


#####First clustering by timepoint:
#plot rtsne
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

A<-x2[x2$Dev_Age_Grouping_new=="0-10 weeks",]
B<-x2[x2$Dev_Age_Grouping_new=="11-40 weeks",]
C<-x2[x2$Dev_Age_Grouping_new=="> 40 weeks",]


A<-A[,c(2:3)]
B<-B[,c(2:3)]
C<-C[,c(2:3)]

brewer.pal(n=8,"Dark2")


tiff("./TSNE_Csort.tiff", res=300,width=3000,height=3000)
par(mai = c(0.8, 1, 0.8,2.1 )) #Bottom,left,top,right
plot(A,col="#1B9E77" ,pch=1,xlim=c(-50,50),ylim=c(-50,50),main="",cex=2.5,cex.lab=2,cex.axis=1.5,cex.main=2.5)
par(new=T)
plot(B,col="#7570B3",pch=2,xlim=c(-50,50),ylim=c(-50,50),cex=2.5,cex.lab=2,cex.axis=1.5)
par(new=T)
plot(C,col="#D95F02",pch=5,xlim=c(-50,50),ylim=c(-50,50),cex=2.5,cex.lab=2,cex.axis=1.5)
par(new=T)


par(fig=c(0,0.6,0,1),new=TRUE)
labels<-c("0-10 weeks","11-40 weeks","> 40 weeks")
add_legend("right",labels,col=c("#1B9E77", "#7570B3","#D95F02"),pch=c(1,2,5),box.col=NULL,bty="n",pt.cex=2,cex=1.5)


dev.off()


#Lets work out the PCA for the expression data

pca.cal<-prcomp(t(expression))
pca.cal<-pca.cal$x
pca.cal<-as.data.frame(pca.cal)


A<-x2[x2$Dev_Age_Grouping_new=="0-10 weeks",]
B<-x2[x2$Dev_Age_Grouping_new=="11-40 weeks",]
C<-x2[x2$Dev_Age_Grouping_new=="> 40 weeks",]

FC1<-subset(sample_data,sample_data$Dev_Age_Grouping_new=="0-10 weeks")
FC1<-merge(FC1,pca.cal,by=0)
row.names(FC1)<-FC1$Row.names
FC1<-FC1[,c(23,79:88)]
FC2<-subset(sample_data,sample_data$Dev_Age_Grouping_new=="11-40 weeks")
FC2<-merge(FC2,pca.cal,by=0)
row.names(FC2)<-FC2$Row.names
FC2<-FC2[,c(23,79:88)]
FC3<-subset(sample_data,sample_data$Dev_Age_Grouping_new=="> 40 weeks")
FC3<-merge(FC3,pca.cal,by=0)
row.names(FC3)<-FC3$Row.names
FC3<-FC3[,c(23,79:88)]

tiff("./PCA_Seq.tiff", res=100,width=1000,height=1000)
par(mai = c(1, 1, 1,2.2 )) 
plot(FC1[,2],FC1[,3],col="#1B9E77",pch=1,cex=2,cex.lab=2,cex.axis=1.5,xlab="PC1 (25.18%)",ylim=c(-100,100),xlim=c(-100,150), ylab = "PC2 (16.55%)", main = "")
par(new=T)
plot(FC2[,2],FC2[,3],col="#7570B3",pch=2,cex=2,cex.axis=1.5,xlab="",ylab="", main = "",ylim=c(-100,100),xlim=c(-100,150))
par(new=T)
plot(FC3[,2],FC3[,3],col="#D95F02",pch=5,cex=2,cex.axis=1.5,xlab="",ylab="", main = "",ylim=c(-100,100),xlim=c(-100,150))
par(new=T)
par(fig=c(0,0.8,0,1),new=TRUE)
labels<-c("0-10 weeks","11-40 weeks","> 40 weeks")
add_legend("right",labels,col=c("#1B9E77", "#7570B3","#D95F02"),pch=c(1,2,5),box.col=NULL,bty="n",pt.cex=2.5,cex=1.4)
dev.off()

corr<-merge(sample_data,pca.cal,by=0)
row.names(corr)<-corr$Row.names
corr<-corr[,c(23,79:88)]

cortt<-t(csort)
cortt<-cortt[,-c(8,10,17,20)]
corr<-merge(corr,cortt,by.x=0,by.y = 0)
row.names(corr)<-corr$Row.names
corr<-corr[,c(2:30)]
corr<-as.matrix(corr)
x<-rcorr(corr,type="spearman")




plot<-x$r
pvalues<-x$P

col5<-viridis(10)

corrplot(plot,type="upper",p.mat = pvalues,col = col5,tl.pos="d",tl.col="black")


tiff("./corr_plot.tiff", res=100,width=1000,height=1000)
corrplot(plot,col = col5,tl.col="black")
dev.off()

corrplot(corr,type="upper",addCoef.col="black",col = col5,tl.pos="d",tl.col="black")



                                

  
  


