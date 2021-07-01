####Code by Vivek Swarup, PhD UC Irvine #############
###contact vswarup@uci.edu################
##please cite doi: 10.1038/s41386-020-00839-w when using the code
#####################################################
#####################################################
rm(list=ls())
library(igraph)
library(impute)
library(qgraph)
library(qgraph)
library(DescTools)

##Processing data file
data=read.csv('cFOS_Intensity.csv') ##a csv file with samples in rows and cFOS intensity in columns.
##Columns 1:3 has info about ratID, sex and treatment group. All other columns correspons to cFOS intensity for that sample in different brain regions

targets=data[,c(1:3)] ##metadata for the mice
data=data[,-c(1:3)] ## cFOS intensity without mice metadata
rownames(data)=targets$RatID

##Imputing data with 10% missing values or less
data.impute<-impute.knn(as.matrix(data))$data
thisdat.HTSC <- t(scale(t(data.impute),scale=F)) ## Centers the mean of all genes - this means the PCA gives us the eigenvectors of the geneXgene covariance matrix, allowing us to assess the proportion of variance each component contributes to the data
orderRegion=c('BNST','CeA','BLA','CAs','DGs','PLC','ILC','Nco','NSh','RVP','CVP','MHb','LHb','Tvta')
thisdat.HTSC=thisdat.HTSC[,orderRegion]

dataVeh=subset(thisdat.HTSC,targets$Group=='Veh')
data0.5THC=subset(thisdat.HTSC,targets$Group=='0.5 THC')
dataTHC=subset(thisdat.HTSC,targets$Group=='5 THC')

cormatVeh=cor(dataVeh,method='s')
cormat0.5THC=cor(data0.5THC,method='s')
cormatTHC=cor(dataTHC,method='s')

##Making correlation plot
g1 <- graph.adjacency(as.matrix(cormatVeh),mode="undirected",weighted=T,diag=FALSE)
layoutVeh <- layout.circle(g1)

g1 <- graph.adjacency(as.matrix(cormat0.5THC),mode="undirected",weighted=T,diag=FALSE)
layout0.5THC <- layout.circle(g1)

g1 <- graph.adjacency(as.matrix(cormatTHC),mode="undirected",weighted=T,diag=FALSE)
layoutTHC <- layout.circle(g1)

pdf("correlationAllSamples.pdf",width=18,height=6)
  par(mfrow=c(1,3))
  qgraph(cormatVeh, shape="circle", posCol="green", negCol="red", layout=layoutVeh, vsize=10)
  qgraph(cormat0.5THC, shape="circle", posCol="green", negCol="red", layout=layout0.5THC, vsize=10)
  qgraph(cormatTHC, shape="circle", posCol="green", negCol="red", layout=layoutTHC, vsize=10)
dev.off()


###Sex Separate

dataVeh.M=subset(thisdat.HTSC,targets$Group=='Veh'&targets$sex=='M')
dataVeh.F=subset(thisdat.HTSC,targets$Group=='Veh'&targets$sex=='F')

data0.5THC.M=subset(thisdat.HTSC,targets$Group=='0.5 THC'&targets$sex=='M')
data0.5THC.F=subset(thisdat.HTSC,targets$Group=='0.5 THC'&targets$sex=='F')

dataTHC.M=subset(thisdat.HTSC,targets$Group=='5 THC'&targets$sex=='M')
dataTHC.F=subset(thisdat.HTSC,targets$Group=='5 THC'&targets$sex=='F')

cormatVeh.M=cor(dataVeh.M,method='s')
cormatVeh.F=cor(dataVeh.F,method='s')

cormat0.5THC.M=cor(data0.5THC.M,method='s')
cormat0.5THC.F=cor(data0.5THC.F,method='s')

cormatTHC.M=cor(dataTHC.M,method='s')
cormatTHC.F=cor(dataTHC.F,method='s')


cormatVeh.M.pos=cormatVeh.M
cormatVeh.M.pos[cormatVeh.M.pos<0]<-0

cormatVeh.F.pos=cormatVeh.F
cormatVeh.F.pos[cormatVeh.F.pos<0]<-0

cormat0.5THC.M.pos=cormat0.5THC.M
cormat0.5THC.M.pos[cormat0.5THC.M.pos<0]<-0

cormat0.5THC.F.pos=cormat0.5THC.F
cormat0.5THC.F.pos[cormat0.5THC.F.pos<0]<-0

cormatTHC.M.pos=cormatTHC.M
cormatTHC.M.pos[cormatTHC.M.pos<0]<-0

cormatTHC.F.pos=cormatTHC.F
cormatTHC.F.pos[cormatTHC.F.pos<0]<-0

#######
cormatVeh.M.neg=cormatVeh.M
cormatVeh.M.neg[cormatVeh.M.neg>0]<-0

cormatVeh.F.neg=cormatVeh.F
cormatVeh.F.neg[cormatVeh.F.neg>0]<-0

cormat0.5THC.M.neg=cormat0.5THC.M
cormat0.5THC.M.neg[cormat0.5THC.M.neg>0]<-0

cormat0.5THC.F.neg=cormat0.5THC.F
cormat0.5THC.F.neg[cormat0.5THC.F.neg>0]<-0

cormatTHC.M.neg=cormatTHC.M
cormatTHC.M.neg[cormatTHC.M.neg>0]<-0

cormatTHC.F.neg=cormatTHC.F
cormatTHC.F.neg[cormatTHC.F.neg>0]<-0


b1<-bartlett.test(cormatVeh.M,cormat0.5THC.M)
b2<-bartlett.test(cormatVeh.M,cormatTHC.M)

b3<-bartlett.test(cormatVeh.F,cormat0.5THC.F)
b4<-bartlett.test(cormatVeh.F,cormatTHC.F)

b5<-bartlett.test(cormatVeh,cormat0.5THC)
b6<-bartlett.test(cormatVeh,cormatTHC)

b7<-bartlett.test(thisdat.HTSC[,1]~targets$Group)



##get the order of nodes using g1[order(degree(g1))] and then rearrage

g1 <- graph.adjacency(as.matrix(cormatVeh.M),mode="undirected",weighted=T,diag=FALSE)
layoutVeh.M <- layout.circle(g1)
g1 <- graph.adjacency(as.matrix(cormatVeh.F),mode="undirected",weighted=T,diag=FALSE)
layoutVeh.F <- layout.circle(g1)

g1 <- graph.adjacency(as.matrix(cormat0.5THC.M),mode="undirected",weighted=T,diag=FALSE)
layout0.5THC.M <- layout.circle(g1)
g1 <- graph.adjacency(as.matrix(cormat0.5THC.F),mode="undirected",weighted=T,diag=FALSE)
layout0.5THC.F <- layout.circle(g1)

g1 <- graph.adjacency(as.matrix(cormatTHC.M),mode="undirected",weighted=T,diag=FALSE)
layoutTHC.M <- layout.circle(g1)
g1 <- graph.adjacency(as.matrix(cormatTHC.F),mode="undirected",weighted=T,diag=FALSE)
layoutTHC.F <- layout.circle(g1)

##########Pos

g1 <- graph.adjacency(as.matrix(cormatVeh.M.pos),mode="undirected",weighted=T,diag=FALSE)
layoutVeh.M.pos <- layout.circle(g1)
g1 <- graph.adjacency(as.matrix(cormatVeh.F.pos),mode="undirected",weighted=T,diag=FALSE)
layoutVeh.F.pos <- layout.circle(g1)

g1 <- graph.adjacency(as.matrix(cormat0.5THC.M.pos),mode="undirected",weighted=T,diag=FALSE)
layout0.5THC.M.pos <- layout.circle(g1)
g1 <- graph.adjacency(as.matrix(cormat0.5THC.F.pos),mode="undirected",weighted=T,diag=FALSE)
layout0.5THC.F.pos <- layout.circle(g1)

g1 <- graph.adjacency(as.matrix(cormatTHC.M.pos),mode="undirected",weighted=T,diag=FALSE)
layoutTHC.M.pos <- layout.circle(g1)
g1 <- graph.adjacency(as.matrix(cormatTHC.F.pos),mode="undirected",weighted=T,diag=FALSE)
layoutTHC.F.pos <- layout.circle(g1)

######Neg

g1 <- graph.adjacency(as.matrix(cormatVeh.M.neg),mode="undirected",weighted=T,diag=FALSE)
layoutVeh.M.neg <- layout.circle(g1)
g1 <- graph.adjacency(as.matrix(cormatVeh.F.neg),mode="undirected",weighted=T,diag=FALSE)
layoutVeh.F.neg <- layout.circle(g1)

g1 <- graph.adjacency(as.matrix(cormat0.5THC.M.neg),mode="undirected",weighted=T,diag=FALSE)
layout0.5THC.M.neg <- layout.circle(g1)
g1 <- graph.adjacency(as.matrix(cormat0.5THC.F.neg),mode="undirected",weighted=T,diag=FALSE)
layout0.5THC.F.neg <- layout.circle(g1)

g1 <- graph.adjacency(as.matrix(cormatTHC.M.neg),mode="undirected",weighted=T,diag=FALSE)
layoutTHC.M.neg <- layout.circle(g1)
g1 <- graph.adjacency(as.matrix(cormatTHC.F.neg),mode="undirected",weighted=T,diag=FALSE)
layoutTHC.F.neg <- layout.circle(g1)


pdf("correlation_SexSeparate.pdf",width=18,height=6)
  par(mfrow=c(1,3))
  qgraph(cormatVeh.M, shape="circle", posCol="green", negCol="red", layout=layoutVeh.M, vsize=10)
  qgraph(cormat0.5THC.M, shape="circle", posCol="green", negCol="red", layout=layout0.5THC.M, vsize=10)
  qgraph(cormatTHC.M, shape="circle", posCol="green", negCol="red", layout=layoutTHC.M, vsize=10)

  qgraph(cormatVeh.F, shape="circle", posCol="green", negCol="red", layout=layoutVeh.F, vsize=10)
  qgraph(cormat0.5THC.F, shape="circle", posCol="green", negCol="red", layout=layout0.5THC.F, vsize=10)
  qgraph(cormatTHC.F, shape="circle", posCol="green", negCol="red", layout=layoutTHC.F, vsize=10)

dev.off()



pdf("Pos_correlation_SexSeparate.pdf",width=18,height=6)
  par(mfrow=c(1,3))
  qgraph(cormatVeh.M.pos, shape="circle", posCol="green", negCol="red", layout=layoutVeh.M.pos, vsize=10)
  qgraph(cormat0.5THC.M.pos, shape="circle", posCol="green", negCol="red", layout=layout0.5THC.M.pos, vsize=10)
  qgraph(cormatTHC.M.pos, shape="circle", posCol="green", negCol="red", layout=layoutTHC.M.pos, vsize=10)

  qgraph(cormatVeh.F.pos, shape="circle", posCol="green", negCol="red", layout=layoutVeh.F.pos, vsize=10)
  qgraph(cormat0.5THC.F.pos, shape="circle", posCol="green", negCol="red", layout=layout0.5THC.F.pos, vsize=10)
  qgraph(cormatTHC.F.pos, shape="circle", posCol="green", negCol="red", layout=layoutTHC.F.pos, vsize=10)

dev.off()


pdf("Neg_correlation_SexSeparate.pdf",width=18,height=6)
  par(mfrow=c(1,3))
  qgraph(cormatVeh.M.neg, shape="circle", posCol="green", negCol="red", layout=layoutVeh.M.neg, vsize=10)
  qgraph(cormat0.5THC.M.neg, shape="circle", posCol="green", negCol="red", layout=layout0.5THC.M.neg, vsize=10)
  qgraph(cormatTHC.M.neg, shape="circle", posCol="green", negCol="red", layout=layoutTHC.M.neg, vsize=10)

  qgraph(cormatVeh.F.neg, shape="circle", posCol="green", negCol="red", layout=layoutVeh.F.neg, vsize=10)
  qgraph(cormat0.5THC.F.neg, shape="circle", posCol="green", negCol="red", layout=layout0.5THC.F.neg, vsize=10)
  qgraph(cormatTHC.F.neg, shape="circle", posCol="green", negCol="red", layout=layoutTHC.F.neg, vsize=10)

dev.off()




a=FisherZ(cormatVeh.M.pos)
b=FisherZ(cormat0.5THC.M.pos)
c=b-a
c[is.na(c)] <- 0

g1 <- graph.adjacency(as.matrix(c),mode="undirected",weighted=T,diag=FALSE)
layoutVeh0.5THC_M <- layout.circle(g1)

a1=FisherZ(cormatVeh.M.pos)
b1=FisherZ(cormatTHC.M.pos)
c1=b1-a1
c1[is.na(c1)] <- 0

g1 <- graph.adjacency(as.matrix(c1),mode="undirected",weighted=T,diag=FALSE)
layoutVehTHC_M <- layout.circle(g1)



a2=FisherZ(cormatVeh.F.pos)
b2=FisherZ(cormat0.5THC.F.pos)
c2=b2-a2
c2[is.na(c2)] <- 0

g1 <- graph.adjacency(as.matrix(c2),mode="undirected",weighted=T,diag=FALSE)
layoutVeh0.5THC_F <- layout.circle(g1)

a3=FisherZ(cormatVeh.F.pos)
b3=FisherZ(cormatTHC.F.pos)
c3=b3-a3
c3[is.na(c3)] <- 0

g1 <- graph.adjacency(as.matrix(c3),mode="undirected",weighted=T,diag=FALSE)
layoutVehTHC_F <- layout.circle(g1)


pdf("FisherZ_Vehicle_Compared.pdf",width=12,height=6)
  par(mfrow=c(1,2))
  qgraph(c, shape="circle", posCol="green", negCol="red", layout=layoutVeh0.5THC_M, vsize=10)
  qgraph(c1, shape="circle", posCol="green", negCol="red", layout=layoutVehTHC_M, vsize=10)

  qgraph(c2, shape="circle", posCol="green", negCol="red", layout=layoutVeh0.5THC_F, vsize=10)
  qgraph(c3, shape="circle", posCol="green", negCol="red", layout=layoutVehTHC_F, vsize=10)
dev.off()

Val_0.5M=sum(rowSums(c))
Val_5M=sum(rowSums(c1))
Val_0.5F=sum(rowSums(c2))
Val_5F=sum(rowSums(c3))

pdf('barplot_SummedZscores.pdf',height=,width=6)
barplot(c(Val_0.5F,Val_5F,Val_0.5M,Val_5M),names.arg=c('0.5F','5F','0.5M','5M'),ylab='Summed Z.scores',col='lightblue')
dev.off()



write.csv(c,'Z.0_5THC_M_allregions.csv')
write.csv(c1,'Z.5THC_M_allregions.csv')
write.csv(c2,'Z.0_5THC_F_allregions.csv')
write.csv(c3,'Z.5THC_F_allregions.csv')

##z-score to pvalue
pvalue_0.5THCM=2*pnorm(-abs(c))
pvalue_5THCM=2*pnorm(-abs(c1))
pvalue_0.5THCF=2*pnorm(-abs(c2))
pvalue_5THCF=2*pnorm(-abs(c3))

write.csv(pvalue_0.5THCM,'Pvalue.0_5THC_M_allregions.csv')
write.csv(pvalue_5THCM,'Pvalue.5THC_M_allregions.csv')
write.csv(pvalue_0.5THCF,'Pvalue.0_5THC_F_allregions.csv')
write.csv(pvalue_5THCF,'Pvalue.5THC_F_allregions.csv')





### STats
a=FisherZ(cormatVeh.F.pos)
b=FisherZ(cormatTHC.F.pos)
d=FisherZ(cormat0.5THC.F.pos)

a[a==Inf] <- 0
b[b==Inf] <- 0
d[d==Inf]<-0

t.test(b,a) ##p-value = 1.138e-05
t.test(d,a) ##p=0.01253

#### in males
a=FisherZ(cormatVeh.M.pos)
b=FisherZ(cormatTHC.M.pos)
e=FisherZ(cormat0.5THC.M.pos)

a[a==Inf] <- 0
b[b==Inf] <- 0
e[e==Inf]<-0

t.test(b,a) ##p-value = 0.01062
t.test(e,a) ##p=0.6867


a1[a1==Inf] <- 0
b1[b1==Inf] <- 0
save(list=ls(),file='Analysis.rda')

sapply(1:ncol(a1),function(x) t.test(a1[,x],b1[,x])$p.value)




group=factor(targets$Group,c('Veh','0.5 THC','5 THC'))

pdf("correlation2.pdf",width=12,height=6)
  par(mfrow=c(1,2))
  qgraph(cormatVeh, shape="circle", posCol="green", negCol="red", layout="spring", vsize=10,threshold=0.2)
  qgraph(cormatTHC, shape="circle", posCol="green", negCol="red", layout="spring", vsize=10,threshold=0.2)
  qgraph(cormatVeh1, shape="circle", posCol="green", negCol="red", layout="spring", vsize=10,threshold=0.2)
  qgraph(cormatTHC1, shape="circle", posCol="green", negCol="red", layout="spring", vsize=10,threshold=0.2)
  #boxplot(as.numeric(data.impute[,"BLA"])~group,col=c("red","blue","green"))
  #boxplot(as.numeric(data.impute[,"CVP"])~group,col=c("red","blue","green"))

dev.off()



###Merge Brain regions
combined_df <- data.frame()
for(i in 1:nrow(data.impute)){
  thisDat=data.impute[i,]
  ExAmy=mean(thisDat[1],thisDat[2],thisDat[3])
  Hippo=mean(thisDat[13],thisDat[14])

  mPFC=mean(thisDat[4],thisDat[5])
  Nac=mean(thisDat[6],thisDat[7])
  VP=mean(thisDat[8],thisDat[9])
  EpThal=mean(thisDat[10],thisDat[11],thisDat[12])
  combMat=c(ExAmy,Hippo,mPFC,Nac,VP,EpThal)
  combined_df <- rbind(combined_df, combMat)
}
rownames(combined_df)=rownames(data.impute)
colnames(combined_df)=c('ExAmy','Hippo','mPFC','Nac','VP','EpThal')



thisdat.HTSC <- t(scale(t(combined_df),scale=F)) ## Centers the mean of all genes - this means the PCA gives us the eigenvectors of the geneXgene covariance matrix, allowing us to assess the proportion of variance each component contributes to the data


###Sex Separate

dataVeh.M=subset(thisdat.HTSC,targets$Group=='Veh'&targets$sex=='M')
dataVeh.F=subset(thisdat.HTSC,targets$Group=='Veh'&targets$sex=='F')

data0.5THC.M=subset(thisdat.HTSC,targets$Group=='0.5 THC'&targets$sex=='M')
data0.5THC.F=subset(thisdat.HTSC,targets$Group=='0.5 THC'&targets$sex=='F')

dataTHC.M=subset(thisdat.HTSC,targets$Group=='5 THC'&targets$sex=='M')
dataTHC.F=subset(thisdat.HTSC,targets$Group=='5 THC'&targets$sex=='F')

cormatVeh.M=cor(dataVeh.M,method='s')
cormatVeh.F=cor(dataVeh.F,method='s')

cormat0.5THC.M=cor(data0.5THC.M,method='s')
cormat0.5THC.F=cor(data0.5THC.F,method='s')

cormatTHC.M=cor(dataTHC.M,method='s')
cormatTHC.F=cor(dataTHC.F,method='s')



cormatVeh.M.pos=cormatVeh.M
cormatVeh.M.pos[cormatVeh.M.pos<0]<-0

cormatVeh.F.pos=cormatVeh.F
cormatVeh.F.pos[cormatVeh.F.pos<0]<-0

cormat0.5THC.M.pos=cormat0.5THC.M
cormat0.5THC.M.pos[cormat0.5THC.M.pos<0]<-0

cormat0.5THC.F.pos=cormat0.5THC.F
cormat0.5THC.F.pos[cormat0.5THC.F.pos<0]<-0

cormatTHC.M.pos=cormatTHC.M
cormatTHC.M.pos[cormatTHC.M.pos<0]<-0

cormatTHC.F.pos=cormatTHC.F
cormatTHC.F.pos[cormatTHC.F.pos<0]<-0






a=FisherZ(cormatVeh.M.pos)
b=FisherZ(cormat0.5THC.M.pos)
c=b-a
c[is.na(c)] <- 0

g1 <- graph.adjacency(as.matrix(c),mode="undirected",weighted=T,diag=FALSE)
layoutVeh0.5THC_M <- layout.circle(g1)

a1=FisherZ(cormatVeh.M.pos)
b1=FisherZ(cormatTHC.M.pos)
c1=b1-a1
c1[is.na(c1)] <- 0

g1 <- graph.adjacency(as.matrix(c1),mode="undirected",weighted=T,diag=FALSE)
layoutVehTHC_M <- layout.circle(g1)



a2=FisherZ(cormatVeh.F.pos)
b2=FisherZ(cormat0.5THC.F.pos)
c2=b2-a2
c2[is.na(c2)] <- 0

g1 <- graph.adjacency(as.matrix(c2),mode="undirected",weighted=T,diag=FALSE)
layoutVeh0.5THC_F <- layout.circle(g1)

a3=FisherZ(cormatVeh.F.pos)
b3=FisherZ(cormatTHC.F.pos)
c3=b3-a3
c3[is.na(c3)] <- 0

g1 <- graph.adjacency(as.matrix(c3),mode="undirected",weighted=T,diag=FALSE)
layoutVehTHC_F <- layout.circle(g1)

write.csv(c,'Z.0_5THC_M.csv')
write.csv(c1,'Z.5THC_M.csv')
write.csv(c2,'Z.0_5THC_F.csv')
write.csv(c3,'Z.5THC_F.csv')

##z-score to pvalue
pvalue_0.5THCM=2*pnorm(-abs(c))
pvalue_5THCM=2*pnorm(-abs(c1))
pvalue_0.5THCF=2*pnorm(-abs(c2))
pvalue_5THCF=2*pnorm(-abs(c3))

write.csv(pvalue_0.5THCM,'Pvalue.0_5THC_M.csv')
write.csv(pvalue_5THCM,'Pvalue.5THC_M.csv')
write.csv(pvalue_0.5THCF,'Pvalue.0_5THC_F.csv')
write.csv(pvalue_5THCF,'Pvalue.5THC_F.csv')


pdf("FisherZ_Vehicle_Compared_MeanRegions.pdf",width=12,height=6)
  par(mfrow=c(1,2))
  qgraph(c, shape="circle", posCol="green", negCol="red", layout=layoutVeh0.5THC_M, vsize=10)
  qgraph(c1, shape="circle", posCol="green", negCol="red", layout=layoutVehTHC_M, vsize=10)

  qgraph(c2, shape="circle", posCol="green", negCol="red", layout=layoutVeh0.5THC_F, vsize=10)
  qgraph(c3, shape="circle", posCol="green", negCol="red", layout=layoutVehTHC_F, vsize=10)
dev.off()
