#####
library(limma)
rt1=read.table("GTExNormalExp.txt",sep="\t",header=T,check.names=F)
rt1=as.matrix(rt1)
rownames(rt1)=rt1[,1]
exp1=rt1[,2:ncol(rt1)]
dimnames1=list(rownames(exp1),colnames(exp1))
data1=matrix(as.numeric(as.matrix(exp1)),nrow=nrow(exp1),dimnames=dimnames1)
data1=avereps(data1)

rt2=read.table("tcgaSymbol.txt",sep="\t",header=T,check.names=F)
rt2=as.matrix(rt2)
rownames(rt2)=rt2[,1]
exp2=rt2[,2:ncol(rt2)]
dimnames2=list(rownames(exp2),colnames(exp2))
data2=matrix(as.numeric(as.matrix(exp2)),nrow=nrow(exp2),dimnames=dimnames2)
data2=avereps(data2)

sameGene=intersect( row.names(data1),row.names(data2) )
data=cbind(data1[sameGene,],data2[sameGene,])

outTab=normalizeBetweenArrays(data)
outTab=rbind(geneNames=colnames(outTab),outTab)
write.table(outTab,file="merge.txt",sep="\t",quote=F,col.names=F)

#######
rt=read.table("merge.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

gene=read.table("gene.txt", header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]

out=rbind(ID=colnames(geneExp),geneExp)
write.table(out,file="MergeGlnExp.txt",sep="\t",quote=F,col.names=F)
############

#引用包
library(limma)
library(reshape2)
library(ggpubr)

logFCfilter=2        
fdrFilter=0.05       
expFile="TCGAGETx.GlnExp.txt"      
geneFile="gene.txt"      

#读取输入文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

gene=read.table(geneFile, header=T, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), row.names(data))
data=data[sameGene,]

gene=read.table("gene.txt", header=F, check.names=F, sep="\t")
sameGene=intersect(as.vector(gene[,1]),rownames(data))
geneExp=data[sameGene,]

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
conNum=length(group[group==1])      
treatNum=length(group[group==0])     
sampleType=c(rep(1,conNum), rep(2,treatNum))

exp=log2(data+1)
exp=as.data.frame(t(exp))
exp=cbind(exp, Type=sampleType)
exp$Type=ifelse(exp$Type==1, "Normal", "Tumor")
data=melt(exp, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

p=ggboxplot(data, x="Gene", y="Expression",fill="Type",
            ylab="Gene expression",
            xlab="",
            legend.title="Type",
            palette = c("CadetBlue3", "LightCoral"),
            width=1)
p=p+rotate_x_text(60)            
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

pdf(file="boxplot.pdf", width=10, height=6)
print(p1)
dev.off()
####################
library(limma)      
expFile="TCGA.GlnExp.txt"     
cliFile="time-OS.txt"            

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
data=data[,group==0]
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(data))
data=t(data)
data=avereps(data)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
out=cbind(cli,data)
out=cbind(id=row.names(out),out)
write.table(out,file="TCGA.expTime.txt",sep="\t",row.names=F,quote=F)

###############
library(survminer)  
library(survival)    
pFilter=0.05           
rt=read.table("TCGAGln.expTime.txt", header=T, sep="\t", check.names=F, row.names=1)  
rt[,3:ncol(rt)]=log2(rt[,3:ncol(rt)]+1)
rt$futime=rt$futime/365

outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
    

    data=rt[,c("futime", "fustat", i)]
    colnames(data)=c("futime", "fustat", "gene")
    res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("gene"))
    res.cat=surv_categorize(res.cut)
    diff=survdiff(Surv(futime, fustat) ~gene,data =res.cat)
    pValue=1-pchisq(diff$chisq, df=1)
    #print(pValue)
    if(pValue<0.001){
      pValue="p<0.001"
    }else{
      pValue=paste0("p=",sprintf("%.03f",pValue))
    }
    fit=survfit(Surv(futime, fustat) ~gene, data = res.cat)
    surPlot=ggsurvplot(fit,
                       data=res.cat,
                       pval=pValue,
                       pval.size=6,
                       legend.title=i,
                       legend.labs=c("high","low"),
                       xlab="Time(years)",
                       palette=c("#E41A1C", "#377EB8"),
                       break.time.by=1,
                       conf.int=T,
                       risk.table=TRUE,
                       risk.table.title="",
                       risk.table.height=.25)
    pdf(file=paste0("survival.",i,".pdf"), onefile=FALSE, width=6, height =5)
    print(surPlot)
    dev.off()
  }
}
write.table(outTab,file="tcga.uniCox.txt",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="tcga.uniSigExp.txt",sep="\t",row.names=F,quote=F)

bioForest=function(coxFile=null, forestFile=null){
  rt <- read.table(coxFile,header=T,sep="\t",check.names=F,row.names=1)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrLow[hrLow<0.001]=0.001
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  height=(nrow(rt)/15)+5.5
  pdf(file=forestFile, width=8, height=height)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,adj=1,)
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  LOGindex=10 
  hrLow = log(as.numeric(hrLow),LOGindex)
  hrHigh = log(as.numeric(hrHigh),LOGindex)
  hr = log(as.numeric(hr),LOGindex)
  xlim = c(floor(min(hrLow,hrHigh)),ceiling(max(hrLow,hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=log(1,LOGindex),col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > log(1,LOGindex), "#E41A1C", "#377EB8")
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
  a1 = axis(1,labels=F,tick=F)
  axis(1,a1,LOGindex^a1)
  dev.off()
}

bioForest(coxFile="TCGA.uniCox.txt", forestFile="forest.pdf")

################
########################
library(limma)
library(pheatmap)
expFile="unicoxgeneExp.txt"    


#??ȡ?????ļ?
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#????????????Ŀ
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
conNum=length(group[group==1])      
treatNum=length(group[group==0])    
sampleType=c(rep(1,conNum), rep(2,treatNum))

sigVec=c()
outTab=data.frame()
for(i in rownames(data)){
  if(sd(data[i,])<0.001){next}
  wilcoxTest=wilcox.test(data[i,] ~ sampleType)
  pvalue=wilcoxTest$p.value
  if(pvalue<0.05){
    Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
    sigVec=c(sigVec, paste0(i, Sig))
    conGeneMeans=mean(data[i,1:conNum])
    treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
    logFC=log2(treatGeneMeans)-log2(conGeneMeans)
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
  }
}

write.table(outTab, file="diff.xls", sep="\t", row.names=F, quote=F)
write.table(outTab, file="diff.txt", sep="\t", row.names=F, quote=F)

exp=data[as.vector(outTab[,1]),]
expOut=rbind(ID=colnames(exp), exp)
write.table(expOut, file="diffGeneExp.txt", sep="\t", col.names=F, quote=F)

exp=log2(exp+0.1)
row.names(exp)=sigVec
Type=c(rep("Normal",conNum),rep("Tumor",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf", width=9, height=6)
pheatmap(exp, 
         annotation=Type, 
         color = colorRampPalette(c(rep("CadetBlue3",10), "white", rep("LightCoral",10)))(50),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames = F,
         show_rownames = T,
         fontsize = 8,
         fontsize_row=8,
         fontsize_col=8)
dev.off()

###########################
library(survival)                
coxFile="tcga.uniCox.txt"        
geneFile="interGenes.txt"        
   
rt=read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
geneRT=read.table(geneFile,header=F,sep="\t",check.names=F)
rt=rt[as.vector(geneRT[,1]),]
gene=rownames(rt)
hr=sprintf("%.3f",rt$"HR")
hrLow=sprintf("%.3f",rt$"HR.95L")
hrHigh=sprintf("%.3f",rt$"HR.95H")
Hazard.ratio=paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal=ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

pdf(file="forest.pdf", width = 6,height = 4.5)
n=nrow(rt)
nRow=n+1
ylim=c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))

xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="Black",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, '#E41A1C', '#377EB8')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()

par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
LOGindex=10 
hrLow = log(as.numeric(hrLow),LOGindex)
hrHigh = log(as.numeric(hrHigh),LOGindex)
hr = log(as.numeric(hr),LOGindex)
xlim = c(floor(min(hrLow,hrHigh)),ceiling(max(hrLow,hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="Black",lwd=2.5)
abline(v=log(1,LOGindex),col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > log(1,LOGindex), "#E41A1C", "#377EB8")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
a1 = axis(1,labels=F,tick=F)
axis(1,a1,LOGindex^a1)
dev.off()

###############

library(venn)                  
outFile="interGenes.txt"
geneList=list()

rt=read.table("TCGAGlndiff-77.txt",header=T,sep="\t",check.names=F)
geneNames=as.vector(rt[,1])                
geneNames=gsub("^ | $","",geneNames)     
uniqGene=unique(geneNames)                
geneList[["DEGs"]]=uniqGene

rt=read.table("tcga.uniCox.txt",header=T,sep="\t",check.names=F)
geneNames=as.vector(rt[,1])                
geneNames=gsub("^ | $","",geneNames)     
uniqGene=unique(geneNames)                 
geneList[["Prognostic genes"]]=uniqGene

mycol=c("LightCoral", "SteelBlue","#029149","#431A3D","#E0367A","#FFD121","#D8D155","#223D6C","#D20A13","#088247","#11AA4D","#7A142C","#5D90BA","#64495D","#7CC767")
pdf(file="venn.pdf", width=5, height=5)
venn(geneList,col=mycol[1:length(geneList)],zcolor=mycol[1:length(geneList)],box=F,ilabels=F)
dev.off()

interGenes=Reduce(intersect,geneList)
write.table(file=outFile,interGenes,sep="\t",quote=F,col.names=F,row.names=F)

#######################

library(corrplot)
library(circlize)

inputFile="unicoxgeneExp.txt"  

data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]
rt=t(data)

cor1=cor(rt)

col = c(rgb(1,0,0,seq(1,0,length=32)),rgb(0,1,0,seq(0,1,length=32)))
cor1[cor1==1]=0
c1 = ifelse(c(cor1)>=0,rgb(1,0,0,abs(cor1)),rgb(0,1,0,abs(cor1)))
col1 = matrix(c1,nc=ncol(rt))

pdf(file="circos.pdf", width=7, height=7)
par(mar=c(2,2,2,4))
circos.par(gap.degree=c(3,rep(2, nrow(cor1)-1)), start.degree = 180)
chordDiagram(cor1, grid.col=rainbow(ncol(rt)), col=col1, transparency = 0.5, symmetric = T)
par(xpd=T)
colorlegend(col, vertical = T,labels=c(1,0,-1),xlim=c(1.1,1.3),ylim=c(-0.4,0.4))
dev.off()
circos.clear()

pdf(file="corrplot.pdf", width=7, height=7)
corrplot(cor1,
         method = "pie",
         order = "hclust",
         type = "upper",
         col=colorRampPalette(c("Turquoise4", "GhostWhite", "Tomato3"))(20),
         tl.col="black")
dev.off()

pdf(file="corrplot.pdf", width=7, height=7)
corrplot(cor1,
         order="hclust",
         method = "circle",
         type = "upper",
         tl.cex=1, pch=T,
         insig = "label_sig",
         pch.cex = 1.6,
         sig.level=0.05,
         number.cex = 1,
         col=colorRampPalette(c("Turquoise4", "GhostWhite", "Tomato3"))(10),
         tl.col="black")
dev.off()

########################
inputFile="cnvMatrix.txt"     
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)   
GAIN=rowSums(rt> 0)       
LOSS=rowSums(rt< 0)       
GAIN=GAIN/ncol(rt)*100     
LOSS=LOSS/ncol(rt)*100      
data=cbind(GAIN, LOSS)
data=data[order(data[,"GAIN"],decreasing = T),]

data.max = apply(data, 1, max)
pdf(file="CNVfreq.pdf", width=9, height=6)
cex=1.3
par(cex.lab=cex, cex.axis=cex, font.axis=2, las=1, xpd=T)
bar=barplot(data.max, col="grey80", border=NA,
            xlab="", ylab="CNV.frequency(%)", space=1.5,
            xaxt="n", ylim=c(0,1.2*max(data.max)))
points(bar,data[,"GAIN"], pch=20, col="#E41A1C", cex=4)
points(bar,data[,"LOSS"], pch=20, col="#377EB8", cex=4)
legend("top", legend=c('GAIN','LOSS'), col=2:3, pch=20, bty="n", cex=4, ncol=2)
par(srt=45)
text(bar, par('usr')[3]-0.2, rownames(data), adj=1)
dev.off()

######################
library("RCircos")  

cytoBandIdeogram=read.table("refer.txt", header=T, sep="\t")
chr.exclude <- NULL
cyto.info <- cytoBandIdeogram
tracks.inside <- 5
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)

rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$text.size=2
rcircos.params$point.size=2
RCircos.Reset.Plot.Parameters(rcircos.params)

pdf(file="RCircos.pdf", width=8, height=8)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

RCircos.Scatter.Data=read.table("Rcircos.scatter.txt", header=T, sep="\t", check.names=F)
data.col <- 4
track.num <- 1
side <- "in"
RCircos.Scatter.Plot(RCircos.Scatter.Data, data.col, track.num, side, by.fold=0.1)

RCircos.Gene.Label.Data=read.table("Rcircos.geneLabel.txt", header=T, sep="\t", check.names=F)
name.col <- 4
side <- "in"
track.num <- 2
RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, track.num, side)
track.num <- 3
RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col, track.num, side)
dev.off()

#################
library(ConsensusClusterPlus)      
expFile="GlnGene-20.txt"          

data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(data)

maxK=9
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=100,
                             pItem=0.8,
                             pFeature=1,
                             title=workDir,
                             clusterAlg="pam",
                             distance="euclidean",
                             seed=123456,
                             plot="png")
calcICL(results, title="consensusScore", plot="png")


clusterNum=3      
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("Cluster")
cluster$Cluster=paste0("C", cluster$Cluster)
outTab=cbind(t(data), cluster)
outTab=rbind(ID=colnames(outTab), outTab)
write.table(outTab, file="Cluster3.txt", sep="\t", quote=F, col.names=F)

##############
library(survival)
library(survminer)

clusterFile="Cluster.txt"      
cliFile="time-OS.txt"          

cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365

sameSample=intersect(row.names(cluster), row.names(cli))
rt=cbind(cli[sameSample,,drop=F], cluster[sameSample,,drop=F])

length=length(levels(factor(rt$Cluster)))
diff=survdiff(Surv(futime, fustat) ~ Cluster, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ Cluster, data = rt)
#print(surv_median(fit))

bioCol=c("CadetBlue3","#FFC266","LightCoral","#7CC767","#6E568C","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="Cluster",
                   legend.labs=levels(factor(rt[,"Cluster"])),
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(years)",
                   ylab="Overall survival",
                   break.time.by = 1,
                   palette = bioCol,
                   #surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.3)
pdf(file="survival.pdf",onefile = FALSE,width=7,height=6)
print(surPlot)
dev.off()

#####################
library(limma)
library(pheatmap)

expFile="gsvaOut-2.txt"   
clusterFile="Cluster3.txt"        
cliFile="clinical.txt"          

data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(data)
data=t(data)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

samSample=intersect(row.names(cluster), row.names(cli))
data=data[samSample,,drop=F]
cli=cli[samSample,,drop=F]
cluster=cluster[samSample,,drop=F]
data=cbind(data, cluster, cli)
data=data[order(data$Cluster),,drop=F]       
Type=data[,((ncol(data)-ncol(cli)):ncol(data))]      
exp=data[,1:(ncol(data)-ncol(cli)-1)]               

sigVec=c("Cluster")
for(clinical in colnames(Type[,2:ncol(Type)])){
  data=Type[c("Cluster", clinical)]
  colnames(data)=c("Cluster", "clinical")
  data=data[(data[,"clinical"]!="unknow"),]
  tableStat=table(data)
  stat=chisq.test(tableStat)
  pvalue=stat$p.value
  Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
  sigVec=c(sigVec, paste0(clinical, Sig))
}
colnames(Type)=c(sigVec)

colorList=list()
#Type=Type[apply(Type,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
bioCol=c("CadetBlue3","LightCoral","#E6D8CF","#FF9900","#246b93", "#cc8e12", "#d561dd", 
         "#6ad157", "#f7aa5d", "#9ed84e", "#39ba30", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
         "#1a918f", "#ddd53e", "#ff66fc", "#2927c4", "#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
         "#4aef7b", "#e86502", "#99db27", "#e07233", "#8249aa","#cebb10", "#03827f", "#931635", "#ff523f",
         "#edd05e", "#6f25e8", "#0dbc21", "#167275", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
         "#dd27ce", "#07a301", "#391c82", "#2baeb5","#925bea", "#09f9f5",  "#63ff4f")
j=0
for(cli in colnames(Type[,1:ncol(Type)])){
  cliLength=length(levels(factor(Type[,cli])))
  cliCol=bioCol[(j+1):(j+cliLength)]
  j=j+cliLength
  names(cliCol)=levels(factor(Type[,cli]))
  if("unknow" %in% levels(factor(Type[,cli]))){
    cliCol["unknow"]="grey75"}
  colorList[[cli]]=cliCol
}

pdf("heatmap-1.pdf", width=9, height=6.5)
pheatmap(t(exp),
         annotation=Type,
         annotation_colors = colorList,
         color = colorRampPalette(c(rep("CadetBlue3",2), "white", rep("LightCoral",2)))(50),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=6,
         fontsize_row=5,
         fontsize_col=6)
dev.off()

#####################
library(limma)
library(GSEABase)
library(GSVA)
library(pheatmap)
expFile="TCGA.TPM.txt"              
clusterFile="Cluster3.txt"    
gmtFile="c5.go.bp.v2022.1.Hs.symbols.gmt"             

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())
gsvaResult=gsva(data, 
                geneSets, 
                min.sz=10, 
                max.sz=500, 
                verbose=TRUE,
                parallel.sz=1)
gsvaOut=rbind(id=colnames(gsvaResult), gsvaResult)
write.table(gsvaOut, file="gsvaOut.txt", sep="\t", quote=F, col.names=F)

cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

gsvaResult=t(gsvaResult)
sameSample=intersect(row.names(gsvaResult), row.names(cluster))
gsvaResult=gsvaResult[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
gsvaCluster=cbind(gsvaResult, cluster)
Project=gsub("(.*?)\\_.*", "\\1", rownames(gsvaCluster))
gsvaCluster=cbind(gsvaCluster, Project)

adj.P.Val.Filter=0.05
allType=as.vector(gsvaCluster$cluster)
comp=combn(levels(factor(allType)), 2)
for(i in 1:ncol(comp)){
  
  treat=gsvaCluster[gsvaCluster$cluster==comp[2,i],]
  con=gsvaCluster[gsvaCluster$cluster==comp[1,i],]
  data=rbind(con, treat)
  
  Type=as.vector(data$m6Acluster)
  ann=data[,c(ncol(data), (ncol(data)-1))]
  data=t(data[,-c((ncol(data)-1), ncol(data))])
  design=model.matrix(~0+factor(Type))
  colnames(design)=levels(factor(Type))
  fit=lmFit(data, design)
  contrast=paste0(comp[2,i], "-", comp[1,i])
  cont.matrix=makeContrasts(contrast, levels=design)
  fit2=contrasts.fit(fit, cont.matrix)
  fit2=eBayes(fit2)
  
  allDiff=topTable(fit2,adjust='fdr',number=200000)
  allDiffOut=rbind(id=colnames(allDiff),allDiff)
  write.table(allDiffOut, file=paste0(contrast, ".all.txt"), sep="\t", quote=F, col.names=F)
  
  diffSig=allDiff[with(allDiff, (abs(logFC)>0.1 & adj.P.Val < adj.P.Val.Filter )), ]
  diffSigOut=rbind(id=colnames(diffSig),diffSig)
  write.table(diffSigOut, file=paste0(contrast, ".diff.txt"), sep="\t", quote=F, col.names=F)
  
  bioCol=c("CadetBlue3","#FFC266","LightCoral","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  ann_colors=list()
  m6aCluCol=bioCol[1:length(levels(factor(allType)))]
  names(CluCol)=levels(factor(allType))
  ann_colors[["cluster"]]=CluCol[c(comp[1,i], comp[2,i])]
  
  termNum=20
  diffTermName=as.vector(rownames(diffSig))
  diffLength=length(diffTermName)
  if(diffLength<termNum){termNum=diffLength}
  hmGene=diffTermName[1:termNum]
  hmExp=data[hmGene,]
  pdf(file=paste0(contrast,".heatmap.pdf"),height=6,width=12)
  pheatmap(hmExp, 
           annotation=ann,
           annotation_colors = ann_colors,
           color = colorRampPalette(c(rep("CadetBlue3",5), "white", rep("LightCoral",5)))(50),
           cluster_cols =F,
           show_colnames = F,
           gaps_col=as.vector(cumsum(table(Type))),
           scale="row",
           fontsize = 15,
           fontsize_row=7,
           fontsize_col=10)
  dev.off()
}

#######################
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)
expFile="TCGA.TPM.txt"         
gmtFile="1-immune.gmt"       
clusterFile="Cluster3.txt"      

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

ssgseaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)

normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)

ssgseaOut=rbind(id=colnames(ssgseaScore), ssgseaScore)
write.table(ssgseaOut,file="1-ssGSEA.result.txt",sep="\t",quote=F,col.names=F)

cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

ssgseaScore=t(ssgseaScore)
sameSample=intersect(row.names(ssgseaScore), row.names(cluster))
ssgseaScore=ssgseaScore[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
scoreCluster=cbind(ssgseaScore, cluster)

data=melt(scoreCluster, id.vars=c("cluster"))
colnames(data)=c("cluster", "Immune", "Fraction")

bioCol=c("CadetBlue3","#FFC266","LightCoral","#7CC767","#6E568C","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"cluster"])))]
p=ggboxplot(data, x="Immune", y="Fraction",fill="cluster",
            ylab="Immune infiltration",
            xlab="",
            legend.title="cluster",
            palette=bioCol)
p=p+rotate_x_text(50)
pdf(file="1-3-boxplot.pdf", width=10, height=6.5)                     
p+stat_compare_means(aes(group=m6Acluster),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()

##############################
library(limma)
library(estimate)
inputFile="symbol.txt"      
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
out=data[rowMeans(data)>0,]
out=rbind(ID=colnames(out),out)

write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)

filterCommonGenes(input.f="uniq.symbol.txt", 
                  output.f="commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "commonGenes.gct",
              output.ds="estimateScore.gct")

scores=read.table("estimateScore.gct", skip=2, header=T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.", "\\-", rownames(scores))
out=rbind(ID=colnames(scores), scores)
write.table(out, file="TMEscores.txt", sep="\t", quote=F, col.names=F)


library(reshape2)
library(ggpubr)
clusterFile="cluster.txt"          #???ͽ????ļ?
estimateFile="TMEscores.txt"       #????΢?????????ļ?

Type=read.table(clusterFile, header=F, sep="\t", check.names=F, row.names=1)
colnames(Type)=c("Subtype")
Type=Type[order(Type[,"Subtype"],decreasing=T),,drop=F]
Type$Subtype=factor(Type$Subtype, levels=unique(Type$Subtype))

score=read.table(estimateFile, header=T, sep="\t", check.names=F, row.names=1)
score=score[,1:3]
score=score[row.names(Type),,drop=F]

rt=cbind(Type, score)

data=melt(rt, id.vars=c("Subtype"))
colnames(data)=c("Subtype", "scoreType", "Score")

p=ggviolin(data, x="scoreType", y="Score", fill = "Subtype",
           xlab="",
           ylab="TME score",
           legend.title="Subtype",
           add = "boxplot", add.params = list(color="white"),
           palette = c("SkyBlue3","LightCoral"), width=1)
p=p+rotate_x_text(45)
p1=p+stat_compare_means(aes(group=Subtype),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

pdf(file="vioplot.pdf", width=6, height=5)
print(p1)
dev.off()

#########################
library("limma")            
inputFile="TCGA.TPM.txt"      

rt=read.table(inputFile, header=T, sep="\t", check.names=F)       
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

v=voom(data, plot=F, save.plot=F)
out=v$E
out=rbind(ID=colnames(out),out)
write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)  

source("ssGSEA18.CIBERSORT.R")
results=CIBERSORT("ref.txt", "uniq.symbol.txt", perm=1000, QN=TRUE)

library(pheatmap)               
ssgseaFile="2-ssGSEA.result.txt"       
clusterFile="Cluster3.txt"        
estimateFile="TMEscores.txt"    

Type=read.table(clusterFile, header=F, sep="\t", check.names=F, row.names=1)
colnames(Type)=c("Subtype")
Type=Type[order(Type[,"Subtype"],decreasing=T),,drop=F]
Type$Subtype=factor(Type$Subtype, levels=unique(Type$Subtype))

rt=read.table(ssgseaFile, header=T, sep="\t", check.names=F, row.names=1)
rt=rt[,row.names(Type)]

score=read.table(estimateFile, header=T, sep="\t", check.names=F, row.names=1)
score=score[row.names(Type),,drop=F]

cluster=cbind(Type, score)

ann_colors=list()
clusterCol=c("CadetBlue3","#FFC266","LightCoral")
names(clusterCol)=levels(factor(Type$Subtype))
ann_colors[["Subtype"]]=clusterCol

pdf("2-estimateHM.pdf", width=9, height=6)
pheatmap(rt, annotation=cluster,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("CadetBlue3",5), "white", rep("LightCoral",5)))(50),
         cluster_cols =F,
         scale="row",
         show_colnames=F,
         fontsize=10,
         fontsize_row=10,
         fontsize_col=5)
dev.off()

#############################
library(limma) 
library(VennDiagram)
expFile="TCGA.TPM.txt"        
cluFile="Cluster3.txt"      
adj.P.Val.Filter=0.05     
s
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

cluster=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)

sameSample=intersect(colnames(data), row.names(cluster))
data=data[,sameSample]
cluster=cluster[sameSample,]

logFCfilter=1
geneList=list()
Type=as.vector(cluster)
design=model.matrix(~0+factor(Type))
colnames(design)=levels(factor(Type))
comp=combn(levels(factor(Type)), 2)
allDiffGenes=c()
for(i in 1:ncol(comp)){
  fit=lmFit(data, design)
  contrast=paste0(comp[2,i], "-", comp[1,i])
  #print(contrast)
  cont.matrix=makeContrasts(contrast, levels=design)
  fit2=contrasts.fit(fit, cont.matrix)
  fit2=eBayes(fit2)
  
  allDiff=topTable(fit2,adjust='fdr',number=200000)
  allDiffOut=rbind(id=colnames(allDiff),allDiff)
  write.table(allDiffOut, file=paste0(contrast, ".all.txt"), sep="\t", quote=F, col.names=F)
  
  diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adj.P.Val.Filter )), ]
  diffSigOut=rbind(id=colnames(diffSig),diffSig)
  write.table(diffSigOut, file=paste0(contrast, ".diff.txt"), sep="\t", quote=F, col.names=F)
  geneList[[contrast]]=row.names(diffSig)
}

venn.plot=venn.diagram(geneList,filename=NULL,
                       fill=c("CadetBlue3","#FFC266","LightCoral"),
                       lwd = 2,
                       col = c("black"))
pdf(file="venn.pdf", width=5, height=5)
grid.draw(venn.plot)
dev.off()

interGenes=Reduce(intersect,geneList)
write.table(file="interGene.txt",interGenes,sep="\t",quote=F,col.names=F,row.names=F)

interGeneExp=data[interGenes,]
interGeneExp=rbind(id=colnames(interGeneExp), interGeneExp)
write.table(interGeneExp, file="interGeneExp.txt", sep="\t", quote=F, col.names=F)

#########################
library(limma)
library(sva)
tcgaExpFile="TCGA-FPKM.txt"       
geoExpFile="GPL570-merge.txt"      
geneFile="diff-GlnUnicox.txt"             

rt=read.table(tcgaExpFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
tcga=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
tcga=avereps(tcga)
tcga=log2(tcga+1)

group=sapply(strsplit(colnames(tcga),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
tcga=tcga[,group==0]
tcga=t(tcga)
rownames(tcga)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tcga))
tcga=t(avereps(tcga))

rt=read.table(geoExpFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
geo=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
geo=avereps(geo)
qx=as.numeric(quantile(geo, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  geo[geo<0]=0
  geo=log2(geo+1)}
geo=normalizeBetweenArrays(geo)

sameGene=intersect(row.names(tcga),row.names(geo))
tcgaOut=tcga[sameGene,]
geoOut=geo[sameGene,]

all=cbind(tcgaOut,geoOut)
batchType=c(rep(1,ncol(tcgaOut)),rep(2,ncol(geoOut)))
outTab=ComBat(all, batchType, par.prior=TRUE)
tcgaOut=outTab[,colnames(tcgaOut)]
tcgaOut[tcgaOut<0]=0
geoOut=outTab[,colnames(geoOut)]
geoOut[geoOut<0]=0

tcgaTab=rbind(ID=colnames(tcgaOut), tcgaOut)
write.table(tcgaTab, file="TCGA.normalize.txt", sep="\t", quote=F, col.names=F)
geoTab=rbind(ID=colnames(geoOut), geoOut)
write.table(geoTab,file="GEO.normalize.txt",sep="\t",quote=F,col.names=F)

gene=read.table(geneFile, header=T, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(tcgaOut))
tcgaShareExp=tcgaOut[sameGene,]
geoShareExp=geoOut[sameGene,]

tcgaShareExp=rbind(ID=colnames(tcgaShareExp),tcgaShareExp)
write.table(tcgaShareExp,file="TCGA.share.txt",sep="\t",quote=F,col.names=F)
geoShareExp=rbind(ID=colnames(geoShareExp),geoShareExp)
write.table(geoShareExp,file="GEO.share.txt",sep="\t",quote=F,col.names=F)

library(limma)              
expFile="TCGA.share.txt"     
cliFile="TCGA-time.txt"          

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data=t(data)

cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)   

sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
out=cbind(cli,data)
out=cbind(id=row.names(out),out)
write.table(out,file="tcga.expTime.txt",sep="\t",row.names=F,quote=F)

library(limma)               
expFile="GEO.share.txt"     
cliFile="GPL570-time.txt"        

rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data=t(data)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)    

sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
out=cbind(cli,data)
out=cbind(id=row.names(out),out)
write.table(out,file="GEO.expTime.txt",sep="\t",row.names=F,quote=F)

####################################
library(survival)                
coxPfilter=0.01                  
inputFile="tcga.expTime.txt"      
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)   
rt$futime=rt$futime/365

outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<coxPfilter){
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}

write.table(outTab,file="tcga.uniCox.txt",sep="\t",row.names=F,quote=F)

uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="tcga.uniSigExp.txt",sep="\t",row.names=F,quote=F)

bioForest=function(coxFile=null,forestFile=null,forestCol=null){
  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  height=nrow(rt)/12.5+5
  pdf(file=forestFile, width = 7,height = height)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))

  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  LOGindex=10 
  hrLow = log(as.numeric(hrLow),LOGindex)
  hrHigh = log(as.numeric(hrHigh),LOGindex)
  hr = log(as.numeric(hr),LOGindex)
  xlim = c(floor(min(hrLow,hrHigh)),ceiling(max(hrLow,hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="grey31",lwd=2.5)
  abline(v=log(1,LOGindex),col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > log(1,LOGindex), "SteelBlue", "Red")
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
  a1 = axis(1,labels=F,tick=F)
  axis(1,a1,LOGindex^a1)
  dev.off()
}
bioForest(coxFile="tcga.uniCox.txt", forestFile="forest.pdf")

################################
library("glmnet")
library("survival")

set.seed(12345)
coxSigFile="tcga.uniSigExp-0.01.txt"      
geoFile="GEO.expTime.txt"            

rt=read.table(coxSigFile, header=T, sep="\t", check.names=F, row.names=1)
geo=read.table(geoFile, header=T, sep="\t", check.names=F, row.names=1)
sameGene=intersect(colnames(rt)[3:ncol(rt)], colnames(geo)[3:ncol(geo)])
rt=rt[,c("futime","fustat",sameGene)]
rt$futime[rt$futime<=0]=0.003

x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime, rt$fustat))
fit=glmnet(x, y, family="cox", maxit=1000)

pdf("lasso.lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

cvfit=cv.glmnet(x, y, family="cox", maxit=1000)
pdf("lasso.cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

coef=coef(fit, s=cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]lassoSigExp=uniSigExp[,c("futime", "fustat", lassoGene)]
lassoSigExpOut=cbind(id=row.names(lassoSigExp), lassoSigExp)
geneCoef=cbind(Gene=lassoGene, Coef=actCoef)
write.table(lassoSigExpOut,file="lasso.SigExp.txt",sep="\t",row.names=F,quote=F)
write.table(geneCoef, file="lasso.geneCoef.txt", sep="\t", quote=F, row.names=F)

trainFinalGeneExp=rt[,lassoGene]
#myFun=function(x){crossprod(as.numeric(x),actCoef)}
#trainScore=apply(trainFinalGeneExp,1,myFun)
trainScore=predict(cvfit, newx=x, s="lambda.min", type="response")
Risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
outCol=c("futime","fustat",lassoGene)
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore),Risk)
write.table(cbind(id=rownames(outTab),outTab),file="risk.TCGA.txt",sep="\t",quote=F,row.names=F)

geo=read.table(geoFile, header=T, sep="\t", check.names=F, row.names=1)
geo$futime=geo$futime/365
geo=geo[,colnames(rt)]
geo[,3:ncol(geo)]=geo[,3:ncol(geo)]*median(x)/median(as.matrix(geo[,3:ncol(geo)]))
#testFinalGeneExp=geo[,lassoGene]
#testScore=apply(testFinalGeneExp,1,myFun)
testScore=predict(cvfit, newx=as.matrix(geo[,c(3:ncol(geo))]), s="lambda.min", type="response")
outCol=c("futime","fustat",lassoGene)
Risk=as.vector(ifelse(testScore>median(trainScore),"high","low"))
outTab=cbind(geo[,outCol],riskScore=as.vector(testScore),Risk)
write.table(cbind(id=rownames(outTab),outTab),file="risk.GEO.txt",sep="\t",quote=F,row.names=F)

library(survival)
library(survminer)

bioSurvival=function(inputFile=null, outFile=null){
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)ֵ
  diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
  pValue=1-pchisq(diff$chisq, df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=T,
                     pval=pValue,
                     pval.size=6,
                     legend.title="Risk",
                     legend.labs=c("High risk", "Low risk"),
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette=c("LightCoral", "CadetBlue3"),
                     risk.table=TRUE,
                     risk.table.title="",
                     risk.table.height=.25)
  pdf(file=outFile, width=6.5, height=5.5, onefile=FALSE)
  print(surPlot)
  dev.off()
}

bioSurvival(inputFile="risk.TCGA.txt", outFile="survival.TCGA.pdf")
bioSurvival(inputFile="risk.GEO.txt", outFile="survival.GEO.pdf")

#####################
library(reshape2)
library(ggpubr)
library(ggExtra)
library(pheatmap)

bioRiskPlot=function(inputFile=null, riskScoreFile=null, survStatFile=null, survCorFile=null){
  rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    #??ȡ?????ļ?
  rt$riskScore[rt$riskScore>quantile(rt$riskScore,0.99)]=quantile(rt$riskScore,0.99)
  rt$risk=factor(rt$risk, levels=c("low", "high"))
  rt=rt[order(rt$riskScore),]     
  riskClass=rt[,"risk"]
  lowLength=length(riskClass[riskClass=="low"])
  highLength=length(riskClass[riskClass=="high"])
  lowMax=max(rt$riskScore[riskClass=="low"])
  line=rt[,"riskScore"]
  pdf(file=riskScoreFile, width=6, height=5)
  plot(line, type="p", pch=20,
       xlab="Patients (increasing risk socre)", ylab="Risk score",
       col=c(rep("CadetBlue3",lowLength),rep("LightCoral",highLength)) )
  abline(h=lowMax,v=lowLength,lty=2)
  legend("topleft", c("High risk", "Low Risk"),bty="n",pch=19,col=c("LightCoral","CadetBlue3"),cex=1.2)
  dev.off()
  
  color=as.vector(rt$fustat)
  color[color==1]="LightCoral"
  color[color==0]="CadetBlue3"
  pdf(file=survStatFile, width=6, height=5)
  plot(rt$futime, pch=19,
       xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
       col=color)
  legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("LightCoral","CadetBlue3"),cex=1.2)
  abline(v=lowLength,lty=2)
  dev.off()
  
  x=as.numeric(rt[,"riskScore"])
  y=as.numeric(rt[,"futime"])
  df1=as.data.frame(cbind(x,y))
  p1=ggplot(df1, aes(x, y)) + 
    xlab("Risk score") + 
    ylab("OS (years)") +
    geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
    stat_cor(method = 'spearman', aes(x =x, y =y))
  p2=ggMarginal(p1, type="density", xparams=list(fill = "LightCoral"), yparams=list(fill = "CadetBlue3"))

  pdf(file=survCorFile, width=5.2, height=5)
  print(p2)
  dev.off()
}

bioRiskPlot(inputFile="trainRisk.txt",
            riskScoreFile="train.riskScore.pdf",
            survStatFile="train.survStat.pdf",
            survCorFile="train.survCor.pdf")

bioRiskPlot(inputFile="testRisk.txt",
            riskScoreFile="test.riskScore.pdf",
            survStatFile="test.survStat.pdf",
            survCorFile="test.survCor.pdf")

###############
library(Rtsne)
library(ggplot2)

bioPCA=function(inputFile=null, pcaFile=null, tsneFile=null){
  rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
  data=rt[c(3:(ncol(rt)-2))]
  risk=rt[,"risk"]
  
  data.pca=prcomp(data, scale. = TRUE)
  pcaPredict=predict(data.pca)
  PCA = data.frame(PC1 = pcaPredict[,1], PC2 = pcaPredict[,2],risk=risk)	
  pdf(file=pcaFile, height=4.5, width=5.5)       #???????????ļ?
  p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = risk)) +
    scale_colour_manual(name="Risk",  values =c("LightCoral", "#50b7c1"))+
    theme_bw()+
    theme(plot.margin=unit(rep(1.5,4),'lines'))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(p)
  dev.off()

  tsneOut=Rtsne(data, dims=2, perplexity=10, verbose=F, max_iter=500,check_duplicates=F)
  tsne=data.frame(tSNE1 = tsneOut$Y[,1], tSNE2 = tsneOut$Y[,2],risk=risk)	

  pdf(file=tsneFile, height=4.5, width=5.5)       #???????????ļ?
  p=ggplot(data = tsne, aes(tSNE1, tSNE2)) + geom_point(aes(color = risk)) +
    scale_colour_manual(name="Risk",  values =c("LightCoral", "#50b7c1"))+
    theme_bw()+
    theme(plot.margin=unit(rep(1.5,4),'lines'))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(p)
  dev.off()
}
bioPCA(inputFile="trainRisk.txt", pcaFile="train.PCA.pdf", tsneFile="train.t-SNE.pdf")
bioPCA(inputFile="testRisk.txt", pcaFile="test.PCA.pdf", tsneFile="test.t-SNE.pdf")

###########################

library(survival)
library(survminer)
library(timeROC)
library(rms)

bioROC=function(inputFile=null, rocFile=null, calFile=null){
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
                 marker=rt$riskScore, cause=1,
                 weighting='aalen',
                 times=c(1,3,5), ROC=TRUE)
  pdf(file=rocFile,width=5,height=5)
  plot(ROC_rt,time=1,col='#FFC266',title=FALSE,lwd=2)
  plot(ROC_rt,time=3,col='CadetBlue3',add=TRUE,title=FALSE,lwd=2)
  plot(ROC_rt,time=5,col='LightCoral',add=TRUE,title=FALSE,lwd=2)
  legend('bottomright',
         c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
           paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
           paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
         col=c("#FFC266","CadetBlue3","LightCoral"),lwd=2,bty = 'n')
  dev.off()
  
  pdf(file=calFile, width=5, height=5)
  f <- cph(Surv(futime, fustat) ~ riskScore, x=T, y=T, surv=T, data=rt, time.inc=1)
  cal <- calibrate(f, cmethod="KM", method="boot", u=1, m=(nrow(rt)/3), B=1000)
  plot(cal, xlim=c(0,1), ylim=c(0,1),
       xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=1.5, col="#FFC266", sub=F)
  f <- cph(Surv(futime, fustat) ~ riskScore, x=T, y=T, surv=T, data=rt, time.inc=3)
  cal <- calibrate(f, cmethod="KM", method="boot", u=3, m=(nrow(rt)/3), B=1000)
  plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=1.5, col="CadetBlue3", sub=F, add=T)
  f <- cph(Surv(futime, fustat) ~ riskScore, x=T, y=T, surv=T, data=rt, time.inc=5)
  cal <- calibrate(f, cmethod="KM", method="boot", u=5, m=(nrow(rt)/3), B=1000)
  plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=1.5, col="LightCoral", sub=F, add=T)
  legend('bottomright', c('1-year', '3-year', '5-year'),
         col=c("#FFC266","CadetBlue3","LightCoral"), lwd=1.5, bty = 'n')
  dev.off()
}

bioROC(inputFile="trainRisk.txt", rocFile="train.ROC.pdf", calFile="train.cal.pdf")
bioROC(inputFile="testRisk.txt", rocFile="test.ROC.pdf", calFile="test.cal.pdf")

###############
library(plyr)
library(ggplot2)
library(ggpubr)
scoreFile="TCGAscore.group.txt"   
cliFile="TCGA-clinical.txt"          
trait="fustat"                    

score=read.table("TCGAscore.group.txt", header=T, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(score), row.names(cli))
rt=cbind(score[sameSample,], cli[sameSample,])
bioCol=c("CadetBlue3","LightCoral","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(unique(rt[,trait]))]
rt1=rt[,c(trait, "group")]
colnames(rt1)=c("trait", "group")
df=as.data.frame(table(rt1))
df=ddply(df, .(group), transform, percent = Freq/sum(Freq) * 100)
df=ddply(df, .(group), transform, pos = (cumsum(Freq) - 0.5 * Freq))
df$label=paste0(sprintf("%.0f", df$percent), "%")
p=ggplot(df, aes(x = factor(group), y = percent, fill = trait)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values=bioCol)+
  xlab("ICI score")+ ylab("Percent weight")+  guides(fill=guide_legend(title=trait))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  #coord_flip()+
  theme_bw()
pdf(file="barplot.pdf", width=4, height=5)
print(p)
dev.off()

rt2=rt[,c(trait, "score")]
colnames(rt2)=c("trait", "score")
type=levels(factor(rt2[,"trait"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

boxplot=ggboxplot(rt2, x="trait", y="score", fill="trait",
                  xlab=trait,
                  ylab="ICI score",
                  legend.title=trait,
                  palette=bioCol
)+ 
  stat_compare_means(comparisons=my_comparisons)
pdf(file="boxplot.pdf",width=4,height=4.5)
print(boxplot)
dev.off()

########################
library(survival)
library(survminer)
library(timeROC)
riskFile="risk.TCGA.txt"      
cliFile="TCGA-clinical.txt"        

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "riskScore")]
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)
bioCol=rainbow(ncol(rt)-1, s=0.9, v=0.9)

ROC_rt=timeROC(T=risk$futime, delta=risk$fustat,
               marker=risk$riskScore, cause=1,
               weighting='aalen',
               times=c(1,3,5), ROC=TRUE)
pdf(file="ROC.pdf", width=5.5, height=5.5)
plot(ROC_rt,time=1,col='#FFC266',title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='CadetBlue3',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col='LightCoral',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=c("#FFC266","CadetBlue3","LightCoral"), lwd=2, bty = 'n')
dev.off()

predictTime=1   
aucText=c()
pdf(file="cliROC.pdf", width=5.5, height=5.5)

i=3
ROC_rt=timeROC(T=risk$futime,
               delta=risk$fustat,
               marker=risk$riskScore, cause=1,
               weighting='aalen',
               times=c(predictTime),ROC=TRUE)
plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=2)
aucText=c(paste0("Risk", ", AUC=", sprintf("%.3f",ROC_rt$AUC[2])))
abline(0,1)

for(i in 4:ncol(rt)){
  ROC_rt=timeROC(T=rt$futime,
                 delta=rt$fustat,
                 marker=rt[,i], cause=1,
                 weighting='aalen',
                 times=c(predictTime),ROC=TRUE)
  plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=2, add=TRUE)
  aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%.3f",ROC_rt$AUC[2])))
}

legend("bottomright", aucText,lwd=2,bty="n",col=bioCol[1:(ncol(rt)-1)])
dev.off()

######################
library(survival)       

bioForest=function(coxFile=null, forestFile=null, forestCol=null){
 
  rt <- read.table(coxFile, header=T, sep="\t", check.names=F, row.names=1)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  pdf(file=forestFile, width=6.5, height=4.5)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="black",lwd=3)
  abline(v=1, col="black", lty=2, lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=2)
  axis(1)
  dev.off()
}

indep=function(riskFile=null,cliFile=null,uniOutFile=null,multiOutFile=null,uniForest=null,multiForest=null){
  risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)   
  cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)    
  sameSample=intersect(row.names(cli),row.names(risk))
  risk=risk[sameSample,]
  cli=cli[sameSample,]
  rt=cbind(futime=risk[,1], fustat=risk[,2], cli, riskScore=risk[,(ncol(risk)-1)])
  uniTab=data.frame()
  for(i in colnames(rt[,3:ncol(rt)])){
    cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
    coxSummary = summary(cox)
    uniTab=rbind(uniTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
  write.table(uniTab,file=uniOutFile,sep="\t",row.names=F,quote=F)
  bioForest(coxFile=uniOutFile, forestFile=uniForest, forestCol="green")
  LightCoral???ض?########
  ab[as.numeric(uniTab[,"pvalue"])<1,]
  rt1=rt[,c("futime", "fustat", as.vector(uniTab[,"id"]))]
  multiCox=coxph(Surv(futime, fustat) ~ ., data = rt1)
  multiCoxSum=summary(multiCox)
  multiTab=data.frame()
  multiTab=cbind(
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  multiTab=cbind(id=row.names(multiTab),multiTab)
  write.table(multiTab,file=multiOutFile,sep="\t",row.names=F,quote=F)
  bioForest(coxFile=multiOutFile, forestFile=multiForest, forestCol="red")
}

indep(riskFile="risk.TCGA.txt",
      cliFile="clinical.txt",
      uniOutFile="uniCox.txt",
      multiOutFile="multiCox.txt",
      uniForest="uniForest.pdf",
      multiForest="multiForest.pdf")

######################
library(survival)
library(survminer)
inputFile="risk.TCGA.txt"      
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
rt=rt[,1:(ncol(rt)-2)]

outTab=data.frame()
for(gene in colnames(rt)[3:ncol(rt)]){
  if(sd(rt[,gene])<0.001){next}
  data=rt[,c("futime", "fustat", gene)]
  colnames(data)=c("futime", "fustat", "gene")
  res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("gene"))
  res.cat=surv_categorize(res.cut)
  fit=survfit(Surv(futime, fustat) ~gene, data = res.cat)
  diff=survdiff(Surv(futime, fustat) ~gene,data =res.cat)
  pValue=1-pchisq(diff$chisq, df=1)
  outVector=cbind(gene, res.cut$cutpoint[1], pValue)
  outTab=rbind(outTab,outVector)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  
  surPlot=ggsurvplot(fit,
                     data=res.cat,
                     pval=pValue,
                     pval.size=6,
                     legend.title=gene,
                     legend.labs=c("high","low"),
                     xlab="Time(years)",
                     break.time.by=1,
                     palette=c("LightCoral", "CadetBlue3"),
                     conf.int=T,
                     risk.table=F,
                     risk.table.title="",
                     risk.table.height=.25)
  pdf(file=paste0("Survival.",gene,".pdf"), width=5.5, height=5.1, onefile=FALSE)
  print(surPlot)
  dev.off()
}
write.table(outTab,file="survival.result.txt",sep="\t",row.names=F,quote=F)

##################################
#???ð?
library(survival)
library(regplot)
library(rms)

#######################
riskFile="risk.TCGA.txt"     
cliFile="clinical-1.txt"       
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
cli$Age=as.numeric(cli$Age)
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1[,c("futime", "fustat", "risk")], cli)
res.cox=coxph(Surv(futime, fustat) ~ . , data = rt)
nom1=regplot(res.cox,
             plots = c("density", "boxes"),
             clickable=F,
             title="",
             points=TRUE,
             droplines=TRUE,
             observation=rt[2,],
             rank="sd",
             failtime = c(1,3,5),
             prfail = F)
nomoRisk=predict(res.cox, data=rt, type="risk")
rt=cbind(risk1, Nomogram=nomoRisk)
outTab=rbind(ID=colnames(rt), rt)
write.table(outTab, file="nomoRisk.txt", sep="\t", col.names=F, quote=F)

ph<-cox.zph(res.cox)

ph$table
install.packages("survminer")
library(survminer)
par(cex=0.5)
pdf(file="pha_test.pdf", width=12, height=8)
ggcoxzph(ph)
dev.off()

pdf(file="calibration.pdf", width=5, height=5)
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=1)
cal <- calibrate(f, cmethod="KM", method="boot", u=1, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=1.5, col="#FFC266", sub=F)
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=3)
cal <- calibrate(f, cmethod="KM", method="boot", u=3, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=1.5, col="LightCoral", sub=F, add=T)
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=5)
cal <- calibrate(f, cmethod="KM", method="boot", u=5, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=1.5, col="LightCoral", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
       col=c("#FFC266","CadetBlue3","LightCoral"), lwd=1.5, bty = 'n')
dev.off()

################################
library(limma)
library(survival)
library(survminer)
library(timeROC)

expFile="TCGA.TPM.txt"            
riskFile="risk.TCGA.txt"     #
geneFiles=c("Chen et al-36330389.txt", "Zhang et al-36424614.txt","Sabatier et al-21654678.txt","Yue et al-31888563.txt","Xiang et al-36185201.txt")     #ģ?ͻ????ļ?

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2", "1", group)
data=data[,group==0]

riskRT=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
riskRT=riskRT[,c("futime","fustat","riskScore")]
colnames(riskRT)=c("futime","fustat","Gln signature")

for(i in geneFiles){
  #??ȡ?????б?
  header=unlist(strsplit(i, "\\."))
  gene=read.table(i, header=F, sep="\t", check.names=F)
  sameGene=intersect(as.vector(gene[,1]), row.names(data))
  data1=data[sameGene,]
  colnames(data1)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(data1))
  data1=t(data1)
  data1=avereps(data1)
  
  cli=riskRT[,c("futime", "fustat")]
  sameSample=intersect(row.names(data1), row.names(cli))
  data1=data1[sameSample,]
  cli=cli[sameSample,]
  data1=cbind(cli,data1)
  
  multiCox=coxph(Surv(futime, fustat) ~ ., data = data1)
  riskScore=predict(multiCox,type="risk", newdata=data1)
  data1=cbind(data1, riskScore)
  data1=data1[row.names(riskRT),]
  riskRT=cbind(riskRT, data1[,"riskScore"])
  colnames(riskRT)[ncol(riskRT)]=header[[1]]
}

riskOut=rbind(ID=colnames(riskRT), riskRT)
write.table(riskOut, file="risk.models.txt", sep="\t", col.names=F, quote=F)

bioSurvival=function(inputFile=null, outFile=null, varName=null){
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  rt$Type=ifelse(rt[,varName]>median(rt[,varName]), "high", "low")
  diff=survdiff(Surv(futime, fustat) ~ Type,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ Type, data = rt)
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     title=varName,
                     legend.title="Risk",
                     legend.labs=c("High risk", "Low risk"),
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette=c("LightCoral", "CadetBlue3"),
                     risk.table=F,
                     risk.table.title="",
                     risk.table.col = "strata",
                     risk.table.height=.25)
  pdf(file=outFile, onefile = FALSE, width=6, height=5)
  print(surPlot)
  dev.off()
}
bioROC=function(inputFile=null, outFile=null, varName=null){
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
                 marker=rt[,varName], cause=1,
                 weighting='aalen',
                 times=c(1,3,5), ROC=TRUE)
  pdf(file=outFile, width=5, height=5)
  plot(ROC_rt,time=1,col='#FFC266',title=FALSE,lwd=2)
  plot(ROC_rt,time=3,col='CadetBlue3',add=TRUE,title=FALSE,lwd=2)
  plot(ROC_rt,time=5,col='LightCoral',add=TRUE,title=FALSE,lwd=2)
  text(0.75, 0.24, varName, cex=1.2)
  legend('bottomright',
         c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
           paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
           paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
         col=c("#FFC266",'CadetBlue3','LightCoral'),lwd=2,bty = 'n')
  dev.off()
}

for(varName in colnames(riskRT)[3:ncol(riskRT)]){
  bioSurvival(inputFile="risk.models.txt", outFile=paste0("sur.",varName,".pdf"), varName=varName)
  bioROC(inputFile="risk.models.txt", outFile=paste0("ROC.",varName,".pdf"), varName=varName)
}

###############################

library(survival)
library(survcomp)
library(ggplot2)
library(ggpubr)

inputFile="risk.models.txt"  
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
df=data.frame()
for(i in colnames(rt)[3:ncol(rt)]){
  cindex=concordance.index(x=rt[,i], surv.time=rt$futime, surv.event=rt$fustat,method="noether")
  df=rbind(df, cbind(i,sprintf("%.03f",cindex$c.index)))
}
colnames(df)=c("signature", "cindex")
df[,"cindex"]=as.numeric(df[,"cindex"])

color=rainbow(nrow(df),alpha=0.75)
p=ggbarplot(df, x="signature", y="cindex", fill="signature",
            xlab="", ylab="C-index", add = "none",
            palette= c( "LightCoral", "Pink", "MistyRose","SkyBlue1", "LightBlue","PaleTurquoise1"),
            label=T, legend="")
p=p+rotate_x_text(50)
p=p+ylim(0,round(max(df[,"cindex"])+0.15,1))
pdf(file="C-index.pdf", width=6, height=5)
print(p)
dev.off()

outdata=list()
legendsname=c()

for(i in 3:ncol(rt)){
  OS=Surv(rt$futime, rt$fustat)
  marker=rt[,i]
  marker.pp<-seq(from=0, to=1, length=100)
  marker.qq<-quantile(marker,marker.pp)
  fitdat.df<-data.frame(marker=marker)
  newdat.df<-data.frame(marker=marker.qq)
  cox.model<-coxph(OS~marker, data=fitdat.df)
  rms.calc <-summary(survfit(cox.model, newdata=newdat.df))
  rms.mean <-rms.calc$table[,"rmean"]
  name=colnames(rt)[i]
  HR=sprintf("%.03f", summary(cox.model)$conf.int[,"exp(coef)"])
  HR.95L=sprintf("%.03f", summary(cox.model)$conf.int[,"lower .95"])
  HR.95H=sprintf("%.03f", summary(cox.model)$conf.int[,"upper .95"])
  pvalue=summary(cox.model)$coefficients[,"Pr(>|z|)"]
  p=ifelse(pvalue<0.001,"p<0.001",paste0("p=",sprintf("%.03f",pvalue)))
  legendsname=c(legendsname,paste0(name,", HR:",HR,"(",HR.95L,"-",HR.95H,"), ",p))     #ͼ??
  outdata[[name]]= data.frame(marker.pp,rms.mean)
}

alldata=do.call("rbind",outdata)
xlim2=max(alldata$rms.mean)
pdf(file="RMS.pdf", width=6, height=6)
par(las=1)
plot(1,xlim=c(0,1),ylim=c(0,xlim2),type="n",xlab="Percentile of scores",ylab="RMS (years)")
names=names(outdata)
for(i in 1:length(outdata)){
  namei=names[i]
  outdatai=outdata[[namei]]
  points(outdatai$marker.pp,outdatai$rms.mean,col=color[i],pch=20,cex=0.8)
}
legend("bottomleft",legend=legendsname,col=color,pch=20,bty="n",cex=1)
dev.off()

######################
library(limma)
library(GSVA)
library(GSEABase)
library(ggpubr)
library(reshape2)

expFile="TCGA.TPM.txt"      
gmtFile="Lipid.gmt"          
riskFile="risk.TCGA.txt"      

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
mat=mat[rowMeans(mat)>0,]

geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())

ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)

normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}

data=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(data), data)
write.table(ssgseaOut, file="Carbohydrate-Score.txt", sep="\t", quote=F, col.names=F)

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=t(data[,group==0])
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

sameSample=intersect(row.names(data),row.names(risk))
data=data[sameSample,,drop=F]
risk=risk[sameSample,"Risk",drop=F]
rt1=cbind(data, risk)
ͼ
data=melt(rt1, id.vars=c("Risk"))
colnames(data)=c("Risk","Type","Score")
data$Risk=factor(data$Risk, levels=c("low","high"))
p=ggboxplot(data, x="Type", y="Score", fill = "Risk",
            xlab="",ylab="Score",add = "none",palette = c("CadetBlue3","LightCoral") )
p=p+rotate_x_text(50)
p=p+stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")

pdf(file="Lipid.pdf", width=8, height=8)
print(p)
dev.off()

##############################################
library(limma)
library(ggpubr)
tideFile="TIDE-OV-2.txt"          #TIDE?ļ?
riskFile="risk.TCGA.txt"     #?????ļ?

tide=read.table(tideFile, header=T, sep="\t", check.names=F, row.names=1)
group=sapply(strsplit(row.names(tide),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
tide=tide[group==0,]
row.names(tide)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(tide))
tide=avereps(tide)

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

sameSample=intersect(row.names(tide), row.names(risk))
tide=tide[sameSample, , drop=F]
risk=risk[sameSample, "Risk", drop=F]
data=cbind(tide, risk)

data$Risk=ifelse(data$Risk=="high", "High-risk", "Low-risk")
group=levels(factor(data$Risk))
data$Risk=factor(data$Risk, levels=c("Low-risk", "High-risk"))
group=levels(factor(data$Risk))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

for(i in colnames(data)[1:(ncol(data)-1)]){
  gg1=ggviolin(data, x="Risk", y=i, fill = "Risk", 
               xlab="", ylab=i,
               palette=c("CadetBlue3", "LightCoral"),
               legend.title="Risk",
               add = "boxplot", add.params = list(fill="white"))+ 
    stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
  
  pdf(file=paste0("violin.", i, ".pdf"), width=6, height=5)
  print(gg1)
  dev.off()
}

###############################

library(ggpubr)              #???ð?
tciaFile="TCIA.txt"          #TCIA?????ļ?
riskFile="risk.TCGA.txt"     #?????ļ?
ips=read.table(tciaFile, header=T, sep="\t", check.names=F, row.names=1)
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

sameSample=intersect(row.names(ips), row.names(risk))
ips=ips[sameSample, , drop=F]
risk=risk[sameSample, "risk", drop=F]
data=cbind(ips, risk)

data$risk=ifelse(data$risk=="high", "High-risk", "Low-risk")
group=levels(factor(data$risk))
data$risk=factor(data$risk, levels=c("Low-risk", "High-risk"))
group=levels(factor(data$risk))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

for(i in colnames(data)[1:(ncol(data)-1)]){
  rt=data[,c(i, "risk")]
  colnames(rt)=c("IPS", "Risk")
  gg1=ggviolin(rt, x="Risk", y="IPS", fill = "Risk", 
               xlab="", ylab=i,
               legend.title="Risk",
               palette=c("CadetBlue3", "LightCoral"),
               add = "boxplot", add.params = list(fill="white"))+ 
    stat_compare_means(comparisons = my_comparisons)
  #stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
  
  pdf(file=paste0(i, ".pdf"), width=6, height=5)
  print(gg1)
  dev.off()
}

########################
#???ð?
library(limma)
library(reshape2)
library(ggplot2)

expFile="TCGA.TPM.txt"   
geneFile="gene.txt"             
riskFile="risk.TCGA.txt"       

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(geneRT[,1]), rownames(data))
data=t(data[sameGene,])

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(risk))
data=data[sameSample,,drop=F]
risk=risk[sameSample,3:(ncol(risk)-1),drop=F]

outTab=data.frame()
for(checkpiont in colnames(data)){
  for(gene in colnames(risk)){
    x=as.numeric(data[,checkpiont])
    y=as.numeric(risk[,gene])
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pvalue=corT$p.value
    text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
    outTab=rbind(outTab,cbind(Gene=gene, checkpiont=checkpiont, cor, text, pvalue))
  }
}

outTab$Gene=factor(outTab$Gene, levels=colnames(risk))
outTab$cor=as.numeric(outTab$cor)
pdf(file="checkpointCor.pdf", width=8, height=7)
ggplot(outTab, aes(Gene, checkpiont)) + 
  geom_tile(aes(fill = cor), colour = "grey", size = 1)+
  scale_fill_gradient2(low = "CadetBlue3", mid = "white", high = "LightCoral")
geom_text(aes(label=text),col ="black",size = 3) +
  theme_minimal() +   
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),   
        axis.text.y = element_text(size = 10, face = "bold")) +      
  labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) +  
  scale_x_discrete(position = "bottom")     
dev.off()

###############################

library(survival)
library(caret)
library(glmnet)
library(survminer)
library(timeROC)

rt=read.table("lasso.SigExp.txt", header=T, sep="\t", check.names=F, row.names=1) 

multiCox <- coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)
score=predict(multiCox, type="risk", newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`", "", coxGene)

vigor=read.table("exp.txt", header=T, sep="\t", check.names=F, row.names=1)
vigor=t(vigor[coxGene,])

cli=read.table("time.txt", header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("futime", "fustat")

sameSample=intersect(row.names(vigor), row.names(cli))
vigorTime=cbind(cli[sameSample,,drop=F], vigor[sameSample,,drop=F])
#vigorTime[,3:ncol(vigorTime)]=vigorTime[,3:ncol(vigorTime)]*median(as.matrix(rt[,coxGene]))/median(as.matrix(vigorTime[,3:ncol(vigorTime)]))

vigorScore=predict(multiCox, type="risk", newdata=vigorTime)
Risk=as.vector(ifelse(vigorScore>median(vigorScore), "high", "low"))
vigorRiskOut=cbind(vigorTime, riskScore=as.vector(vigorScore), Risk)
vigorRiskOut=cbind(id=rownames(vigorRiskOut), vigorRiskOut)
write.table(vigorRiskOut,file="risk.IMvigor.txt",sep="\t",quote=F,row.names=F)

bioSurvival=function(inputFile=null, outFile=null){
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  diff=survdiff(Surv(futime, fustat) ~Risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt)
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     legend.title="Risk",
                     legend.labs=c("High risk", "Low risk"),
                     xlab="Time(years)",
                     break.time.by = 2,
                     palette=c("LightCoral", "CadetBlue3"),
                     risk.table=F,
                     risk.table.title="",
                     risk.table.col = "strata",
                     risk.table.height=.25)
  pdf(file=outFile,onefile = FALSE,width = 6,height =5)
  print(surPlot)
  dev.off()
}

bioROC=function(inputFile=null, outFile=null){
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,
                 marker=rt$riskScore,cause=1,
                 weighting='aalen',
                 times=c(1,3,5),ROC=TRUE)
  pdf(file=outFile, width=5, height=5)
  plot(ROC_rt,time=1,col='Orange',title=FALSE,lwd=2)
  plot(ROC_rt,time=3,col='CadetBlue3',add=TRUE,title=FALSE,lwd=2)
  plot(ROC_rt,time=5,col='LightCoral',add=TRUE,title=FALSE,lwd=2)
  legend('bottomright',
         c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
           paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
           paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
         col=c("green",'blue','red'),lwd=2,bty = 'n')
  dev.off()
}
bioSurvival(inputFile="risk.IMvigor.txt", outFile="sur.IMvigor.pdf")
bioROC(inputFile="risk.IMvigor.txt", outFile="ROC.IMvigor.pdf")

#####################

library(limma)
library(ggpubr)
riskFile="risk.IMvigor.txt"      #?????ļ?
cliFile="clinical-TC.txt"           #?ٴ??????ļ?
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk$riskScore[risk$riskScore>quantile(risk$riskScore,0.99)]=quantile(risk$riskScore,0.99)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

samSample=intersect(row.names(risk), row.names(cli))
risk=risk[samSample,"riskScore",drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk, cli)
for(clinical in colnames(rt)[2:ncol(rt)]){
  data=rt[c("riskScore", clinical)]
  colnames(data)=c("riskScore", "clinical")
  data=data[(data[,"clinical"]!="unknow"),]
  #???ñȽ???
  group=levels(factor(data$clinical))
  data$clinical=factor(data$clinical, levels=group)
  comp=combn(group,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  #????????ͼ
  boxplot=ggboxplot(data, x="clinical", y="riskScore", color="clinical",
                    xlab="",
                    ylab="Risk score",
                    legend.title=clinical,
                    add = "jitter")+ 
    stat_compare_means(comparisons = my_comparisons)
  #stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
  #????ͼƬ
  pdf(file=paste0("cliCor.", clinical, ".pdf"), width=6, height=5)
  print(boxplot)
  dev.off()
}


#####################################################
library(plyr)
library(ggplot2)
library(ggpubr)
scoreFile="risk.IMvigor-2.txt"    #m6A?????ļ?
cliFile="clinical-response.txt"            #?ٴ??????ļ?
trait="Response"                    #?ٴ???״
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(score), row.names(cli))
rt=cbind(score[sameSample,,drop=F], cli[sameSample,,drop=F])

bioCol=c("LightCoral","CadetBlue3","#FFB395","SkyBlue3","LightSeaGreen","#6E568C","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(unique(rt[,trait]))]
rt1=rt[,c(trait, "group")]
colnames(rt1)=c("trait", "group")
df=as.data.frame(table(rt1))
df=ddply(df, .(group), transform, percent = Freq/sum(Freq) * 100)
df=ddply(df, .(group), transform, pos = (cumsum(Freq) - 0.5 * Freq))
df$label=paste0(sprintf("%.0f", df$percent), "%")
df$group=factor(df$group, levels=c("low", "high"))
p=ggplot(df, aes(x = factor(group), y = percent, fill = trait)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values=bioCol)+
  xlab("m6Ascore")+ ylab("Percent weight")+  guides(fill=guide_legend(title=trait))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  #coord_flip()+
  theme_bw()
pdf(file="barplot.pdf", width=4, height=5)
print(p)
dev.off()


###############################3
library(survival)
library(survminer)
tmbFile="clinical-Neo.txt"                  
scoreFile="risk.IMvigor-2.txt"     

score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)   
tmb=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)     

sameSample=intersect(row.names(tmb), row.names(score))
tmb=tmb[sameSample,,drop=F]
score=score[sameSample,,drop=F]
data=cbind(score, tmb)

res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("Neo"))
cutoff=as.numeric(res.cut$cutpoint[1])
tmbType=ifelse(data[,"Neo"]<=cutoff, "L-Neo", "H-Neo")
scoreType=ifelse(data$group=="low", "L-m6Ascore", "H-m6Ascore")
mergeType=paste0(tmbType, "+", scoreType)

bioSurvival=function(surData=null, outFile=null){
  diff=survdiff(Surv(futime, fustat) ~ group, data=surData)
  length=length(levels(factor(surData[,"group"])))
  pValue=1-pchisq(diff$chisq, df=length-1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ group, data = surData)
  #print(surv_median(fit))
  
  width=7
  height=7
  if(length(levels(factor(surData[,"group"])))>2){
    width=10
    height=10
  }
  bioCol=c( "LightCoral","CadetBlue3","LightSalmon1","SkyBlue3","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  bioCol=bioCol[1:length]
  surPlot=ggsurvplot(fit, 
                     data=surData,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     legend.title="",
                     legend.labs=levels(factor(surData[,"group"])),
                     font.legend=10,
                     legend = c(0.8, 0.8),
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette = bioCol,
                     surv.median.line = "hv",
                     risk.table=T,
                     cumevents=F,
                     risk.table.height=.25)
  pdf(file=outFile, onefile = FALSE, width=width, height=height)
  print(surPlot)
  dev.off()
}

data$group=tmbType
bioSurvival(surData=data, outFile="Neo.survival.pdf")

data$group=mergeType
bioSurvival(surData=data, outFile="Neo-score.survival.pdf")

####################################################
library(limma)
library(ggpubr)
library(pRRophetic)
library(ggplot2)
library(ridge)
library(preprocessCore)
library(car)
library(genefilter)
library(sva)

Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
display.progress = function (index, totalN, breakN=20) {
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
}  

expr <- read.table("TCGAGTExmerge-TPM.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
normsam <- colnames(expr[,which(substr(colnames(expr),16,17) == "11")])
tumosam <- colnames(expr[,which(substr(colnames(expr),14,15) == "01")])

expr <- expr[,c(tumosam,normsam)]
tp <- c(tumosam)
names(tp) <- tumosam
tp.mutsam <- names(tp) 

########################
mydata <- tp.mutsam
write.csv(mydata, "tp.mutsam.csv", quote = F)

signature <- read.table("7gene.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
dim(signature)

pps <- as.numeric(apply(t(log2(expr[rownames(signature),tp.mutsam] + 1)), 1, function(x) {x %*% signature$Coefficient}))

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
npps <- range01(pps)

mydata <- npps
write.csv(mydata, "npps.csv", quote = F)

Sinfo <- data.frame(PPS = npps,
                    TP = 1,
                    row.names = tp.mutsam,
                    stringsAsFactors = F)
head(Sinfo)

write.csv(Sinfo, "output_PPS.csv", quote = F)

normexpr <- as.matrix(expr[,normsam])
tumoexpr <- as.matrix(expr[,tp53.mutsam])

runpure <- F 
if(runpure) {
  set.seed(123)
  ISOpureS1model <- ISOpure.step1.CPE(tumoexpr, normexpr)
  set.seed(456);
  ISOpureS2model <- ISOpure.step2.PPE(tumoexpr,normexpr,ISOpureS1model)
  pure.tumoexpr <- ISOpureS2model$cc_cancerprofiles
}

if(!runpure) {
  pure.tumoexpr <- tumoexpr
}

auc <- read.table("CTRP_AUC_raw.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T) # Supplementary Data Set 3
auc$comb <- paste(auc$master_cpd_id,auc$master_ccl_id,sep = "-")
auc <- apply(auc[,"area_under_curve",drop = F], 2, function(x) tapply(x, INDEX=factor(auc$comb), FUN=max, na.rm=TRUE)) # 重复项取最大AUC
auc <- as.data.frame(auc)
auc$master_cpd_id <- sapply(strsplit(rownames(auc),"-",fixed = T),"[",1)
auc$master_ccl_id <- sapply(strsplit(rownames(auc),"-",fixed = T),"[",2)
auc <- reshape(auc, 
               direction = "wide",
               timevar = "master_cpd_id",
               idvar = "master_ccl_id")
colnames(auc) <- gsub("area_under_curve.","",colnames(auc),fixed = T)
ctrp.ccl.anno <- read.table("CTRP_ccl_anno.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T) # Supplementary Data Set 1
ctrp.cpd.anno <- read.table("CTRP_cpd_anno.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T) # Supplementary Data Set 2
write.table(auc,"CTRP_AUC.txt",sep = "\t",row.names = F,col.names = T,quote = F)

ctrp.auc <- read.table("CTRP_AUC.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
prism.auc <- read.delim("PRISM_AUC.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T) # 数据来自https://depmap.org/portal/download/ Drug sensitivity AUC (PRISM Repurposing Secondary Screen) 19Q4
prism.ccl.anno <- prism.auc[,1:5] 
prism.auc <- prism.auc[,-c(1:5)]

ctrp.auc <- ctrp.auc[,apply(ctrp.auc,2,function(x) sum(is.na(x))) < 0.5*nrow(ctrp.auc)]
write.table(ctrp.auc,"ctrp.auc.txt",sep = "\t",row.names = T,col.names = T,quote = F)

prism.auc <- prism.auc[,apply(prism.auc,2,function(x) sum(is.na(x))) < 0.5*nrow(prism.auc)]
write.table(prism.auc ,"prism.auc.txt",sep = "\t",row.names = T,col.names = T,quote = F)

rmccl <- paste0("CCL",na.omit(ctrp.ccl.anno[which(ctrp.ccl.anno$ccle_primary_site == "haematopoietic_and_lymphoid_tissue"),"master_ccl_id"]))
rownames(ctrp.auc) <- paste0("CCL",rownames(ctrp.auc))
ctrp.auc <- ctrp.auc[setdiff(rownames(ctrp.auc),rmccl),]

ctrp.auc.knn <- impute.knn(as.matrix(ctrp.auc))$data
prism.auc.knn <- impute.knn(as.matrix(prism.auc))$data

ctrp.auc.knn <- ctrp.auc.knn/ceiling(max(ctrp.auc.knn)) # 参考Expression Levels of Therapeutic Targets as Indicators of Sensitivity to Targeted Therapeutics (2019, Molecular Cancer Therapeutics)
prism.auc.knn <- prism.auc.knn/ceiling(max(prism.auc.knn))

ccl.expr <- read.table("CCLE_RNAseq_rsem_genes_tpm_20180929.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T) 
Ginfo <- read.table("overlapTable27.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T) 
ccl.expr <- ccl.expr[,-1]; rownames(ccl.expr) <- sapply(strsplit(rownames(ccl.expr),".",fixed = T),"[",1)
comgene <- intersect(rownames(ccl.expr),rownames(Ginfo))
ccl.expr <- ccl.expr[comgene,]
ccl.expr$gene <- Ginfo[comgene,"genename"]; ccl.expr <- ccl.expr[!duplicated(ccl.expr$gene),]; rownames(ccl.expr) <- ccl.expr$gene; ccl.expr <- ccl.expr[,-ncol(ccl.expr)]
keepgene <- apply(ccl.expr, 1, mad) > 0.15
trainExpr <- log2(ccl.expr[keepgene,] + 1)
colnames(trainExpr) <- sapply(strsplit(colnames(trainExpr),"_",fixed = T),"[",1) # 重置细胞系名
trainPtype <- as.data.frame(ctrp.auc.knn)
ccl.name <- ccl.miss <- c() 
for (i in rownames(trainPtype)) {
  if(!is.element(gsub("CCL","",i),ctrp.ccl.anno$master_ccl_id)) {
    cat(i,"\n")
    ccl.miss <- c(ccl.miss, i) 
    ccl.name <- c(ccl.name, i) 
  } else {
    ccl.name <- c(ccl.name,  ctrp.ccl.anno[which(ctrp.ccl.anno$master_ccl_id == gsub("CCL","",i)),"ccl_name"]) # 插入匹配的细胞系
  }
}

cpd.name <- cpd.miss <- c() 
for (i in colnames(trainPtype)) {
  if(!is.element(i,ctrp.cpd.anno$master_cpd_id)) {
    cat(i,"\n")
    cpd.miss <- c(cpd.miss, i) 
    cpd.name <- c(cpd.name, i) 
  } else {
    cpd.name <- c(cpd.name,  ctrp.cpd.anno[which(ctrp.cpd.anno$master_cpd_id == i),"cpd_name"]) # 插入匹配的药???
  }
}

rownames(trainPtype) <- ccl.name
trainPtype <- trainPtype[setdiff(rownames(trainPtype),ccl.miss),] # 去除未匹配的细胞???
colnames(trainPtype) <- cpd.name
trainPtype <- trainPtype[,setdiff(colnames(trainPtype),cpd.miss)] # 去除未匹配的药物
comccl <- intersect(rownames(trainPtype),colnames(trainExpr)) # 提取有表达且有药敏的细胞???
trainExpr <- trainExpr[,comccl]
trainPtype <- trainPtype[comccl,]

keepgene <- apply(pure.tumoexpr, 1, mad) > 0.15 
testExpr <- log2(pure.tumoexpr[keepgene,] + 1) 
comgene <- intersect(rownames(trainExpr),rownames(testExpr)) 
write.table(comgene,"comgene.txt",sep = "\t",row.names = T,col.names = T,quote = F)
trainExpr <- as.matrix(trainExpr[comgene,])
testExpr <- testExpr[comgene,]
outTab <- NULL
for (i in 1:ncol(trainPtype)) { 
  display.progress(index = i,totalN = ncol(trainPtype))
  d <- colnames(trainPtype)[i]
  tmp <- log2(as.vector(trainPtype[,d]) + 0.00001) # 由于CTRP的AUC可能???0值，因此加一个较小的数值防止报???

  ptypeOut <- quiet(calcPhenotype(trainingExprData = trainExpr,
                                  trainingPtype = tmp,
                                  testExprData = testExpr,
                                  powerTransformPhenotype = F,
                                  selection = 1))
  ptypeOut <- 2^ptypeOut - 0.00001 # 反对???
  outTab <- rbind.data.frame(outTab,ptypeOut)
}
dimnames(outTab) <- list(colnames(trainPtype),colnames(testExpr))
ctrp.pred.auc <- outTab

keepgene <- apply(ccl.expr, 1, mad) > 0.25
trainExpr <- log2(ccl.expr[keepgene,] + 1)
colnames(trainExpr) <- sapply(strsplit(colnames(trainExpr),"_",fixed = T),"[",1)
trainPtype <- as.data.frame(prism.auc.knn)
rownames(trainPtype) <- prism.ccl.anno[rownames(trainPtype),"cell_line_display_name"]
#colnames(trainPtype) <- sapply(strsplit(colnames(trainPtype)," (",fixed = T), "[",1)
comccl <- intersect(rownames(trainPtype),colnames(trainExpr))
trainExpr <- trainExpr[,comccl]
trainPtype <- trainPtype[comccl,]

keepgene <- apply(pure.tumoexpr, 1, mad) > 0.25
testExpr <- log2(pure.tumoexpr[keepgene,] + 1)
comgene <- intersect(rownames(trainExpr),rownames(testExpr))
trainExpr <- as.matrix(trainExpr[comgene,])
testExpr <- testExpr[comgene,]

outTab <- NULL
for (i in 1:ncol(trainPtype)) { 
  display.progress(index = i,totalN = ncol(trainPtype))
  d <- colnames(trainPtype)[i]
  tmp <- log2(as.vector(trainPtype[,d]) + 0.00001) # 由于PRISM的AUC可能???0值，因此加一个较小的数值防止报???
  ptypeOut <- quiet(calcPhenotype(trainingExprData = trainExpr,
                                  trainingPtype = tmp,
                                  testExprData = testExpr,
                                  powerTransformPhenotype = F,
                                  selection = 1))
  ptypeOut <- 2^ptypeOut - 0.00001 # 反对???
  outTab <- rbind.data.frame(outTab,ptypeOut)
}
dimnames(outTab) <- list(colnames(trainPtype),colnames(testExpr))
prism.pred.auc <- outTab

top.pps <- Sinfo[Sinfo$PPS >= quantile(Sinfo$PPS,probs = seq(0,1,0.1))[10],] 
write.table(top.pps ,"top.pps.txt",sep = "\t",row.names = T,col.names = T,quote = F)
bot.pps <- Sinfo[Sinfo$PPS <= quantile(Sinfo$PPS,probs = seq(0,1,0.1))[2],] 
write.table(bot.pps,"bot.pps.txt",sep = "\t",row.names = T,col.names = T,quote = F)
ctrp.log2fc <- c()
for (i in 1:nrow(ctrp.pred.auc)) {
  display.progress(index = i,totalN = nrow(ctrp.pred.auc))
  d <- rownames(ctrp.pred.auc)[i]
  a <- mean(as.numeric(ctrp.pred.auc[d,rownames(top.pps)])) 
  b <- mean(as.numeric(ctrp.pred.auc[d,rownames(bot.pps)])) 
  fc <- b/a
  log2fc <- log2(fc); names(log2fc) <- d
  ctrp.log2fc <- c(ctrp.log2fc,log2fc)
}
candidate.ctrp <- ctrp.log2fc[ctrp.log2fc > 0.05] 
prism.log2fc <- c()
for (i in 1:nrow(prism.pred.auc)) {
  display.progress(index = i,totalN = nrow(prism.pred.auc))
  d <- rownames(prism.pred.auc)[i]
  a <- mean(as.numeric(prism.pred.auc[d,rownames(top.pps)])) 
  b <- mean(as.numeric(prism.pred.auc[d,rownames(bot.pps)])) 
  fc <- b/a
  log2fc <- log2(fc); names(log2fc) <- d
  prism.log2fc <- c(prism.log2fc,log2fc)
}
candidate.prism <- prism.log2fc[prism.log2fc > 0.05] 
ctrp.cor <- ctrp.cor.p <- c()
for (i in 1:nrow(ctrp.pred.auc)) {
  display.progress(index = i,totalN = nrow(ctrp.pred.auc))
  d <- rownames(ctrp.pred.auc)[i]
  a <- as.numeric(ctrp.pred.auc[d,rownames(Sinfo)]) 
  b <- as.numeric(Sinfo$PPS)
  r <- cor.test(a,b,method = "spearman")$estimate; names(r) <- d
  p <- cor.test(a,b,method = "spearman")$p.value; names(p) <- d
  ctrp.cor <- c(ctrp.cor,r)
  ctrp.cor.p <- c(ctrp.cor.p,p)
}
candidate.ctrp2 <- ctrp.cor[ctrp.cor < -0.05]  #Spearman相关系数
ctrp.candidate <- intersect(names(candidate.ctrp),names(candidate.ctrp2))

prism.cor <- prism.cor.p <- c()
for (i in 1:nrow(prism.pred.auc)) {
  display.progress(index = i,totalN = nrow(prism.pred.auc))
  d <- rownames(prism.pred.auc)[i]
  a <- as.numeric(prism.pred.auc[d,rownames(Sinfo)]) 
  b <- as.numeric(Sinfo$PPS)
  r <- cor.test(a,b,method = "spearman")$estimate; names(r) <- d
  p <- cor.test(a,b,method = "spearman")$p.value; names(p) <- d
  prism.cor <- c(prism.cor,r)
  prism.cor.p <- c(prism.cor.p,p)
}
candidate.prism2 <- prism.cor[prism.cor < -0.15]  #Spearman相关系数
prism.candidate <- intersect(names(candidate.prism),names(candidate.prism2))

darkblue <- "#0772B9"
lightblue <- "#48C8EF"

cor.data <- data.frame(drug = ctrp.candidate,
                       r = ctrp.cor[ctrp.candidate],
                       p = -log10(ctrp.cor.p[ctrp.candidate]))
p1 <- ggplot(data = cor.data,aes(r,forcats::fct_reorder(drug,r,.desc = T))) +
  geom_segment(aes(xend=0,yend=drug),linetype = 2) +
  geom_point(aes(size=p),col = darkblue) +
  scale_size_continuous(range =c(2,8)) +
  scale_x_reverse(breaks = c(0, -0.3, -0.5),
                  expand = expansion(mult = c(0.01,.1))) + #左右留空
  theme_classic() +
  labs(x = "Correlation coefficient", y = "", size = bquote("-log"[10]~"("~italic(P)~"-value)")) + 
  theme(legend.position = "bottom", 
        axis.line.y = element_blank())

cor.data <- data.frame(drug = prism.candidate,
                       r = prism.cor[prism.candidate],
                       p = -log10(prism.cor.p[prism.candidate]))
cor.data$drug <- sapply(strsplit(cor.data$drug," (",fixed = T), "[",1)

p2 <- ggplot(data = cor.data,aes(r,forcats::fct_reorder(drug,r,.desc = T))) +
  geom_segment(aes(xend=0,yend=drug),linetype = 2) +
  geom_point(aes(size=p),col = darkblue) +
  scale_size_continuous(range =c(2,8)) +
  scale_x_reverse(breaks = c(0, -0.3, -0.5),
                  expand = expansion(mult = c(0.01,.1))) + #左右留空
  theme_classic() +
  labs(x = "Correlation coefficient", y = "", size = bquote("-log"[10]~"("~italic(P)~"-value)")) + 
  theme(legend.position = "bottom", 
        axis.line.y = element_blank())

ctrp.boxdata <- NULL
for (d in ctrp.candidate) {
  a <- as.numeric(ctrp.pred.auc[d,rownames(top.pps)]) 
  b <- as.numeric(ctrp.pred.auc[d,rownames(bot.pps)])
  p <- wilcox.test(a,b)$p.value
  s <- as.character(cut(p,c(0,0.001,0.01,0.05,1),labels = c("***","**","*","")))
  ctrp.boxdata <- rbind.data.frame(ctrp.boxdata,
                                   data.frame(drug = d,
                                              auc = c(a,b),
                                              p = p,
                                              s = s,
                                              group = rep(c("High PPS","Low PPS"),c(nrow(top.pps),nrow(bot.pps))),
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
}

drug = d

p3 <- ggplot(ctrp.boxdata, aes(drug, auc, fill=group)) + 
  geom_boxplot(aes(col = group),outlier.shape = NA) + 
  # geom_text(aes(drug, y=min(auc) * 1.1, 
  #               label=paste("p=",formatC(p,format = "e",digits = 1))),
  #           data=ctrp.boxdata, 
  #           inherit.aes=F) + 
  geom_text(aes(drug, y=max(auc)), 
            label=ctrp.boxdata$s,
            data=ctrp.boxdata, 
            inherit.aes=F) + 
  scale_fill_manual(values = c(darkblue, lightblue)) + 
  scale_color_manual(values = c(darkblue, lightblue)) + 
  xlab(NULL) + ylab("Estimated AUC value") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5,vjust = 0.5,size = 10),
        legend.position = "bottom",
        legend.title = element_blank()) 

dat <- ggplot_build(p3)$data[[1]]

p3 <- p3 + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)

prism.boxdata <- NULL
for (d in prism.candidate) {
  a <- as.numeric(prism.pred.auc[d,rownames(top.pps)]) 
  b <- as.numeric(prism.pred.auc[d,rownames(bot.pps)])
  p <- wilcox.test(a,b)$p.value
  s <- as.character(cut(p,c(0,0.001,0.01,0.05,1),labels = c("***","**","*","")))
  prism.boxdata <- rbind.data.frame(prism.boxdata,
                                    data.frame(drug = d,
                                               auc = c(a,b),
                                               p = p,
                                               s = s,
                                               group = rep(c("High PPS","Low PPS"),c(nrow(top.pps),nrow(bot.pps))),
                                               stringsAsFactors = F),
                                    stringsAsFactors = F)
}
prism.boxdata$drug <- sapply(strsplit(prism.boxdata$drug," (",fixed = T), "[",1)

p4 <- ggplot(prism.boxdata, aes(drug, auc, fill=group)) + 
  geom_boxplot(aes(col = group),outlier.shape = NA) + 
  # geom_text(aes(drug, y=min(auc) * 1.1, 
  #               label=paste("p=",formatC(p,format = "e",digits = 1))),
  #           data=prism.boxdata, 
  #           inherit.aes=F) + 
  geom_text(aes(drug, y=max(auc)), 
            label=prism.boxdata$s,
            data=prism.boxdata, 
            inherit.aes=F) + 
  scale_fill_manual(values = c(darkblue, lightblue)) + 
  scale_color_manual(values = c(darkblue, lightblue)) + 
  xlab(NULL) + ylab("Estimated AUC value") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5,vjust = 0.5,size = 10),
        legend.position = "bottom",
        legend.title = element_blank())
dat <- ggplot_build(p4)$data[[1]]

p4 <- p4 + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)

plot_grid(p1, p3, p2, p4, labels=c("A", "", "B", ""), 
          ncol=2, 
          rel_widths = c(2, 2)) 
ggsave(filename = "drug target.pdf",width = 8,height = 8)

sessionInfo()

######################################
rm(list=ls())
gc()
library(Seurat) 
library(tidyverse)
library(dplyr)
library(magrittr)
library(ggplot2)
library(Matrix)
library(clusterProfiler)
library(org.Hs.eg.db)
library(RColorBrewer)
library(gplots)
library(patchwork)

samples=list.files("data/")
samples
dir <- file.path('data',samples)
dir

#合并矩阵
counts <- lapply(1:length(dir),function(i){read.csv2(dir[i],header = T,sep = ",")})
sample <- counts[1]
sample <- data.frame(sample,row.names = "X")
sample <- as.matrix(sample,"dgCMatrix")
colnames(sample) <- paste0(str_split(samples[1],'_')[[1]][2],'_',colnames(sample))
sample_all <- sample
sample_all[1:10,1:10]
for (i in 2:5) {
  sample <- counts[i]
  sample <- data.frame(sample,row.names = "X")
  sample <- as.matrix(sample,"dgCMatrix")
  colnames(sample) <- paste0(str_split(samples[i],'_')[[1]][2],'_',colnames(sample))
  sample_all <- cbind(sample_all,sample)
}
sample_all[1:10,8000:8010]
class(sample_all)
sample_all <- as(sample_all,"sparseMatrix")#变为稀疏矩阵

rownames(sample_all)
which(duplicated(colnames(sample_all)))#查看有无相同的细胞

#基因名转换
dim(sample_all)
which(duplicated(row.names(sample_all)))
gene.df <- bitr(row.names(sample_all), fromType = "ENSEMBL", # 从SYMBOL到ENSEMBL和ENTREZID,得到新命名的gene.df
                toType = c("SYMBOL"),
                OrgDb = org.Hs.eg.db) #小鼠换成org.Mm.eg.db
head(gene.df)
dim(gene.df)
anyNA(gene.df$ENSEMBL)
anyNA(gene.df$SYMBOL)
gene.df <- na.omit(gene.df)
dim(gene.df)
which(duplicated(gene.df$SYMBOL))
which(duplicated(gene.df$ENSEMBL))
gene.df <- gene.df[!duplicated(gene.df$SYMBOL),]
gene.df <- gene.df[!duplicated(gene.df$ENSEMBL),]
dim(gene.df)
pos=match(gene.df$ENSEMBL,rownames(sample_all) )
sample_all <- sample_all[pos,]
dim(sample_all)
rownames(sample_all)=gene.df$SYMBOL
sample_all[1:10,1:10]

saveRDS(sample_all,file="scRNA_all_STEP1.1.rds")

### STEP 1.2 创建seurat对象 ----
#整理metadata
sample_name <- str_split(colnames(sample_all), "_", simplify = T)[,1]
# sample <- str_match(colnames(counts),".+_\\d+d")[,1]#'.表示任意字符+表示匹配不止1次，\\d+表示多个数字
table(sample_name)

#创建seurat对象
scRNA_all = CreateSeuratObject(sample_all)

#写入metadata
scRNA_all@meta.data$sample_name<- sample_name
current.cluster.ids <- c("Normal","OM1",'OM2','PS1','PS2')
new.cluster.ids <- c('Normal','Metastasis','Metastasis','Primary','Primary')
scRNA_all@meta.data$sample_site <- plyr::mapvalues(x = as.character(scRNA_all@meta.data$sample_name), from = current.cluster.ids, to = new.cluster.ids)
table(scRNA_all$sample_name)
table(scRNA_all$orig.ident)
table(scRNA_all$sample_site)
dim(scRNA_all)   #查看基因数和细胞总数
tail(rownames(scRNA_all),3)
head(rownames(scRNA_all),3)
head(scRNA_all@meta.data) 

saveRDS(scRNA_all,file="scRNA_all_STEP1.2.rds")

#### STEP 2 数据质控 ----
### STEP 2.1 质控前处理 ----
rm(list=ls())
library(limma)
library(tidyverse)
library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)
library(patchwork)
library(Matrix)

dir.create("QC")
raw_sce <- readRDS("scRNA_all_STEP1.2.rds")
tail(rownames(raw_sce),3)
head(rownames(raw_sce),3)
head(raw_sce@meta.data)

#检测线粒体基因
raw_sce[["percent.MT"]] <- PercentageFeatureSet(object = raw_sce, pattern = "^MT-")
fivenum(raw_sce[["percent.MT"]][,1])

##检测核糖体体基因
raw_sce[["percent.RP"]] <- PercentageFeatureSet(object = raw_sce, pattern = "^RP[SL]")
fivenum(raw_sce[["percent.RP"]][,1])

#计算红细胞比例（即血红蛋白的编码基因）
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(raw_sce@assays$RNA)) 
HB.genes <- rownames(raw_sce@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
raw_sce[["percent.HB"]] <- PercentageFeatureSet(object = raw_sce, features = HB.genes) #鼠的血红蛋白基因以Hba,Hbb,Hbq开头
fivenum(raw_sce[["percent.HB"]][,1])

violin <- VlnPlot(raw_sce,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.MT","percent.HB","percent.RP"), 
                  group.by = "sample_name", 
                  pt.size = 0.0001, #不需要显示点，可以设置pt.size = 0
                  ncol = 5)
ggsave("QC/vlnplot_before_qc.pdf", plot = violin, width = 25, height = 5) 
ggsave("QC/vlnplot_before_qc.png", plot = violin, width = 25, height = 5)

#查看各指标的细胞分布
#转录本数量
p1=raw_sce@meta.data %>% 
  ggplot(aes(color=sample_name, x=nCount_RNA, fill= sample_name)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 1000)
#基因数
p2=raw_sce@meta.data %>% 
  ggplot(aes(color=sample_name, x=nFeature_RNA, fill= sample_name)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
#线粒体基因占比
p3=raw_sce@meta.data %>% 
  ggplot(aes(color=sample_name, x=percent.MT, fill= sample_name)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 20)
#核糖体基因占比
p4=raw_sce@meta.data %>% 
  ggplot(aes(color=sample_name, x=percent.RP, fill= sample_name)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 20)
#查看基因数和转录本数之间的关系，看是否存在基因/转录本特别低的细胞
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
p5=raw_sce@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.MT)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 1000) +
  geom_hline(yintercept = 500) +
  facet_wrap(~sample_name)

pc <- p1+p2+p3+p4+plot_layout(nrow=2,guides = 'collect') 
ggsave("QC/densityplot_before_qc.pdf", plot = pc, width = 10, height = 10) 
ggsave("QC/densityplot_before_qc.png", plot = pc, width = 10, height = 10)
ggsave("QC/FeaturetoCount_before_qc.pdf", plot = p5, width = 10, height = 10) 
ggsave("QC/FeaturetoCount_before_qc.png", plot = p5, width = 10, height = 10)

plot1 <- FeatureScatter(raw_sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="sample_name",raster=F)+
  theme(plot.title  = element_blank(),legend.position="none")
plot2 <- FeatureScatter(raw_sce, feature1 = "nCount_RNA", feature2 = "percent.MT", group.by="sample_name", raster=F)+
  theme(plot.title  = element_blank(),legend.position="none")
plot3 <- FeatureScatter(raw_sce, feature1 = "nCount_RNA", feature2 = "percent.HB", group.by="sample_name", raster=F)+
  theme(plot.title  = element_blank(),legend.position="none")
plot4 <- FeatureScatter(raw_sce, feature1 = "nCount_RNA", feature2 = "percent.RP", group.by="sample_name", raster=F)+
  theme(plot.title  = element_blank(),legend.position="none")
pearplot <- plot1+plot2+plot3+plot4+plot_layout(nrow=2) 
ggsave("QC/pearplot_before_qc.pdf", plot = pearplot, width = 10, height = 10) 
ggsave("QC/pearplot_before_qc.png", plot = pearplot, width = 10, height = 10)

saveRDS(raw_sce, file = "scRNA_all_STEP2.1.rds")
head(raw_sce@meta.data)

### STEP 2.2质控后处理 ----
rm(list=ls())
raw_sce <- readRDS("scRNA_all_STEP2.1.rds")
raw_sce
head(raw_sce@meta.data)
#按照总基因数目的百分比来取质控值
nfeature_q <-quantile(raw_sce@meta.data$nFeature_RNA, c(0.05, 0.95))

#设置固定值过滤线粒体基因并考虑线粒体基因表达分布情况
Filter_high_mito_genes <- 
  function(mito_quantile, max_ratio){
    for (i in 1:length(mito_quantile)){
      while (mito_quantile[i] >= max_ratio){
        i = i - 1
        print(mito_quantile[i])
        return(m_raw[i])
        break
      }
    }
  }
m_raw <- quantile(raw_sce@meta.data$percent.MT, probs=seq(0,1,0.01)) 
mt_q <- Filter_high_mito_genes(m_raw, 15)#最高的线粒体基因占比设置为15%

#线粒体基因、核糖体基因占比的百分比95%的值，要去除占比过高的细胞，所以只有上限
# mt_q <- quantile(raw_sce@meta.data$percent.mt, 0.95)#直接设定阈值
rp_q <- quantile(raw_sce@meta.data$percent.rp, 0.95)

# QC
#去除表达基因过多或过少的细胞,可以根据分布来提取，也可以指定具体的数值
#根据分布
dim(raw_sce)
table(raw_sce$sample_name)
scRNA1 <- subset(raw_sce, 
                 subset = nFeature_RNA > as.numeric(nfeature_q)[1] & 
                   nFeature_RNA < as.numeric(nfeature_q)[2] )
dim(scRNA1)
table(scRNA1$sample_name)
#指定数值
scRNA2 <- subset(raw_sce, 
                 subset = nFeature_RNA > 500 & 
                   nFeature_RNA < 5000 )
raw_sce
scRNA1
scRNA2
table(raw_sce@meta.data$sample_name)
table(scRNA1@meta.data$sample_name)
table(scRNA2@meta.data$sample_name)

scRNA <- scRNA1
rm(raw_sce,scRNA1,scRNA2)
#去除线粒体基因表达过高的细胞
scRNA <- subset(scRNA, 
                subset = percent.mt < mt_q )
scRNA
# #去除核糖体基因表达过高的细胞
# #（核糖体基因表达过高提示RNA降解，核糖体基因对分群没有帮助，但是很少有依据某一个阈值过滤核糖体基因，要么不过滤，要么就把核糖体基因全部删除）
# scRNA <- subset(scRNA, 
#                    subset = percent.rp < rp_q )
# scRNA
#去除红细胞基因表达过高的细胞
dim(scRNA)
table(scRNA$sample_name)
scRNA1 <- subset(scRNA, 
                 subset = percent.HB < 0.1)#红细胞基因一般表达均比较低，故阈值根据上面的分布选择
dim(scRNA1)
table(scRNA1$sample_name)
scRNA <- scRNA1
rm(scRNA1)
#重新查看各指标的细胞分布
#转录本数量
p1=scRNA@meta.data %>% 
  ggplot(aes(color=sample_name, x=nCount_RNA, fill= sample_name)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 1000)
#基因数
p2=scRNA@meta.data %>% 
  ggplot(aes(color=sample_name, x=nFeature_RNA, fill= sample_name)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
#线粒体基因占比
p3=scRNA@meta.data %>% 
  ggplot(aes(color=sample_name, x=percent.MT, fill= sample_name)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 20)
#核糖体基因占比
p4=scRNA@meta.data %>% 
  ggplot(aes(color=sample_name, x=percent.RP, fill= sample_name)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 20)
#查看基因数和转录本数之间的关系，看是否存在基因/转录本特别低的细胞
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
p5=scRNA@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.MT)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 1000) +
  geom_hline(yintercept = 500) +
  facet_wrap(~sample_name)

pc <- p1+p2+p3+p4+plot_layout(nrow=2,guides = 'collect') 
ggsave("QC/densityplot_after_qc.pdf", plot = pc, width = 10, height = 10) 
ggsave("QC/densityplot_after_qc.png", plot = pc, width = 10, height = 10)
ggsave("QC/FeaturetoCount_after_qc.pdf", plot = p5, width = 10, height = 10) 
ggsave("QC/FeaturetoCount_after_qc.png", plot = p5, width = 10, height = 10)

violin <-VlnPlot(scRNA,
                 features = c("nFeature_RNA", "nCount_RNA", "percent.MT","percent.HB","percent.RP"), 
                 group.by = "sample_name", 
                 pt.size = 0.0001, 
                 ncol = 5)
ggsave("QC/vlnplot_after_qc.pdf", plot = violin, width = 25, height = 5) 
ggsave("QC/vlnplot_after_qc.png", plot = violin, width = 25, height = 5)

plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="sample_name",raster=F)+
  theme(plot.title  = element_blank(),legend.position="none")
plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.MT", group.by="sample_name", raster=F)+
  theme(plot.title  = element_blank(),legend.position="none")
plot3 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.HB", group.by="sample_name", raster=F)+
  theme(plot.title  = element_blank(),legend.position="none")
plot4 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.RP", group.by="sample_name", raster=F)+
  theme(plot.title  = element_blank(),legend.position="none")
pearplot <- plot1+plot2+plot3+plot4+plot_layout(nrow=2) 
ggsave("QC/pearplot_after_qc.pdf", plot = pearplot, width = 10, height = 10) 
ggsave("QC/pearplot_after_qc.png", plot = pearplot, width = 10, height = 10)

saveRDS(scRNA, file = "scRNA_all_STEP2.2.rds")

#### STEP 3 标准化、降维和聚类----
rm(list=ls())
options(stringsAsFactors = F)
gc()
library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)
library(magrittr)
library(ggplot2)
library(gplots)
library(ggsci)
library(reshape2)
library(glmGamPoi)
dir.create("Cluster")
dir.create("Marker")
dir.create("SCT")
scRNA <- readRDS("scRNA_all_STEP2.2.rds")
d1="SCT" #或者"SCT"
d2="" #或者"_Scaleall"
head(scRNA@meta.data)

### STEP 3.2 SCTransform ----
#包括了标准化，中心化寻找高变基因3个步骤
#先根据SCT过程看mt、rp和细胞周期基因相关蛋白对分群的影响，再决定SCT过程是否排除这几个基因
scRNA <- SCTransform(scRNA, method= 'glmGamPoi', verbose = F)
top10 <- head(VariableFeatures(scRNA), 10) 
plot1 <- VariableFeaturePlot(scRNA)+
  theme(legend.position="none")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge =0, ynudge=0,size=2.5) +
  theme(legend.position="none")
plot <- plot1|plot2 
ggsave("SCT/3.2_VariableFeatures_before.pdf", plot = plot, width = 10, height = 5) 
ggsave("SCT/3.2_VariableFeatures_before.png", plot = plot, width = 10, height = 5)
#查看我们选择的高变基因中有哪些细胞周期相关基因
# CaseMatch(c('Mcm5','Pcna','Tyms','Fen1','Mcm2','Mcm4'),VariableFeatures(scRNAsub))
#在meta.data中添加S.Score、G2M.Score和Phase三列有关细胞周期的信息。
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(scRNA))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(scRNA))
scRNA <- CellCycleScoring(object=scRNA,  g2m.features=g2m_genes,  s.features=s_genes)
head(scRNA@meta.data)
scRNAa <- RunPCA(scRNA, features = c(s_genes, g2m_genes))
p <- DimPlot(scRNAa, reduction = "pca", group.by = "Phase")
ggsave("SCT/3.2 cellcycle_pca_before.pdf", p, width = 5, height = 5)
ggsave("SCT/3.2 cellcycle_pca_before.png", p, width = 5, height = 5)
saveRDS(scRNA, file="scRNA_all_STEP3.2.1.rds")

# 去除线粒体基因、核糖体基因和细胞周期基因的影响
scRNA <- SCTransform(scRNA, method= 'glmGamPoi', verbose = F, 
                     vars.to.regress = c("percent.RP","S.Score", "G2M.Score"))
top10 <- head(VariableFeatures(scRNA), 10) 
plot1 <- VariableFeaturePlot(scRNA)+
  theme(legend.position="none")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge =0, ynudge=0,size=2.5) +
  theme(legend.position="none")
plot <- plot1|plot2 
ggsave("SCT/3.2_VariableFeatures_after_SCTregress.pdf", plot = plot, width = 10, height = 5) 
ggsave("SCT/3.2_VariableFeatures_after_SCTregress.png", plot = plot, width = 10, height = 5)
scRNAa <- RunPCA(scRNA, features = c(s_genes, g2m_genes))
p <- DimPlot(scRNAa, reduction = "pca", group.by = "Phase")
ggsave("SCT/3.2 cellcycle_pca_after_SCTregress.pdf", p, width = 5, height = 5)
ggsave("SCT/3.2 cellcycle_pca_after_SCTregress.png", p, width = 5, height = 5)
saveRDS(scRNA, file="scRNA_all_STEP3.2.2.rds")

### STEP 3.3 PCA降维 ----
rm(list=ls())
gc()
scRNA <- readRDS(file="scRNA_all_STEP3.2.2.rds")
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))
# #每个细胞在PC轴上的坐标
# head(scRNAsub@reductions$pca@cell.embeddings)
# #每个基因对每个PC轴的贡献度（loading值）
# head(scRNAsub@reductions$pca@feature.loadings)

# # Get the feature loadings for a given DimReduc
# t(Loadings(object = scRNAsub[["pca"]])[1:5,1:5])
# # Get the feature loadings for a specified DimReduc in a Seurat object
# t(Loadings(object = scRNAsub, reduction = "pca")[1:5,1:5])
# Set the feature loadings for a given DimReduc
new.loadings <- Loadings(object = scRNA[["pca"]])
new.loadings <- new.loadings + 0.01
Loadings(object = scRNA[["pca"]]) <- new.loadings
VizDimLoadings(scRNA)

# Examine and visualize PCA results a few different ways
# print(scRNAsub[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(scRNA, reduction = "pca",group.by = 'sample_name')
DimHeatmap(scRNA, dims = 1:30, cells = 500, balanced = TRUE)
# DimHeatmap(scRNAsub, dims = 20:40, cells = 500, balanced = TRUE)

#用SCT则不用Jackstraw和ScoreJackstraw
# scRNAsub <- JackStraw(scRNAsub, num.replicate = 100)
# scRNAsub <- ScoreJackStraw(scRNAsub, dims = 1:20)
# plot3 <-JackStrawPlot(scRNAsub, dims = 1:20)
plot1 <- DimPlot(scRNA, reduction = "pca", group.by="sample_name") 
plot2 <- ElbowPlot(scRNA, ndims=40, reduction="pca") 
plotc <- plot1|plot2
ggsave("Cluster/3.3 pca_afterregress.pdf", plot = plotc, width = 10, height = 5) 
ggsave("Cluster/3.3 pca_afterregress.png", plot = plotc, width = 10, height = 5)
ElbowPlot(scRNA, ndims=40, reduction="pca")

saveRDS(scRNA, file="scRNA_all_STEP3.3.rds")

### STEP 3.4 细胞聚类 ----
rm(list=ls())
options(stringsAsFactors = F)
gc()
scRNAsub <- readRDS("scRNA_all_STEP3.3.rds")
pc.num=1:40
# 运行harmony ----
library(harmony)
scRNA_Harmony <- RunHarmony(scRNAsub,"sample_name",plot_convergence=T, assay.use = "SCT")

saveRDS(scRNA_Harmony,file = paste0("scRNA_all_STEP3.4.1_Harmony.rds"))

scRNAsub <- scRNA_Harmony
rm(scRNA_Harmony)
scRNAsub <- FindNeighbors(scRNAsub,reduction = "harmony",dims = pc.num) 

for (res in c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
              0.6, 0.7, 0.8,0.9, 1, 1.2)) {
  print(res)
  scRNAsub <- FindClusters(scRNAsub, reduction = "harmony", resolution = res)
}
#查看分群的变化
library(clustree)
pdf(file = "Cluster/3.4 FindClusters_clustree.pdf", width = 15, height = 10)
clustree(scRNAsub) +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Paired") +
  scale_edge_color_continuous(low = "grey80", high = "red")
dev.off()

d1="Cluster"
re=0.2
tr="afterharmony_"
# Idents(object = scRNAsub) <- "SCT_snn_res.0.05"
scRNAsub <- FindClusters(scRNAsub,reduction = "harmony", resolution = re)
table(scRNAsub@meta.data$seurat_clusters)
metadata <- scRNAsub@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,paste0(d1,'/3.4_cell_cluster_',tr,re,'.csv'),row.names = F)
##非线性降维 ----
#tSNE
scRNAsub = RunTSNE(scRNAsub, reduction = "harmony", dims = pc.num)
embed_tsne <- Embeddings(scRNAsub, 'tsne')
write.csv(embed_tsne,paste0(d1,'/3.4_embed_tsne_',tr,re,'.csv'))
#UMAP
scRNAsub <- RunUMAP(scRNAsub, reduction = "harmony", dims = pc.num)
embed_umap <- Embeddings(scRNAsub, 'umap')
write.csv(embed_umap,paste0(d1,'/3.4_embed_umap_',tr,re,'.csv'))
#可视化
plot1 = DimPlot(scRNAsub, reduction = "tsne", label = T) 
ggsave(paste0(d1,"/3.4_tSNE_",tr,re,".pdf"), plot = plot1, width = 8, height = 7)
ggsave(paste0(d1,"/3.4_tSNE_",tr,re,".png"), plot = plot1, width = 8, height = 7)
plot2 = DimPlot(scRNAsub, reduction = "umap", label = T) 
ggsave(paste0(d1,"/3.4_UMAP_",tr,re,".pdf"), plot = plot2, width = 8, height = 7)
ggsave(paste0(d1,"/3.4_UMAP_",tr,re,".png"), plot = plot2, width = 8, height = 7)
#合并tSNE与UMAP
plotc <- plot1+plot2+ plot_layout(guides = 'collect')
ggsave(paste0(d1,"/3.4_tSNE_UMAP_",tr,re,".pdf"), plot = plotc, width = 15, height = 7)
ggsave(paste0(d1,"/3.4_tSNE_UMAP_",tr,re,".png"), plot = plotc, width = 15, height = 7)

p1=DimPlot(scRNAsub, reduction = "tsne",label = F, group.by  = 'sample_name')
ggsave(paste0(d1,"/3.4_tSNE_groupby_sample_name_",tr,re,".pdf"), plot = p1, width = 8, height = 7)
ggsave(paste0(d1,"/3.4_tSNE_groupby_sample_name_",tr,re,".png"), plot = p1, width = 8, height = 7)
p2=DimPlot(scRNAsub, reduction = "umap", label = F,group.by  = 'sample_name')
ggsave(paste0(d1,"/3.4_UMAP_groupby_sample_name_",tr,re,".pdf"), plot = p2, width = 8, height = 7)
ggsave(paste0(d1,"/3.4_UMAP_groupby_sample_name_",tr,re,".png"), plot = p2, width = 8, height = 7)
pc <- p1|p2+ plot_layout(guides = 'collect')
ggsave(paste0(d1,"/3.4_tSNE_UMAP_groupby_v_",tr,re,".pdf"), plot = pc, width = 15, height = 7)
ggsave(paste0(d1,"/3.4_tSNE_UMAP_groupby_sample_name_",tr,re,".png"), plot = pc, width = 15, height = 7)

p3=DimPlot(scRNAsub, reduction = "tsne", split.by = 'sample_name')
ggsave(paste0(d1,"/3.4_tSNE_splitby_sample_name_",tr,re,".pdf"), plot = p3, width = 29, height = 7)
ggsave(paste0(d1,"/3.4_tSNE_splitby_sample_name_",tr,re,".png"), plot = p3, width = 29, height = 7)
p4=DimPlot(scRNAsub, reduction = "umap", split.by = 'sample_name')
ggsave(paste0(d1,"/3.4_UMAP_splitby_sample_name_",tr,re,".pdf"), plot = p4, width = 29, height = 7)
ggsave(paste0(d1,"/3.4_UMAP_splitby_sample_name_",tr,re,".png"), plot = p4, width = 29, height = 7)
pc2 <- p3/p4+ plot_layout(guides = 'collect')
ggsave(paste0(d1,"/3.4_tSNE_UMAP_splitby_sample_name_",tr,re,".pdf"), plot = pc2, width = 29, height = 14)
ggsave(paste0(d1,"/3.4_tSNE_UMAP_splitby_sample_name_",tr,re,".png"), plot = pc2, width = 29, height = 14)

p5=DimPlot(scRNAsub, reduction = "tsne",label = F, group.by  = 'sample_site')
ggsave(paste0(d1,"/3.4_tSNE_groupby_sample_site_",tr,re,".pdf"), plot = p1, width = 8, height = 7)
ggsave(paste0(d1,"/3.4_tSNE_groupby_sample_site_",tr,re,".png"), plot = p1, width = 8, height = 7)
p6=DimPlot(scRNAsub, reduction = "umap", label = F,group.by  = 'sample_site')
ggsave(paste0(d1,"/3.4_UMAP_groupby_sample_site_",tr,re,".pdf"), plot = p2, width = 8, height = 7)
ggsave(paste0(d1,"/3.4_UMAP_groupby_sample_site_",tr,re,".png"), plot = p2, width = 8, height = 7)
pc <- p5|p6+ plot_layout(guides = 'collect')
ggsave(paste0(d1,"/3.4_tSNE_UMAP_groupby_sample_site_",tr,re,".pdf"), plot = pc, width = 15, height = 7)
ggsave(paste0(d1,"/3.4_tSNE_UMAP_groupby_sample_site_",tr,re,".png"), plot = pc, width = 15, height = 7)

p7=DimPlot(scRNAsub, reduction = "tsne", split.by = 'sample_site')
ggsave(paste0(d1,"/3.4_tSNE_splitby_sample_id_",tr,re,".pdf"), plot = p3, width = 29, height = 7)
ggsave(paste0(d1,"/3.4_tSNE_splitby_sample_id_",tr,re,".png"), plot = p3, width = 29, height = 7)
p8=DimPlot(scRNAsub, reduction = "umap", split.by = 'sample_site')
ggsave(paste0(d1,"/3.4_UMAP_splitby_sample_site_",tr,re,".pdf"), plot = p4, width = 29, height = 7)
ggsave(paste0(d1,"/3.4_UMAP_splitby_sample_site_",tr,re,".png"), plot = p4, width = 29, height = 7)
pc2 <- p7/p8+ plot_layout(guides = 'collect')
ggsave(paste0(d1,"/3.4_tSNE_UMAP_splitby_sample_site_",tr,re,".pdf"), plot = pc2, width = 29, height = 14)
ggsave(paste0(d1,"/3.4_tSNE_UMAP_splitby_sample_site_",tr,re,".png"), plot = pc2, width = 29, height = 14)


p9=DimPlot(scRNAsub, reduction = "tsne", label = F, group.by = 'Phase')
p10=DimPlot(scRNAsub, reduction = "umap", label = F, group.by = 'Phase')
pc1=p9|p10
ggsave(paste0(d1,"/3.4_tSNE_UMAP_groupby_CellcylePhase_",tr,re,".pdf"), plot = pc1, width = 15, height = 7)
ggsave(paste0(d1,"/3.4_tSNE_UMAP_groupby_CellcylePhase_",tr,re,".pdf"), plot = pc1, width = 15, height = 7)

#查看每一个cluster分别来自哪几个sample
pdf(file = paste0(d1,"/3.4_cluster_to_samples_ballonplot_",tr,re,".pdf"), width = 10, height = 8)
tab.1=table(scRNAsub@meta.data$seurat_clusters,scRNAsub@meta.data$sample_name) 
balloonplot(tab.1)
dev.off()

#分析各个群所占百分比
tab.2=table(scRNAsub@meta.data$sample_name, scRNAsub@meta.data$seurat_clusters) 
tab.2
statistics = apply(tab.2, 1, sum)   # 得到每个样本的观测值总和
plot(statistics)
data_percent = data.frame()   # 建立空数据框
for (n in 1:5) {
  data_percent = rbind(data_percent, tab.2[n,] / statistics[n] )  
} 
# 再来看下，每个样本总和都等于1，现在符合要求了
statistics = apply(data_percent, 1, sum)   
plot(statistics)
# 再加上样本的命名信息，方便看图，已有命名的请忽略
colnames(data_percent)=c(0:13)
data_percent$names = c('Normal',"OM1",'OM2',"PS1","PS2")
data_percent 
data_plot = melt(data_percent)
colnames(data_plot) = c('sample_id','cluster','percent')
data_plot
dev.off()

pdf (file = paste0(d1,"/3.4_clusterpercent_to_samples_barplot_",tr,re,".pdf"), width = 5, height = 5)
ggplot( data_plot, aes( x = sample_id, weight = percent, fill = cluster))+
  geom_bar( position = "stack")+theme_classic()
dev.off()

pdf (file = paste0(d1,"/3.4_clusterpercent_to_samples_lineplot_",tr,re,".pdf"), width = 5, height = 5)
ggplot(data_plot, aes(x =sample_id, y = percent, group =cluster, color = cluster))+
  geom_line()+theme_classic()
dev.off()


saveRDS(scRNAsub, file=paste0("scRNA_all_STEP3.4.2_",tr,re,".rds"))

# Cluster差异分析 ----
library(Seurat)
library(patchwork)
library(magrittr)
library(ggplot2)
library(gplots)
library(Matrix)
library(harmony)
library(tidyverse)
rm(list=ls())
gc()
options(stringsAsFactors = F)
d1="Cluster"
d2="Marker"
re=0.2
tr='afterharmony_'
scRNAsub <- readRDS(paste0("scRNA_all_STEP3.4.2_",tr,re,".rds"))
head(scRNAsub@meta.data)

diff.wilcox = FindAllMarkers(scRNAsub)
all.markers = diff.wilcox %>% dplyr::select(.,gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top20 = all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(all.markers, paste0(d2,"/3.4_diff_genes_wilcox_",tr,re,".csv"), row.names = F)
write.csv(top10, paste0(d2,"/3.4_top10_diff_genes_wilcox_",tr,re,".csv"), row.names = F)
write.csv(top20, paste0(d2,"/3.4_top20_diff_genes_wilcox_",tr,re,".csv"), row.names = F)

# all.markers <- read.csv("subsample_CDDP/subcluster_CAFs/diff_genes_wilcox_0.3.csv")
# head(all.markers)
# all.markers$gene <- str_sub(all.markers$gene,3)
# write.csv(all.markers, paste0("subsample_CDDP/subcluster_CAFs/diff_genes_wilcox_0.3_2.csv"), row.names = F)


##top10基因绘制热图 ----
top10_genes <- read.csv(paste0(d2,"/3.4_top10_diff_genes_wilcox_",tr,re,".csv"))
top10_genes = CaseMatch(search = as.vector(top10_genes$gene), match = rownames(scRNAsub)) 
plot1 = DoHeatmap(scRNAsub, features = top10_genes, group.by = "seurat_clusters", group.bar = T, size = 4)
# scale_fill_gradientn(colors = c("#74add1", "#FFFFBF", "#f46d43"))#改颜色
ggsave(paste0(d2,"/3.4_top10_markers_Heatmap_",tr,re,".pdf"), plot=plot1, width=15, height=20) 
ggsave(paste0(d2,"/3.4_top10_markers_Heatmap_",tr,re,".png"), plot=plot1, width=15, height=20)

top20_genes <- read.csv(paste0(d2,"/3.4_top20_diff_genes_wilcox_",tr,re,".csv"))
top20_genes = CaseMatch(search = as.vector(top20_genes$gene), match = rownames(scRNAsub)) 
plot1 = DoHeatmap(scRNAsub, features = top20_genes, group.by = "seurat_clusters", group.bar = T, size = 4)
ggsave(paste0(d2,"/3.4_top20_markers_Heatmap_",tr,re,".pdf"), plot=plot1, width=20, height=30) 
ggsave(paste0(d2,"/3.4_top20_markers_Heatmap_",tr,re,".png"), plot=plot1, width=20, height=30)

which(duplicated(top10_genes))
top10_genes <- top10_genes[-which(duplicated(top10_genes))]
plot2 = DotPlot(scRNAsub, features = top10_genes)+
  theme(axis.text.x =element_text(angle = 90),axis.text.y =element_text(angle = 90),
        strip.text=element_text(size=0))
ggsave(paste0(d2,"/3.4_top10_markers_",tr,re,".pdf"), plot=plot2, width=49, height=10) 
ggsave(paste0(d2,"/3.4_top10_markers_",tr,re,".png"), plot=plot2, width=49, height=10)


## Step 4 基因可视化 ----
rm(list=ls())
gc()
options(stringsAsFactors = F)
d1="Cluster"
d2="Marker"
re=0.2
tr='afterharmony_'
scRNA <- readRDS(paste0("scRNA_all_STEP3.4.2_",tr,re,".rds"))

# 挑选部分感兴趣基因可视化----
# 参考PMID: 32505533，Single-cell transcriptomic architecture and intercellular crosstalk of human intrahepatic cholangiocarcinoma
celltype_marker=c(
  "PTPRC",#免疫细胞 immune'
  "EPCAM",#上皮细胞 epithelial
  "PECAM1","CDH5",#内皮细胞 endothelial
  "COL1A1","VIM",#成纤维细胞 fibroblasts
  "P2RX5",#脂肪细胞 adipocyte
  "CXCL13","GPM6A",#间皮细胞 Mesothelial cell
  "RGS5")#周皮细胞 Pericyte
pdf(file = "Marker/04.markerFeature_afterregress_tSNE_0.2.pdf", width = 20, height = 8)
FeaturePlot(object = scRNA, 
            features = celltype_marker,
            order = TRUE,
            reduction = 'tsne',
            ncol = 5,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)
dev.off()
pdf(file = "Marker/04.markerFeature_afterregress_UMAP_0.2.pdf", width = 20, height = 8)
FeaturePlot(object = scRNA, 
            features = celltype_marker,
            order = TRUE,
            ncol = 5,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)
dev.off()
pdf(file = "Marker/04.markerBubble_afterregress_0.2.pdf", width = 12, height = 4)
DotPlot(object = scRNA, features = celltype_marker)
dev.off()
pdf(file = "Marker/04.markerViolin_afterregress_0.2.pdf", width = 12, height = 4)
VlnPlot(scRNA,  features = celltype_marker ,pt.size = 0,stack=T)
dev.off()


# 散点图和气泡图 ----
# 大类鉴定 ----
immunecell.Marker=c('PTPRC','CD33','CD14','ITGAM','CD3E','CD79B','ITGAX','CD74')
pdf(file = "Marker/05.Immune_markerFeature.pdf", width = 12, height = 12)
FeaturePlot(object = scRNA, 
            features = immunecell.Marker,
            order = TRUE,
            reduction = 'tsne',
            ncol = 3,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)
dev.off()
pdf(file = "Marker/05.Immune_markerBubble.pdf", width = 12, height = 4)
DotPlot(object = scRNA, features = immunecell.Marker)
dev.off()

epi.cell.Marker=c('EPCAM')
pdf(file = "Marker/05.Epi_markerFeature.pdf", width = 4, height = 4)
FeaturePlot(object = scRNA, 
            features = epi.cell.Marker,
            order = TRUE,
            reduction = 'tsne',
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)
dev.off()
pdf(file = "Marker/05.Epi_markerBubble.pdf", width = 4, height = 4)
DotPlot(object = scRNA, features = epi.cell.Marker)
dev.off()

Fibroblasts.Marker=c("COL1A2","COL1A1",'PDGFRB',"GSN","COL6A1","VIM")
pdf(file = "Marker/05.Fibroblasts_markerFeature.pdf", width = 12, height = 8)
FeaturePlot(object = scRNA, 
            features = Fibroblasts.Marker,
            order = TRUE,
            reduction = 'tsne',
            ncol = 3,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)
dev.off()
pdf(file = "Marker/05.Fibroblasts_markerBubble.pdf", width = 12, height = 4)
DotPlot(object = scRNA, features = Fibroblasts.Marker)
dev.off()

endothelial.Marker=c("PECAM1","CDH5","VWF")
pdf(file = "Marker/06.Endothelium_markerFeature.pdf", width = 12, height = 4)
FeaturePlot(object = scRNA, 
            features = endothelial.Marker,
            order = TRUE,
            reduction = 'tsne',
            ncol = 3,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)
dev.off()
pdf(file = "Marker/05.Endothelial_markerBubble.pdf", width = 6, height = 4)
DotPlot(object = scRNA, features = endothelial.Marker)
dev.off()

OC.Marker=c('WT1', 'MUC16', 'PAX8', 'MUC1','BRCA1','BRCA2','FAS','TP53','EPCAM')
pdf(file = "Marker/05.OC_markerfeatureplot.pdf", width = 12, height = 12)
FeaturePlot(object = scRNA,
            features = OC.Marker,
            order = TRUE,
            reduction = 'tsne',
            ncol = 3,
            min.cutoff = 'q10',
            label = TRUE,
            repel = TRUE)
dev.off()
pdf(file = "Marker/05.OC_markerBubble.pdf", width = 10, height = 4)
DotPlot(object = scRNA, features = OC.Marker)
dev.off()

OCstem.Marker=c('CD44', 'PROM1', 'CD24','ALDH1A1','ALDH1A2','ALDH1A3','TP53','EPCAM')
pdf(file = "Marker/05.OC_Stemcell_markerfeatureplot.pdf", width = 16, height = 8)
FeaturePlot(object = scRNA,
            features = OCstem.Marker,
            order = TRUE,
            reduction = 'tsne',
            ncol = 4,
            min.cutoff = 'q10',
            label = TRUE,
            repel = TRUE)
dev.off()
pdf(file = "Marker/05.OC_Stemcell_markerBubble.pdf", width = 10, height = 4)
DotPlot(object = scRNA, features = OCstem.Marker)
dev.off()

# 免疫细胞鉴定 ----
macrophage.Marker=c('CD14','CD68','CD163','MERTK','MRC1')
pdf(file = "Marker/06.macrophage_markerFeature.pdf", width = 20, height = 4)
FeaturePlot(object = scRNA,
            features = macrophage.Marker,
            order  = TRUE,
            reduction = 'tsne',
            ncol = 5,
            min.cutoff = 'q10',
            label = TRUE)
dev.off()
pdf(file = "Marker/06.macrophage_markerBubble.pdf", width = 10, height = 4)
DotPlot(object = scRNA, features = macrophage.Marker)
dev.off()

Bcell.Marker=c('CD74','CD69','CD79B','LTB','BCL11A',
               'CD79A','FCRL1','CD52')
pdf(file = "Marker/06.Bcell_markerFeature.pdf", width = 16, height = 8)
FeaturePlot(object = scRNA, 
            features = Bcell.Marker,
            order = TRUE,
            reduction = 'tsne',
            ncol = 4,
            min.cutoff = 'q10', 
            label = TRUE)
dev.off()
pdf(file = "Marker/06.Bcell_markerBubble.pdf", width = 10, height = 4)
DotPlot(object = scRNA,features = Bcell.Marker)
dev.off()

Tcell.Marker=c('CD3D','CD8A','CD3E','CD4','CD2','CD8B')
pdf(file = "Marker/06.Tcell_markerFeature.pdf", width = 12, height = 8)
FeaturePlot(object = scRNA, 
            features = Tcell.Marker,
            order = TRUE,
            reduction = 'tsne',
            ncol = 3,
            min.cutoff = 'q10', 
            label = TRUE)
dev.off()
pdf(file = "Marker/06.Tcell_markerBubble.pdf", width = 12, height = 4)
DotPlot(object = scRNA, features = Tcell.Marker)
dev.off()

NK.Marker=c('NKG7','GZMB','PRF1','IFNG','KLRD1',
            'GNLY','GZMA','IL32','SH2D2A','IL2RG','GZMM')
pdf(file = "Marker/06.NK_markerBubble.pdf", width = 10, height = 4)
DotPlot(object = scRNA, features = NK.Marker)
dev.off()

pdf(file = "Marker/06.NKT_markerBubble.pdf", width = 10, height = 4)
NKT.Marker=c('CD3D','FCGR3A','PTPRC','NCAM1',
             'GZMB','GZMA','CD3E','GNLY','IL2RG','GZMM')
DotPlot(object = scRNA, features = NKT.Marker)
dev.off()

DC.Marker=c('ITGAX','CD74','CD83','CD86')
pdf(file = "Marker/06.DC_markerFeature.pdf", width = 16, height = 4)
FeaturePlot(object = scRNA, 
            features = DC.Marker,
            order = TRUE,
            reduction = 'tsne',
            ncol = 4,
            min.cutoff = 'q10', 
            label = TRUE)
dev.off()
pdf(file = "Marker/06.DC_markerBubble.pdf", width = 10, height = 4)
DotPlot(object = scRNA, features = DC.Marker)
dev.off()

mast.cell.Marker=c('TPSAB1','CD274','MGST2','CPA3','CMA1',
                   'IL1RL1','RGS1','TPSB2','MS4A2')
pdf(file = "Marker/06.mast.cell_markerBubble.pdf", width = 10, height = 4)
DotPlot(object = scRNA, features = mast.cell.Marker)
dev.off()


Neutrophil.Marker=c('G0S2','S100A8','S100A9','SORL1')
pdf(file = "Marker/06.Neutrophil_markerFeature.pdf", width = 16, height = 4)
FeaturePlot(object = scRNA, 
            features = Neutrophil.Marker,
            order = TRUE,
            reduction = 'tsne',
            ncol = 4,
            min.cutoff = 'q10', 
            label = TRUE)
dev.off()
pdf(file = "Marker/06.Neutrophil_markerBubble.pdf", width = 10, height = 4)
DotPlot(object = scRNA, features = Neutrophil.Marker)
dev.off()

Treg.Marker=c('CTLA4','TIGIT','FOXP3','IKZF2','TNFRSF18',
              'IL10','ITGB8','MAF','CNGB1','TNFRSF4')
pdf(file = "Marker/06.Treg_markerBubble.pdf", width = 10, height = 4)
DotPlot(object = scRNA, features = Treg.Marker)
dev.off()

Tcytotoxic.Marker=c('CD2','CD8A','CD5','GZMB',
                    'PRF1','IFNG','CD69','LAG3','EOMES')
pdf(file = "Marker/06.Tcytotoxic_markerBubble.pdf", width = 10, height = 4)
DotPlot(object = scRNA, features = Tcytotoxic.Marker)
dev.off()

plasma.cell.Marker=c('SPAG4','PDK1','IGLL5','TNFRSF17',
                     'SSR4','CYCS','PRDM1','EAF2','MZB1')
pdf(file = "Marker/06.plasma.cell_markerBubble.pdf", width = 10, height = 4)
DotPlot(object = scRNA, features = plasma.cell.Marker)
dev.off()

# 其他细胞鉴定 ----
Smooth_muscle_cells_Marker=c("ACTA2","MYLK","MYL9")
pdf(file = "Marker/06.Smooth_muscle_cells_markerFeature.pdf", width = 12, height = 4)
FeaturePlot(object = scRNA, 
            features = Smooth_muscle_cells_Marker,
            order = TRUE,
            reduction = 'tsne',
            ncol = 3,
            min.cutoff = 'q10', 
            label = TRUE)
dev.off()
pdf(file = "Marker/06.Smooth_muscle_cells_markerBubble.pdf", width = 10, height = 4)
DotPlot(object = scRNA,features = Smooth_muscle_cells_Marker)
dev.off()

Adipocyte_Marker=c("APOC1","P2RX5","PNPLA2")
pdf(file = "Marker/06.Adipocyte_markerFeature.pdf", width = 12, height = 4)
FeaturePlot(object = scRNA, 
            features = Adipocyte_Marker,
            order = TRUE,
            reduction = 'tsne',
            ncol = 3,
            min.cutoff = 'q10', 
            label = TRUE)
dev.off()
pdf(file = "Marker/06.Adipocyte_markerBubble.pdf", width = 10, height = 4)
DotPlot(object = scRNA,features = Adipocyte_Marker)
dev.off()

Mesothelial_Marker=c("CALB2","CXCL13","GPM6A")
pdf(file = "Marker/06.Mesothelial_markerFeature.pdf", width = 12, height = 4)
FeaturePlot(object = scRNA, 
            features = Mesothelial_Marker,
            order = TRUE,
            reduction = 'tsne',
            ncol = 3,
            min.cutoff = 'q10', 
            label = TRUE)
dev.off()
pdf(file = "Marker/06.Mesothelial_markerBubble.pdf", width = 10, height = 4)
DotPlot(object = scRNA,features = Mesothelial_Marker)
dev.off()

Pericyte_Marker=c("PDGFRB","RGS5","MCAM","DLK1")
pdf(file = "Marker/06.Pericyte_markerFeature.pdf", width = 16, height = 4)
FeaturePlot(object = scRNA, 
            features = Pericyte_Marker,
            order = TRUE,
            reduction = 'tsne',
            ncol = 4,
            min.cutoff = 'q10', 
            label = TRUE)
dev.off()
pdf(file = "Marker/06.Pericyte_markerBubble.pdf", width = 10, height = 4)
DotPlot(object = scRNA,features = Pericyte_Marker)
dev.off()

# 特定基因可视化 ----
cellmarker=c('CCR1','CCR2','CCR3','CCR4','CCR5','CCR6','CCR7','CCR8','CCR9','CCR10','ACKR1','ACKR2','ACKR4')
cellmarker=c('CXCR1','CXCR2','CXCR3','CXCR4','CXCR5','CXCR6','ACKR1','ACKR3')
cellmarker='GPR85'
pdf(file = "Marker/06.GPR85_markerFeature.pdf", width = 4, height = 4)
FeaturePlot(object = scRNA, 
            features = cellmarker,
            order = TRUE,
            reduction = 'tsne',
            min.cutoff = 'q10', 
            label = TRUE)
dev.off()
pdf(file = "Marker/06.GPR85_markerBubble.pdf", width = 4, height = 4)
DotPlot(object = scRNA,features = cellmarker)
dev.off()

#### STEP 4 注释重命名 ----
### 4.1 将注释信息添加到metadata
Macro = c(1,11)
OC = c(3,6,10)
Tcell = c(0,4,9)
NK = c(2,5)
Bcell = c(7)
Normalcell =c(12)
deadcell=c(8)

current.cluster <- c(Macro,OC,Tcell,NK,Bcell,Normalcell,deadcell)

new.cluster <- c(rep("Macrophage",length(Macro)),
                 rep("OC_cell",length(OC)),
                 rep("T_cell",length(Tcell)),
                 rep("NK_cell",length(NK)),
                 rep("B_cell",length(Bcell)),
                 rep("Normal_cell",length(Normalcell)),
                 rep("Dead_cell",length(deadcell)))
cell.id <- c( "Tcell_1",
              "Macro_1",
              "NKcell_1",
              "OC_1",
              "Tcell_2",
              "NKcell_2",
              "OC_2",
              "Bcell_1",
              "Deadcell",
              "Treg",
              "OC_3",
              "Macro_2",
              "NE_1")

scRNA@meta.data$cell_type <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster, to = new.cluster)
table(scRNA@meta.data$cell_type)
scRNA@meta.data$cell_id <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = c(0:12), to = cell.id)
table(scRNA@meta.data$cell_id)

plot1 = DimPlot(scRNA, reduction = "tsne", group.by='cell_type', label = F, label.size = 3, repel = T) 
ggsave("Cluster/tSNE_groupby_cell_type_afterregress_0.2.pdf", plot = plot1, width = 8, height = 7)
ggsave("Cluster/tSNE_groupby_cell_type_afterregress_0.2.png", plot = plot1, width = 8, height = 7)
plot2 = DimPlot(scRNA, reduction = "umap", group.by='cell_type', label = F, label.size = 3, repel = T)
ggsave("Cluster/UMAP_groupby_cell_type_afterregress_0.2.pdf", plot = plot2, width = 8, height = 7)
ggsave("Cluster/UMAP_groupby_cell_type_afterregress_0.2.png", plot = plot2, width = 8, height = 7)
pc=plot1|plot2
ggsave("Cluster/tSNE_UMAP_groupby_cell_type_afterregress_0.2.pdf", plot = pc, width = 15, height = 7)
ggsave("Cluster/tSNE_UMAP_groupby_cell_type_afterregress_0.2.png", plot = pc, width = 15, height = 7)

Idents(scRNA) <- 'cell_id'
plot3 = DimPlot(scRNA, reduction = "tsne", label = T, label.size = 3, repel = T) 
ggsave("Cluster/tSNE_renamed_afterregress_0.2.pdf", plot = plot3, width = 8, height = 7)
ggsave("Cluster/tSNE_renamed_afterregress_0.2.png", plot = plot3, width = 8, height = 7)
plot4 = DimPlot(scRNA, reduction = "umap", label = T, label.size = 3, repel = T)
ggsave("Cluster/UMAP_renamed_afterregress_0.2.pdf", plot = plot4, width = 8, height = 7)
ggsave("Cluster/UMAP_renamed_afterregress_0.2.png", plot = plot4, width = 8, height = 7)
pc2=plot1|plot2
ggsave("Cluster/tSNE_UMAP_renamed_afterregress_0.2.pdf", plot = pc2, width = 15, height = 7)
ggsave("Cluster/tSNE_UMAP_renamed_afterregress_0.2.png", plot = pc2, width = 15, height = 7)
Idents(scRNA) <- 'seurat_clusters'

saveRDS(scRNA, file = "scRNA_all_STEP4.1.rds")


# 根据细胞亚群是否跨越病人存在来区分上皮细胞的恶性与否
# 0,2,3,6,7,9 是跨越病人的聚类情况
# Cancer cells were identified in 16 of the original 19 tumor biopsy samples
# All cancer cells (n = 8,980) were re-clustered, resulting in 9 unique clusters
# The non-cancer epithelial cell clusters (n = 3) were further annotated into cell subtype

#查看每个sample里面各个细胞的变化
scRNA <- readRDS(file = "scRNA_all_STEP4.1.rds")

library(gplots)
pdf (file = "Cluster/celltype_to_samples_ballonplot_afterregress_0.2.pdf", width = 12, height = 8)
tab.1=table(scRNA@meta.data$cell_type,scRNA@meta.data$sample_name) 
balloonplot(tab.1)
dev.off()


#分析各种细胞所占百分比
tab.2=table(scRNA@meta.data$sample_name, scRNA@meta.data$cell_type) 
tab.2
statistics = apply(tab.2, 1, sum)   # 得到每个样本的观测值总和
plot(statistics)
data_percent = data.frame()   # 建立空数据框
for (n in 1:5) {
  data_percent = rbind(data_percent, tab.2[n,] / statistics[n] )  
} 
# 再来看下，每个样本总和都等于1，现在符合要求了
statistics = apply(data_percent, 1, sum)   
plot(statistics)
# 再加上样本的命名信息，方便看图，已有命名的请忽略
colnames(data_percent)=colnames(tab.2)
data_percent$names = c('Normal',"OM1",'OM2',"PS1","PS2")
data_percent 
library(reshape2)
library(ggsci)
data_plot = melt(data_percent)
colnames(data_plot) = c('sample_id','celltype','percent')
data_plot

pdf (file = "Cluster/celltypepercent_to_samples_barplot_afterregress_0.2.pdf", width = 8, height = 7)
ggplot( data_plot, aes( x = sample_id, weight = percent, fill = celltype))+
  geom_bar( position = "stack")  +theme_classic()#只留下坐标轴
dev.off()

pdf (file = "Cluster/celltypepercent_to_samples_lineplot_afterregress_0.2.pdf", width = 8, height = 7)
ggplot(data_plot, aes(x =sample_id, y = percent, group =celltype, color = celltype))+
  geom_line()+theme_classic()
dev.off()

#引用包
library(limma)
library(NMF)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(CellChat)

expFile="04.expMatirx.txt"     #表达数据文件
annFile="04.cellAnn.txt"       #细胞注释文件
geneFile="hubGenes.csv"        #基因列表文件

#读取表达数据文件，并对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#读取细胞注释文件
meta=read.table(annFile, header=T, sep="\t", check.names=F, row.names=1)

#创建cellchat对象
cellchat <- createCellChat(object = data, meta = meta, group.by = "labels")
cellchat <- setIdent(cellchat, ident.use="labels")
groupSize <- as.numeric(table(cellchat@idents))      #每个细胞类型的细胞数目

#导入配体和受体数据库
CellChatDB <- CellChatDB.human       #如果物种是小鼠改成CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

#查看配体和受体对的类型
pdf(file="COMM01.DatabaseCategory.pdf", width=7, height=5)
showDatabaseCategory(CellChatDB)
dev.off()

#对表达数据进行预处理
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)      #识别每种细胞中过表达的基因
cellchat <- identifyOverExpressedInteractions(cellchat)      #识别基因的相互作用
cellchat <- projectData(cellchat, PPI.human)  

#计算细胞通讯的概率
cellchat <- computeCommunProb(cellchat)
#过滤掉小于10个细胞的细胞通讯
cellchat <- filterCommunication(cellchat, min.cells = 10)
#输出细胞间的通讯关系
df.net=subsetCommunication(cellchat)
write.table(file="COMM02.Comm.network.xls", df.net, sep="\t", row.names=F, quote=F)

#在信号通路的水平进一步推测胞间的通讯, 推断通路水平的互作网络
cellchat <- computeCommunProbPathway(cellchat)

#对计算结果汇总整合，展示整体细胞通讯状态
cellchat <- aggregateNet(cellchat)
#输出细胞通讯的图形(互作数量的图形)
pdf(file="COMM03.cellNetworkCount.pdf", width=7, height=6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
#输出细胞通讯的图形(互作强度的图形)
pdf(file="COMM04.cellNetworkWeight.pdf", width=7, height=6)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction strength")
dev.off()

#分细胞类型展示(把单个细胞类型提取出来,观察这个细胞与其他细胞的通讯)
pdf(file="COMM05.singleCell.pdf", width=8, height=6)
weight_mat <- cellchat@net$weight
par(mfrow = c(2,3), mgp=c(0,0,0), xpd=TRUE)
for (cel in unique(cellchat@idents)){
  cir_mat <- matrix(0, nrow = nrow(weight_mat), ncol = ncol(weight_mat), dimnames = dimnames(weight_mat))
  cir_mat[cel, ] <- weight_mat[cel, ]
  netVisual_circle( cir_mat, vertex.weight= groupSize, weight.scale= T,edge.weight.max = max(weight_mat), vertex.label.cex=0.8,title.name=cel)
}
dev.off()

#绘制受体配体对的气泡图
pdf(file=paste0("COMM06.bubble.pdf"), width=8, height=5.5)
netVisual_bubble(cellchat, remove.isolate = FALSE, angle.x = 45)
dev.off()

#读取基因列表文件,找到核心基因参与的细胞通讯,获取核心基因的通路
geneRT=read.csv(geneFile, header=T, sep=",", check.names=F)
hubGenes=as.vector(geneRT[,1])
def.hub=df.net[((df.net$ligand %in% hubGenes) | (df.net$receptor %in% hubGenes)),]
write.table(file="COMM07.Comm.hubNetwork.xls", def.hub, sep="\t", row.names=F, quote=F)

#通路水平的可视化
cellchat@netP$pathways     #展示所有相关通路的名称
pathways.show="SPP1"       #选择需要展示的通路(可修改)
#通路的细胞通讯图
pdf(file=paste0("COMM08.", pathways.show , ".circle.pdf"), width=8, height=6)
circle=netVisual_aggregate(cellchat, signaling=pathways.show, layout="circle", vertex.size = groupSize)
print(circle)
dev.off()
#使用层次图展示通路的细胞通讯
pdf(file=paste0("COMM09.", pathways.show , ".hierarchy.pdf"), width=12, height=6)
hierarchy=netVisual_aggregate(cellchat, signaling=pathways.show, layout="hierarchy",  vertex.receiver=seq(1,4), vertex.size = groupSize)
print(hierarchy)
dev.off()
#细胞通讯的热图
pdf(file=paste0("COMM10.", pathways.show , ".heatmap.pdf"), width=8, height=6)
heatmap=netVisual_heatmap(cellchat, signaling=pathways.show, color.heatmap = "Reds", measure= 'weight')	
print(heatmap)
dev.off()
#细胞作用类型分析
pdf(file=paste0("COMM11.", pathways.show , ".netAnalysis.pdf"), width=6, height=5)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
netAnalysis=netAnalysis_signalingRole_network(cellchat, signaling =pathways.show, width = 8, height = 5, font.size = 12)
print(netAnalysis)
dev.off()

#查看哪些配体和受体对在通路中起作用(配体受体对的贡献程度)
pdf(file=paste0("COMM12.", pathways.show , ".contribution.pdf"), width=8, height=6)
contribution=netAnalysis_contribution(cellchat, signaling= pathways.show)
print(contribution)
dev.off()
#查看通路中互作基因的表达水平
pdf(file=paste0("COMM13.", pathways.show , ".geneExp.pdf"), width=8, height=6)
geneExp=plotGeneExpression(cellchat, signaling=pathways.show)
print(geneExp)
dev.off()

#配体受体对水平的细胞通讯展示
pairLR <- extractEnrichedLR(cellchat, signaling=pathways.show, geneLR.return=FALSE)
pdf(file=paste0("COMM14.", pathways.show , ".pairLR.pdf"), width=9, height=8)
pairCircos=netVisual_individual(cellchat, signaling=pathways.show, pairLR.use=pairLR[1] , layout="circle" )
print(pairCircos)
dev.off()
#对通路中的受体配体对进行循环, 以和弦图形式展示
for(i in 1:nrow(pairLR)){
  pdf(file=paste0("COMM15.", pairLR[i,], ".pairLR.pdf"), width=8, height=6)
  pairChord=netVisual_individual(cellchat, signaling=pathways.show, pairLR.use=pairLR[i,], layout="chord" )
  print(pairChord)
}


library(monocle)
library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
data <- as.matrix(OC_cell@assays$RNA@counts)
data <- as(data, 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = OC_cell@meta.data[,c("risk_cellTypes",
                                                            "Seurat_riskscore",
                                                            "cell_type")])
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
## 以下代码一律不得修改
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores = 4, relative_expr = TRUE)
##使用monocle选择的高变基因，不修改
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
pdf("monocle2_1.plot_ordering_genes.pdf")
plot_ordering_genes(mycds)
dev.off()
#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#排序
mycds <- orderCells(mycds)
#State轨迹分布图
plot1 <- plot_cell_trajectory(mycds, color_by = "State")
pdf("monocle2_2.cell_trajectory_state.pdf")
plot1
dev.off()
##Pseudotime轨迹图
plot2 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
pdf("monocle2_3.cell_trajectory_Pseudotime.pdf")
plot2
dev.off()
plot4 <- plot_cell_trajectory(mycds, color_by = "risk_cellTypes")
pdf("monocle2_4.cell_trajectory_risk_cellTypes.pdf")
plot4
dev.off()
saveRDS(scRNA, "scRNA_all_STEP4.1.rds")



#################################################################3
setwd("../Ucell")
# UCell score
scRNA <- AddModuleScore_UCell(scRNA, features = list(genes))
colnames(scRNA@meta.data)[32] <- "UCell_riskscore"
genes %in% rownames(scRNA)
Idents(scRNA) <- scRNA$cell_type
scRNA <- subset(scRNA, idents = c("NK_T cell","CD8_T cell","naive_T cell","B_cell",
                                  "OC_cell","Treg_cell","Monocyte","Macrophage","Fibroblast",
                                  "CD4_T cell"))
DefaultAssay(scRNA) <- "RNA"
pdf("1.dotplot.pdf", width = 6, height = 6)
DotPlot(scRNA, features = c(genes,"UCell_riskscore")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()
pdf("2.featureplot.pdf", width = 11, height = 9)
FeaturePlot(scRNA, features = c(genes,"UCell_riskscore"),
            ncol = 3, order = T)
dev.off()

pdf("2.Vinplot.pdf", width = 11, height = 9)
VinPlot(scRNA, features = c(genes,"UCell_riskscore"),
        ncol = 3, order = T)
dev.off()

scRNA$risk_cellTypes <- scRNA$cell_type
median_score <- median(scRNA$UCell_riskscore[which(scRNA$cell_type %in% c("OC_cell"))])
scRNA$risk_cellTypes <- ifelse(scRNA$UCell_riskscore > median_score,
                               "OC_cell_highrisk", "OC_cell_lowrisk")
scRNA$risk_cellTypes[which(scRNA$cell_type %in% c("NK_T cell","CD8_T cell","naive_T cell","B_cell",
                                                  "Treg_cell","Monocyte","Macrophage","Fibroblast",
                                                  "CD4_T cell"))] <- scRNA$cell_type[which(scRNA$cell_type %in% c("NK_T cell","CD8_T cell","naive_T cell","B_cell",
                                                                                                                  "Treg_cell","Monocyte","Macrophage","Fibroblast",
                                                                                                                  "CD4_T cell"))]
library(CellChat)
meta <- scRNA@meta.data # a dataframe with rownames containing cell mata data
data_input <- as.matrix(scRNA@assays$RNA@data)
#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "risk_cellTypes")
CellChatDB <- CellChatDB.human
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use 
dplyr::glimpse(CellChatDB$interaction)
##配体-受体分析
# 提取数据库支持的数据子集
cellchat <- subsetData(cellchat)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.human)
unique(cellchat@idents)
# 等待
cellchat <- computeCommunProb(cellchat)  ##默认计算方式为type = "truncatedMean",
##默认cutoff的值为20%，即表达比例在25%以下的基因会被认为是0， trim = 0.1可以调整比例阈值
cellchat <- filterCommunication(cellchat, min.cells = 10)#去掉通讯数量很少的细胞
cellchat <- computeCommunProbPathway(cellchat) ##每对配受体的预测结果存在net中，每条通路的预测结果存在netp中
df.net<- subsetCommunication(cellchat) #将细胞通讯预测结果以数据框的形式取出
DT::datatable(df.net)
write.csv(df.net,'01.df.net.csv')
cellchat <- aggregateNet(cellchat) ##计算联路数与通讯概率，可用sources.use and targets.use指定来源与去向
#dev.off()

###可视化
groupSize <- as.numeric(table(cellchat@idents))
pdf("./cellchat_1.netVisual.pdf", width = 12, height = 6)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

##############################################################
#互作数量与重要性图
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {#
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
}#循环绘出每种与其他细胞之间互作的关系

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
}#循环绘出每种与其他细胞之间互作的关系

pathways.show <- df.net$pathway_name#计算到的所有通路
# Hierarchy plot
par(mfrow = c(1,2), xpd=TRUE)
vertex.receiver = seq(1,4) 
netVisual_aggregate(cellchat, signaling =pathways.show[1],  
                    vertex.receiver = vertex.receiver,layout = 'hierarchy')





table(cellchat@meta$risk_cellTypes)
p_bubble = netVisual_bubble(cellchat,
                            targets.use = c("OC_cell_highrisk", "OC_cell_lowrisk"),
                            sources.use = c("NK_T cell","CD8_T cell","naive_T cell","B_cell",
                                            "Treg_cell","Monocyte","Macrophage","Fibroblast",
                                            "CD4_T cell"),
                            remove.isolate = FALSE) + coord_flip()
pdf("./cellchat_2_netVisual_bubble.pdf", width = 7, height = 7)
p_bubble
dev.off()

pathways <- unique(cellchat@netP[["pathways"]])
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
for (i in pathways) {
  pdf(paste0("./pathways/cellchat_3_", i, "_netVisual_aggregate.pdf"),
      width = 7, height = 7)
  p <- netVisual_aggregate(cellchat, signaling = i)
  print(p)
  dev.off()
  
  pdf(paste0("./pathways/cellchat_3_", i, "_signalingRole_network.pdf"),
      width = 7, height = 7)
  p <- netAnalysis_signalingRole_network(cellchat, signaling = i, width = 8,
                                         height = 2.5, font.size = 10)
  print(p)
  dev.off()
}

Idents(scRNA) <- scRNA$cell_type
OC_cell <- subset(scRNA, idents = "OC_cell")
pdf("3.featureplot_OC_cell.pdf", width = 11, height = 9)
FeaturePlot(OC_cell, features = c(genes,"Seurat_riskscore"),
            ncol = 3, order = T)
dev.off()


#############################################
library(monocle)
library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
data <- as.matrix(OC_cell@assays$RNA@counts)
data <- as(data, 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = OC_cell@meta.data[,c("risk_cellTypes",
                                                            "Seurat_riskscore",
                                                            "cell_type")])
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
## 以下代码一律不得修改
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores = 4, relative_expr = TRUE)
##使用monocle选择的高变基因，不修改
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
pdf("monocle2_1.plot_ordering_genes.pdf")
plot_ordering_genes(mycds)
dev.off()
#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#排序
mycds <- orderCells(mycds)
#State轨迹分布图
plot1 <- plot_cell_trajectory(mycds, color_by = "State")
pdf("monocle2_2.cell_trajectory_state.pdf")
plot1
dev.off()
##Pseudotime轨迹图
plot2 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
pdf("monocle2_3.cell_trajectory_Pseudotime.pdf")
plot2
dev.off()
plot4 <- plot_cell_trajectory(mycds, color_by = "risk_cellTypes")
pdf("monocle2_4.cell_trajectory_risk_cellTypes.pdf")
plot4
dev.off()
saveRDS(scRNA, "scRNA_all_STEP4.1.rds")

library(GSVA)
library(GSEABase)
setwd("../ssGSEA")
# UCell score
data_matrix <- scRNA@assays$RNA@data
gs.exp <- gsva(data_matrix, list(genes), method = "ssgsea")
gs.exp <- as.matrix(gs.exp)
scRNA$ssGSEA_riskscore <- t(gs.exp)[,1]
Idents(scRNA) <- scRNA$cell_type
scRNA <- subset(scRNA, idents = c("NK_T cell","CD8_T cell","naive_T cell","B_cell",
                                  "OC_cell","Treg_cell","Monocyte","Macrophage","Fibroblast",
                                  "CD4_T cell"))
DefaultAssay(scRNA) <- "RNA"
pdf("1.dotplot.pdf", width = 6, height = 6)
DotPlot(scRNA, features = c(genes,"ssGSEA_riskscore")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()
pdf("2.featureplot.pdf", width = 11, height = 9)
FeaturePlot(scRNA, features = c(genes,"ssGSEA_riskscore"),
            ncol = 3, order = T)
dev.off()

scRNA$risk_cellTypes <- scRNA$cell_type
median_score <- median(scRNA$ssGSEA_riskscore[which(scRNA$cell_type %in% c("OC_cell"))])
scRNA$risk_cellTypes <- ifelse(scRNA$ssGSEA_riskscore > median_score,
                               "OC_cell_highrisk", "OC_cell_lowrisk")
scRNA$risk_cellTypes[which(scRNA$cell_type %in% c("NK_T cell","CD8_T cell","naive_T cell","B_cell",
                                                  "Treg_cell","Monocyte","Macrophage","Fibroblast",
                                                  "CD4_T cell"))] <- scRNA$cell_type[which(scRNA$cell_type %in% c("NK_T cell","CD8_T cell","naive_T cell","B_cell",
                                                                                                                  "Treg_cell","Monocyte","Macrophage","Fibroblast",
                                                                                                                  "CD4_T cell"))]
library(CellChat)
meta <- scRNA@meta.data # a dataframe with rownames containing cell mata data
data_input <- as.matrix(scRNA@assays$RNA@data)
#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "risk_cellTypes")
CellChatDB <- CellChatDB.human
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use 
dplyr::glimpse(CellChatDB$interaction)
##配体-受体分析
# 提取数据库支持的数据子集
cellchat <- subsetData(cellchat)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.human)
unique(cellchat@idents)
# 等待
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
df.net<- subsetCommunication(cellchat)
cellchat <- aggregateNet(cellchat)
#dev.off()
groupSize <- as.numeric(table(cellchat@idents))
pdf("./cellchat_1.netVisual.pdf", width = 12, height = 6)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
table(cellchat@meta$risk_cellTypes)
p_bubble = netVisual_bubble(cellchat,
                            targets.use = c("OC_cell_highrisk", "OC_cell_lowrisk"),
                            sources.use = c("NK_T cell","CD8_T cell","naive_T cell","B_cell",
                                            "Treg_cell","Monocyte","Macrophage","Fibroblast",
                                            "CD4_T cell"),
                            remove.isolate = FALSE) + coord_flip()
pdf("./cellchat_2_netVisual_bubble.pdf", width = 7, height = 7)
p_bubble
dev.off()

pathways <- unique(cellchat@netP[["pathways"]])
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
for (i in pathways) {
  pdf(paste0("./pathways/cellchat_3_", i, "_netVisual_aggregate.pdf"),
      width = 7, height = 7)
  p <- netVisual_aggregate(cellchat, signaling = i)
  print(p)
  dev.off()
  
  pdf(paste0("./pathways/cellchat_3_", i, "_signalingRole_network.pdf"),
      width = 7, height = 7)
  p <- netAnalysis_signalingRole_network(cellchat, signaling = i, width = 8,
                                         height = 2.5, font.size = 10)
  print(p)
  dev.off()
}

Idents(scRNA) <- scRNA$cell_type
OC_cell <- subset(scRNA, idents = "OC_cell")
pdf("3.featureplot_OC_cell.pdf", width = 11, height = 9)
FeaturePlot(OC_cell, features = c(genes,"Seurat_riskscore"),
            ncol = 3, order = T)
dev.off()
library(monocle)
library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
data <- as.matrix(OC_cell@assays$RNA@counts)
data <- as(data, 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = OC_cell@meta.data[,c("risk_cellTypes",
                                                            "Seurat_riskscore",
                                                            "cell_type")])
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
## 以下代码一律不得修改
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores = 4, relative_expr = TRUE)
##使用monocle选择的高变基因，不修改
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
pdf("monocle2_1.plot_ordering_genes.pdf")
plot_ordering_genes(mycds)
dev.off()
#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#排序
mycds <- orderCells(mycds)
#State轨迹分布图
plot1 <- plot_cell_trajectory(mycds, color_by = "State")
pdf("monocle2_2.cell_trajectory_state.pdf")
plot1
dev.off()
##Pseudotime轨迹图
plot2 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
pdf("monocle2_3.cell_trajectory_Pseudotime.pdf")
plot2
dev.off()
plot4 <- plot_cell_trajectory(mycds, color_by = "risk_cellTypes")
pdf("monocle2_4.cell_trajectory_risk_cellTypes.pdf")
plot4
dev.off()
saveRDS(scRNA, "scRNA_all_STEP4.1.rds")


























