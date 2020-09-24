###Video source: http://study.163.com/provider/1026136977/index.htm?share=2&shareId=1026136977
######Video source: http://www.biowolf.cn/shop/
######生信自学网: http://www.biowolf.cn/
######合作邮箱：2749657388@qq.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma", version = "3.8")

library("limma")

setwd("C:\\Users\\26024\\Desktop\\m6a\\WHJ\\12.StromalDiff")     #设置工作目录
inputFile="symbol.txt"                                           #输入文件
scoreFile="StromalScore.txt"                                           #score文件
fdrFilter=0.05                                                   #fdr临界值
logFCfilter=1                                                    #logFC临界值

#读取score文件
score=read.table(scoreFile,sep="\t",header=T,check.names=F)
#score=score[order(score[,2]),]
med=median(score[,2])
conTab=score[score[,2]<med,]
treatTab=score[score[,2]>=med,]
con=as.vector(conTab[,1])
treat=as.vector(treatTab[,1])
conNum=length(con)
treatNum=length(treat)

#读取输入文件
outTab=data.frame()
grade=c(rep(1,conNum),rep(2,treatNum))
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.1,]
data=cbind(data[,con],data[,treat])

#差异分析
for(i in row.names(data)){
  geneName=unlist(strsplit(i,"\\|",))[1]
  geneName=gsub("\\/", "_", geneName)
  rt=rbind(expression=data[i,],grade=grade)
  rt=as.matrix(t(rt))
  wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)  
  pvalue=wilcoxTest$p.value
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
	if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
		  outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
	 }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)

#输出所有基因的差异情况
write.table(outTab,file="all.xls",sep="\t",row.names=F,quote=F)

#输出差异表格
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="StromalDiff.xls",sep="\t",row.names=F,quote=F)
up=outTab[( as.numeric(as.vector(outTab$logFC))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(up,file="StromalUp.txt",sep="\t",row.names=F,quote=F)
down=outTab[( as.numeric(as.vector(outTab$logFC))< -logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(down,file="StromalDown.txt",sep="\t",row.names=F,quote=F)


#绘制热图需要的文件
heatmap=rbind(ID=colnames(data[as.vector(outDiff[,1]),]),data[as.vector(outDiff[,1]),])
write.table(heatmap,file="heatmap.txt",sep="\t",col.names=F,quote=F)
Type=c(rep("low",conNum),rep("high",treatNum))    #修改正常和癌症样品数目
names(Type)=colnames(data)
Type=as.data.frame(Type)
Type=cbind(ID=rownames(Type),Type)
write.table(Type,file="type.txt",sep="\t",row.names=F,quote=F)

###Video source: http://study.163.com/provider/1026136977/index.htm?share=2&shareId=1026136977
######Video source: http://www.biowolf.cn/shop/
######生信自学网: http://www.biowolf.cn/
######合作邮箱：2749657388@qq.com
######答疑微信: 18520221056