###Video source: http://study.163.com/provider/1026136977/index.htm?share=2&shareId=1026136977
######Video source: http://www.biowolf.cn/shop/
######生信自学网: http://www.biowolf.cn/
######合作邮箱：2749657388@qq.com
######答疑微信: 18520221056
inputFile="scoersClinical4.txt"                                         #输入文件
setwd("C:\\Users\\26024\\Desktop\\m6a\\WHJ\\11.clinicalCor")         #修改工作目录
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
clinical="T"                                                       #定义临床类型

xlabel=vector()
tab1=table(rt[,clinical])
labelNum=length(tab1)
dotCol=c(2,3)
if(labelNum==3){
	dotCol=c(2,3,4)
}
if(labelNum==4){
	dotCol=c(2,3,4,5)
}
if(labelNum>4){
	dotCol=rainbow(labelNum)
}
for(i in 1:labelNum){
  xlabel=c(xlabel,names(tab1[i]) )
}

outTab=data.frame()

for(i in colnames(rt[,3:ncol(rt)])){
  rt1=rbind(expression=rt[,i],clinical=rt[,clinical])
  rt1=as.matrix(t(rt1))
  if(labelNum==2){
    wilcoxTest<-wilcox.test(expression ~ clinical, data=rt1)
  }else{
    wilcoxTest<-kruskal.test(expression ~ clinical, data = rt1)}
  pValue=wilcoxTest$p.value
  outTab=rbind(outTab,cbind(gene=i,pVal=pValue))
  pval=0
  if(pValue<0.001){
  pval=signif(pValue,4)
    pval=format(pval, scientific = TRUE)
     }else{
    pval=round(pValue,3)
  }
  
   b = boxplot(expression ~ clinical, data = rt1,outline = FALSE, plot=F) 
   yMin=min(b$stats)
   yMax = max(b$stats/5+b$stats)
   ySeg = max(b$stats/10+b$stats)
   ySeg2 = max(b$stats/12+b$stats)
   n = ncol(b$stats)

   outPdf=paste0(i,".",clinical,".pdf")
   pdf(file=outPdf,
       width=9,
       height=6,)
   par(mar = c(4,7,3,3))
   boxplot(expression ~ clinical, data = rt1,names=xlabel,
     ylab = i,col=dotCol,
     cex.main=1.6, cex.lab=1.4, cex.axis=1.3,ylim=c(yMin,yMax),outline = FALSE)
   segments(1,ySeg, n,ySeg);
   segments(1,ySeg, 1,ySeg2)
   segments(n,ySeg, n,ySeg2)
   text((1+n)/2,ySeg,labels=paste0("p=",pval),cex=1.5,pos=3)
   dev.off()

}
write.table(outTab,file=paste0(clinical,".xls"),sep="\t",row.names=F,quote=F)

###Video source: http://study.163.com/provider/1026136977/index.htm?share=2&shareId=1026136977
######Video source: http://www.biowolf.cn/shop/
######生信自学网: http://www.biowolf.cn/
######合作邮箱：2749657388@qq.com
######答疑微信: 18520221056