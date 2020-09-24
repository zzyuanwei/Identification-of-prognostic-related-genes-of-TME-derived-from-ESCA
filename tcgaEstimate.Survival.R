###Video source: http://study.163.com/provider/1026136977/index.htm?share=2&shareId=1026136977
######Video source: http://www.biowolf.cn/shop/
######生信自学网: http://www.biowolf.cn/
######合作邮箱：2749657388@qq.com
######答疑微信: 18520221056

#install.packages("survival")

setwd("C:\\Users\\lexb4\\Desktop\\tcgaEstimate\\23.geneSurvival")         #工作目录（需修改）

library(survival)
rt=read.table("clinicalExp.txt",header=T,sep="\t",check.names=F)      #读取文件
rt$futime=rt$futime/365                                               #如果以月为单位，除以30；以年为单位，除以365
outTab=data.frame()

for(gene in colnames(rt[,4:ncol(rt)])){
  a=rt[,gene]<=median(rt[,gene])
  diff=survdiff(Surv(futime, fustat) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  outTab=rbind(outTab,cbind(gene=gene,pvalue=pValue))

  fit <- survfit(Surv(futime, fustat) ~ a, data = rt)
  summary(fit)

  if(pValue<0.05){
    if(pValue<0.001){
      pValue="<0.001"
    }else{
      pValue=round(pValue,3)
      pValue=paste0("=",pValue)
    }
	  pdf(file=paste(gene,".survival.pdf",sep=""),
	      width=6,
	      height=6)
	  plot(fit, 
	     lwd=2,
	     col=c("red","blue"),
	     xlab="Time (year)",
	     #mark.time=T,
	     ylab="Survival rate",
	     main=paste(gene,"(p", pValue ,")",sep=""） )
	  legend("topright", 
	       c("High","Low"), 
	       lwd=2, 
	       col=c("red","blue"))
	  dev.off()
	}
}
write.table(outTab,file="survival.xls",sep="\t",row.names=F,quote=F)

###Video source: http://study.163.com/provider/1026136977/index.htm?share=2&shareId=1026136977
######Video source: http://www.biowolf.cn/shop/
######生信自学网: http://www.biowolf.cn/
######合作邮箱：2749657388@qq.com
######答疑微信: 18520221056
library(survival)

library(survminer)
library(ggplot2)
library(ggpubr)
group <- ifelse(rt$TRAV16>median(rt$TRAV16),'high','low')
sfit <- survfit(Surv(futime, fustat)~group, data=rt)

sfit

summary(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE,
           xlab="Time (year)",
           ylab="Survival rate",
           legend.title = "TRAV16", # 设置图例标题
           legend.labs = c("high", "low"),
           pval.size =6  ,
           font.title= 8
        )
