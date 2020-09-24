###Video source: http://study.163.com/provider/1026136977/index.htm?share=2&shareId=1026136977
######Video source: http://www.biowolf.cn/shop/
######生信自学网: http://www.biowolf.cn/
######合作邮箱：2749657388@qq.com
######答疑微信: 18520221056

#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DOSE", version = "3.8")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("clusterProfiler", version = "3.8")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("enrichplot", version = "3.8")


library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

setwd("C:\\Users\\lexb4\\Desktop\\tcgaEstimate\\19.KEGG")
rt=read.table("id.txt",sep="\t",header=T,check.names=F)
rt=rt[is.na(rt[,"entrezID"])==F,]
gene=rt$entrezID

#kegg富集分析
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)

#柱状图
tiff(file="barplot.tiff",width = 20,height = 12,units ="cm",compression="lzw",bg="white",res=600)
barplot(kk, drop = TRUE, showCategory = 30)
dev.off()

#气泡图
tiff(file="dotplot.tiff",width = 20,height = 12,units ="cm",compression="lzw",bg="white",res=600)
dotplot(kk, showCategory = 30)
dev.off()

###Video source: http://study.163.com/provider/1026136977/index.htm?share=2&shareId=1026136977
######Video source: http://www.biowolf.cn/shop/
######生信自学网: http://www.biowolf.cn/
######合作邮箱：2749657388@qq.com
######答疑微信: 18520221056