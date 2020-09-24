######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!require("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install("maftools")


library(maftools)           #引用包
setwd("D:\\biowolf\\TMBimmune\\08.maftools")      #设置工作目录
maf = read.maf(maf = 'input.maf')          #读取输入文件

#summary图
pdf(file="summary.pdf",width=7,height=6)
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

#瀑布图
pdf(file="waterfall.pdf",width=7,height=6)
oncoplot(maf = maf, top = 30, fontSize = 0.8 ,showTumorSampleBarcodes = F )
dev.off()

#相关性图
pdf(file="interaction.pdf",width=7,height=6)
somaticInteractions(maf = maf, top = 20, pvalue = c(0.05, 0.001),showCounts = FALSE, fontSize = 0.6)
dev.off()

#基因云图
pdf(file="Genecloud.pdf",width=7,height=6)
geneCloud(input = maf, minMut = 30)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
