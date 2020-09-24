######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056

#if (!require("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install("maftools")


library(maftools)           #���ð�
setwd("D:\\biowolf\\TMBimmune\\08.maftools")      #���ù���Ŀ¼
maf = read.maf(maf = 'input.maf')          #��ȡ�����ļ�

#summaryͼ
pdf(file="summary.pdf",width=7,height=6)
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

#�ٲ�ͼ
pdf(file="waterfall.pdf",width=7,height=6)
oncoplot(maf = maf, top = 30, fontSize = 0.8 ,showTumorSampleBarcodes = F )
dev.off()

#�����ͼ
pdf(file="interaction.pdf",width=7,height=6)
somaticInteractions(maf = maf, top = 20, pvalue = c(0.05, 0.001),showCounts = FALSE, fontSize = 0.6)
dev.off()

#������ͼ
pdf(file="Genecloud.pdf",width=7,height=6)
geneCloud(input = maf, minMut = 30)
dev.off()


######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056