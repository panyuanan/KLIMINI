
#调用R包
library(tidyverse)
library(microeco)
library(magrittr)
library(patchwork)
#读取数据
otutab <- read.table(file="otutab.txt",sep="\t",header=T,check.names=FALSE ,row.names=1)
metadata <- read.table("metadata.tsv", sep='\t', header=T,check.names=FALSE )
rownames(metadata)<-metadata$SampleID

taxonomy<-read.table("taxonomy.tsv",sep='\t', header=T,check.names=FALSE)
rownames(taxonomy)<-taxonomy$OTUID
taxonomy %<>% tidy_taxonomy


df <- microtable$new(sample_table = metadata,
                     otu_table = otutab,
                     tax_table = taxonomy,
                     auto_tidy = F)
df$sample_table$Treatment %<>% factor(., levels = c("C","PS","PV_1","PV_3","PV_10","PI_1","PI_3","PI_10","NS","NE"))
df
#去除不属于非古菌和细菌的OTU
df$tax_table %<>% subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")
df$tax_table %<>% .[grepl("Bacteria|Archaea", .$Kingdom), ]
df
# 
# ##去除“线粒体”和“叶绿体”污染
# df$filter_pollution(taxa = c("mitochondria", "chloroplast"))
# df
# 
#统一各数据的样本和OTU信息
df$tidy_dataset()
df
# 
# ##检查序列号
df$sample_sums() %>% range
# 
#重采样以减少测序深度对多样性测量的影响，使每个样本的序列号相等。
df$rarefy_samples(sample.size = 1890)
df$sample_sums() %>% range

#开始LEfse分析
lefse <- trans_diff$new(dataset = df, 
                        method = "lefse", 
                        group = "Treatment", 
                        alpha = 0.01, 
                        lefse_subgroup = NULL)
# 查看分析结果
head(lefse$res_diff)


# 绘制前30个具有最高LDA（log10）的分类单元的差异特征柱状图
lefse$plot_diff_bar(use_number = 1:30, 
                    width = 0.8, 
                    )

#fungiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
otutab_f <- read.table(file="otutab.txt",sep="\t",header=T,check.names=FALSE ,row.names=1)
metadata_f <- read.table("metadata.tsv", sep='\t', header=T,check.names=FALSE )
rownames(metadata_f)<-metadata$SampleID

taxonomy_f<-read.table("taxonomy.txt",sep='\t', header=T,check.names=FALSE)
rownames(taxonomy_f)<-taxonomy_f$`Feature ID`
taxonomy_f %<>% tidy_taxonomy


df_f <- microtable$new(sample_table = metadata_f,
                       otu_table = otutab_f,
                       tax_table = taxonomy_f,
                       auto_tidy = F)
df_f$sample_table$Treatment %<>% factor(., levels = c("C","PS","PV_1","PV_3","PV_10","PI_1","PI_3","PI_10","NS","NE"))
df_f
#去除不属于fungi的OTU
df_f$tax_table %<>% subset(Kingdom == "k__Fungi")

df_f
# 

#统一各数据的样本和OTU信息
df_f$tidy_dataset()
df_f
# 
# ##检查序列号
df_f$sample_sums() %>% range
# 
#重采样以减少测序深度对多样性测量的影响，使每个样本的序列号相等。
df_f$rarefy_samples(sample.size = 3256)
df_f$sample_sums() %>% range
#lefse
lefse_f <- trans_diff$new(dataset = df_f, 
                        method = "lefse", 
                        group = "Treatment", 
                        alpha = 0.01, 
                        lefse_subgroup = NULL)
# 查看分析结果
head(lefse_f$res_diff)


# 绘制前30个具有最高LDA（log10）的分类单元的差异特征柱状图
lefse$plot_diff_bar(use_number = 1:30, 
                    width = 0.8, 
)
