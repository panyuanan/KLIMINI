library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggprism)
library(vegan)
library(picante)
library(dplyr)
library(RColorBrewer)
library(rstatix)
library(gridExtra)

otutab <- read.table(file="otutab.txt",sep="\t",header=T,check.names=FALSE ,row.names=1)
metadata <- read.table("metadata.tsv", sep='\t', header=T,check.names=FALSE )
##导入数据，所需是数据行名为样本名、列名为OTUxxx的数据表
df <- otutab
#使用vegan包计算多样性指数
Shannon <- diversity(df, index = "shannon", MARGIN = 2, base = exp(1))
Simpson <- diversity(df, index = "simpson", MARGIN = 2, base =  exp(1))
Richness <- specnumber(df, MARGIN = 2)#spe.rich =sobs
###将以上多样性指数统计成表格
index <- as.data.frame(cbind(Shannon, Simpson, Richness))
tdf <- t(df)#转置表格
tdf<-ceiling(as.data.frame(t(df)))
#计算obs，chao，ace指数
obs_chao_ace <- t(estimateR(tdf))
obs_chao_ace <- obs_chao_ace[rownames(index),]#统一行名
#将obs，chao，ace指数与前面指数计算结果进行合并
index$Chao <- obs_chao_ace[,2]
index$Ace <- obs_chao_ace[,4]
index$obs <- obs_chao_ace[,1]
#计算Pielou及覆盖度
index$Pielou <- Shannon / log(Richness, 2)
index$Goods_coverage <- 1 - colSums(df ==1) / colSums(df)

index$samples <- rownames(index)#将样本名写到文件中
#读入分组文件
groups <- metadata
colnames(groups)[1:2] <- c('samples','treatment')#改列名
#合并分组信息与多样性指数
df2 <- merge(index,groups,by = 'samples')
df2<-df2[-12,]
row.names(df2)<-1:29
df2$treatment<-factor(df2$treatment,levels = c("C","PS","PV_1","PV_3","PV_10","PI_1","PI_3","PI_10","NS","NE"))
pig1<-df2[c(10:11,18:23),]
pig3<-df2[c(12:14,18:20,24:26),]
pig10<-df2[c(15:20,27:29),]
chemi<-df2[c(4:9),]

#Shannon Pig1
P1col = c( "white","darkolivegreen1","khaki1")

P1 <- ggbarplot(pig1, x = "treatment", y = "Shannon", fill = "treatment",
               add = "mean_se", legend="none",width=0.7,title="Dosage 1",
               palette = P1col) +  theme(plot.title = element_text(hjust = 0.5))+labs(y="Bacteria Shannon Diversity",x="")
              
P1<-P1 + stat_compare_means(method = "kruskal.test", #统计方法
                     aes(label =paste0("p=",after_stat(p.format))), #显示方式
                           label.x = 0.8, label.y = 6,#位置
                           size = 5)+ylim(0,6.5) #大小
P1
#Shannon Pig3
P3col = c( "white","darkolivegreen3","goldenrod2")

P3 <- ggbarplot(pig3, x = "treatment", y = "Shannon", fill = "treatment",
                add = "mean_se", legend="none",width=0.7,title="Dosage 3",
                palette = P3col) +  theme(plot.title = element_text(hjust = 0.5))+ylab("")+xlab("")

P3<-P3 + stat_compare_means(method = "kruskal.test", #统计方法
                            aes(label =paste0("p=",after_stat(p.format))), #显示方式
                            label.x = 0.8, label.y = 6,#位置
                            size = 5)+ylim(0,6.5) #大小
P3

#Shannon Pig10
P10col = c( "white","darkolivegreen","goldenrod4")

P10 <- ggbarplot(pig10, x = "treatment", y = "Shannon", fill = "treatment",
                add = "mean_se", legend="none",width=0.7,title="Dosage 10",
                palette = P10col) +  theme(plot.title = element_text(hjust = 0.5))+ylab("")+xlab("")

P10<-P10 + stat_compare_means(method = "kruskal.test", #统计方法
                            aes(label =paste0("p=",after_stat(p.format))), #显示方式
                            label.x = 0.8, label.y = 6,#位置
                            size = 5)+ylim(0,6.5) #大小
P10
#Shannon Chemical
Checol <- c('lightblue', 'cadetblue')

Chemi <- ggbarplot(chemi, x = "treatment", y = "Shannon", fill = "treatment",
                 add = "mean_se", legend="none",width=0.5,title="Dosage 1",
                 palette = Checol) +  theme(plot.title = element_text(hjust = 0.5))+ylab("")+xlab("")

Chemi<-Chemi + stat_compare_means(method = "kruskal.test", #统计方法
                              aes(label =paste0("p=",after_stat(p.format))), #显示方式
                              label.x = 0.8, label.y = 6,#位置
                              size = 5)+ylim(0,6.5) #大小
Chemi

#pairwise
stat.test <- pig1 %>%
  pairwise_wilcox_test (Shannon ~ group)
stat.test <- stat.test %>% add_y_position()

P + stat_pvalue_manual(stat.test,label = "p.adj",tip.length = 0.01)

P + stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0.01,                       
                       hide.ns = T)
#Merge
Bacteria_Shannow <- grid.arrange(P1, P3, P10,Chemi, ncol = 4)

#Fungi

otutab_f <- read.table(file="otutab.txt",sep="\t",header=T,check.names=FALSE ,row.names=1)
metadata_f <- read.table("metadata.tsv", sep='\t', header=T,check.names=FALSE )
##导入数据，所需是数据行名为样本名、列名为OTUxxx的数据表

#使用vegan包计算多样性指数
Shannon <- diversity(otutab_f, index = "shannon", MARGIN = 2, base = exp(1))
Simpson <- diversity(otutab_f, index = "simpson", MARGIN = 2, base =  exp(1))
Richness <- specnumber(otutab_f, MARGIN = 2)#spe.rich =sobs
###将以上多样性指数统计成表格
index <- as.data.frame(cbind(Shannon, Simpson, Richness))
t_f <- t(otutab_f)#转置表格
t_f<-ceiling(as.data.frame(t(otutab_f)))
#计算obs，chao，ace指数
obs_chao_ace <- t(estimateR(t_f))
obs_chao_ace <- obs_chao_ace[rownames(index),]#统一行名
#将obs，chao，ace指数与前面指数计算结果进行合并
index$Chao <- obs_chao_ace[,2]
index$Ace <- obs_chao_ace[,4]
index$obs <- obs_chao_ace[,1]
#计算Pielou及覆盖度
index$Pielou <- Shannon / log(Richness, 2)
index$Goods_coverage <- 1 - colSums(otutab_f ==1) / colSums(otutab_f)

index$samples <- rownames(index)#将样本名写到文件中
#读入分组文件
groups_f <- metadata_f
colnames(groups_f)[1:2] <- c('samples','treatment')#改列名
#合并分组信息与多样性指数
df_f <- merge(index,groups,by = 'samples')


df_f$Treatment<-factor(df_f$Treatment,levels = c("C","PS","PV_1","PV_3","PV_10","PI_1","PI_3","PI_10","NS","NE"))
pig1_f<-df_f[c(9:11,18:23),]
pig3_f<-df_f[c(12:14,18:20,24:26),]
pig10_f<-df_f[c(15:20,27:29),]
chemi_f<-df_f[c(3:8),]

#Shannon Pig1
P1col = c( "white","darkolivegreen1","khaki1")

P1_f <- ggbarplot(pig1_f, x = "Treatment", y = "Shannon", fill = "Treatment",
                add = "mean_se", legend="none",width=0.7,
                palette = P1col) +  theme(plot.title = element_text(hjust = 0.5))+labs(y="Fungi Shannon Diversity",x="")

P1_f<-P1_f + stat_compare_means(method = "kruskal.test", #统计方法
                            aes(label =paste0("p=",after_stat(p.format))), #显示方式
                            label.x = 0.8, label.y = 6,#位置
                            size = 5)+ylim(0,6.5) #大小
P1_f
#Shannon Pig3
P3col = c( "white","darkolivegreen3","goldenrod2")

P3_f <- ggbarplot(pig3_f, x = "Treatment", y = "Shannon", fill = "Treatment",
                add = "mean_se", legend="none",width=0.7,
                palette = P3col) +  theme(plot.title = element_text(hjust = 0.5))+ylab("")+xlab("")

P3_f<-P3_f + stat_compare_means(method = "kruskal.test", #统计方法
                            aes(label =paste0("p=",after_stat(p.format))), #显示方式
                            label.x = 0.8, label.y = 6,#位置
                            size = 5)+ylim(0,6.5) #大小
P3_f

#Shannon Pig10
P10col = c( "white","darkolivegreen","goldenrod4")

P10_f <- ggbarplot(pig10_f, x = "Treatment", y = "Shannon", fill = "Treatment",
                 add = "mean_se", legend="none",width=0.7,
                 palette = P10col) +  theme(plot.title = element_text(hjust = 0.5))+ylab("")+xlab("")

P10_f<-P10_f + stat_compare_means(method = "kruskal.test", #统计方法
                              aes(label =paste0("p=",after_stat(p.format))), #显示方式
                              label.x = 0.8, label.y = 6,#位置
                              size = 5)+ylim(0,6.5) #大小
P10_f
#Shannon Chemical
Checol <- c('lightblue', 'cadetblue')

che_f <- ggbarplot(chemi_f, x = "Treatment", y = "Shannon", fill = "Treatment",
                   add = "mean_se", legend="none",width=0.5,
                   palette = Checol) +  theme(plot.title = element_text(hjust = 0.5))+ylab("")+xlab("")

che_f<-che_f + stat_compare_means(method = "kruskal.test", #统计方法
                                  aes(label =paste0("p=",after_stat(p.format))), #显示方式
                                  label.x = 0.8, label.y = 6,#位置
                                  size = 5)+ylim(0,6.5) #大小
che_f

#pairwise
stat.test <- pig1 %>%
  pairwise_wilcox_test (Shannon ~ group)
stat.test <- stat.test %>% add_y_position()

P + stat_pvalue_manual(stat.test,label = "p.adj",tip.length = 0.01)

P + stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0.01,                       
                       hide.ns = T)
#Merge
Shannon_p <- grid.arrange(P1, P3, P10,Chemi, P1_f,P3_f,P10_f,che_f,ncol = 4)
Shannon_p
