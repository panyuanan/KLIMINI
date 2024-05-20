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

> setwd("C:/Users/win/Desktop/数据/KLIMILI")
> library(readxl)
> chemical <- read_excel("chemical.xlsx")
> View(chemical)

chemical$Treatment<-factor(chemical$Treatment,levels = c("C","PS","PV_1","PV_3","PV_10","PI_1","PI_3","PI_10","NS","NE"))
col<-c("beige","white","darkolivegreen1","darkolivegreen3","darkolivegreen","khaki1","goldenrod2","goldenrod4",'lightblue', 'cadetblue')

N<- ggbarplot(chemical, x = "Treatment", y = "N", fill = "Treatment",
                   add = "mean_se", legend="none",width=0.7,
                   palette = col) +  theme(plot.title = element_text(hjust = 0.5))+ylab("")+xlab("")+ggtitle('N %')

N<-N + stat_compare_means(method = "kruskal.test", #统计方法
                                  aes(label =paste0("p=",after_stat(p.format))), #显示方式
                                  label.x = 0.8, label.y = 0.23,#位置
                                  size = 5)+ylim(0,0.23) #大小
N

C<- ggbarplot(chemical, x = "Treatment", y = "C", fill = "Treatment",
              add = "mean_se", legend="none",width=0.7,
              palette = col) +  theme(plot.title = element_text(hjust = 0.5))+ylab("")+xlab("")+ggtitle("C %")

C<-C + stat_compare_means(method = "kruskal.test", #统计方法
                          aes(label =paste0("p=",after_stat(p.format))), #显示方式
                          label.x = 0.8, label.y = 2.1,#位置
                          size = 5)+ylim(0,2.1) #大小
C

W<- ggbarplot(chemical, x = "Treatment", y = "Water Content", fill = "Treatment",
              add = "mean_se", legend="none",width=0.7,
              palette = col) +  theme(plot.title = element_text(hjust = 0.5))+ylab("")+xlab("")+ggtitle("Water Content %")

W<-W + stat_compare_means(method = "kruskal.test", #统计方法
                          aes(label =paste0("p=",after_stat(p.format))), #显示方式
                          label.x = 0.8, label.y = 0.25,#位置
                          size = 5)+ylim(0,0.257) #大小
W


O<- ggbarplot(chemical, x = "Treatment", y = "Organic Matter", fill = "Treatment",
              add = "mean_se", legend="none",width=0.7,
              palette = col) +  theme(plot.title = element_text(hjust = 0.5))+ylab("")+xlab("")+ggtitle("Organic Matter %")

O<-O + stat_compare_means(method = "kruskal.test", #统计方法
                          aes(label =paste0("p=",after_stat(p.format))), #显示方式
                          label.x = 0.8, label.y = 5.35,#位置
                          size = 5)+ylim(0,5.5) #大小
O

Che<-grid.arrange(C,N,W,O,ncol=2)

stat.test <- chemical %>%
  pairwise_wilcox_test (C ~ Treatment)
stat.test <- stat.test %>% add_y_position()

C + stat_pvalue_manual(stat.test,label = "p.adj",tip.length = 0.01)

C + stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0.01,                       
                       hide.ns = T)

group_A <- subset(chemical, Treatment == "PS")$"Organic Matter"
group_B<-subset(chemical, Treatment == "C")$"Organic Matter"
t.test(group_A,group_B)
