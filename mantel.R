#mantel
devtools::install_github("Hy4m/linkET", force = TRUE)
library(linkET)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(vegan)

fungi<-t(df_f$otu_table)%>%as.data.frame()
bac<-t(df$otu_table)%>%as.data.frame()
bac<-bac[-2,]
collem<- read_excel("Yuan Pan Nov 2023 KLIMINI.xlsx", sheet = "soilbiostore-project-data-2 (2)")
rownames(collem)<-collem$`Sample Name`
collem<-collem[,-1]

#chemi
chemi<-chemical[c(-2,-18),-6]
rownames(chemi)<-chemi$SampleID
chemi<-chemi[,-1]

#pig slurry
f_p<-fungi[3:5,]
b_p<-bac[3:5,]
c_p<-collem[3:5,]
chemi_p<-chemi[3:5,]
spe_p<-cbind(f_p,b_p,c_p)
# 土壤理化性质之间的相关性分析
mdata_p <- correlate(chemi_p)
# 把varechem拆分为2个类别的矩阵，分别与spec的每列进行mantel_test
mantel_p <- mantel_test(spe_p, chemi_p,  
                      spec_select = list(Bacteria = 1:2728,Fungi= 2729:4216,Collembola=4217:4243))
# 把连续型变量划分为区间，转换为因子型变量
mantel_p <- mantel_p %>% mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                                     labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
                            pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                                     labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))


PS<-qcorrplot(mdata_p,
          type = "lower", # 热图展示下半部分
          diag = FALSE) + # 不展示对角线
  geom_square() +
  geom_couple(data = mantel_p, 
              aes(colour = pd, # 网络线的颜色映射mantel_test结果的相关性
                  size = rd), # 网络线的粗细映射mantel_test结果的显著性
              curvature = nice_curvature()) +ggtitle("PS")+ # 基于起点或终点绘制最合适的曲线 
  scale_fill_gradientn(colours =brewer.pal(11, "RdYlBu")) + # 设置热图填充颜色,也可选"RdYlBu""RdYlGn"
  scale_size_manual(values = c(0.5, 1, 2)) + theme(legend.position = "none")+ # 设置网络线粗细范围
  scale_colour_manual(values =c('#CFCECC','#2166ac','#b2182b')) + # 设置网络线颜色
  guides(size = guide_legend(title = "Mantel's r",  # 调整图例次序、标题、样式
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))
PS
#pig1
f_p1<-fungi[c(6:8,15:16),]
b_p1<-bac[c(6:8,15:16),]
c_p1<-collem[c(6:8,15:16),]
chemi_p1<-chemi[c(6:8,15:16),]
spe_p1<-cbind(f_p1,b_p1,c_p1)
mdata_1 <- correlate(chemi_p1)
# 把varechem拆分为2个类别的矩阵，分别与spec的每列进行mantel_test
mantel_p1 <- mantel_test(spe_p1, chemi_p1,  
                        spec_select = list(Bacteria = 1:2728,Fungi= 2729:4216,Collembola=4217:4243))
# 把连续型变量划分为区间，转换为因子型变量
mantel_p1 <- mantel_p1 %>% mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                                         labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
                                pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                                         labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))


PS1<-qcorrplot(mdata_1,
              type = "lower", # 热图展示下半部分
              diag = FALSE) + # 不展示对角线
  geom_square() +
  geom_couple(data = mantel_p1, 
              aes(colour = pd, # 网络线的颜色映射mantel_test结果的显著性
                  size = rd), # 网络线的粗细映射mantel_test结果的相关性
              curvature = nice_curvature()) +ggtitle("PS_1")+ # 基于起点或终点绘制最合适的曲线 
  scale_fill_gradientn(colours =brewer.pal(11, "RdYlBu")) + # 设置热图填充颜色,也可选"RdYlBu""RdYlGn"
  scale_size_manual(values = c(0.5, 1, 2)) + theme(legend.position = "none")+ # 设置网络线粗细范围
  scale_colour_manual(values =c('#CFCECC','#2166ac','#b2182b')) + # 设置网络线颜色
  guides(size = guide_legend(title = "Mantel's r",  # 调整图例次序、标题、样式
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))

PS1
#pig3
f_p3<-fungi[c(9:11,17:19),]
b_p3<-bac[c(9:11,17:19),]
c_p3<-collem[c(9:11,17:19),]
chemi_p3<-chemi[c(9:11,17:19),]
spe_p3<-cbind(f_p3,b_p3,c_p3)
mdata <- correlate(chemi_p3)
# 把varechem拆分为2个类别的矩阵，分别与spec的每列进行mantel_test
mantel_p3 <- mantel_test(spe_p3, chemi_p3,  
                         spec_select = list(Bacteria = 1:2728,Fungi= 2729:4216,Collembola=4217:4243))
# 把连续型变量划分为区间，转换为因子型变量
mantel_p3 <- mantel_p3 %>% mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                                           labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
                                  pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                                           labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))


PS3<-qcorrplot(mdata,
               type = "lower", # 热图展示下半部分
               diag = FALSE) + # 不展示对角线
  geom_square() +
  geom_couple(data = mantel_p3, 
              aes(colour = pd, # 网络线的颜色映射mantel_test结果的显著性
                  size = rd), # 网络线的粗细映射mantel_test结果的相关性
              curvature = nice_curvature()) +ggtitle("PS_3")+ # 基于起点或终点绘制最合适的曲线 
  scale_fill_gradientn(colours =brewer.pal(11, "RdYlBu")) + # 设置热图填充颜色,也可选"RdYlBu""RdYlGn"
  scale_size_manual(values = c(0.5, 1, 2)) + # 设置网络线粗细范围
  scale_colour_manual(values =c('#b2182b','#CFCECC','#2166ac')) + # 设置网络线颜色
  guides(size = guide_legend(title = "Mantel's r",  # 调整图例次序、标题、样式
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))

PS3




#pig10+legend
f_p10<-fungi[c(12:14,20:22),]
b_p10<-bac[c(12:14,20:22),]
c_p10<-collem[c(12:14,20:22),]
chemi_p10<-chemi[c(12:14,20:22),]
spe_p10<-cbind(f_p10,b_p10,c_p10)
mdata_p10 <- correlate(chemi_p10)
# 把varechem拆分为2个类别的矩阵，分别与spec的每列进行mantel_test
mantel_p10 <- mantel_test(spe_p10, chemi_p10,  
                         spec_select = list(Bacteria = 1:2728,Fungi= 2729:4216,Collembola=4217:4243))
# 把连续型变量划分为区间，转换为因子型变量
mantel_p10 <- mantel_p10 %>% mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                                           labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
                                  pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                                           labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))


PS10<-qcorrplot(mdata_p10,
               type = "lower", # 热图展示下半部分
               diag = FALSE) + # 不展示对角线
  geom_square() +
  geom_couple(data = mantel_p10, 
              aes(colour = pd, # 网络线的颜色映射mantel_test结果的显著性
                  size = rd), # 网络线的粗细映射mantel_test结果的相关性
              curvature = nice_curvature()) +ggtitle("PS_10")+ # 基于起点或终点绘制最合适的曲线 
  scale_fill_gradientn(colours =brewer.pal(11, "RdYlBu")) + # 设置热图填充颜色,也可选"RdYlBu""RdYlGn"
  scale_size_manual(values = c(0.5, 1, 2)) +# 设置网络线粗细范围
  scale_colour_manual(values =c('#b2182b','#CFCECC','#2166ac')) + # 设置网络线颜色
  guides(size = guide_legend(title = "Mantel's r",  # 调整图例次序、标题、样式
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))

PS10

legend<-get_legend(PS10)
PS10<-PS10+theme(legend.position = "none")
#ns
f_ns<-fungi[23:25,]
b_ns<-bac[23:25,]
c_ns<-collem[23:25,]
chemi_ns<-chemi[23:25,]
spe_ns<-cbind(f_ns,b_ns,c_ns)
mdata_ns <- correlate(chemi_ns)
# 把varechem拆分为2个类别的矩阵，分别与spec的每列进行mantel_test
mantel_ns <- mantel_test(spe_ns, chemi_ns,  
                          spec_select = list(Bacteria = 1:2728,Fungi= 2729:4216,Collembola=4217:4243))
# 把连续型变量划分为区间，转换为因子型变量
mantel_ns <- mantel_ns %>% mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                                             labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
                                    pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                                             labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))


NS<-qcorrplot(mdata_ns,
                type = "lower", # 热图展示下半部分
                diag = FALSE) + # 不展示对角线
  geom_square() +
  geom_couple(data = mantel_ns, 
              aes(colour = pd, # 网络线的颜色映射mantel_test结果的显著性
                  size = rd), # 网络线的粗细映射mantel_test结果的相关性
              curvature = nice_curvature()) +ggtitle("NS")+ # 基于起点或终点绘制最合适的曲线 
  scale_fill_gradientn(colours =brewer.pal(11, "RdYlBu")) + # 设置热图填充颜色,也可选"RdYlBu""RdYlGn"
  scale_size_manual(values = c(0.5, 1, 2)) +theme(legend.position = "none")+ # 设置网络线粗细范围
  scale_colour_manual(values =c('#CFCECC','#2166ac','#b2182b')) + # 设置网络线颜色
  guides(size = guide_legend(title = "Mantel's r",  # 调整图例次序、标题、样式
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))
NS
#ne
f_ne<-fungi[26:28,]
b_ne<-bac[26:28,]
c_ne<-collem[26:28,]
chemi_ne<-chemi[26:28,]
spe_ne<-cbind(f_ne,b_ne,c_ne)
mdata_ne <- correlate(chemi_ne)
# 把varechem拆分为2个类别的矩阵，分别与spec的每列进行mantel_test
mantel_ne <- mantel_test(spe_ne, chemi_ne,  
                         spec_select = list(Bacteria = 1:2728,Fungi= 2729:4216,Collembola=4217:4243))
# 把连续型变量划分为区间，转换为因子型变量
mantel_ne <- mantel_ne %>% mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                                           labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
                                  pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                                           labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))


NE<-qcorrplot(mdata_ne,
              type = "lower", # 热图展示下半部分
              diag = FALSE) + # 不展示对角线
  geom_square() +
  geom_couple(data = mantel_ne, 
              aes(colour = pd, # 网络线的颜色映射mantel_test结果的显著性
                  size = rd), # 网络线的粗细映射mantel_test结果的相关性
              curvature = nice_curvature()) +ggtitle("NE")+ # 基于起点或终点绘制最合适的曲线 
  scale_fill_gradientn(colours =brewer.pal(11, "RdYlBu")) + # 设置热图填充颜色,也可选"RdYlBu""RdYlGn"
  scale_size_manual(values = c(0.5, 1, 2)) +theme(legend.position = "none")+ # 设置网络线粗细范围
  scale_colour_manual(values =c('#CFCECC','#2166ac','#b2182b')) + # 设置网络线颜色
  guides(size = guide_legend(title = "Mantel's r",  # 调整图例次序、标题、样式
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))
NE

#Merge
MAN<-grid.arrange(PS,PS1,PS3,PS10,NS,NE,ncol=2)
MAN$widths <- unit(c(0.5, 0.5), "npc")
grid.arrange(MAN)
# Add legends back
MAN_legend <- arrangeGrob(MAN,right=legend)
MAN_legend 
# Show the arranged plots with legends
grid.arrange(MAN_legend)

