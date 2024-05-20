otu <- t(otutab)
otu<-otu[-18,]
otu.distance <- vegdist(otu)
PCoA <- cmdscale (otu.distance,eig=TRUE)
pc12 <- PCoA$points[,1:2]
pc <- round(PCoA$eig/sum(PCoA$eig)*100,digits=2)#解释度
pc12 <- as.data.frame(pc12)
pc12$samples <- row.names(pc12)
groups <- metadata
colnames(groups)[1:2] <- c('samples','Treatment')#改列名
df <- merge(pc12,groups,by="samples")

head(df)
df$Treatment<-factor(df$Treatment,levels = c("C","PS","PV_1","PV_3","PV_10","PI_1","PI_3","PI_10","NS","NE"))
#绘图
p1<-ggplot(df,aes(x=V1, y=V2,color=Treatment,shape=Treatment,))+#指定数据、X轴、Y轴
  geom_point(size=4)+scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9))+
  theme_bw()+labs(x=paste0("PC1 ",pc[1],"%"), y=paste0("PC2 ",pc[2],"%"))+ggtitle("PCoA Bacteria Community")+
  theme(plot.title = element_text(hjust = 0.5))

p1<-p1+annotate("text", x = 0.05, y = 0.4, label = "Adonis test: R² = 0.456***",  fontface = "bold")+
  theme(legend.position = "none")

p1+ theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

p1


Adonis <- adonis2(otu.distance~Treatment,data=df,
                  distance = "bray",
                  permutations = 999)
Adonis

#Fungi
otu_f <- t(otutab_f)
otu.distance_f <- vegdist(otu_f)
PCoA_f <- cmdscale (otu.distance_f,eig=TRUE)
pc12_f <- PCoA_f$points[,1:2]
pc_f <- round(PCoA_f$eig/sum(PCoA_f$eig)*100,digits=2)#解释度
pc12_f <- as.data.frame(pc12_f)
pc12_f$samples <- row.names(pc12_f)

colnames(metadata_f)[1:2] <- c('samples','Treatment')#改列名
df_f <- merge(pc12_f,metadata_f,by="samples")

head(df_f)
df_f$Treatment<-factor(df$Treatment,levels = c("C","PS","PV_1","PV_3","PV_10","PI_1","PI_3","PI_10","NS","NE"))
#绘图
p1_f<-ggplot(df_f,aes(x=V1, y=V2,color=Treatment,shape=Treatment,))+#指定数据、X轴、Y轴
  geom_point(size=4)+scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9))+
  theme_bw()+labs(x=paste0("PC1 ",pc[1],"%"), y=paste0("PC2 ",pc[2],"%"))+ggtitle("PCoA Fungi Community")+
  theme(plot.title = element_text(hjust = 0.5))

p1_f<-p1_f+annotate("text", x = 0.3, y = 0.3, label = "Adonis test: R² = 0.382**",  fontface = "bold")

p1_f+ theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))


Adonis_f <- adonis2(otu.distance_f~Treatment,data=df_f,
                  distance = "bray",
                  permutations = 999)
Adonis_f

#Merge
beta_p <- grid.arrange(p1,p1_f,ncol = 2)
beta_p