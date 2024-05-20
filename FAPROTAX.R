library(ggsci)
library(aplot)
library(reshape2)
library(tidyverse)
library(ggsci)

fap <- trans_func$new(df)
fap$cal_spe_func(prok_database = "FAPROTAX")
fap$res_spe_func[1:5, 1:2]


#PER SAMPLE
fap$res_spe_func_perc %>% 
  select(1:30) %>% 
  scale() %>% 
  as.data.frame() %>% 
  rownames_to_column('sample') %>% 
  melt(id.vars = 'sample') %>% 
  ggplot(aes(sample, variable)) +
  labs(x="",y="") + 
  geom_point(aes(size=abs(value),color=value)) +
  scale_size_area(max_size = 2) +
  scale_color_gsea() +
  theme_classic(base_size = 4.5) +
  theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))

#PER TREATMENT
data<-fap$res_spe_func_perc
data$Treatment<-metadata$Treatment
data$Treatment %<>% factor(., levels = c("C","PS","PV_1","PV_3","PV_10","PI_1","PI_3","PI_10","NS","NE"))
metadata$Treatment%<>% factor(., levels = c("C","PS","PV_1","PV_3","PV_10","PI_1","PI_3","PI_10","NS","NE"))

mean_abundance <- aggregate(. ~ Treatment, data = data, FUN = mean)


mean<-mean_abundance %>% select(c(9:11,15:18,21:22,30,40:42))

  mean<-mean%>%scale()%>%as.data.frame()
  
 mean$Treatment<-unique(metadata$Treatment)
 
me<- melt(mean,id.vars = 'Treatment') 
  
bac<-ggplot(me,aes(Treatment, variable)) +
  labs(x="",y="") + labs(title="Bacteria function prediction")+
  geom_point(aes(size=abs(value),color=value)) +
  scale_size_area(max_size = 2) +
  scale_color_gsea() +theme(plot.title = element_text(hjust = 0.5))+
  theme_classic(base_size = 4.5) +
  theme(axis.text.y=element_text(size = 10),axis.text.x = element_text(size=7))

bac<-bac+theme(plot.title = element_text(size = 16, hjust = 0.5),legend.position = "none")
bac
  
  