devtools::install_github("brendanf/FUNGuildR")
library(FUNGuildR)
library(R.filesets)
library(xlsx)
# load ITS data
otutab_f <- read.table(file="otutab.txt",sep="\t",header=T,check.names=FALSE ,row.names=1)
metadata_f <- read.table("metadata.tsv", sep='\t', header=T,check.names=FALSE )
rownames(metadata_f)<-metadata$SampleID

taxonomy_f<-read.table("taxonomy.tsv",sep='\t', header=T,check.names=FALSE)
rownames(taxonomy_f)<-taxonomy_f$`Feature ID`
taxonomy_f %<>% tidy_taxonomy
# create microtable object
meco_fungi<- microtable$new( otu_table = otutab_f,
                      tax_table = taxonomy_f,
                         auto_tidy = F)
# remove the taxa not assigned in the Kingdom "k__Fungi"
meco_fungi$tax_table %<>% base::subset(Kingdom == "k__Fungi")
# use tidy_dataset() to make OTUs and samples information consistent across files
meco_fungi$tidy_dataset()



#fungilde
otutab_f$Feature_ID<-rownames(otutab_f)
guild<-cbind(taxonomy_f,otutab_f)
guild<-guild[,-1]
guild<-guild[,-2]


fung <- get_funguild_db()
saveRDS(fung, 'funguild.rds')
fung <- loadRDS('funguild.rds')

#读取 FUNGuild 自带的测试数据运行

fung_guilds <- funguild_assign(guild, db = fung, tax_col = 'Taxon')  #tax_col 是 OTU 表中的物种注释列

gi<-fung_guilds[fung_guilds$Taxon!="Unassigned",]
gi <- fung_guilds[!is.na(fung_guilds$taxon), ]
fun<-gi[,c(2:30,36)]
fun_ma<-fun[,-30]


#raw to scale
fun_ma$Function<-fun$trophicMode%>%as.factor()
abundance <- aggregate(. ~ Function, data = fun_ma, FUN = sum)
funguilde <- read_excel("funguilde.xls")
num<-funguilde[,-1]
num<-scale(num)
num<-as.data.frame(num)
num<-t(num)
colnames(num)<-funguilde$Function
num<-as.data.frame(num)
num$Treatment<-rownames(num)
num$Treatment%<>% factor(., levels = c("C","PS","PV_1","PV_3","PV_10","PI_1","PI_3","PI_10","NS","NE"))
me_f<- melt(num,id.vars = 'Treatment') 

fungi<-ggplot(me_f,aes(Treatment, variable)) +
  labs(x="",y="") + labs(title="Fungi function prediction")+
  geom_point(aes(size=abs(value),color=value)) +
  scale_size_area(max_size = 2) +
  scale_color_gsea() +theme(plot.title = element_text(hjust = 0.5))+
  theme_classic(base_size = 4.5) +
  theme(axis.text.y=element_text(size = 8),axis.text.x = element_text(size=7))

fungi<-fungi+theme(plot.title = element_text(size = 16, hjust = 0.5),legend.text = element_text(size=8),
                  legend.title = element_text(size=10) )
fungi

legend<-get_legend(fungi)
fungi<-fungi+theme(legend.position = "none")

merge<-grid.arrange(bac,fungi,ncol=2)

# Add legends back
merge_legend <- arrangeGrob(merge,right=legend)
merge_legend 
# Show the arranged plots with legends
grid.arrange(merge_legend)

