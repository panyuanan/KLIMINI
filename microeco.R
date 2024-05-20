library(magrittr)
library(microeco)
library(gridExtra)

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


Phylum <- trans_abund$new(dataset = df, taxrank = "Phylum", ntaxa = 10)
Phylum$plot_bar(others_color = "grey70",#剩余分类的填充色
                facet = "Group", #根据组进行分面
                xtext_keep = T, #是否显示样本名称
                legend_text_italic = F)#设置图例物种名称是否斜体显示

Phylum2 <- trans_abund$new(dataset = df, taxrank = "Phylum", ntaxa = 10, groupmean = "Group")
Phylum2$plot_bar(others_color = "grey70",#剩余分类的填充色
                 legend_text_italic = FALSE)+
  theme_bw() + theme(axis.title.y = element_text(size = 18),
                     panel.grid = element_blank(),
                     axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
                     axis.text.y = element_text(size = 10))

Phylum_nitro <- trans_abund$new(dataset = df, taxrank = "Phylum",input_taxaname = "Nitrospirae")
Phylum_nitro$plot_box(group = "Group", xtext_angle = 30)


# 2.1.1 构建物种组成丰度表-根据自己要求设置细节
a1 <- trans_abund$new(
  dataset = df, 
  taxrank = "Phylum", 
  show = 0, # 根据相对丰度进行过滤
  ntaxa = 10, # 展示top10丰度门
  groupmean = "Treatment", # 根据某类分类因子计算均值。
  use_percentage = TRUE, # 展示相对丰度。
)
a1$data_abund # 存储着门组成数据

# 2.1.2 绘制每个treats的门水平物种组成饼图
pie_plot <- a1$plot_pie(
  # 颜色顺序与丰度顺序正好相反，设置颜色是按着相反顺序放置颜色即可。
  color_values = rev(c(ggsci::pal_d3("category10")(10)[-8],"grey50")),
  facet_nrow = 4, # 分面设置4行
  strip_text = 12, # treats分组标题字体大小
  legend_text_italic = FALSE # 图例字体非斜体
)

pie_plot

t1 <- trans_abund$new(dataset = df, taxrank = "Phylum", ntaxa = 7, groupmean = "Treatment")
dou_plot<-t1$plot_donut(facet_nrow = 4,strip_text=10,label=FALSE)
dou_plot

# use "Type" column in sample_table
huaban <- df$merge_samples(use_group = "Treatment")
h1 <- trans_venn$new(huaban)
ven<-h1$plot_venn(petal_plot = TRUE, petal_center_size = 25, petal_r = 1.5, petal_a = 3, petal_move_xy = 3.8, petal_color_center = "#BEBADA")
ven
#FungiIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
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


Phylum <- trans_abund$new(dataset = df, taxrank = "Phylum", ntaxa = 10)
Phylum$plot_bar(others_color = "grey70",#剩余分类的填充色
                facet = "Group", #根据组进行分面
                xtext_keep = T, #是否显示样本名称
                legend_text_italic = F)#设置图例物种名称是否斜体显示

Phylum2 <- trans_abund$new(dataset = df, taxrank = "Phylum", ntaxa = 10, groupmean = "Group")
Phylum2$plot_bar(others_color = "grey70",#剩余分类的填充色
                 legend_text_italic = FALSE)+
  theme_bw() + theme(axis.title.y = element_text(size = 18),
                     panel.grid = element_blank(),
                     axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
                     axis.text.y = element_text(size = 10))

Phylum_nitro <- trans_abund$new(dataset = df, taxrank = "Phylum",input_taxaname = "Nitrospirae")
Phylum_nitro$plot_box(group = "Group", xtext_angle = 30)


# 2.1.1 构建物种组成丰度表-根据自己要求设置细节
a1 <- trans_abund$new(
  dataset = df, 
  taxrank = "Phylum", 
  show = 0, # 根据相对丰度进行过滤
  ntaxa = 10, # 展示top10丰度门
  groupmean = "Treatment", # 根据某类分类因子计算均值。
  use_percentage = TRUE, # 展示相对丰度。
)
a1$data_abund # 存储着门组成数据

# 2.1.2 绘制每个treats的门水平物种组成饼图
pie_plot <- a1$plot_pie(
  # 颜色顺序与丰度顺序正好相反，设置颜色是按着相反顺序放置颜色即可。
  color_values = rev(c(ggsci::pal_d3("category10")(10)[-8],"grey50")),
  facet_nrow = 4, # 分面设置4行
  strip_text = 12, # treats分组标题字体大小
  legend_text_italic = FALSE # 图例字体非斜体
)

pie_plot
#dounought plot
t1_f <- trans_abund$new(dataset = df_f, taxrank = "Order", ntaxa = 7, groupmean = "Treatment")
dou_plot_f<-t1_f$plot_donut(facet_nrow = 4,strip_text=10,label=FALSE)
dou_plot_f

# use "Type" column in sample_table
huaban_f <- df_f$merge_samples(use_group = "Treatment")
h1_f <- trans_venn$new(huaban_f)
ven_f<-h1_f$plot_venn(petal_plot = TRUE, petal_center_size = 25, petal_r = 1.5, petal_a = 3, petal_move_xy = 3.8, petal_color_center = "#BEBADA")
ven_f
#merge
ven_p <- grid.arrange(ven,ven_f,ncol = 2)

h1$data_summary
