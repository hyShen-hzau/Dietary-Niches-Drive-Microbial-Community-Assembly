library(ggplot2)
library(readr)
library(readxl)
library(ggsignif)
library(RColorBrewer)
library(tidyverse)
library(corrplot)
library(psych)
library(purrr)
library(permute)
library(lattice)
library(vegan)
library(phyloseq)#parse_taxonomy_qiime
library(ggh4x)
library(ggsci)
library(patchwork)
library(cowplot)

#--------------------------------color
cols3 <- c("#FF7802","#81D8D0","#5F559BE5","#2CA02CFF")


#Bac_color----------------------------------
library(cols4all)
c4a_gui()  # 弹出Shiny界面预览
palette <- c("#FF9898FF", "#D9636CFF","#1BB6AFFF","#D72000FF",  "#9093A2FF", 
             "#E64B35E5","#4DBBD5E5","#00A087E5","#3C5488E5",
             "#F39B7FE5","#8491B4E5","#91D1C2E5","#DC0000E5",
             "#7E6148E5","#B09C85E5","#3B4992E5", "#EE0000E5",
             "#FF7802","#002EA6","#FFE76F","#01847F",
             "#008B45E5", "#631879E5", "#008280E5" ,"#BB0021E5" ,
             "#5F559BE5", "#FF0000","#FFDC91E5", "#91D1C2FF",
             "#3C5488FF","#D3D3D3","#9EC8FFFF","#000000","#BC3C29E5" ,"#0072B5E5" ,"#E18727E5",
             "#20854EE5" , "#7876B1E5" ,"#6F99ADE5" ,"#FFDC91E5", "#EE4C97E5")


pal <- c("#F56455FF", "#00A087FF", "#87C785FF", "#FFDC91E5",
         "#FD0C81FF", "#FFED4DFF", "#C34582FF","#EBA49EFF", "#73D8FEFF",
         "#FF7802","#81D8D0","#5F559BE5","#2CA02CFF","#4DBBD5FF","#00A087FF",
         "#FF7802","#81D8D0","#5F559BE5","#2CA02CFF","#4DBBD5FF","#00A087FF",
         "#FF7802","#81D8D0","#5F559BE5","#2CA02CFF","#4DBBD5FF","#00A087FF",
         "#FF7802","#81D8D0","#5F559BE5","#2CA02CFF","#4DBBD5FF","#00A087FF",
         "#FF7802","#81D8D0","#5F559BE5","#2CA02CFF","#4DBBD5FF","#00A087FF"
         )

#Bac_plot--------------------------------
taxa <- read_csv("../data/myPlan.csv")
x1 <- read.csv("../data/Tax/Tax_32/Bac_30%.csv", row.names = 1)
taxa <- read_csv("../data/myPlan.csv")
colnames(taxa)[1] <- "Sampleid"
x1 <- x1 %>% as.data.frame() %>% mutate(Taxa=taxa$Sampleid) 
x1 <- x1 %>% remove_rownames()
#%>%   column_to_rownames("Taxa")

# x1 <- x1 %>% column_to_rownames("Taxa")
# rownames(x1) <- taxa$Sampleid

x1 <- x1 %>% pivot_longer(cols = !Taxa,names_to = "Samples",
             values_to = "number") 
colnames(x1)[1] <- "Sampleid"
#taxa <- taxa[-5]

meta_taxa <- taxa %>% 
  inner_join(.,x1,by="Sampleid")

title <- taxa$Sampleid
meta_taxa$Sampleid <- factor(meta_taxa$Sampleid , levels = title , ordered = T)
meta_taxa$Plan2 <- factor(meta_taxa$Plan2)
meta_taxa$Repeat<- factor(meta_taxa$Repeat ,
                                      levels = c(1:10) , ordered = T)

its1 <- read.csv("../data/Tax/Tax_32/Bac_30%.csv", row.names = 1)
its1 <- t(its1) %>% as.data.frame()
sum2 <- apply(its1, 1, sum) %>% as.data.frame()
orsum <- rownames(sum2)
meta_taxa$Samples <- factor(meta_taxa$Samples,levels = orsum,ordered = T)



p2 <- ggplot(meta_taxa,aes(Sampleid,number,fill=Samples))+
  geom_col(position="stack") +
  facet_nested(.~Plan2+Plan3+Repeat,drop=T,scale="free",
               space="free_x",switch="x",
               strip =strip_nested(
                 background_x = elem_list_rect(fill =pal),
                 by_layer_x = F
               ))+
  scale_fill_manual(values=palette,limits = orsum)+
  labs(x=NULL, y="Percent Genus Abundance")+
  scale_y_continuous(expand = c(0,0),labels=scales::percent)+
  theme(strip.background = element_rect(fill="white",color="black"),
        panel.spacing = unit(0,"lines"),
        strip.text.x = element_text(size=7,color="black"),
        axis.text.y=element_text(size=8,color="black"),
        axis.title.y = element_text(size=10,color="black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key=element_blank(),   # 图例键为空
        #legend.position = "top",
        legend.text = element_text(color="black",size=9,face="italic"), # 定义图例文本
        legend.spacing.x=unit(0.1,'cm'), # 定义文本书平距离
        legend.key.width=unit(0.3,'cm'), # 定义图例水平大小
        legend.key.height=unit(0.3,'cm'), # 定义图例垂直大小
        legend.background=element_blank(),
        panel.grid.major=element_blank(), # 移除主网格线
        panel.grid.minor=element_blank())+
  # guides(fill = guide_legend( ncol = 1, byrow = TRUE))+#图例设为i一列
  guides(fill=guide_legend(fill=guide_legend(reverse = TRUE),
                           title = NULL,ncol = 1,
                           keywidth = unit(8,units = "points"),
                           keyheight = unit(8,units = "points")))
  
p2

unique(meta_taxa$Plan3)
meta_taxa$Plan3 <- factor(meta_taxa$Plan3 ,
                          levels=c("M.salmoides","Siniperca chuatsi","H.molitrix",
                                   "C.idellus","P.fulvidraco"))

colnames(meta_taxa)
p1 <- ggplot(meta_taxa,aes(Sampleid,Bacteria,fill=Plan3))+
  geom_col(width = 0.9)+theme_grey()+
  scale_fill_manual(values=c("#FD0C81FF", "#FFED4DFF", "#C34582FF","#EBA49EFF", "#73D8FEFF"))+
  facet_nested(.~Plan3,drop=TRUE,scale="free",space="free")+
  scale_x_discrete(expand=c(0,0)) +
  scale_y_continuous(expand = c(0,0),labels=scales::scientific_format(digits=1))+
  theme(strip.text = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill='white'),
        panel.spacing = unit(0.01,"lines"),
        axis.text.x = element_blank(),
        legend.key=element_blank(),   # 图例键为空
        # legend.text = element_text(color="black",size=10), # 定义图例文本
        # legend.spacing.x=unit(0.1,'cm'), # 定义文本书平距离
        # legend.key.width=unit(0.5,'cm'), # 定义图例水平大小
        # legend.key.height=unit(0.5,'cm'), # 定义图例垂直大小
        # legend.background=element_blank(),
        panel.grid.major=element_blank(), # 移除主网格线
        panel.grid.minor=element_blank())+
  guides(fill="none")+
  labs(x=NULL, y="Clean Reads")
p1

p3 <- p1/p2 + plot_layout(ncol = 1, heights = c(1,3))

p3

ggsave(p2,device = "pdf",file="../result/Tax/Bac_genus.pdf",width = 21, height = 16, units = "cm",onefile=F)   
ggsave(p2,device = "png",file="../result/Tax/Bac_genus.png",
       dpi = 300,width = 20, height = 16, units = "cm") 
ggsave(p3,device = "pdf",file="../result/Tax/Bac_genus_reads.pdf",width = 21, height = 16, units = "cm",onefile=F)   
ggsave(p3,device = "png",file="../result/Tax/Bac_genus_reads.png",
       dpi = 300,width = 20, height = 16, units = "cm") 


its <- read.csv("../data/Tax/Tax_32/Bac_30%.csv")
group <- read.csv("../data/myPlan.csv")
colnames(group)[1] <- "Sampleid"

colnames(its)[1] <- "Sampleid"
its1 <- its %>% left_join(group,by = "Sampleid")

its2 <- its1 %>% select(orsum,Plan3) %>% group_by(Plan3) %>% 
  summarise(across(.fns = mean))

colnames(its2)[1] <- "Sampleid"
x1 <- its2 %>% pivot_longer(cols = !Sampleid,names_to = "Samples",
                            values_to = "number") 
colnames(x1)[1] <- "Sampleid"

x1$Samples <- factor(x1$Samples,levels = orsum,ordered = T)



q2 <- ggplot(x1,aes(Sampleid,number,fill=Samples))+
  geom_bar(stat = "identity",position = "stack",
           width = 0.4) +
  scale_fill_manual(values=palette,limits = orsum)+
  labs(x=NULL, y="Average Percent Genus Abundance")+
  scale_y_continuous(expand = c(0,0),labels=scales::percent)+
  theme(strip.background = element_rect(fill="white",color="black"),
        panel.spacing = unit(0,"lines"),
        strip.text.x = element_text(size=10,color="black"),
        axis.text.x=element_text(angle = 45, size=8,color="black",hjust = .5,vjust = .5),
        axis.text.y=element_text(size=8,color="black"),
        axis.title.y = element_text(size=10,color="black"),
        legend.key=element_blank(),   # 图例键为空
        #legend.position = "top",
        legend.text = element_text(color="black",size=9,face="italic"), # 定义图例文本
        legend.spacing.x=unit(0.1,'cm'), # 定义文本书平距离
        legend.key.width=unit(0.3,'cm'), # 定义图例水平大小
        legend.key.height=unit(0.3,'cm'), # 定义图例垂直大小
        legend.background=element_blank(),
        panel.grid.major=element_blank(), # 移除主网格线
        panel.grid.minor=element_blank())+
  # guides(fill = guide_legend( ncol = 1, byrow = TRUE))+#图例设为i一列
  guides(fill=guide_legend(fill=guide_legend(reverse = TRUE),
                           title = NULL,ncol = 1,
                           keywidth = unit(8,units = "points"),
                           keyheight = unit(8,units = "points")))

q2

ggsave(q2,device = "pdf",file="../result/Tax/Bac_genus_aev.pdf",width = 15, height = 15, units = "cm",onefile=F)   
ggsave(q2,device = "png",file="../result/Tax/Bac_genus_aev.png",
       dpi = 300,width = 12, height = 12, units = "cm") 






#Euk-----------------------------------------------------------------------
Euk <- read.csv("../data/Tax/Euk.csv",row.names = 1)
group1 <- read.csv("../data/myPlan.csv")
var <- c(t(group1[1]))

euk <- Euk %>% select(var, Genus)
#euk[1:32] <- apply(euk[1:32],1,as.numeric)
length(unique(euk$Genus))
euk <- euk %>%
  group_by(Genus) %>%
  summarise(across(.cols = everything(), .fns = ~sum(.x, na.rm = TRUE)))
  #summarise(across(.cols = everything(), .fns = sum, na.rm = TRUE))
euk <- euk %>%
  column_to_rownames("Genus") %>%
  decostand(method = "total", MARGIN = 2)

write.csv(euk,"../data/Tax/Tax_32/euk_%.csv")

euk #<- read.csv("../data/Tax/Tax_32/euk_%.csv" , row.names = 1)

#----------------------------------- major microorganisms
eukdata <- euk %>% 
  rownames_to_column("sampleid") %>% 
  rowwise() %>% 
  mutate(total=sum(across(-1)),.after=1) %>% 
  ungroup() %>% 
  mutate(prop=total/sum(total),.after=1) %>% 
  arrange(desc(prop)) %>% 
  filter(row_number() <=20) %>% 
  select(-c("prop","total")) %>% 
  column_to_rownames("sampleid")

x <- t(eukdata)
sum <- apply(x,1,sum)
x <- as.data.frame(x)
title <- rownames(x)
x1 <- cbind( x , sum)

x1 <-  x1 %>% 
  mutate(Other=1-sum) %>% 
  select(-c("sum"))

write.csv( x1,"../data/Tax/Tax_32/euk_20%.csv")


#euk_color----------------------------------

palette <- c(
             "#F39B7FE5","#8491B4E5","#91D1C2E5","#DC0000E5",
             "#4DBBD5E5","#00A087E5","#3C5488E5",
             "#7E6148E5","#B09C85E5","#3B4992E5", "#EE0000E5",
             "#FF7802","#002EA6","#000000","#FFE76F","#01847F",
             "#008B45E5", "#631879E5", "#008280E5" ,"#BB0021E5" ,"#D3D3D3",
             "#5F559BE5", "#FF0000","#00A087FF", "#91D1C2FF",
             "#3C5488FF","#BC3C29E5" ,"#0072B5E5" ,"#E18727E5",
             "#20854EE5" , "#7876B1E5" ,"#6F99ADE5" ,"#FFDC91E5", "#EE4C97E5")


pal <- c("#F56455FF", "#00A087FF", "#87C785FF", "#FFDC91E5",
         "#FD0C81FF", "#FFED4DFF", "#C34582FF","#EBA49EFF", "#73D8FEFF",
         "#FF7802","#81D8D0","#5F559BE5","#2CA02CFF","#4DBBD5FF","#00A087FF",
         "#FF7802","#81D8D0","#5F559BE5","#2CA02CFF","#4DBBD5FF","#00A087FF",
         "#FF7802","#81D8D0","#5F559BE5","#2CA02CFF","#4DBBD5FF","#00A087FF",
         "#FF7802","#81D8D0","#5F559BE5","#2CA02CFF","#4DBBD5FF","#00A087FF",
         "#FF7802","#81D8D0","#5F559BE5","#2CA02CFF","#4DBBD5FF","#00A087FF"
)


#euk_plot--------------------------------
taxa <- read_csv("../data/myPlan.csv")
x1 <- read.csv("../data/Tax/Tax_32/euk_20%.csv", row.names = 1)
taxa <- read_csv("../data/myPlan.csv")
colnames(taxa)[1] <- "Sampleid"
x1 <- x1 %>% as.data.frame() %>% mutate(Taxa=taxa$Sampleid)
x1 <- x1 %>% remove_rownames()
#%>%   column_to_rownames("Taxa")

# x1 <- x1 %>% column_to_rownames("Taxa")
# rownames(x1) <- taxa$Sampleid

x1 <- x1 %>% pivot_longer(cols = !Taxa,names_to = "Samples",
                          values_to = "number") 
colnames(x1)[1] <- "Sampleid"
#taxa <- taxa[-5]

meta_taxa <- taxa %>% 
  inner_join(.,x1,by="Sampleid")

title <- taxa$Sampleid
meta_taxa$Sampleid <- factor(meta_taxa$Sampleid , levels = title , ordered = T)
meta_taxa$Plan2 <- factor(meta_taxa$Plan2)
meta_taxa$Repeat<- factor(meta_taxa$Repeat ,
                          levels = c(1:10) , ordered = T)

its1 <- read.csv("../data/Tax/Tax_32/euk_20%.csv", row.names = 1)
its1 <- t(its1) %>% as.data.frame()
sum2 <- apply(its1, 1, sum) %>% as.data.frame()
orsum <- rownames(sum2)
meta_taxa$Samples <- factor(meta_taxa$Samples,levels = orsum,ordered = T)



p5 <- ggplot(meta_taxa,aes(Sampleid,number,fill=Samples))+
  geom_col(position="stack") +
  facet_nested(.~Plan2+Plan3+Repeat,drop=T,scale="free",
               space="free_x",switch="x",
               strip =strip_nested(
                 background_x = elem_list_rect(fill =pal),
                 by_layer_x = F
               ))+
  scale_fill_manual(values=palette,limits = orsum)+
  labs(x=NULL, y="Percent Genus Abundance")+
  scale_y_continuous(expand = c(0,0),labels=scales::percent)+
  theme(strip.background = element_rect(fill="white",color="black"),
        panel.spacing = unit(0,"lines"),
        strip.text.x = element_text(size=7,color="black"),
        axis.text.y=element_text(size=8,color="black"),
        axis.title.y = element_text(size=10,color="black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key=element_blank(),   # 图例键为空
        #legend.position = "top",
        legend.text = element_text(color="black",size=9,face="italic"), # 定义图例文本
        legend.spacing.x=unit(0.1,'cm'), # 定义文本书平距离
        legend.key.width=unit(0.3,'cm'), # 定义图例水平大小
        legend.key.height=unit(0.3,'cm'), # 定义图例垂直大小
        legend.background=element_blank(),
        panel.grid.major=element_blank(), # 移除主网格线
        panel.grid.minor=element_blank())+
  # guides(fill = guide_legend( ncol = 1, byrow = TRUE))+#图例设为i一列
  guides(fill=guide_legend(fill=guide_legend(reverse = TRUE),
                           title = NULL,ncol = 1,
                           keywidth = unit(8,units = "points"),
                           keyheight = unit(8,units = "points")))

p5

#unique(meta_taxa$Plan3)
meta_taxa$Plan3 <- factor(meta_taxa$Plan3 ,
                          levels=c("M.salmoides","Siniperca chuatsi","H.molitrix",
                                   "C.idellus","P.fulvidraco"))

colnames(meta_taxa)
p4 <- ggplot(meta_taxa,aes(Sampleid,Eukaryota,fill=Plan3))+
  geom_col(width = 0.9)+theme_grey()+
  scale_fill_manual(values=c("#FD0C81FF", "#FFED4DFF", "#C34582FF","#EBA49EFF", "#73D8FEFF"))+
  facet_nested(.~Plan3,drop=TRUE,scale="free",space="free")+
  scale_x_discrete(expand=c(0,0)) +
  scale_y_continuous(expand = c(0,0),labels=scales::scientific_format(digits=1))+
  theme(strip.text = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill='white'),
        panel.spacing = unit(0.01,"lines"),
        axis.text.x = element_blank(),
        legend.key=element_blank(),   # 图例键为空
        # legend.text = element_text(color="black",size=10), # 定义图例文本
        # legend.spacing.x=unit(0.1,'cm'), # 定义文本书平距离
        # legend.key.width=unit(0.5,'cm'), # 定义图例水平大小
        # legend.key.height=unit(0.5,'cm'), # 定义图例垂直大小
        # legend.background=element_blank(),
        panel.grid.major=element_blank(), # 移除主网格线
        panel.grid.minor=element_blank())+
  guides(fill="none")+
  labs(x=NULL, y="Clean Reads")
p4


p6 <- p4/p5 + plot_layout(ncol = 1, heights = c(1,3))

p6

ggsave(p5,device = "pdf",file="../result/Tax/Euk_genus.pdf",width = 21, height = 16, units = "cm",onefile=F)   
ggsave(p5,device = "png",file="../result/Tax/Euk_genus.png",
       dpi = 300,width = 20, height = 16, units = "cm") 
ggsave(p6,device = "pdf",file="../result/Tax/Euk_genus_reads.pdf",width = 21, height = 16, units = "cm",onefile=F)   
ggsave(p6,device = "png",file="../result/Tax/Euk_genus_reads.png",
       dpi = 300,width = 20, height = 16, units = "cm") 


its <- read.csv("../data/Tax/Tax_32/euk_20%.csv")
group <- read.csv("../data/myPlan.csv")
colnames(group)[1] <- "Sampleid"

colnames(its)[1] <- "Sampleid"
its1 <- its %>% left_join(group,by = "Sampleid")

its2 <- its1 %>% select(orsum,Plan3) %>% group_by(Plan3) %>% 
  summarise(across(.cols = where(is.numeric), .fns = mean, na.rm = TRUE))

colnames(its2)[1] <- "Sampleid"
x1 <- its2 %>% pivot_longer(cols = !Sampleid,names_to = "Samples",
                            values_to = "number") 
colnames(x1)[1] <- "Sampleid"

x1$Samples <- factor(x1$Samples,levels = orsum,ordered = T)



q4 <- ggplot(x1,aes(Sampleid,number,fill=Samples))+
  geom_bar(stat = "identity",position = "stack",
           width = 0.4) +
  scale_fill_manual(values=palette,limits = orsum)+
  labs(x=NULL, y="Average Percent Genus Abundance")+
  scale_y_continuous(expand = c(0,0),labels=scales::percent)+
  theme(strip.background = element_rect(fill="white",color="black"),
        panel.spacing = unit(0,"lines"),
        strip.text.x = element_text(size=10,color="black"),
        axis.text.x=element_text(angle =45, size=8,color="black",hjust = .5,vjust = .5),
        axis.text.y=element_text(size=8,color="black"),
        axis.title.y = element_text(size=10,color="black"),
        legend.key=element_blank(),   # 图例键为空
        #legend.position = "top",
        legend.text = element_text(color="black",size=9,face="italic"), # 定义图例文本
        legend.spacing.x=unit(0.1,'cm'), # 定义文本书平距离
        legend.key.width=unit(0.3,'cm'), # 定义图例水平大小
        legend.key.height=unit(0.3,'cm'), # 定义图例垂直大小
        legend.background=element_blank(),
        panel.grid.major=element_blank(), # 移除主网格线
        panel.grid.minor=element_blank())+
  # guides(fill = guide_legend( ncol = 1, byrow = TRUE))+#图例设为i一列
  guides(fill=guide_legend(fill=guide_legend(reverse = TRUE),
                           title = NULL,ncol = 1,
                           keywidth = unit(8,units = "points"),
                           keyheight = unit(8,units = "points")))

q4

ggsave(q4,device = "pdf",file="../result/Tax/Euk_genus_aev.pdf",width = 15, height = 15, units = "cm",onefile=F)   
ggsave(q4,device = "png",file="../result/Tax/Euk_genus_aev.png",
       dpi = 300,width = 12, height = 12, units = "cm") 







