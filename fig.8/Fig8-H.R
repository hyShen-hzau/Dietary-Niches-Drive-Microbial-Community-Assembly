library(readxl)
library(ggplot2)
library(ggalt)
library(tidyverse)
library(dplyr)
library(tibble)
library(ggsci)

mytheme <-  theme_bw() +
  theme(axis.text.x = element_text(size = 8),axis.text.y = element_text(size = 8))+
  theme(axis.title.y= element_text(size=6))+theme(axis.title.x = element_text(size = 12))+
  theme(legend.title=element_text(size=0),legend.text = element_text(size=6))+theme(legend.position = "right")


#sheet9-------------------------
library(readxl)
data9 <- read_xlsx("../ARGs.xlsx",sheet = 2)

group <- data9[c(3,4)]
data9 <- data9[c(3,6:10)]
  
mypal <- pal_aaas(alpha = 0.7)(10)
library("scales")
show_col(mypal)

# library(circlize)
# col_fun = colorRamp2(c("R","I","S"), c("#DC0000B2", "#91D1C2B2", "#4DBBD5B2"))


hdata <- data9 %>% left_join(group,by="name")
hdata$Group <- factor(hdata$Group)
hdata <- arrange(hdata,Group) %>% column_to_rownames("name") 

hdata$MC9 <- factor(hdata$MC9,levels = c("R","I","S"))
hdata$MF6 <- factor(hdata$MF6,levels = c("R","I","S"))
hdata$PF9 <- factor(hdata$PF9,levels = c("R","I","S"))
hdata$SB1 <- factor(hdata$SB1,levels = c("R","I","S"))
hdata$SF1 <- factor(hdata$SF1,levels = c("R","I","S"))



library(ComplexHeatmap)
col = c("R"="#EABB77","I"="#91D1C2B2","S"= "#4DBBD5B2")
group1 <- hdata[5] %>% remove_rownames()
ha <-  rowAnnotation(df=group1,
                         annotation_name_side = "left"
  #gp = gpar(col = "black")
)
plot(ha)
ha1 <-  rowAnnotation(bar = sample(letters[1:7], 19, replace = TRUE))
p10 <- Heatmap(hdata[6],column_names_rot = 0,name="Types of antibiotics",
        
        border_gp = gpar(col = "black", lwd =2 ))+
        Heatmap(hdata[1:5],col=col,cluster_rows = F,
        column_names_rot = 0,name="Antibiotics resistance",
        border_gp = gpar(col = "black", lwd = 2), #设置热图的边框与粗度
        rect_gp = gpar(col = "white", lwd = 1),
        heatmap_legend_param =  list(
        at = c("R", "I", "S"),
        labels = c("R = tolerance", "I = median", "S = intolerance"))
        )
p10   

pdf("../result/heatmap_anti.pdf",width = 12 , height = 8)
draw(p10)
dev.off()

png(file="../result/heatmap_anti.png", width = 2400, height = 1600, res = 300)
draw(p10)
dev.off()
cat("已保存为: heatmap.png (300 DPI)\n")



library(patchwork)
q1 <- p4+p5 + plot_layout(guides = 'collect')
q2 <- p7+p8 + plot_layout(guides = 'collect')

(p2 + p3 + p4 + plot_layout(widths = c(2,1,1)))/(p5+p6+plot_spacer()+plot_spacer()+plot_layout(widths = c(1,1,1,1)))/(p7+p8+plot_spacer()+plot_spacer()+plot_layout(widths = c(1,1,1,1)))

a1 <- (p2+p3+p6+ plot_layout(widths = c(2,1,1)))/(q1|plot_spacer())/(q2|plot_spacer())
ggsave(a1 ,device = pdf,file="../result/resisitance/cowplot1.pdf")
ggsave(a1 ,device = png,file="../result/resisitance/cowplot.png",dpi = 300)
