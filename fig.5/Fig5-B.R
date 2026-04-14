library(readxl)
library(ggplot2)
library(ggsci)
library(tidyverse)


#color+theme-------------------------------------------
mytheme <-  theme_bw() +
  theme(axis.text.x = element_text(size = 8),axis.text.y = element_text(size = 8))+
  theme(axis.title.y= element_text(size=6))+theme(axis.title.x = element_text(size = 12))+
  theme(legend.title=element_text(size=0),legend.text = element_text(size=6))+theme(legend.position = "right")


theme_custom <- 
  theme(
    # 背景设置
    panel.background = element_rect(fill = "white", color = NA),  # 白色背景，无边框
    plot.background = element_rect(fill = "white", color = NA),    # 整个图表背景白色
    
    # 坐标轴设置
    axis.line.x = element_blank(),        # 隐藏x轴线
    axis.line.y = element_line(color = "black", size = 0.5),  # 显示y轴线
    axis.ticks.x = element_blank(),      # 隐藏x轴刻度
    axis.ticks.y = element_line(color = "black"),  # 显示y轴刻度
    axis.text.x = element_blank(),        # 隐藏x轴标签
    axis.text.y = element_text(color = "black", size = rel(0.8)),  # 显示y轴标签
    
    # 边框设置
    panel.border = element_blank(),       # 隐藏面板边框
    panel.grid.major = element_blank(),   # 隐藏主要网格线
    panel.grid.minor = element_blank(),   # 隐藏次要网格线
    
    # 标题和标签设置
    axis.title.x = element_blank(),       # 隐藏x轴标题
    axis.title.y = element_text(color = "black", size = rel(1.1), 
                                margin = margin(r = 10)),
    
    # 图例设置（可选）
    legend.position = "none"              # 隐藏图例（根据需要调整）
  )

phy.cols <- c("#85BA8F", "#A3C8DC",
              "#349839","#EA5D2D",
              "#EABB77","#F09594","#2072A8","black")



#acid-------------------------
data1 <- read_xls("../data/factor/acid.xls",sheet = 1)

data1$name<-factor(
  data1$name,levels = data1$name)
colnames(data1)
#"name"           "Acetic_acid"    "sd_ace"         "Propionic_acid" "sd_pro"         "Butyric_acid"   "sd_but"     
data1$type <- "1"

p1<- ggplot(data=data1,
            aes(x=name,y= Butyric_acid))+
  geom_bar(data=data1,
           aes(x=name,y= Butyric_acid),
           stat = 'identity',#position = position_dodge2(),
           width = 0.7)+
  geom_errorbar(aes(ymin = Butyric_acid - sd_but, ymax = Butyric_acid + sd_but),
                position = position_dodge(width = 0.7),
                # 误差线粗细和宽度
                size = 0.5,width = 0.2)+
  labs(x="",y="")+
  theme(axis.text.x = element_text(size = 12),
        legend.position = 'none')+ 
  theme_custom
p1

p <- ggplot()+
  geom_point(data = data1,
             aes(x=name,y=Acetic_acid,
             stroke = 0.8) , size = 2)+
  #geom_text(data = data1_1,aes(label = sig, y = 14),color = "black", size = 5) +
  geom_line(data = data1, aes(x=name,y=Acetic_acid,group = type),
            color = "#EA5D2D", linewidth = 1.6)+
  geom_errorbar(data = data1, aes(x = name,
                                  ymin = Acetic_acid - sd_ace, ymax = Acetic_acid + sd_ace),
                position = position_dodge(width = 0.7), color ="black",
                size = 0.5, linewidth = 0.8, width = 0.2) +
  geom_point(data = data1,
             aes(x=name,y=Propionic_acid,
                 color = "black",stroke = 0.8), size = 2, fill = "#3C548899")+
  geom_line(data = data1, aes(x=name,y=Propionic_acid,group = type), linewidth = 1.6)+
  geom_errorbar(data = data1, aes(x = name,
                                  ymin = Propionic_acid - sd_pro, ymax = Propionic_acid + sd_pro),
                position = position_dodge(width = 0.7), color ="black",
                size = 0.5, linewidth = 0.8, width = 0.2) +
  scale_y_break(c(150, 799),  # 在201-799之间添加断点
                scales = 2,  # 上下占比2：1
                ticklabels = c(seq(100, 150, by = 50), 
                               seq(800, 1700, by = 200))) +
  labs(y="Organic acid production (mg/L)")+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.y = element_line(colour = "black", linewidth = 0.8),
    axis.line.x = element_line(colour = "black", linewidth = 0.8),
    axis.text.x = element_text(
      angle = 90,     # 旋转90度
      hjust = 1,      # 右对齐（通常效果最好）
      vjust = 0.5     # 垂直居中
    ),
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )


p
p_1 <- p %>% 
  insert_top(p1, height=.3)

p_1
ggsave(p_1 ,file="../result/factor/acid.pdf",width = 20 , height = 10,
       units = "cm", onefile=F)


