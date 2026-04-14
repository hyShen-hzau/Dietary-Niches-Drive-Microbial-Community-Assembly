library(ReporterScore)
library(dplyr)
library(readr)
library(tidyverse)
library(ggplot2)
#update_KEGG(download_dir = "../result/Fun_enrich/")


#fish---------------------------------------------------------------------------
dat1 <- read.csv("../data/function/KEGG_NORM.csv", row.names = 1)
Kfun <- dat1[-31]

rownames(Kfun) <- as.character(Kdata$KEGG_KO)

groupdat <- read.csv("../data/myPlan.csv")
rownames(groupdat) <- groupdat$name
group <- groupdat %>% select(c("Plan2","Plan3"))


all(rownames(group) %in% colnames(Kfun))
intersect(rownames(group) , colnames(Kfun))>0


cat("Comparison: ", levels(factor(group$Plan3)))
#Comparison:  C.idellus H.molitrix M.salmoides P.fulvidraco S.chuatsi


rsa_cm_res <- GRSA(Kfun, "Plan3", group,
                        mode = "directed",
                        method = "spearman", perm = 999
)

p1 <- plot_KOs_in_pathway(rsa_cm_res,box_color = c("#FD0C81FF", "#FFED4DFF", "#C34582FF","#EBA49EFF", "#73D8FEFF"), 
                    map_id = "map00220") + scale_y_log10()

p2 <- plot_KOs_in_pathway(rsa_cm_res,box_color = c("#FD0C81FF", "#FFED4DFF", "#C34582FF","#EBA49EFF", "#73D8FEFF"), 
                          map_id = "map00010") + scale_y_log10()

p3 <- plot_KOs_in_pathway(rsa_cm_res,box_color = c("#FD0C81FF", "#FFED4DFF", "#C34582FF","#EBA49EFF", "#73D8FEFF"), 
                          map_id = "map00330") + scale_y_log10()

p4 <- plot_KOs_in_pathway(rsa_cm_res,box_color = c("#FD0C81FF", "#FFED4DFF", "#C34582FF","#EBA49EFF", "#73D8FEFF"), 
                          map_id = "map00620") + scale_y_log10()


p5 <- plot_KOs_in_pathway(rsa_cm_res,box_color = c("#FD0C81FF", "#FFED4DFF", "#C34582FF","#EBA49EFF", "#73D8FEFF"), 
                          map_id = "map00121") + scale_y_log10()
#上调，胆汁酸代谢通路

p6 <- plot_KOs_in_pathway(rsa_cm_res,box_color = c("#FD0C81FF", "#FFED4DFF", "#C34582FF","#EBA49EFF", "#73D8FEFF"), 
                    map_id = "map00564") + scale_y_log10()
#上调，甘油磷脂代谢
p7 <- plot_KOs_in_pathway(rsa_cm_res,box_color = c("#FD0C81FF", "#FFED4DFF", "#C34582FF","#EBA49EFF", "#73D8FEFF"), 
                    map_id = "map04974") + scale_y_log10()
#蛋白质降解

p8 <- plot_KOs_in_pathway(rsa_cm_res,box_color = c("#FD0C81FF", "#FFED4DFF", "#C34582FF","#EBA49EFF", "#73D8FEFF"), 
                          map_id = "map00500") + scale_y_log10()
#淀粉代谢

plot_KOs_in_pathway(rsa_cm_res,box_color = c("#FD0C81FF", "#FFED4DFF", "#C34582FF","#EBA49EFF", "#73D8FEFF"), 
                          map_id = "map00730") + scale_y_log10()

p9 <- plot_KOs_in_pathway(rsa_cm_res,box_color = c("#FD0C81FF", "#FFED4DFF", "#C34582FF","#EBA49EFF", "#73D8FEFF"), 
                          map_id = "map00740") + scale_y_log10()
#核黄素代谢，重要

p10 <- plot_KOs_in_pathway(rsa_cm_res,box_color = c("#FD0C81FF", "#FFED4DFF", "#C34582FF","#EBA49EFF", "#73D8FEFF"), 
                    map_id = "map00750") + scale_y_log10()

plot_KOs_in_pathway(rsa_cm_res,box_color = c("#FD0C81FF", "#FFED4DFF", "#C34582FF","#EBA49EFF", "#73D8FEFF"), 
                    map_id = "map00770") + scale_y_log10()

plot_KOs_in_pathway(rsa_cm_res,box_color = c("#FD0C81FF", "#FFED4DFF", "#C34582FF","#EBA49EFF", "#73D8FEFF"), 
                    map_id = "map00780") + scale_y_log10()

plot_KOs_in_pathway(rsa_cm_res,box_color = c("#FD0C81FF", "#FFED4DFF", "#C34582FF","#EBA49EFF", "#73D8FEFF"), 
                    map_id = "map00790") + scale_y_log10()

plot_KOs_in_pathway(rsa_cm_res,box_color = c("#FD0C81FF", "#FFED4DFF", "#C34582FF","#EBA49EFF", "#73D8FEFF"), 
                    map_id = "map00830") + scale_y_log10()

plot_KOs_in_pathway(rsa_cm_res,box_color = c("#FD0C81FF", "#FFED4DFF", "#C34582FF","#EBA49EFF", "#73D8FEFF"), 
                    map_id = "map00380") + scale_y_log10()


plot_list <- paste0("p", 1:10)

for(plot_name in plot_list) {
  filename <- paste0(plot_name, ".pdf")
  
  ggsave(
  filename =  paste0("../result/Fun_enrich/map-5/",filename),
  plot = get(plot_name),
  device = "pdf",
  width = 8,   
  height = 6,  
  dpi = 300    
)

message("xxx: ", filename)
}



#pattern------------------------------------------------------------------------
rsa_cm_res2 <- RSA_by_cm(Kfun, "Plan3", group,
                        method = "spearman",
                        k_num = 5, perm = 99
)

# show the patterns
plot_c_means(rsa_cm_res2, filter_membership = 0.7) + theme_bw()

plot_report_bar(rsa_cm_res2, rs_threshold = 4.5, y_text_size = 10)



#step-step----------------------------------------------------------------------
ko_pvalue <- ko.test(Kfun, "Plan3", group, method = "kruskal.test")
ko_stat <- pvalue2zs(ko_pvalue, mode = "mixed")
reporter_s1 <- get_reporter_score(ko_stat, perm = 499)
load_KOlist()

rownames(KO_abundance)



#ci vs sc-------------------------------------------------------------------
ci_sc <- Kfun[c(1:6,25:30)]
g_cisc <- group[c(1:6,25:30),]

all(rownames(g_cisc) %in% colnames(ci_sc))
intersect(rownames(g_cisc) , colnames(ci_sc))>0


cat("Comparison: ", levels(factor(g_cisc$Plan3)))
#Comparison:  C.idellus S.chuatsi

rsa_cm_res_ci <- GRSA(ci_sc, "Plan3", g_cisc,
                   mode = "directed",
                   method = "wilcox.test", perm = 999
)

plot_report_bar(rsa_cm_res_ci, rs_threshold = 4)
plot_report_bar(rsa_cm_res_ci, rs_threshold = c(-4, 4), facet_level = TRUE,
                mode = 2)


plot_KOs_in_pathway(rsa_cm_res_ci,box_color = c("#FD0C81FF",  "#73D8FEFF"), 
                          map_id = "map00220") + scale_y_log10()

plot_KOs_in_pathway(rsa_cm_res_ci,box_color = c("#FD0C81FF", "#73D8FEFF"), 
                          map_id = "map00010") + scale_y_log10()
 
plot_KOs_in_pathway(rsa_cm_res_ci,box_color = c("#FD0C81FF",  "#73D8FEFF"), 
                          map_id = "map00330") + scale_y_log10()

plot_KOs_in_pathway(rsa_cm_res_ci,box_color = c("#FD0C81FF", "#73D8FEFF"), 
                          map_id = "map00620") + scale_y_log10()

plot_KOs_in_pathway(rsa_cm_res,box_color = c("#FD0C81FF", "#FFED4DFF", "#C34582FF","#EBA49EFF", "#73D8FEFF"), 
                    map_id = "map00740") + scale_y_log10()

plot_KOs_in_pathway(rsa_cm_res_ci,box_color = c("#FD0C81FF", "#FFED4DFF", "#C34582FF","#EBA49EFF", "#73D8FEFF"), 
                    map_id = "map00740") + scale_y_log10()
#sc中富集*****************

plot_KOs_in_pathway(rsa_cm_res_ci,box_color = c("#FD0C81FF", "#73D8FEFF"), 
                    map_id = "map00780") + scale_y_log10()

plot_KOs_in_pathway(rsa_cm_res_ci,box_color = c("#FD0C81FF", "#73D8FEFF"), 
                    map_id = "map00830") + scale_y_log10()



plot_KOs_network(rsa_cm_res,
                 map_id = c("map0121", "map00071", "map00564"),
                 near_pathway = TRUE, main = "", mark_module = TRUE, 
                 str_width = 70)

plot_KEGG_map(rsa_cm_res_ci, map_id = "map00121", color_var = "Z_score",
              save_dir = "../result/Fun_enrich/map/")


#hm vs sc-------------------------------------------------------------------
hm_sc <- Kfun[c(7:12,25:30)]
g_hmsc <- group[c(7:12,25:30),]

all(rownames(g_hmsc) %in% colnames(hm_sc))
intersect(rownames(g_hmsc) , colnames(hm_sc))>0


cat("Comparison: ", levels(factor(g_hmsc$Plan3)))
#Comparison:  C.idellus S.chuatsi

rsa_cm_res_hm <- GRSA(hm_sc, "Plan3", g_hmsc,
                      mode = "directed",
                      method = "wilcox.test", perm = 999
)

plot_report_bar(rsa_cm_res_hm, rs_threshold = 4)
plot_report_bar(rsa_cm_res_hm, rs_threshold = c(-4, 4), facet_level = TRUE,
                mode = 2)


p1_2 <- plot_KOs_in_pathway(rsa_cm_res_hm,box_color = c("#FFED4DFF",  "#73D8FEFF"), 
                            map_id = "map00220") + scale_y_log10()

p2_2 <- plot_KOs_in_pathway(rsa_cm_res_hm,box_color = c("#FFED4DFF", "#73D8FEFF"), 
                            map_id = "map00010") + scale_y_log10()

p3_2 <- plot_KOs_in_pathway(rsa_cm_res_hm,box_color = c("#FFED4DFF",  "#73D8FEFF"), 
                            map_id = "map00330") + scale_y_log10()

p4_2 <- plot_KOs_in_pathway(rsa_cm_res_hm,box_color = c("#FFED4DFF", "#73D8FEFF"), 
                            map_id = "map00620") + scale_y_log10()

plot_KOs_in_pathway(rsa_cm_res_hm,box_color = c("#FFED4DFF", "#73D8FEFF"), 
                    map_id = "map04360") + scale_y_log10()

plot_KOs_network(rsa_cm_res_hm,
                 map_id = c("map00220", "map00010", "map00330", "map00620"),
                 near_pathway = FALSE, main = "", mark_module = TRUE, 
                 str_width = 70)

# plot_KEGG_map(rsa_cm_res_hm, map_id = "map00220", color_var = "Z_score",
#               save_dir = "../result/Fun_enrich/map/")

#ms vs sc-------------------------------------------------------------------
ms_sc <- Kfun[c(13:18,25:30)]
g_mssc <- group[c(13:18,25:30),]

all(rownames(g_mssc) %in% colnames(ms_sc))
intersect(rownames(g_mssc) , colnames(ms_sc))>0


cat("Comparison: ", levels(factor(g_mssc$Plan3)))
#Comparison:  C.idellus S.chuatsi

rsa_cm_res_ms <- GRSA(ms_sc, "Plan3", g_mssc,
                      mode = "directed",
                      method = "wilcox.test", perm = 999
)

plot_report_bar(rsa_cm_res_ms, rs_threshold = 3)
plot_report_bar(rsa_cm_res_ms, rs_threshold = c(-3, 3), facet_level = TRUE,
                mode = 2)


p1_2 <- plot_KOs_in_pathway(rsa_cm_res_ms,box_color = c("#FFED4DFF",  "#73D8FEFF"), 
                            map_id = "map00220") + scale_y_log10()

p2_2 <- plot_KOs_in_pathway(rsa_cm_res_ms,box_color = c("#FFED4DFF", "#73D8FEFF"), 
                            map_id = "map00010") + scale_y_log10()

p3_2 <- plot_KOs_in_pathway(rsa_cm_res_ms,box_color = c("#FFED4DFF",  "#73D8FEFF"), 
                            map_id = "map00330") + scale_y_log10()

p4_2 <- plot_KOs_in_pathway(rsa_cm_res_ms,box_color = c("#FFED4DFF", "#73D8FEFF"), 
                            map_id = "map00620") + scale_y_log10()

plot_KOs_in_pathway(rsa_cm_res_ms,box_color = c("#FFED4DFF", "#73D8FEFF"), 
                    map_id = "map04360") + scale_y_log10()

plot_KOs_network(rsa_cm_res_ms,
                 map_id = c("map00220", "map00010", "map00330", "map00620"),
                 near_pathway = FALSE, main = "", mark_module = TRUE, 
                 str_width = 70)

# plot_KEGG_map(rsa_cm_res_hm, map_id = "map00220", color_var = "Z_score",
#               save_dir = "../result/Fun_enrich/map/")


#pf vs sc-------------------------------------------------------------------
pf_sc <- Kfun[c(19:24,25:30)]
g_pfsc <- group[c(19:24,25:30),]

all(rownames(g_pfsc) %in% colnames(pf_sc))
intersect(rownames(g_pfsc) , colnames(pf_sc))>0


cat("Comparison: ", levels(factor(g_pfsc$Plan3)))
#Comparison:  C.idellus S.chuatsi

rsa_cm_res_pf <- GRSA(pf_sc, "Plan3", g_pfsc,
                      mode = "directed",
                      method = "wilcox.test", perm = 999
)

plot_report_bar(rsa_cm_res_pf, rs_threshold = 3)
plot_report_bar(rsa_cm_res_pf, rs_threshold = c(-3, 3), facet_level = TRUE,
                mode = 2)


p1_3 <- plot_KOs_in_pathway(rsa_cm_res_pf,box_color = c("#FFED4DFF",  "#73D8FEFF"), 
                            map_id = "map00220") + scale_y_log10()

p2_3 <- plot_KOs_in_pathway(rsa_cm_res_pf,box_color = c("#FFED4DFF", "#73D8FEFF"), 
                            map_id = "map00010") + scale_y_log10()

p3_3 <- plot_KOs_in_pathway(rsa_cm_res_pf,box_color = c("#FFED4DFF",  "#73D8FEFF"), 
                            map_id = "map00330") + scale_y_log10()

p4_3 <- plot_KOs_in_pathway(rsa_cm_res_pf,box_color = c("#FFED4DFF", "#73D8FEFF"), 
                            map_id = "map00620") + scale_y_log10()

plot_KOs_in_pathway(rsa_cm_res_pf,box_color = c("#FFED4DFF", "#73D8FEFF"), 
                    map_id = "map04360") + scale_y_log10()

plot_KOs_network(rsa_cm_res_pf,
                 map_id = c("map00220", "map00010", "map00330", "map00620"),
                 near_pathway = FALSE, main = "", mark_module = TRUE, 
                 str_width = 70)

# plot_KEGG_map(rsa_cm_res_hm, map_id = "map00220", color_var = "Z_score",
#               save_dir = "../result/Fun_enrich/map/")