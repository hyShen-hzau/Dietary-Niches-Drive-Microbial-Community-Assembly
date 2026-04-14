library(ggplot2)
library(readr)
library(RColorBrewer)
library(tidyverse)
library(corrplot)
library(psych)
library(purrr)
library(permute)
library(lattice)
library(vegan)
library(ggh4x)
library(ggsci)
library(patchwork)
library(cowplot) 
library(ggNetView)
library(ggnewscale)
library(igraph)
library(agricolae)


#fish-----------------------------------------------------------
Bac <- read.csv("../data/tax.csv",row.names = 1)
group1 <- read.csv("../data/myPlan.csv")

otu_group <- group1 %>% 
  select(Sample = name, Group = Plan3)

var <- colnames(Bac)[1:30]

rownames(Bac) <- paste0("asv", 1:nrow(Bac))

tax <- Bac %>% rownames_to_column("OTUID") %>% 
  select(-var)

dat_all <- Bac[1:30] %>%
  as.data.frame() %>%
  filter(rowSums(. > 0) >= 6) %>% 
  rownames_to_column("ASV") %>%
  mutate(total_abundance = rowSums(select(., -ASV))) %>%
  arrange(desc(total_abundance)) %>%
  # 3. 计算累积占比
  mutate(cumulative_prop = cumsum(total_abundance) / sum(total_abundance)) %>%
  # 4. 筛选 Top 80% 或保留前 N 个 (为了网络稀疏性，通常保留前 500-1000 个 ASV)
  filter(cumulative_prop <= 0.80) %>%
  select(-total_abundance, -cumulative_prop) %>%
  column_to_rownames("ASV")

cat("筛选后保留", nrow(dat_all), "个主要ASV\n")

tax_tab <- dat_all %>% rownames_to_column("OTUID") %>% 
              select(OTUID) %>% 
              left_join(tax)
  

graph_obj <- build_graph_from_mat(
  mat = dat_all,
  transfrom.method = "none",
  method = "WGCNA",
  proc = "Bonferroni",
  r.threshold = 0.8,
  p.threshold = 0.01,
  node_annotation = tax_tab,
  top_modules = 15,
  seed = 1115
)

p <- ggNetView_multi(
  mat = dat_all,
  group_info = otu_group,
  transfrom.method = "none",
  r.threshold = 0.8,
  p.threshold = 0.01,
  method = "WGCNA",
  pointsize = c(0.5,4),
  cor.method = "pearson",
  proc = "BH",
  module.method = "Fast_greedy",
  node_annotation = NULL,
  top_modules = 15,
  center = T,
  shrink = 0.5,
  scale = F,
  jitter = T,
  jitter_sd = 0.3,
  layout = "gephi",
  layout.module = "adjacent",
  seed = 1115,
  nrow = 2 
)

# [1] "C.idellus"
# [1] "H.molitrix"
# [1] "M.salmoides"
# [1] "P.fulvidraco"
# [1] "S.chuatsi"

p

ggsave(filename = "../result/NET_WGCNA/p1.pdf",
       plot = p,
       height = 20,
       width = 21
)


bac4 <- dat_all %>% 
  rownames_to_column("OTUID") %>% 
  left_join(tax[c("OTUID", "Phylum")], by = "OTUID") %>%
  # 移除没有匹配ASV号的行
  filter(!is.na(OTUID)) %>%
  # 使用ASV号作为新行名，移除原genus列
  column_to_rownames("OTUID")


obj1 <- build_graph_from_mat(
  mat = bac4[1:30],
  transfrom.method = "none",
  method = "WGCNA",
  cor.method = "pearson",
  proc = "BH",
  r.threshold = 0.7,
  p.threshold = 0.05,
  node_annotation = tax_final
)

p5 <- ggNetView(
  graph_obj = obj1,
  layout = "gephi",
  layout.module = "adjacent",
  group.by = "Modularity",
  center = F,
  shrink = 0.9,
  linealpha = 0.6,
  linecolor = "#d9d9d9",
  # add_outer = T,
  # label = F
) + 
  scale_size(range = c(1,5)) + 
  scale_color_discrete(breaks = c(1:16, "Others")) 


p6 <- ggNetView(
  graph_obj = obj,
  layout = "gephi",
  layout.module = "adjacent",
  group.by = "Phylum_10",
  center = F,
  shrink = 0.9,
  linealpha = 0.6,
  linecolor = "#d9d9d9"
)

p1 <- ggNetView_multi(
  mat = bac4,
  group_info = otu_group,
  transfrom.method = "none",
  r.threshold = 0.8,
  p.threshold = 0.01,
  method = "WGCNA",
  pointsize = c(0.5,4),
  cor.method = "pearson",
  proc = "BH",
  module.method = "Fast_greedy",
  group.by = "Phylum",
  fill.by = "Phylum",
  node_annotation = NULL,
  top_modules = 15,
  center = T,
  shrink = 0.5,
  scale = F,
  jitter = T,
  jitter_sd = 0.3,
  layout = "gephi",
  layout.module = "adjacent",
  seed = 1115,
  nrow = 2 
)+
  scale_fill_manual(values = rep(c("#8dd3c7", 
                                   "#ffffb3", "#bebada", "#fb8072", "#80b1d3", 
                                   "#fdb462", "#b3de69", "#fccde5", "#cab2d6", 
                                   "#bc80bd", "#ccebc5", "#ffed6f", "#a6cee3", 
                                   "#b2df8a", "#fb9a99", "#bdbdbd", "#a6cee3", 
                                   "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", 
                                   "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", 
                                   "#6a3d9a", "#ffff99", "#b15928"), times = 3))

p1

ggsave(filename = "../result/NET_WGCNA/p2.pdf",
       plot = p,
       height = 20,
       width = 21
)



# ==============================================================================
# 3. 定义全指标计算函数
# ==============================================================================
# 辅助：自然连通度 (Robustness)
get_natural_connectivity <- function(g, weighted = FALSE) {
  if (vcount(g) == 0) return(NA)
  if (weighted && "weight" %in% edge_attr_names(g)) {
    adj <- as_adjacency_matrix(g, attr = "weight", sparse = FALSE)
  } else {
    adj <- as_adjacency_matrix(g, sparse = FALSE)
    adj[adj != 0] <- 1
  }
  evals <- tryCatch(eigen(adj, only.values = TRUE)$values, error = function(e) NA)
  if(all(is.na(evals))) return(NA)
  nc <- log(mean(exp(evals)))
  return(nc)
}
# 辅助：全局效率
get_global_efficiency <- function(g) {
  if (vcount(g) < 2) return(0)
  d_mat <- distances(g)
  d_mat[d_mat == Inf] <- NA
  diag(d_mat) <- NA
  eff <- mean(1/d_mat, na.rm = TRUE)
  return(eff)
}
# 辅助：易损性 (Vulnerability)
get_vulnerability <- function(g) {
  nodes_num <- vcount(g)
  if (nodes_num < 3) return(0)
  E_global <- get_global_efficiency(g)
  if (E_global == 0) return(0)
  top_nodes <- names(sort(degree(g), decreasing = TRUE)[1:min(nodes_num, 3)])
  max_vul <- 0
  for (node in top_nodes) {
    g_rem <- delete_vertices(g, node)
    E_rem <- get_global_efficiency(g_rem)
    v <- (E_global - E_rem) / E_global
    if(v > max_vul) max_vul <- v
  }
  return(max_vul)
}
# --- 主计算函数 ---
calc_topology <- function(sub_g, sample_name) {
  num_nodes <- vcount(sub_g)
  num_edges <- ecount(sub_g)
  if (num_nodes < 3 || num_edges == 0) return(NULL) 
  if ("weight" %in% edge_attr_names(sub_g)) {
    weights <- E(sub_g)$weight
  } else {
    weights <- rep(1, num_edges)
  }
  # --- 1. Network size & geometry ---
  density_val <- edge_density(sub_g)
  avg_path_len <- mean_distance(sub_g, directed = FALSE)
  diam_val <- diameter(sub_g, directed = FALSE)
  # --- 2. Clustering ---
  trans_global <- transitivity(sub_g, type = "global") 
  trans_local <- mean(transitivity(sub_g, type = "local"), na.rm = TRUE) 
  # --- 3. Centrality ---
  avg_degree <- mean(degree(sub_g), na.rm = TRUE)
  avg_betweenness <- mean(betweenness(sub_g, normalized = TRUE), na.rm = TRUE)
  avg_edge_betweenness <- mean(edge_betweenness(sub_g, directed = FALSE), na.rm = TRUE) 
  avg_closeness <- mean(closeness(sub_g, normalized = TRUE), na.rm = TRUE) 
  eig <- tryCatch(mean(eigen_centrality(sub_g)$vector, na.rm=TRUE), error=function(e) NA)
  info_centrality_sum <- sum(closeness(sub_g, normalized = TRUE), na.rm = TRUE)
  # --- 4. Modularity ---
  mod_obj <- cluster_fast_greedy(as.undirected(sub_g))
  modularity_val <- tryCatch(modularity(mod_obj), error=function(e) NA)
  cores <- coreness(sub_g)
  k_core_mean <- mean(cores)
  k_core_max <- max(cores) 
  k_core_min <- min(cores) 
  # --- 5. Efficiency ---
  efficiency_val <- get_global_efficiency(sub_g)
  cohesion_pos <- sum(weights[weights > 0], na.rm = TRUE)
  cohesion_neg <- abs(sum(weights[weights < 0], na.rm = TRUE)) 
  # --- 6. Robustness ---
  robust_unweight <- get_natural_connectivity(sub_g, weighted = FALSE)
  robust_weight <- get_natural_connectivity(sub_g, weighted = TRUE)    
  vulnerability_val <- get_vulnerability(sub_g)
  stability_val <- (1 - vulnerability_val)
  res <- data.frame(
    Sample = sample_name,
    Node = num_nodes,
    Edge = num_edges,
    Density = density_val,
    Distance = avg_path_len,
    Diameter = diam_val,
    Transitivity_global = trans_global,
    Transitivity_local = trans_local,
    Degree = avg_degree,
    Betweenness = avg_betweenness,
    Betweenness_edge = avg_edge_betweenness,
    Closeness = avg_closeness,
    Eigen_centrality = eig,
    Network_info_centrality = info_centrality_sum,
    Modularity = modularity_val,
    K_core_mean = k_core_mean,
    K_core_max = k_core_max,
    K_core_min = k_core_min,
    Network_efficiency = efficiency_val,
    Cohesion_Positive = cohesion_pos,
    Cohesion_Negative = cohesion_neg,
    Robustness_unweight = robust_unweight,
    Robustness_weight = robust_weight,
    Vulnerability = vulnerability_val,
    Stability = stability_val
  )
  return(res)
}


# ==============================================================================
# 4. 遍历样本提取子网络
# ==============================================================================
cat("开始遍历样本计算子网络属性 (指标较多，请耐心等待)...\n")
sample_names <- colnames(dat_all)
otu_binary <- dat_all
otu_binary[otu_binary > 0] <- 1
otu_binary[otu_binary <= 0] <- 0
results_list <- list()
pb <- txtProgressBar(min = 0, max = length(sample_names), style = 3)
for (i in seq_along(sample_names)) {
  s_name <- sample_names[i]
  present_otus <- rownames(otu_binary)[otu_binary[, i] == 1]
  valid_nodes <- intersect(present_otus, V(graph_obj)$name)
  if (length(valid_nodes) > 5) {
    sub_g <- induced_subgraph(graph_obj, valid_nodes)
    stats <- calc_topology(sub_g, s_name)
    if (!is.null(stats)) results_list[[s_name]] <- stats
  }
  setTxtProgressBar(pb, i)
}
close(pb)
final_topology_df <- do.call(rbind, results_list)
final_topology_df[is.na(final_topology_df)] <- 0
final_topology_df <- do.call(data.frame, lapply(final_topology_df, function(x) replace(x, is.infinite(x), 0)))
write.csv(final_topology_df, "../result/NET_WGCNA/SubNetwork_Topology_Properties_Full.csv", row.names = FALSE)
cat("\n子网络属性计算完成。\n")



# ==============================================================================
# 5. 统计分析与绘图 (绘制所有指标)
# ==============================================================================

group_use <- otu_group

cat("开始进行统计分析与绘图...\n")
if(exists("final_topology_df") && nrow(final_topology_df) > 0) {
  # 5.1 合并分组
  plot_data <- final_topology_df %>%
    inner_join(group_use, by = "Sample")
  # 5.2 转换长格式
  metrics_to_test <- colnames(final_topology_df)[-1] 
  plot_data_long <- plot_data %>%
    pivot_longer(cols = all_of(metrics_to_test), names_to = "Metric", values_to = "Value")
  # 5.3 循环统计检验
  stat_res_df <- data.frame()
  group_letters_list <- list()
  for(met in metrics_to_test) {
    sub_dat <- subset(plot_data_long, Metric == met)
    sub_dat <- sub_dat %>% filter(!is.na(Value) & !is.infinite(Value))
    # 统计检验跳过条件
    if(nrow(sub_dat) < 3 || length(unique(sub_dat$Group)) < 2 || var(sub_dat$Value) == 0) {
      next
    }
    letters <- NULL
    method <- NULL
    p_global <- NA
    tryCatch({
      # 检查正态性与方差齐性
      normality_p <- tryCatch(
        min(tapply(sub_dat$Value, sub_dat$Group, function(x) shapiro.test(x)$p.value), na.rm=TRUE),
        error = function(e) 0
      )
      homogeneity_p <- tryCatch(
        leveneTest(Value ~ Group, data = sub_dat)$`Pr(>F)`[1],
        error = function(e) 0
      )
      if(!is.na(normality_p) && !is.na(homogeneity_p) && 
         normality_p > 0.05 && homogeneity_p > 0.05) {
        # ANOVA
        model <- aov(Value ~ Group, data = sub_dat)
        tukey <- HSD.test(model, "Group", console = FALSE)
        if(!is.null(tukey$groups)) {
          letters <- tukey$groups
          method <- "ANOVA"
          p_global <- summary(model)[[1]]$`Pr(>F)`[1]
        }
      } else {
        # Kruskal-Wallis
        kw_test <- kruskal.test(Value ~ Group, data = sub_dat)
        p_global <- kw_test$p.value
        kw_letters <- kruskal(sub_dat$Value, sub_dat$Group, console = FALSE)
        if(!is.null(kw_letters$groups)) {
          letters <- kw_letters$groups
          method <- "Kruskal-Wallis"
        }
      }
      if(!is.null(letters)) {
        letters$Group <- rownames(letters)
        if(is.null(p_global)) p_global <- NA
        stat_res_df <- rbind(stat_res_df, data.frame(Metric = met, Method = method, P_value = p_global))
        group_letters_list[[met]] <- letters
      }
    }, error = function(e) {
      cat("  指标", met, "统计出错:", e$message, "\n")
    })
  }
  if(nrow(stat_res_df) > 0) write.csv(stat_res_df, "SubNetwork_Stats_Result.csv", row.names = FALSE)
  # --- 5.4 批量绘图 (绘制表格中所有指标) ---
  cat("正在绘制箱线图 (所有指标)...\n")
  # 这里不筛选，直接使用全部指标
  plot_list <- list()
  for(met in metrics_to_test) {
    # 提取数据并清洗
    sub_dat <- subset(plot_data_long, Metric == met) %>% 
      filter(!is.na(Value) & !is.infinite(Value))
    # 如果该指标全是NA或没有数据，跳过
    if(nrow(sub_dat) == 0) next
    # 尝试获取显著性字母 (可能为 NULL)
    letters_dat <- group_letters_list[[met]]
    # 计算Y轴范围以便放置字母
    y_max <- max(sub_dat$Value, na.rm = TRUE)
    y_min <- min(sub_dat$Value, na.rm = TRUE)
    range_val <- y_max - y_min
    #if(range_val == 0) range_val <- 1
    p <- ggplot(sub_dat, aes(x = Group, y = Value, fill = Group)) +
      geom_boxplot(alpha = 0.8, outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5, size = 1, color = "grey30") +
      theme_bw() +
      scale_y_continuous(
        limits = c(y_min, y_max * 1.2),  # 增加10%的空间
        expand = expansion(mult = c(0.05, 0.1))  # 调整expand参数
      ) +
      labs(title = met, x = NULL, y = NULL) +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      scale_fill_brewer(palette = "Set2") 
    # 只有当统计成功且有字母时才添加文字
    if(!is.null(letters_dat)) {
      p <- p + geom_text(data = letters_dat,
                         aes(x = Group, y = y_max * 1.1, label = groups),
                         inherit.aes = FALSE, size = 3.5,
                         vjust = -1.2)
    }
    plot_list[[met]] <- p
  }
  # 拼图输出
  if(length(plot_list) > 0){
    # 自动计算布局
    n_plots <- length(plot_list)
    n_cols <- 5 # 固定4列
    n_rows <- ceiling(n_plots / n_cols)
    # 动态计算高度: 每行给 3.5 英寸高度
    total_height <- max(8, n_rows * 3.5)
    final_plot <- wrap_plots(plot_list, ncol = n_cols) + 
      plot_annotation(
        title = "Sub-network Topology Comparison (All Metrics)",
        subtitle = "Significant differences marked by letters (p < 0.05)",
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
      )
    print(final_plot)
    ggsave("../result/NET_WGCNA/SubNetwork_Stats_Boxplot_All_1.pdf", final_plot, width = 16, height = total_height, limitsize = FALSE)
    cat("绘图完成！结果已保存为 'SubNetwork_Stats_Boxplot_All.pdf'。\n")
  } else {
    cat("没有可用的数据用于绘图。\n")
  }
} else {
  cat("错误：未生成有效的子网络数据。\n")
}

if(exists("final_topology_df") && nrow(final_topology_df) > 0) {
  # 5.1 合并分组
  plot_data <- final_topology_df %>%
    inner_join(group_use, by = "Sample")
  # 5.2 转换长格式
  metrics_to_test <- colnames(final_topology_df)[-1] 
  plot_data_long <- plot_data %>%
    pivot_longer(cols = all_of(metrics_to_test), names_to = "Metric", values_to = "Value")
  # 5.3 循环统计检验
  stat_res_df <- data.frame()
  group_letters_list <- list()
  for(met in metrics_to_test) {
    sub_dat <- subset(plot_data_long, Metric == met)
    sub_dat <- sub_dat %>% filter(!is.na(Value) & !is.infinite(Value))
    # 统计检验跳过条件
    if(nrow(sub_dat) < 3 || length(unique(sub_dat$Group)) < 2 || var(sub_dat$Value) == 0) {
      next
    }
    letters <- NULL
    method <- NULL
    p_global <- NA
    tryCatch({
      # 检查正态性与方差齐性
      normality_p <- tryCatch(
        min(tapply(sub_dat$Value, sub_dat$Group, function(x) shapiro.test(x)$p.value), na.rm=TRUE),
        error = function(e) 0
      )
      homogeneity_p <- tryCatch(
        leveneTest(Value ~ Group, data = sub_dat)$`Pr(>F)`[1],
        error = function(e) 0
      )
      if(!is.na(normality_p) && !is.na(homogeneity_p) && 
         normality_p > 0.05 && homogeneity_p > 0.05) {
        # ANOVA
        model <- aov(Value ~ Group, data = sub_dat)
        tukey <- HSD.test(model, "Group", console = FALSE)
        if(!is.null(tukey$groups)) {
          letters <- tukey$groups
          method <- "ANOVA"
          p_global <- summary(model)[[1]]$`Pr(>F)`[1]
        }
      } else {
        # Kruskal-Wallis
        kw_test <- kruskal.test(Value ~ Group, data = sub_dat)
        p_global <- kw_test$p.value
        kw_letters <- kruskal(sub_dat$Value, sub_dat$Group, console = FALSE)
        if(!is.null(kw_letters$groups)) {
          letters <- kw_letters$groups
          method <- "Kruskal-Wallis"
        }
      }
      if(!is.null(letters)) {
        letters$Group <- rownames(letters)
        if(is.null(p_global)) p_global <- NA
        stat_res_df <- rbind(stat_res_df, data.frame(Metric = met, Method = method, P_value = p_global))
        group_letters_list[[met]] <- letters
      }
    }, error = function(e) {
      cat("  指标", met, "统计出错:", e$message, "\n")
    })
  }
  if(nrow(stat_res_df) > 0) write.csv(stat_res_df, "SubNetwork_Stats_Result.csv", row.names = FALSE)
  # --- 5.4 批量绘图 (绘制表格中所有指标) ---
  cat("正在绘制箱线图 (所有指标)...\n")
  # 这里不筛选，直接使用全部指标
  plot_list <- list()
  for(met in metrics_to_test) {
    # 提取数据并清洗
    sub_dat <- subset(plot_data_long, Metric == met) %>% 
      filter(!is.na(Value) & !is.infinite(Value))
    
    # 如果该指标全是NA或没有数据，跳过
    if(nrow(sub_dat) == 0) next
    
    # 尝试获取显著性字母 (可能为 NULL)
    letters_dat <- group_letters_list[[met]]
    
    # 计算Y轴范围
    y_max <- max(sub_dat$Value, na.rm = TRUE)
    y_min <- min(sub_dat$Value, na.rm = TRUE)
    range_val <- y_max - y_min
    
    # 如果有显著性字母，计算每个组的数据最大值
    if(!is.null(letters_dat)) {
      # 计算每个组的最大值
      group_max <- sub_dat %>%
        group_by(Group) %>%
        summarise(max_val = max(Value, na.rm = TRUE)) %>%
        ungroup()
      
      # 将每个组的最大值合并到letters_dat
      letters_dat <- letters_dat %>%
        left_join(group_max, by = "Group")
      
      # 计算字母标记的位置：每个组的最大值 + 偏移量
      # 偏移量设为总范围的5%
      offset <- range_val * 0.05
      # 如果数据范围很小，使用固定偏移
      if(offset < 0.01 * y_max) offset <- 0.1 * y_max
      if(offset == 0) offset <- 0.1  # 防止range_val为0的情况
      
      letters_dat$y_pos <- letters_dat$max_val + offset
    }
    
    # 计算绘图Y轴的上限
    # 如果有显著性标记，确保Y轴上限包含标记
    if(!is.null(letters_dat)) {
      y_upper <- max(c(sub_dat$Value, letters_dat$y_pos), na.rm = TRUE)
    } else {
      y_upper <- y_max
    }
    
    # 创建绘图
    p <- ggplot(sub_dat, aes(x = Group, y = Value, fill = Group)) +
      geom_boxplot(alpha = 0.8, outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5, size = 1, color = "grey30") +
      theme_bw() +
      # 设置Y轴范围，为显著性标记留出空间
      scale_y_continuous(
        limits = c(y_min, y_upper * 1.1),  # 增加10%的空间
        expand = expansion(mult = c(0.05, 0.1))  # 调整expand参数
      ) +
      labs(title = met, x = NULL, y = NULL) +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      scale_fill_brewer(palette = "Set2")
    
    # 只有当统计成功且有字母时才添加文字
    if(!is.null(letters_dat)) {
      p <- p + geom_text(
        data = letters_dat,
        aes(x = Group, y = y_pos, label = groups),
        inherit.aes = FALSE,
        size = 3.5,
        vjust = 0  # 改为0，因为y_pos已经是计算好的位置
      )
    }
    
    plot_list[[met]] <- p
  }
  # 拼图输出
  if(length(plot_list) > 0){
    # 自动计算布局
    n_plots <- length(plot_list)
    n_cols <- 5 # 固定4列
    n_rows <- ceiling(n_plots / n_cols)
    # 动态计算高度: 每行给 3.5 英寸高度
    total_height <- max(8, n_rows * 3.5)
    final_plot <- wrap_plots(plot_list, ncol = n_cols) + 
      plot_annotation(
        title = "Sub-network Topology Comparison (All Metrics)",
        subtitle = "Significant differences marked by letters (p < 0.05)",
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
      )
    print(final_plot)
    ggsave("../result/NET_WGCNA/SubNetwork_Stats_Boxplot_All_1.pdf", final_plot, width = 16, height = total_height, limitsize = FALSE)
    cat("绘图完成！结果已保存为 'SubNetwork_Stats_Boxplot_All.pdf'。\n")
  } else {
    cat("没有可用的数据用于绘图。\n")
  }
} else {
  cat("错误：未生成有效的子网络数据。\n")
}




