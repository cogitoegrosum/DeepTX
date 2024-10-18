# import R package
options(warn = -1)
library("ggplot2")
library("reshape2")
library("ggsci")
library("ggpubr")
library(Matrix)
library(SingleCellExperiment)
library(ggthemes)
library(stringr)
library(scales)

save_figure_dir = "D:/academic_relate_code_two/Nessie-main/DeepTX/plotsAndAnalysis/figure/"
setwd("D:/academic_relate_code_two/Nessie-main/DeepGTM/analyseScRNAdata/Fifty5FU")
source("utils.R")
gene_estimated_matrix_dmso = read.csv("result/gene_estimated_stats_matrix.csv",row.names = 1,
                                      header = T, 
                                      sep = ',')
gene_estimated_matrix_dmso$CV = gene_estimated_matrix_dmso$var_true/(gene_estimated_matrix_dmso$mean_true^2)
gene_estimated_matrix_dmso$CV_idu = gene_estimated_matrix_dmso$var_true_idu/(gene_estimated_matrix_dmso$mean_true_idu^2)

dash.color = '#63B8FF'
axis.title.size = 6
point.size = 0.5

# mean
figureS9A = ggplot(gene_estimated_matrix_dmso[filter_gene,], aes(x = mean_true, y =mean_true_idu)) +
  geom_point(color = point.color,shape = 16,size = point.size,alpha = 0.5) +
  geom_abline(slope = 1,color =dash.color,linetype = "dashed") +
  stat_density_2d(
    aes(fill = ..level..,),
    linewidth = 0.15,
    colour = "black",
    bins = 5
  ) +
  scale_size_continuous(range = c(0.1, 1)) +
  labs(x = "Untreatments mean",y = expression(paste("50 ", mu,"M 5FU treatments mean"))) +
  theme_bw() +
  #   xlim(-1, 2) +
  #   ylim(-1, 2) +
  scale_y_continuous(breaks = 10^(-100:100),
                     labels = trans_format("log10", math_format(10^.x)),
                     trans="log10",limits = c(0.1, 100) )+
  scale_x_continuous(breaks = 10^(-100:100),
                     labels = trans_format("log10", math_format(10^.x)),
                     trans="log10",limits = c(0.1, 100))+
  theme(
    title = element_text(colour = 'black', size = 6),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 6, color = "black"),
    legend.key.size = unit(15, "pt"),
    axis.title = element_text(colour = 'black', size = 6),
    axis.text = element_text(colour = 'black', size = 6),
    axis.ticks = element_line(size = 0.25, lineend = 10),
    panel.grid = element_blank()
  )


# variance
figureS9B = ggplot(gene_estimated_matrix_dmso[filter_gene,], aes(x = var_true, y =var_true_idu)) +
  geom_point(color = point.color,shape = 16,size = point.size,alpha = 0.5) +
  geom_abline(slope = 1,color = dash.color,linetype = "dashed") +
  stat_density_2d(aes(fill = ..level..,),linewidth = 0.15,colour = "black",bins = 5) +
  scale_size_continuous(range = c(0.1, 1)) +
  labs(x = "Untreatments variance",y = expression(paste("50 ", mu,"M 5FU treatments variance"))) +
  theme_bw() +
  #   xlim(-1, 2) +
  #   ylim(-1, 2) +
  scale_y_continuous(breaks = 10^(-100:100),
                     labels = trans_format("log10", math_format(10^.x)),
                     trans="log10",limits = c(0.1, 100) )+
  scale_x_continuous(breaks = 10^(-100:100),
                     labels = trans_format("log10", math_format(10^.x)),
                     trans="log10",limits = c(0.1, 100))+
  theme(
    title = element_text(colour = 'black', size = 6),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 6, color = "black"),
    legend.key.size = unit(15, "pt"),
    axis.title = element_text(colour = 'black', size = 6),
    axis.text = element_text(colour = 'black', size = 6),
    axis.ticks = element_line(size = 0.25, lineend = 10),
    panel.grid = element_blank()
  )
# PF_E

# CV
figureS9C = ggplot(gene_estimated_matrix_dmso[filter_gene,], aes(x = CV, y =CV_idu)) +
  geom_point(color = point.color,shape = 16,size = point.size,alpha = 0.5) +
  geom_abline(slope = 1,color = dash.color,linetype = "dashed") +
  stat_density_2d(aes(fill = ..level..,),linewidth = 0.15,colour = "black",bins = 5) +
  scale_size_continuous(range = c(0.1, 1)) +
  labs(x = "Untreatments CV",y = expression(paste("50 ", mu,"M 5FU treatments CV"))) +
  theme_bw() +
  #   xlim(-1, 2) +
  #   ylim(-1, 2) +
  scale_y_continuous(breaks = 10^(-100:100),
                     labels = trans_format("log10", math_format(10^.x)),
                     trans="log10",limits = c(0.01, 10) )+
  scale_x_continuous(breaks = 10^(-100:100),
                     labels = trans_format("log10", math_format(10^.x)),
                     trans="log10",limits = c(0.01, 10))+
  theme(
    title = element_text(colour = 'black', size = 6),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 6, color = "black"),
    legend.key.size = unit(15, "pt"),
    axis.title = element_text(colour = 'black', size = 6),
    axis.text = element_text(colour = 'black', size = 6),
    axis.ticks = element_line(size = 0.25, lineend = 10),
    panel.grid = element_blank()
  )
figureS9C

figureS9A_C <-
  ggarrange(figureS9A,figureS9B,figureS9C,ncol = 3,nrow = 1,widths = c(5.6),heights = c(1.8),align = "v",
            labels = c("a","b","c"), # 添加标签
            font.label = list(size = 12, face = "plain")
  )
#plot busrt density 

library(RColorBrewer)
colormap<- rev(brewer.pal(9,"Blues")[1:6])

theme_pubr_new = function() {
  theme_pubr(base_family = "Helvetica",
             base_size = 14,
             legend = "right") %+replace%
    theme(
      axis.title = element_text(face = "bold", size = 18),
      legend.text = element_text(face = "bold", size = 12),
      axis.text = element_text(size = 16)
    )
}
color.feedback <- c("#3AA438", "#547CBE", "#E71419","#8abd44","#db2f2c")
point.color = color.feedback[2]
data_type = "50"
burst_type = "negative"
burst_item = "bf_distance"

burst_items =c("bf_distance","bs_distance")
burst_types = c("negative","positive")
gene_estimated_matrix_dmso = read.csv("result/gene_estimated_stats_matrix.csv",row.names = 1,
                                      header = T, 
                                      sep = ',')
gene_estimated_matrix_dmso$CV = gene_estimated_matrix_dmso$var_true/(gene_estimated_matrix_dmso$mean_true^2)
gene_estimated_matrix_dmso$CV_idu = gene_estimated_matrix_dmso$var_true_idu/(gene_estimated_matrix_dmso$mean_true_idu^2)

gene_estimated_matrix_dmso_sub = filterBurstDistance(burst_items[1], gene_estimated_matrix_dmso, burst_types[1])
gene_estimated_matrix_dmso_sub_bs = filterBurstDistance(burst_items[2], gene_estimated_matrix_dmso, burst_types[2])
gene_estimated_matrix_dmso_sub_bf_pos = filterBurstDistance(burst_items[1], gene_estimated_matrix_dmso, burst_types[2])
gene_estimated_matrix_dmso_sub_bs_neg = filterBurstDistance(burst_items[2], gene_estimated_matrix_dmso, burst_types[1])

filter_gene = row.names(gene_estimated_matrix_dmso)
dash.color = '#63B8FF'
axis.title.size = 6
point.size = 0.5

figureS9E = ggplot(gene_estimated_matrix_dmso[filter_gene,], aes(x = bf, y = bf_idu)) +                                                    
  
  stat_density2d(
    geom ="raster",
    aes(fill = ..density..,),
    contour = F
  ) +
  stat_density2d(linewidth = 0.15,colour = "grey")+
  scale_fill_gradientn(colours=rev(colormap))+
  geom_point(data = gene_estimated_matrix_dmso_sub,aes(x =bf, y = bf_idu),
             color = color.feedback[4],shape = 16,size = 0.3,alpha = 1) +
  geom_point(data = gene_estimated_matrix_dmso_sub_bf_pos,aes(x =bf, y = bf_idu),
             color = color.feedback[5],shape = 16,size = 0.3,alpha = 1) +
  
  geom_abline(slope = 1,color = dash.color,linetype = "dashed") +
  labs(x = " BF(Untreatment)",y = expression(paste("BF(50 ", mu,"M 5FU treatment)")),) +
  scale_y_continuous(breaks = 10^(-100:100),
                     labels = trans_format("log10", math_format(10^.x)),
                     trans="log10",limits = c(0.001, 10) )+
  scale_x_continuous(breaks = 10^(-100:100),
                     labels = trans_format("log10", math_format(10^.x)),
                     trans="log10",limits = c(0.001, 10))+
  theme_bw() +
  theme(
    title = element_text(colour = 'black', size = 6),
    # legend.position = c(0.9,0.3),
    # legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 6, color = "black"),
    legend.key.size = unit(6, "pt"),
    axis.title = element_text(colour = 'black', size = 6),
    axis.text = element_text(colour = 'black', size = 6),
    axis.ticks = element_line(size = 0.25, lineend = 10),
    panel.grid = element_blank()
  )
figureS9E

# BS
figureS9F = ggplot(gene_estimated_matrix_dmso[filter_gene,], aes(x = bs, y =bs_idu)) +
  stat_density2d(
    geom ="raster",
    aes(fill = ..density..,),
    contour = F
  ) +
  stat_density2d(linewidth = 0.15,colour = "grey")+
  scale_fill_gradientn(colours=rev(colormap))+
  geom_point(data = gene_estimated_matrix_dmso_sub_bs,aes(x =bs, y = bs_idu),color = color.feedback[5],shape = 16,size = 0.3,alpha = 1) +
  geom_point(data = gene_estimated_matrix_dmso_sub_bs_neg,aes(x =bs, y = bs_idu),color = color.feedback[4],shape = 16,size = 0.3,alpha = 1) +
  geom_abline(slope = 1,color = dash.color,linetype = "dashed") +
  labs(x = " BS(Untreatemnt)",y = expression(paste("BS(50 ", mu,"M 5FU treatment)")),) +
  theme_bw() +
  scale_y_continuous(breaks = 10^(-100:100),
                     labels = trans_format("log10", math_format(10^.x)),
                     trans="log10",limits = c(0.1, 1000) )+
  scale_x_continuous(breaks = 10^(-100:100),
                     labels = trans_format("log10", math_format(10^.x)),
                     trans="log10",limits = c(0.1, 1000))+
  theme(
    title = element_text(colour = 'black', size = 6),
    # legend.position = "none",
    # legend.position = c(0.9,0.3),
    legend.title = element_blank(),
    legend.text = element_text(size = 6, color = "black"),
    legend.key.size = unit(6, "pt"),
    axis.title = element_text(colour = 'black', size = 6),
    axis.text = element_text(colour = 'black', size = 6),
    axis.ticks = element_line(size = 0.25, lineend = 10),
    panel.grid = element_blank()
  )
figureS9F

figureS9_EF <-
  ggarrange(figureS9E,figureS9F,ncol = 2,nrow = 1,widths = c(5.6),heights = c(1.8),align = "v",
            labels = c("e","f"), 
            font.label = list(size = 12, face = "plain")
  )
figureS9_EF

figureS9 <-
  ggarrange(figureS9A_C,figureS9_EF,ncol = 1,nrow = 2,widths = c(4.4),heights = c(1.4),align = "v")
figureS9
ggsave(sprintf('%sfigureS9.pdf',save_figure_dir),width = 4.4,height =2.8)
