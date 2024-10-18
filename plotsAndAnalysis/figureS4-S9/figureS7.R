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
setwd("D:/academic_relate_code_two/Nessie-main/DeepTX/inferObservedData")
source("utils.R")
gene_estimated_matrix_dmso = read.csv("result/Ten5FU/gene_estimated_stats_matrix.csv",row.names = 1,
                                      header = T, 
                                      sep = ',')
gene_estimated_matrix_dmso$CV = gene_estimated_matrix_dmso$var_true/(gene_estimated_matrix_dmso$mean_true^2)
gene_estimated_matrix_dmso$CV_idu = gene_estimated_matrix_dmso$var_true_idu/(gene_estimated_matrix_dmso$mean_true_idu^2)
dash.color = '#63B8FF'
axis.title.size = 6
point.size = 0.5

# mean
figureS7A = ggplot(gene_estimated_matrix_dmso[filter_gene,], aes(x = mean_true, y =mean_true_idu)) +
  geom_point(color = point.color,shape = 16,size = point.size,alpha = 0.5) +
  geom_abline(slope = 1,color = dash.color,linetype = "dashed") +
  stat_density_2d(
    aes(fill = ..level..,),
    linewidth = 0.15,
    colour = "black",
    bins = 5
  ) +
  scale_size_continuous(range = c(0.1, 1)) +
  labs(x = "Untreatemnts mean",y = expression(paste("10 ", mu,"M 5FU treatments mean")),) +
  theme_bw() +
  scale_y_continuous(breaks = 10^(-100:100),
                     labels = trans_format("log10", math_format(10^.x)),
                     trans="log10",limits = c(0.1, 100) )+
  scale_x_continuous(breaks = 10^(-100:100),
                     labels = trans_format("log10", math_format(10^.x)),
                     trans="log10",limits = c(0.1, 100))+
  theme(
    title = element_text(colour = 'black', size = 6),
    # legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 6, color = "black"),
    legend.key.size = unit(15, "pt"),
    axis.title = element_text(colour = 'black', size = axis.title.size),
    axis.text = element_text(colour = 'black', size = 6),
    axis.ticks = element_line(size = 0.25, lineend = 10),
    panel.grid = element_blank()
  )

figureS7A

# variance
figureS7B = ggplot(gene_estimated_matrix_dmso[filter_gene,], aes(x = var_true, y =var_true_idu)) +
  geom_point(color = point.color,shape = 16,size = point.size,alpha = 0.5) +
  geom_abline(slope = 1,color = dash.color,linetype = "dashed") +
  stat_density_2d(
    aes(fill = ..level..,),
    linewidth = 0.15,
    colour = "black",
    bins = 5
  ) +
  scale_size_continuous(range = c(0.1, 1)) +
  labs(x = "Untreatemnts variance",y = expression(paste("10 ", mu,"M 5FU treatments variance")),) +
  theme_bw() +
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
    axis.title = element_text(colour = 'black', size = axis.title.size),
    axis.text = element_text(colour = 'black', size = 6),
    axis.ticks = element_line(size = 0.25, lineend = 10),
    panel.grid = element_blank()
  )

# CV
figureS7C = ggplot(gene_estimated_matrix_dmso[filter_gene,], aes(x = CV, y =CV_idu)) +
  geom_point(color = point.color,shape = 16,size = point.size,alpha = 0.5) +
  geom_abline(slope = 1,color = dash.color,linetype = "dashed") +
  stat_density_2d(aes(fill = ..level..,),linewidth = 0.15,colour = "black",bins = 5) +
  scale_size_continuous(range = c(0.1, 1)) +
  labs(x = "Untreatments CV",y = expression(paste("10 ", mu,"M 5FU treatments CV"))) +
  theme_bw() +
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
figureS7C

figureS7A_C <-
  ggarrange(figureS7A,figureS7B,figureS7C,ncol = 3,nrow = 1,widths = c(5.6),heights = c(1.8),align = "v",
            labels = c("a","b","c"), 
            font.label = list(size = 12, face = "plain")
  )

#plot busrt density 
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
data_type = "10"
burst_type = "negative"
burst_item = "bf_distance"
burst_items =c("bf_distance","bs_distance")
# negative means bf/bs_fu > bf/bs_untreated
burst_types = c("negative","positive")

gene_estimated_matrix_dmso_sub = filterBurstDistance(burst_items[1], gene_estimated_matrix_dmso, burst_types[1])
gene_estimated_matrix_dmso_sub_bs = filterBurstDistance(burst_items[2], gene_estimated_matrix_dmso, burst_types[2])
gene_estimated_matrix_dmso_sub_bf_pos = filterBurstDistance(burst_items[1], gene_estimated_matrix_dmso, burst_types[2])
gene_estimated_matrix_dmso_sub_bs_neg = filterBurstDistance(burst_items[2], gene_estimated_matrix_dmso, burst_types[1])

filter_gene = row.names(gene_estimated_matrix_dmso)

figureS7E = ggplot(gene_estimated_matrix_dmso[filter_gene,], aes(x = bf, y = bf_idu)) +                                                    
  
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
  labs(x = "Untreatemnts BF",y = expression(paste("10 ", mu,"M 5FU treatments BF")),) +
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
figureS7E

# BS
figureS7F = ggplot(gene_estimated_matrix_dmso[filter_gene,], aes(x = bs, y =bs_idu)) +
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
  labs(x = "Untreatments BS",y =expression(paste("10 ", mu,"M 5FU treatments BS"))) +
  scale_y_continuous(breaks = 10^(-100:100),
                     labels = trans_format("log10", math_format(10^.x)),
                     trans="log10",limits = c(0.1, 1000) )+
  scale_x_continuous(breaks = 10^(-100:100),
                     labels = trans_format("log10", math_format(10^.x)),
                     trans="log10",limits = c(0.1, 1000))+
  theme_bw() +
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
figureS7F

figureS7_EF <-
  ggarrange(figureS7E,figureS7F,ncol = 2,nrow = 1,widths = c(5.6),heights = c(1.8),align = "v",
            labels = c("e","f"), 
            font.label = list(size = 12, face = "plain")
  )
figureS7_EF

figureS7 <-
  ggarrange(figureS7A_C,figureS7_EF,ncol = 1,nrow = 2,widths = c(4.4),heights = c(1.4),align = "v")
  
figureS7
ggsave(sprintf('%sfigureS7.pdf',save_figure_dir),width = 4.4,height =2.8)
