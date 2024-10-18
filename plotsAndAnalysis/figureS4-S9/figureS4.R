# import R package
library("ggplot2")
library("reshape2")
library("ggsci")
library("ggpubr")
library(Matrix)
library(SingleCellExperiment)
library(ggthemes)
library(stringr)
# library(philentropy)
library(scales)
options(warn=-1)

save_figure_dir = "D:/academic_relate_code_two/Nessie-main/DeepTX/plotsAndAnalysis/figure/"
setwd("D:/academic_relate_code_two/Nessie-main/DeepTX/inferObservedData")
source("utils.R")

axis.text.fontSize =  6
axis.ticks.size = 0.25
axis.title.fontSize = 6
plot_distribution <- function(RNA_histogram_DATA,RNA_prob_DATA,col_name, x_label = "mRNA count"){
  color.histogram <- "#D6D8DB"
  color.PF <- "#E84622"
  color.OnOff <- "#00837E"
  p_A1 <- ggplot(RNA_histogram_DATA) + 
    geom_histogram(aes(x = RNA_histogram_DATA[,col_name],y = after_stat(density)), bins = length(RNA_prob_DATA[,"x1"]), na.rm = TRUE, fill = color.histogram) +
    
    geom_line(data = RNA_prob_DATA, aes(x = 0:(length(RNA_prob_DATA[,"x1"])-1), y =x1), colour = color.PF, size = 0.5, linetype = 1) +
    labs(x = x_label, y = "Probability") + 
    theme_bw() + 
    theme(axis.title = element_text(size = 6, color = "black"),
          axis.text = element_text(size = 6, color = "black"),
          axis.ticks = element_line(size = 0.25, lineend = 10),
          panel.grid = element_blank())
  return(p_A1)
}
idu_kl_df = read.csv( "result/IdUDMSO/IdU_kl.csv", row.names =2) 
gene_estimated_matrix_idu = read.csv(file = "result/IdUDMSO/idu_estimated_model_stats_prob.csv", row.names =6)

kl_gene_estimated_matrix_idu = merge(idu_kl_df,gene_estimated_matrix_idu,by.x = 0, by.y = 0)
kl_gene_estimated_matrix_idu[which(kl_gene_estimated_matrix_idu$kl_list <= 0.1),"kl_list"] <-0
kl_gene_estimated_matrix_idu[which(kl_gene_estimated_matrix_idu$kl_list > 0.1),"kl_list"] <-1

point.size = 0.1

# mean
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
color.feedback <- c("#3AA438", "#547CBE", "#E71419")
point.color = color.feedback[2]
P8_A = ggplot(kl_gene_estimated_matrix_idu, aes(x = mean_val, y = mean_true)) +                                               
  geom_point( color = point.color,shape = 16,size = point.size,alpha = 1,) +
  geom_abline(slope = 1,color = 'blue',linetype = "dashed") +
  labs(x = "Mean (DeepTX)",y = "Mean (Data)") +
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
    axis.title = element_text(colour = 'black', size = axis.text.fontSize),
    axis.text = element_text(colour = 'black', size = axis.text.fontSize),
    axis.ticks = element_line(size =  axis.ticks.size, lineend = 10),
    panel.grid = element_blank()
  )
P8_A

# variance
P8_B = ggplot(kl_gene_estimated_matrix_idu, aes(x = var_val, y = var_true)) +
  geom_point(color = point.color,shape = 16,size = point.size, alpha = 1,) +
  geom_abline(slope = 1,color = 'blue', linetype = "dashed") +
  labs(x = "Variance (DeepTX)",y = "Variance (Data)") +
  theme_bw() +
  #   xlim(-1, 2) +
  #   ylim(-1, 2) +
  scale_y_continuous(breaks = 10^(-100:100),
                     labels = trans_format("log10", math_format(10^.x)),
                     trans="log10",limits = c(0.1, 10000) )+
  scale_x_continuous(breaks = 10^(-100:100),
                     labels = trans_format("log10", math_format(10^.x)),
                     trans="log10",limits = c(0.1, 10000))+
  theme(
    title = element_text(colour = 'black', size = 6),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 6, color = "black"),
    legend.key.size = unit(15, "pt"),
    axis.title = element_text(colour = 'black', size = axis.title.fontSize),
    axis.text = element_text(colour = 'black', size = axis.text.fontSize),
    axis.ticks = element_line(size = 0.25, lineend = 10),
    panel.grid = element_blank()
  )
P8_B

P8_C <- ggplot(kl_gene_estimated_matrix_idu, aes(x = mean_true, fill = as.factor(kl_list))) +
  geom_histogram(bins = 30) + 
  scale_fill_manual(values = c( "#71A1C6","#E26463")) +
  
  labs(x = "Mean (Data)", y = "Numbers of gene") +
  scale_x_continuous(breaks = 10^(-100:100),
                     labels = trans_format("log10", math_format(10^.x)),
                     trans="log10",limits = c(0.5, 80))+
  theme_bw()+
  theme(legend.position = "none",
        axis.title = element_text(size = axis.text.fontSize, color = "black"),
        axis.text = element_text(size = axis.text.fontSize, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())

P8_C

RNA_histogram_DATA = read.table("data/IdUDMSO/idu_norm_filter.csv",header=T, sep = ',')
genes_name = names(RNA_histogram_DATA)
i=280
gene_estimated_matrix_idu = read.csv(file = sprintf("result/IdUDMSO/iduSSADistribution/distribution_%s.csv",i))
P8_D = plot_distribution(RNA_histogram_DATA,gene_estimated_matrix_idu,genes_name[i])
P8_D
gene_estimated_matrix_idu = read.csv(file = "result/IdUDMSO/iduSSADistribution/distribution_2.csv")
P8_E = plot_distribution(RNA_histogram_DATA,gene_estimated_matrix_idu,genes_name[2])
gene_estimated_matrix_idu = read.csv(file = "result/IdUDMSO/iduSSADistribution/distribution_8.csv")
P8_F = plot_distribution(RNA_histogram_DATA,gene_estimated_matrix_idu,genes_name[8])
P8_F


dmso_kl_df = read.csv( "result/IdUDMSO/dmso_kl.csv", row.names =2) 
gene_estimated_matrix_dmso = read.csv(file = "result/IdUDMSO/dmso_estimated_model_stats_prob.csv", row.names =6)
kl_gene_estimated_matrix_dmso=merge(dmso_kl_df,gene_estimated_matrix_dmso,by.x = 0, by.y = 0)
kl_gene_estimated_matrix_dmso[which(kl_gene_estimated_matrix_dmso$kl_list <= 0.1),"kl_list"] <-0
kl_gene_estimated_matrix_dmso[which(kl_gene_estimated_matrix_dmso$kl_list > 0.1),"kl_list"] <-1

P9_A = ggplot(kl_gene_estimated_matrix_dmso[kl_gene_estimated_matrix_dmso$kl_list==0,], aes(x = mean_val, y = mean_true)) +                                               
  geom_point(color = point.color,shape = 16,size = point.size, alpha = 1,) +
  geom_abline(slope = 1,color = 'blue', linetype = "dashed") +
  labs(x = "Mean (DeepTX)",y = "Mean (Data)") +
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
    axis.title = element_text(colour = 'black', size = axis.title.fontSize),
    axis.text = element_text(colour = 'black', size = axis.text.fontSize),
    axis.ticks = element_line(size = 0.25, lineend = 10),
    panel.grid = element_blank()
  )
P9_A

# variance
P9_B = ggplot(kl_gene_estimated_matrix_dmso[kl_gene_estimated_matrix_dmso$kl_list==0,], aes(x = var_val, y = var_true)) +
  geom_point(color = point.color,shape = 16,size = point.size, alpha = 1,) +
  geom_abline(slope = 1,color = 'blue', linetype = "dashed") +
  labs(x = "Variance (DeepTX)",y = "Variance (Data)") +
  theme_bw() +
  #   xlim(-1, 2) +
  #   ylim(-1, 2) +
  scale_y_continuous(breaks = 10^(-100:100),
                     labels = trans_format("log10", math_format(10^.x)),
                     trans="log10",limits = c(0.1, 10000) )+
  scale_x_continuous(breaks = 10^(-100:100),
                     labels = trans_format("log10", math_format(10^.x)),
                     trans="log10",limits = c(0.1, 10000))+
  theme(
    title = element_text(colour = 'black', size = 6),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 6, color = "black"),
    legend.key.size = unit(15, "pt"),
    axis.title = element_text(colour = 'black', size = axis.title.fontSize),
    axis.text = element_text(colour = 'black', size = axis.text.fontSize),
    axis.ticks = element_line(size = 0.25, lineend = 10),
    panel.grid = element_blank()
  )
P9_B


P9_C <- ggplot(kl_gene_estimated_matrix_dmso, aes(x = mean_true, fill = as.factor(kl_list))) +
  geom_histogram(bins = 30) + 
  scale_fill_manual(values = c( "#71A1C6","#E26463")) +
  
  labs(x = "Mean (Data)", y = "Numbers of gene") +
  scale_x_continuous(breaks = 10^(-100:100),
                     labels = trans_format("log10", math_format(10^.x)),
                     trans="log10",limits = c(0.5, 80))+
  theme_bw()+
  theme(legend.position = "none",
        axis.title = element_text(size = axis.text.fontSize, color = "black"),
        axis.text = element_text(size = axis.text.fontSize, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())

P9_C

RNA_histogram_DATA = read.table("data/IdUDMSO/dmso_norm_filter.csv",header=T, sep = ',')
genes_name = names(RNA_histogram_DATA)
gene_estimated_matrix_dmso = read.csv(file = "result/IdUDMSO/dmsoSSADistribution/distribution_22.csv")
P9_D = plot_distribution(RNA_histogram_DATA,gene_estimated_matrix_dmso,genes_name[22])
P9_D
gene_estimated_matrix_dmso = read.csv(file = "result/IdUDMSO/dmsoSSADistribution/distribution_35.csv")
P9_E = plot_distribution(RNA_histogram_DATA,gene_estimated_matrix_dmso,genes_name[35])
P9_E
gene_estimated_matrix_dmso = read.csv(file = "result/IdUDMSO/dmsoSSADistribution/distribution_55.csv")
P9_F = plot_distribution(RNA_histogram_DATA,gene_estimated_matrix_dmso,genes_name[55])
P9_F
figure9 <-ggarrange(P8_A,P8_B,P8_C,P8_D,P8_E,P8_F,P9_A,P9_B,P9_C,P9_D,P9_E,P9_F,ncol = 3,nrow = 4,widths = 8,heights = 9,align = "v",
                    labels =c("A","B","C","D","E","F","G","H","I","J","K","L"), font.label = list(size = 12, face = "plain"))
figure9
ggsave(sprintf( '%sfigureS4.pdf',save_figure_dir),width = 6,height = 7.2,useDingbats = FALSE)
