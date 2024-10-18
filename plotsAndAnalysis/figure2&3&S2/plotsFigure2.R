# import R package
options(warn=-1)
library("ggplot2")
library("reshape2")
library("ggsci")
library("ggpubr")
library(Matrix)
library(SingleCellExperiment)
library(ggthemes)
library(stringr)
library(scales)
setwd("D:/academic_relate_code_two/Nessie-main/DeepTX/trainInferSSA/result/")
figure_path = "D:/academic_relate_code_two/Nessie-main/DeepTX/plotsAndAnalysis/figure/"

theo_NN_stats <- read.table("synthetic/theo_NN_stats.csv", sep = ',',header=T)

# color.feedback <- c("#3AA438", "#547CBE", "#E71419")
color.feedback <- c("#3AA438", "#547CBE", "#E54551")

point_color_sta = color.feedback[3]
color_abline = "#77C5E0"
point_size_sta = 0.001
axis.title.size=7
line.width = 0.2

# mean
max(theo_NN_stats$m_true)
min(theo_NN_stats$m_true)
p_B <- ggplot(theo_NN_stats, aes(x = m_NN,y = m_true))+
  geom_point(size = point_size_sta,color =point_color_sta) +
  geom_abline(slope = 1, color =color_abline,linetype = "dashed",lwd=line.width) +
  stat_cor(method = "pearson", size = 2) + 
  # labs(x = "DeepTX mean", y = "Theor. mean") + 
  labs(x = "DeepTX", y = "Theoretical") + 
  theme_bw() +
  scale_color_manual(values = color.feedback) +
  xlim(0,400) +
  ylim(0,400) +
  theme(title = element_text(colour = 'black', size = 6),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(colour = 'black', size = axis.title.size),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()) 
p_B

# fano_factor
color.feedback <- c("#3AA438", "#547CBE", "#E71419")
max(log10(theo_NN_stats$fano_NN))
min(log10(theo_NN_stats$fano_NN))

p_C <- ggplot(theo_NN_stats, aes(x = log10(fano_NN),y = log10(fano_true)))+
# p_C <- ggplot(theo_NN_stats, aes(x = fano_NN,y = fano_true))+
  geom_point(size = point_size_sta,color =point_color_sta) +
  geom_abline(slope = 1, color = color_abline,linetype = "dashed",lwd=line.width) +
  stat_cor(method = "pearson", size = 2) + 
#   labs(x = expression("Log"["10"]*"(DeepTx Fano factor)"), y =expression("Log"["10"]*"(Theo Fano factor)")) + 
  # labs(x = expression("Log"["10"]*"(DeepTX FF)"), y =expression("Log"["10"]*"(Theor. FF)")) + 
  labs(x = "DeepTX", y = "Theoretical") + 
  theme_bw() +
  scale_color_manual(values = color.feedback) +
  xlim(0,2.5) +
  ylim(0,2.5) +
  theme(title = element_text(colour = 'black', size = 6),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(colour = 'black', size = axis.title.size),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank())
p_C

# bimodalcoeff
color.feedback <- c("#3AA438", "#547CBE", "#E71419")
max(log10(theo_NN_stats$bimodcoeff_NN))
min(log10(theo_NN_stats$bimodcoeff_NN))
p_D <- ggplot(theo_NN_stats, aes(x = log10(bimodcoeff_NN),y = log10(bimodcoeff_true)))+
  geom_point(size = point_size_sta,color =point_color_sta) +
  geom_abline(slope = 1, color = color_abline,linetype = "dashed",lwd=line.width) +
  stat_cor(method = "pearson", size = 2) + 
#   labs(x = expression(Log[10](DeepTx~bimodal_coeff)), y = expression(Log[10](Theo~bimodal_coeff))) + 
  # labs(x = expression(Log[10](DeepTX~BC)), y = expression(Log[10](Theor.~BC))) + 
  labs(x = "DeepTX", y = "Theoretical") + 
  theme_bw() +
  scale_color_manual(values = color.feedback) +
  xlim(-0.7,0.1) +
  ylim(-0.7,0.1) +
  theme(title = element_text(colour = 'black', size = 6),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(colour = 'black', size = axis.title.size),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()) 
p_D
picB_r_u <- ggarrange(p_B, p_C, p_D, ncol = 3, nrow = 1, widths=c(5.1), heights=c(3), align = "v")
picB_r_u


# Compare the consistency of the simulated distribution of the test data and the model predicted distribution
RNA_histogram_DATA = read.table("synthetic/distributionResultFig/ssa_result.csv",header=T, sep = ',')
RNA_prob_DATA=read.table("synthetic/distributionResultFig/NN_predicted_10.csv",header=T, sep = ',')
gene_names= names(RNA_histogram_DATA)
plot_distribution <- function(RNA_histogram_DATA,RNA_prob_DATA,col_name, x_label = "mRNA number"){
    color.histogram <- "#D6D8DB"
    color.PF <- "#E84622"
    color.OnOff <- "#00837E"
    p_A1 <- ggplot(RNA_histogram_DATA) + 
      geom_histogram(aes(x = RNA_histogram_DATA[,col_name],y = after_stat(density)), bins = length(RNA_prob_DATA[,"x1"]), na.rm = TRUE, fill = color.histogram) +

      geom_line(data = RNA_prob_DATA, aes(x = 0:(length(RNA_prob_DATA[,"x1"])-1), y =x1), colour = color.PF, size = 0.5, linetype = 1) +
      labs(x = x_label, y = "Probability") + 
      theme_bw() + 
      theme(axis.title = element_text(size = 7, color = "black"),
            axis.text = element_text(size = 6, color = "black"),
            axis.ticks = element_line(size = 0.25, lineend = 10),
            panel.grid = element_blank())
    return(p_A1)
}

# valid_param_list = [10,21,72,106,114,291]
i=21
NN_predicted_file_name = str_glue("synthetic/distributionResultFig/","NN_predicted_{i}.csv")
RNA_prob_DATA=read.table(NN_predicted_file_name,header=T, sep = ',')
p_A1 = plot_distribution(RNA_histogram_DATA,RNA_prob_DATA,gene_names[2])
p_A1

i=106
NN_predicted_file_name = str_glue("synthetic/distributionResultFig/","NN_predicted_{i}.csv")
RNA_prob_DATA=read.table(NN_predicted_file_name,header=T, sep = ',')
p_A2 = plot_distribution(RNA_histogram_DATA,RNA_prob_DATA,gene_names[4])
p_A2

i=291
NN_predicted_file_name = str_glue("synthetic/distributionResultFig/","NN_predicted_{i}.csv")
RNA_prob_DATA=read.table(NN_predicted_file_name,header=T, sep = ',')
p_A3 = plot_distribution(RNA_histogram_DATA,RNA_prob_DATA,gene_names[6])
p_A3

i=72
NN_predicted_file_name = str_glue("synthetic/distributionResultFig/","NN_predicted_{i}.csv")
RNA_prob_DATA=read.table(NN_predicted_file_name,header=T, sep = ',')
p_A4 = plot_distribution(RNA_histogram_DATA,RNA_prob_DATA,gene_names[3])
p_A4

picB_left <- ggarrange(p_A1, p_A2, p_A3, p_A4, ncol = 2, nrow = 2, widths=c(5), heights=c(4), align = "v")

# draw the loss Figure
loss_only_prob = read.csv(file = "lossVal/train_only_prob.csv")
loss_stats_prob = read.csv(file = "lossVal/train_stats_prob.csv")
loss_only_prob$x_value =  1:dim(loss_only_prob)[1]
loss_stats_prob$x_value =  1:dim(loss_stats_prob)[1]
loss_only_prob$x2 = loss_stats_prob[1:nrow(loss_only_prob),c("x1")]

# loss_only_prob$x3 = loss_only_stats[1:nrow(loss_only_prob),c("x1")]
PB_3 = ggplot(data = loss_only_prob[2:nrow(loss_only_prob),], mapping = aes(x = x_value, y = x1)) +
geom_line()+
# geom_text(aes(label = "label"), hjust = -0.2, vjust = 0.5)+
geom_line(aes(x = x_value, y = x2),colour="red")+
xlim(2,200)+
labs(x = "Epoch", y = "Loss") +
  theme_bw() + 
  theme(title = element_text(colour = 'black', size = 6),
      legend.position = "none",
      legend.title = element_blank(),
      legend.text = element_text(size = 6, color="black"),
      legend.key.size = unit(15, "pt"),
      axis.title = element_text(colour = 'black', size = axis.title.size),
      axis.text = element_text(colour = 'black', size = 6),
      axis.ticks = element_line(size=0.25, lineend = 10),
      panel.grid = element_blank())
PB_3

picB_right <- ggarrange(picB_r_u,PB_3, ncol = 1, nrow = 2, widths=c(6.3), heights=c(1.9,2.1), align = "v")
picB_right

picB = ggarrange(picB_left,picB_right, ncol = 2, nrow =1,widths=c(4.5,5.5), heights=c(4,4), align = "v")
picB
ggsave(str_glue(figure_path,'figure2.pdf'), width = 8, height = 2.5, useDingbats = FALSE)
