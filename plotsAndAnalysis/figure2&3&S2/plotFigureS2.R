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


setwd("D:/academic_relate_code_two/Nessie-main/DeepTX/trainInferSSA")
figure_path = "D:/academic_relate_code_two/Nessie-main/DeepTX/plotsAndAnalysis/figure/"


load_data = function(i) {
    cellNumbers = c(100,500,1000,1500,2000)
    estimated_params_one <- read.table(sprintf("result/synthetic/validation/theo_NN_estimated_param_es2_%d-%d.csv",i,cellNumbers[1]), sep = ',',header=T)
    estimated_params_one$cellNumber = cellNumbers[1]

    estimated_params_two <- read.table(sprintf("result/synthetic/validation/theo_NN_estimated_param_es2_%d-%d.csv",i,cellNumbers[2]), sep = ',',header=T)
    estimated_params_two$cellNumber = cellNumbers[2]

    estimated_params_three <- read.table(sprintf("result/synthetic/validation/theo_NN_estimated_param_es2_%d-%d.csv",i,cellNumbers[3]), sep = ',',header=T)
    estimated_params_three$cellNumber = cellNumbers[3]

    estimated_params_four <- read.table(sprintf("result/synthetic/validation/theo_NN_estimated_param_es2_%d-%d.csv",i,cellNumbers[4]), sep = ',',header=T)
    estimated_params_four$cellNumber = cellNumbers[4]

    estimated_params_five <- read.table(sprintf("result/synthetic/validation/theo_NN_estimated_param_es2_%d-%d.csv",i,cellNumbers[5]), sep = ',',header=T)
    estimated_params_five$cellNumber = cellNumbers[5]
    estimated_params_all <- rbind(estimated_params_one,estimated_params_two,estimated_params_three,estimated_params_four,estimated_params_five)
    return(estimated_params_all)
}


cellNumbers = c(100,500,1000,1500,2000)
theo_NN_estimated_param <- read.table("result/synthetic/theo_NN_estimated_param_GTM.csv", sep = ',',header=T)
valid_loss <- read.table("result/lossVal/valid_loss.csv", sep = ',',header=T)

# BF
legend.text.fontSize = 6
axis.title.fontSize = 6
axis.text.fontSize = 6
color.feedback <- c( "#547CBE", "#3AA438","#E71419")
point.color = color.feedback[2]
PB <- ggplot(theo_NN_estimated_param, aes(x = log10(bf_true),y = log10(bf_es)))+
  geom_point(size = 0.05,color =point.color,fill=point.color) +
  geom_abline(slope = 1, color = '#717070',linetype = "dashed") +
  stat_cor(method = "pearson", size =2) + 
  labs(x = expression(log[10](Theoretical~BF)), y = expression(log[10](Estimated~BF) )) + 
  theme_bw() +
  scale_color_manual(values = color.feedback) +
  xlim(-1,1) +
  ylim(-1,1) +
  theme(title = element_text(colour = 'black', size = 6),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = legend.text.fontSize, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(colour = 'black', size = axis.title.fontSize),
        axis.text = element_text(colour = 'black', size =axis.text.fontSize),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()) 
PB

# BS
PC <- ggplot(theo_NN_estimated_param, aes(x = log10(bs_true),y = log10(bs_es)))+
  geom_point(size = 0.05,color =point.color,fill=point.color) +
  geom_abline(slope = 1, color = '#717070',linetype = "dashed") +
  stat_cor(method = "pearson", size = 2) + 
  labs(x = expression(log[10](Theoretical~BS)), y = expression(log[10](Estimated~BS) )) + 
  theme_bw() +
  scale_color_manual(values = color.feedback) +
  xlim(-1,4) +
  ylim(-1,4) +
  theme(title = element_text(colour = 'black', size = 6),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = legend.text.fontSize, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(colour = 'black', size = axis.title.fontSize),
        axis.text = element_text(colour = 'black', size =axis.text.fontSize),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()) 
PC

PA = ggplot(valid_loss, aes(x="Valid set", y=loss_hel_sp)) + 
  geom_boxplot(utlier.colour = '#D8D8D8', outlier.size = 0, size = 0.15,)+
  labs(x = "Dataset", y = "Helling distance") + 

  theme_bw() +
  scale_color_manual(values = color.feedback) +
  
  ylim(0,0.5) +
  theme(title = element_text(colour = 'black', size = 6),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = legend.text.fontSize, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(colour = 'black', size = axis.title.fontSize),
        axis.text = element_text(colour = 'black', size =axis.text.fontSize),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()) 
PA

i=9
estimated_params_all = load_data(i)
true_params <- read.table("result/synthetic/validation/theo_NN_estimated_param_true.csv", sep = ',',header=T)

PD <- ggplot(estimated_params_all, aes(x = factor(cellNumber, level = cellNumbers), y = log10(bs_es))) +
geom_boxplot(outlier.colour = '#D8D8D8', width=0.3,outlier.size = 0, size = 0.15, colour = '#3E3938', fill = '#EEEEEE', fatten = 1) +
geom_hline(aes(yintercept = log10(true_params$bs_true[i])), colour = '#D95319', linetype = 'dashed', size = 0.2) + 
labs(x = 'Cell number', y = expression(log[10](BS))) +
theme_bw() +
theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
PD

PE <- ggplot(estimated_params_all, aes(x = factor(cellNumber, level = cellNumbers), y = log10(bf_es))) +
geom_boxplot(outlier.colour = '#D8D8D8',  width=0.3,outlier.size = 0, size = 0.15, colour = '#3E3938', fill = '#EEEEEE', fatten = 1) +
geom_hline(aes(yintercept = log10(true_params$bf_true[i])), colour = '#D95319', linetype = 'dashed', size = 0.2) + 
labs(x = 'Cell number', y = expression(log[10](BF))) +
theme_bw() +
theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
PE

FigureS3_1 = ggarrange(PA,PB,PC, ncol = 3, nrow =1,widths=c(1,1.6,1.6),heights=c(2), labels = c("A","B","C"),align = "hv",
                       font.label = list(size = 12,face = "plain"))
FigureS3_1

FigureS3_2 = ggarrange(PD,PE, ncol = 2, nrow =1, align = "hv",labels = c("D","E"), # 添加标签
         font.label = list(size = 12,face = "plain"))
FigureS3_2

FigureS3 = ggarrange(FigureS3_1,FigureS3_2, ncol = 1, nrow =2, align = "v")
FigureS3

ggsave(str_glue(figure_path,'figureS3.pdf'), width = 6.3, height =4, useDingbats = FALSE)


# plot figure2E
df1 <- data.frame(loss_kl_sp = valid_loss$loss_kl_sp)
df1$label = "Distribution + statistics"
df2 <- data.frame(loss_kl_sp = valid_loss$loss_kl_op)
df2$label = "Distribution"
df3 <- data.frame(loss_hel_sp = valid_loss$loss_hel_sp)
df3$label = "Distribution + statistics"
df4 <- data.frame(loss_hel_sp = valid_loss$loss_hel_op)
df4$label =  "Distribution"
loss_valid_kl = rbind(df1,df2)
loss_valid_hel = rbind(df3,df4)

figure2E <- ggplot(loss_valid_hel, aes(x = label,y = loss_hel_sp))+
geom_boxplot(outlier.colour = '#D8D8D8', outlier.size = 0, size = 0.15, fill = c('#EEEEEE','#db2f2c'), fatten = 1,width=0.3) +
  labs(x = "Loss type", y = "Helling distance") + 
  theme_bw() +
  scale_color_manual(values = color.feedback) +
  theme(title = element_text(colour = 'black', size = 6),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = legend.text.fontSize, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(colour = 'black', size = axis.title.fontSize),
        axis.text = element_text(colour = 'black', size =axis.text.fontSize),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()) 
figure2E

ggsave(str_glue(figure_path,'figure2E.jpg'), width = 1.6, height = 1.2)
ggsave(str_glue(figure_path,'figure2E.jpg'), width =1.2, height = 1.2)
