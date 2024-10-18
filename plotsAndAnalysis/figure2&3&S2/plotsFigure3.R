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

library(RColorBrewer)
colormap<- brewer.pal(9,"Blues")[2:6]
colors = c("#D2E0FF","#7BA4FF","#D2E0FF")

# colors = c("#89A1CF","#bac7e4","#F2F3F9")
# load the estimated data
theo_NN_estimated_param <- read.table("result/synthetic/theo_NN_estimated_param_GTM.csv", sep = ',',header=T)
color.feedback <- c( "#547CBE", "#3AA438","#E71419")
point.color = color.feedback[2]
theo_NN_estimated_param$bs_bf_dist = abs(log(theo_NN_estimated_param$bs_true /theo_NN_estimated_param$bs_es)) + abs(log(theo_NN_estimated_param$bf_true / theo_NN_estimated_param$bf_es))

PC_A = ggplot(theo_NN_estimated_param, aes(x=bs_true, y=bf_true) ) +

geom_point(aes(color = log10(bs_bf_dist ) ),size =0.01,alpha=1)+

  scale_colour_gradient2( 
                          # low = "#F902FF" ,mid = "#00DBDE",high = muted("red"),
                          # low = colors[1] ,mid = colors[3],high = colors[5],
                          low = colors[1] ,mid = colors[2],high = colors[3],

                         midpoint = mean(theo_NN_estimated_param$bs_bf_dist))+
stat_density_2d(aes(fill =..level..,), linewidth = 0.15,bins=5)+

scale_size_continuous(range = c(0.1,1))+
labs(x ="Theor. BS", y = "Theor. BF") + 

 theme_bw()+
  scale_y_continuous(breaks = 10^(-100:100),
                       labels = trans_format("log10", math_format(10^.x)),
                       trans="log10")+
  scale_x_continuous(breaks = 10^(-100:100),
                       labels = trans_format("log10", math_format(10^.x)),
                       trans="log10")+
  theme(
        legend.position = "none",
        # legend.position = c(.6,.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(6, "pt"),
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()
       ) 

PC_A

ggsave('figure/PC_A.pdf', width = 2.1, height = 2, useDingbats = FALSE)

theo_NN_stats <- read.table("result/theo_NN_stats.csv", sep = ',',header=T)
density.color = "darkblue"
vline.color = "green"
trueSolution = read.csv(file = "result/synthetic/posteriorDist/trueSolution.csv")
i=35
estimated_file_path = sprintf("result/synthetic/posteriorDist/AllSolution_%s.csv", i )
estimatedSolution = read.csv(file =estimated_file_path)
density.color = "#AFCCA0"
vline.color = "#AAB7D9"
P1 = ggplot(estimatedSolution, aes(x=x1)) +
geom_density(aes(x = x1, y = after_stat(density),weight=loss),colour=density.color ) +
geom_vline(xintercept = trueSolution[i,c("x1")],colour=vline.color, linetype = "longdash")+
  labs(x = expression(r[on]), y = "Posteriori") +

xlim(0.6,15)+
theme_bw() +
  theme(title = element_text(colour = 'black', size = 6),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()
       )
P2 = ggplot(estimatedSolution) +
geom_density(aes(x = x2, y = after_stat(density),weight=loss),colour=density.color ) +
geom_vline(xintercept = trueSolution[i,c("x2")],colour=vline.color, linetype = "longdash")+
  labs(x = expression(k[on]), y =  NULL) +
xlim(0.6,15)+
theme_bw() +
  theme(title = element_text(colour = 'black', size = 6),
  # plot.margin =margin(t = 0, r = 0, b = 0, l = -2,  unit = "pt"),
#         legend.position = "none",
#         legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title.x = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()
       )
P2
P3 = ggplot(estimatedSolution) +
geom_density(aes(x = x3, y = after_stat(density),weight=loss),colour=density.color ) +
geom_vline(xintercept = trueSolution[i,c("x3")],colour=vline.color, linetype = "longdash")+
  labs(x = expression(r[off]), y = NULL) +

xlim(6,15)+
theme_bw() +
  theme(title = element_text(colour = 'black', size = 6),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(colour = 'black', size = 6),
        axis.title.y = element_text(colour = 'black', size = 1),

        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()
       )

P4 = ggplot(estimatedSolution) +
geom_density(aes(x = x4, y = after_stat(density),weight=loss),colour=density.color ) +
geom_vline(xintercept = trueSolution[i,c("x4")],colour=vline.color, linetype = "longdash")+
  labs(x = expression(k[off]), y =NULL) +

xlim(6,11)+
theme_bw() +
  theme(title = element_text(colour = 'black', size = 6),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(colour = 'black', size = 6),
        axis.title.y = element_text(colour = 'black', size = 1),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()
       )

P5 = ggplot(estimatedSolution) +
geom_density(aes(x = x5, y = after_stat(density),weight=loss),colour=density.color ) +
geom_vline(xintercept = trueSolution[i,c("x5")],colour=vline.color, linetype = "longdash")+
  labs(x = expression(mu), y =NULL) +
xlim(125,175)+
theme_bw() +
  theme(title = element_text(colour = 'black', size = 6),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(colour = 'black', size = 6),
        axis.title.y = element_text(colour = 'black', size = 1),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()
       )
P5
P6 = ggplot(estimatedSolution) +
geom_density(aes(x = bs, y = after_stat(density),weight=loss),colour=density.color ) +
geom_vline(xintercept = trueSolution[i,c("bs")],colour=vline.color, linetype = "longdash")+
  labs(x = expression(mu), y = "Posteriori") +

xlim(200,260)+
theme_bw() +
  theme(title = element_text(colour = 'black', size = 6),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()
       )
P6

P7 = ggplot(estimatedSolution) +
geom_density(aes(x = bf, y = after_stat(density),weight=loss),colour=density.color ) +
geom_vline(xintercept = trueSolution[i,c("bf")],colour=vline.color, linetype = "longdash")+
  labs(x = expression(mu), y = "Posteriori") +

xlim(0.3,0.4)+
theme_bw() +
  theme(title = element_text(colour = 'black', size = 6),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()
       )
P7
picD_sub_one <- ggarrange(P1,P2,P3,P4,P5, ncol = 5, nrow = 1, widths=c(10), heights=c(2), align = "v")
picD_sub_one

picD_sub_two<- ggarrange(P6,P7, ncol = 2, nrow = 1, widths=c(6.3), heights=c(2), align = "v")
picD <- ggarrange(picD_sub_one,picD_sub_two, ncol = 1, nrow = 2, widths=c(6.3), heights=c(1.4,2), align = "v")
# picB_right <- ggarrange(picB_r_u,PB_3, ncol = 1, nrow = 2, widths=c(6.3), heights=c(4), align = "v")
picD

figure_C = ggarrange(PC_A,picD,ncol = 2, nrow = 1, widths=c(2.1,4.2), heights=c(3.6), align = "v")
figure_C
ggsave(str_glue(figure_path,'figure3.pdf'), width = 6.3, height = 2, useDingbats = FALSE)
