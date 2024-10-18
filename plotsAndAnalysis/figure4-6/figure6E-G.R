library(Seurat)
library(dplyr)
library(cowplot)
library(GSVA)
library(data.table)
library(ggplot2)
library(umap)
library(ggpubr)
library(ggsci)

setwd("D:/academic_relate_code_two/Nessie-main/DeepTX/inferObservedData")
save_figure_dir = "D:/academic_relate_code_two/Nessie-main/DeepTX/plotsAndAnalysis/figure/"

anno <-read.csv("data/5FUcommon/GSE149224_meta.information.csv", header = TRUE)
# The data set can be downloaded from link https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149224
expr.data <- read.table("data/5FUcommon/GSE149224_RSH.all.txt", header = TRUE,sep = " ") 
load("data/5FUcommon/upregulated_common.RData")
load("data/5FUcommon/downregulated_common.RData")
font_theme = theme(title = element_text(colour = 'black', size = 6),
                   legend.position = "none",
                   legend.title = element_blank(),
                   legend.text = element_text(size = 6, color="black"),
                   legend.key.size = unit(15, "pt"),
                   axis.title = element_text(colour = 'black', size = 6),
                   axis.text = element_text(colour = 'black', size = 6),
                   axis.ticks = element_line(size=0.25, lineend = 10),
                   panel.grid = element_blank())

#raw counts of expr
counts_data_zero = read.csv("data/5FUcommon/rkoDoseZero.csv", header = TRUE,row.names = 1)
counts_data_ten = read.csv("data/5FUcommon/rkoDoseTen.csv", header = TRUE,row.names = 1)
counts_data_fifty = read.csv("data/5FUcommon/rkoDoseFifty.csv", header = TRUE,row.names = 1)

gene_lists <- list(upregulated_common, downregulated_common)
expr.data <- as.matrix(expr.data)
es.dif <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
es.dif <- t(es.dif)
es.dif <- data.frame(es.dif)
es.dif$Common_score <- es.dif$X1 - es.dif$X2
z_scores <- es.dif
z_scores$X1 <- NULL
z_scores$X2 <- NULL

#Merge with annotation data:
z_scores$SampleID <- rownames(z_scores)
z_scores$SampleID <- gsub('\\.', '-', z_scores$SampleID)
z_scores$SampleID <- gsub('\\X', '', z_scores$SampleID)
all(z_scores$SampleID == anno$X)
anno$QS <- z_scores$Common_score
rownames(anno) <- rownames(z_scores)
anno$SampleID <- rownames(anno)

save_cell_state_counts_data <- function(counts_data,anno_data,cell_state,file_name){
  anno_RKO_state = anno_data[anno_data$quiescence_group %in% cell_state,]
  counts_data_state <- counts_data[,colnames(counts_data) %in% as.character(rownames(anno_RKO_state))]
  write.table(counts_data_state,file_name,row.names=TRUE,col.names=TRUE,sep=",")
}

#Select RKO cells
anno_RKO <- anno[anno$df.gid %in% "RKO",]
anno_RKO <- as.character(rownames(anno_RKO))
expr.data_RKO <- expr.data[,colnames(expr.data) %in% anno_RKO]

#Dose 0
anno_RKO <- anno[anno$df.gid %in% "RKO",]
anno_RKO <- anno_RKO[anno_RKO$dose %in% 0,]
# save the data by cell state
save_cell_state_counts_data(counts_data_zero,anno_RKO,"Quiescent","data/5FUcommon/rkoDoseZeroQuiescent.csv")
save_cell_state_counts_data(counts_data_zero,anno_RKO,"Proliferating","data/5FUcommon/rkoDoseZeroProliferating.csv")

anno_RKO <- as.character(rownames(anno_RKO))
umap.expr <- expr.data_RKO[,colnames(expr.data_RKO) %in% anno_RKO]
umap.expr <- umap.expr[rownames(umap.expr) %in% c(downregulated_common, upregulated_common),]
umap.expr <- as.matrix(t(umap.expr))
common.umap = umap(umap.expr, random_state=123)
common.umap$layout 
UMAP_coordinates <- data.frame(common.umap$layout)
colnames(UMAP_coordinates) <- c("UMAP1","UMAP2")
UMAP_coordinates$Sample <- rownames(UMAP_coordinates)
UMAP_coordinates_D0 <- merge(UMAP_coordinates, anno,
                             by.x = "Sample", by.y = "SampleID")

#Dose 10
anno_RKO <- anno[anno$df.gid %in% "RKO",]
anno_RKO <- anno_RKO[anno_RKO$dose %in% 10,]

# save the data by cell state
save_cell_state_counts_data(counts_data_ten,anno_RKO,"Quiescent","data/5FUcommon/rkoDoseTenQuiescent.csv")
save_cell_state_counts_data(counts_data_ten,anno_RKO,"Proliferating","data/5FUcommon/rkoDoseTenProliferating.csv")

anno_RKO <- as.character(rownames(anno_RKO))
umap.expr <- expr.data_RKO[,colnames(expr.data_RKO) %in% anno_RKO]
umap.expr <- umap.expr[rownames(umap.expr) %in% c(downregulated_common, upregulated_common),]
umap.expr <- as.matrix(t(umap.expr))
common.umap = umap(umap.expr, random_state=123)
common.umap$layout 
UMAP_coordinates <- data.frame(common.umap$layout)
colnames(UMAP_coordinates) <- c("UMAP1","UMAP2")
UMAP_coordinates$Sample <- rownames(UMAP_coordinates)
UMAP_coordinates_D10 <- merge(UMAP_coordinates, anno,
                              by.x = "Sample", by.y = "SampleID")

#Dose 50
anno_RKO <- anno[anno$df.gid %in% "RKO",]
anno_RKO <- anno_RKO[anno_RKO$dose %in% 50,]

# save the data by cell state
save_cell_state_counts_data(counts_data_fifty,anno_RKO,"Quiescent","data/5FUcommon/rkoDoseFiftyQuiescent.csv")
save_cell_state_counts_data(counts_data_fifty,anno_RKO,"Proliferating","data/5FUcommon/rkoDoseFiftyProliferating.csv")

anno_RKO <- as.character(rownames(anno_RKO))
umap.expr <- expr.data_RKO[,colnames(expr.data_RKO) %in% anno_RKO]
umap.expr <- umap.expr[rownames(umap.expr) %in% c(downregulated_common, upregulated_common),]
umap.expr <- as.matrix(t(umap.expr))
common.umap = umap(umap.expr, random_state=123)
common.umap$layout 
UMAP_coordinates <- data.frame(common.umap$layout)
colnames(UMAP_coordinates) <- c("UMAP1","UMAP2")
UMAP_coordinates$Sample <- rownames(UMAP_coordinates)
UMAP_coordinates_D50 <- merge(UMAP_coordinates, anno,
                              by.x = "Sample", by.y = "SampleID")

colors = c("#00c3ff","#7de190","#ffff1c")

# colors <- pal_jco("default")(3)
###Plot combined UMAP plots:
UMAP_coordinates <- rbind(UMAP_coordinates_D0, UMAP_coordinates_D50)

p1 = ggplot(UMAP_coordinates, aes(x=UMAP1, y=UMAP2, colour = QS)) +
  geom_point(size = 0.2) +
  scale_color_gradient2(low = colors[1], mid = colors[2], high = colors[3], midpoint = median(UMAP_coordinates$QS)) + theme_classic() + facet_wrap(~dose,nrow = 1)
p1

anno$CellStatus <- sapply(anno$QS, function(x)
  ifelse(x < 0, "Proliferating","Quiescent"))
Dose <- c(0,0,10,10,50,50)
CellStatus <- c("Proliferating","Quiescent","Proliferating","Quiescent","Proliferating","Quiescent")
N <- NULL
for (i in c(0,10,50)) {
  
  print(i)
  test <- anno[anno$df.gid %in% "RKO",]
  test <- test[test$dose %in% i,]
  test <- table(test$CellStatus)
  n <- test[1]
  N <- c(N,n)
  n <- test[2]
  N <- c(N,n)
  
}
Summary <- data.frame(Dose, CellStatus, N)
Summary$Dose <- factor(Summary$Dose, levels = c(0,10,50))
# take off 0 and 50 dose
Summary <- Summary[!(Summary$Dose %in% "10"),] 
p2 <- ggplot(Summary, aes(fill=CellStatus, y=N, x=Dose, width = 0.75)) +
  geom_bar(stat = "identity", position = "fill",width=0.1)+
  labs(x="Dose",y = "Percentage of cells")+ 
  theme_classic()+font_theme
p2 + rotate_x_text(45) + scale_fill_manual(values = c("#666666", "#D95F02"))
p2
picB = ggarrange(p1,p2, ncol = 2, nrow =1,widths=c(9,3), heights=c(4,4), align = "hv")
picB

#########################################################
#Investigation of downregulated_quiescence burst of different dose  treatment cell. 
##########################################################
calculateStats_cell <- function(data,cell_type) {
  bs.data.frame=data.frame(Value =data$bs,sta.type="BS",CellState=cell_type)
  bf.data.frame=data.frame(Value = data$bf,sta.type="BF",CellState=cell_type)
  mean.data.frame=data.frame(Value = data$mean_true,sta.type="Mean",CellState=cell_type)
  var.data.frame=data.frame(Value =data$var_true,sta.type="Variance",CellState=cell_type)
  result = rbind(bs.data.frame,bf.data.frame,mean.data.frame,var.data.frame)
  return(result)
}

bs_zero_dose_quiescent <-read.csv("data/5FUcommon/new_rkoDoseZeroQuiescent_estimated_model_only_prob.csv", header = TRUE)
bs_ten_dose_proliferating <-read.csv("data/5FUcommon/new_rkoDoseZeroProliferating_estimated_model_only_prob.csv", header = TRUE)
bs_one_dose_proliferating <-read.csv("data/5FUcommon/new_rkoDoseTwoProliferating_estimated_model_only_prob.csv", header = TRUE)
bs_fifty_dose_quiescent <-read.csv("data/5FUcommon/new_rkoDoseFiftyQuiescent_estimated_model_only_prob.csv", header = TRUE)
load("data/5FUcommon/upregulated_common.RData")
load("data/5FUcommon/downregulated_common.RData")
apoptosis.genes <- read.table("data/5FUcommon/Apoptosis_hallmarks.txt", header = FALSE,sep = "\t") #downloaded from MSigDb (HALLMARK_APOPTOSIS)
apoptosis.genes <- as.character(apoptosis.genes$V1)

bs_zero_dose = bs_zero_dose_quiescent
bs_ten_dose = bs_ten_dose_proliferating
bs_one_dose = bs_one_dose_proliferating
bs_fifty_dose = bs_fifty_dose_quiescent

gene_type_li = c("downregulated_quiescence","apoptosis","all")
gene_type =gene_type_li[1]
if(gene_type == "downregulated_quiescence"){
  quiescence.genes = downregulated_common
} else if(gene_type == "apoptosis"){
  quiescence.genes = apoptosis.genes
} else if(gene_type == "all"){
  quiescence.genes = Reduce(intersect, list(bs_zero_dose$gene_name,bs_fifty_dose$gene_name, bs_ten_dose$gene_name))
}

quiescence.genes = setdiff(quiescence.genes,"SRRM1")
bs_zero_dose.quiescence <- bs_zero_dose[bs_zero_dose$gene_name %in% quiescence.genes,]
bs_ten_dose.quiescence <- bs_ten_dose[bs_ten_dose$gene_name %in% quiescence.genes,]
bs_one_dose.quiescence <- bs_one_dose[bs_one_dose$gene_name %in% quiescence.genes,]
bs_fifty_dose.quiescence <- bs_fifty_dose[bs_fifty_dose$gene_name %in% quiescence.genes,]

bs_zero_dose.stats = calculateStats_cell(bs_zero_dose.quiescence,"Quie (0 D)")
bs_ten_dose.stats = calculateStats_cell(bs_ten_dose.quiescence,"prol (0 D)")
bs_one_dose.stats = calculateStats_cell(bs_one_dose.quiescence,"Prol (0 D)")
bs_fifty_dose.stats = calculateStats_cell(bs_fifty_dose.quiescence,"Quie (50 D)")

stas.all.data.frame = rbind(bs_one_dose.stats,bs_fifty_dose.stats)
stas.all.data.frame$Value = log(stas.all.data.frame$Value)

figure_file_name = sprintf ("figure/%s_genes_bs_compare_cell_state_0_50.pdf",gene_type)
stas.all.data.frame = stas.all.data.frame[!(stas.all.data.frame$sta.type %in% "BF"),]
p3 = ggpaired(stas.all.data.frame, x = "CellState", y = "Value",
              outlier.size = 0.2,
              color = "CellState", palette = "jco", 
              line.color = "gray", line.size = 0.2,
              width = 0.6,point.size=0.25,
              short.panel.labs = FALSE,)+
               stat_compare_means(label = "p.format",method = "kruskal.test") + facet_wrap(~sta.type, nrow = 1)+
              # No T-test was performed because the sample size was too small
 # facet_wrap(~sta.type, nrow = 1)+ 
  theme_classic()+font_theme
p3
picC = ggarrange(picB,p3, ncol = 2, nrow =1,widths=c(4,4), heights=c(4,4), align = "hv")
picC
ggsave(figure_file_name, width = 7, height = 1.8, useDingbats = FALSE)

figure_file_name = sprintf ("Figures/%s_genes_bs_bf_compare_cell_state_0_0.pdf",gene_type)
# stas.all.data.frame = rbind(bs_zero_dose.stats,bs_ten_dose.stats,bs_fifty_dose.stats)
stas.all.data.frame = rbind(bs_zero_dose.stats,bs_ten_dose.stats)

stas.all.data.frame$Value = log(stas.all.data.frame$Value)
p <- ggboxplot(stas.all.data.frame, x = "CellState", y = "Value",
               color = "CellState", palette = "jco", short.panel.labs = FALSE)

p
# Use only p.format as label. Remove method name.
p + stat_compare_means(label = "p.format") + facet_wrap(~sta.type, nrow = 1)
ggsave(figure_file_name, width = 9, height = 3, useDingbats = FALSE)
ggpaired(stas.all.data.frame, x = "CellState", y = "Value",
         color = "CellState", palette = "jco", 
         line.color = "gray", line.size = 0.4,
         short.panel.labs = FALSE,
)+ stat_compare_means(label = "p.format") + facet_wrap(~sta.type, nrow = 1)


figure_file_name = sprintf ("Figures/%s_genes_bs_bf_compare_cell_state_0_50.pdf",gene_type)
stas.all.data.frame = rbind(bs_one_dose.stats,bs_fifty_dose.stats)
stas.all.data.frame$Value = log(stas.all.data.frame$Value)
p <- ggboxplot(stas.all.data.frame, x = "CellState", y = "Value",width=0.5,size=0.1,
               color = "CellState", palette = "jco", short.panel.labs = FALSE,outlier.size=0.1)

# p <- ggboxplot(stas.all.data.frame, x = "Dose", y = "Value",width=0.5,size=0.1,
#                color = "Dose", palette = "jco", short.panel.labs = FALSE,outlier.size=0.1)

p + stat_compare_means(label = "p.format") + facet_wrap(~sta.type, nrow = 1)+font_theme
ggsave(figure_file_name, width = 6, height = 2.5, useDingbats = FALSE)

figure_file_name = sprintf ("D:/academic_relate_code_two/Nessie-main/DeepGTM/analyseScRNAdata/Fifty5FU/%s_genes_bs_compare_cell_state_0_50.pdf",gene_type)
stas.all.data.frame = stas.all.data.frame[!(stas.all.data.frame$sta.type %in% "BF"),]
p3 = ggpaired(stas.all.data.frame, x = "CellState", y = "Value",
              outlier.size = 0.2,
              color = "CellState", palette = "jco", 
              line.color = "gray", line.size = 0.2,
              width = 0.6,point.size=0.25,
              short.panel.labs = FALSE,
              )+ stat_compare_means(label = "p.format",method = "kruskal.test") + facet_wrap(~sta.type, nrow = 1)+
              # No T-test was performed because the sample size was too small
# )+ facet_wrap(~sta.type, nrow = 1)+ 
  theme_classic()+font_theme
p3
picC = ggarrange(picB,p3, ncol = 2, nrow =1,widths=c(4,4), heights=c(4,4), align = "hv")
picC
ggsave(sprintf( "%sfigure6E-G.pdf",save_figure_dir),width = 7, height = 1.8,useDingbats = FALSE)

























