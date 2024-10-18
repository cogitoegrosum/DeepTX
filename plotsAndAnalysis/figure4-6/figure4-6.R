library(Seurat)
library(dplyr)
library(cowplot)
library(GSVA)
library(data.table)
library(ggplot2)
library(ggpubr)

font_theme = theme(
  title = element_text(colour = 'black', size = 6),
  legend.position = "none",
  legend.title = element_blank(),
  legend.text = element_text(size = 8, color = "black"),
  legend.key.size = unit(15, "pt"),
  axis.title = element_text(colour = 'black', size = 8),
  axis.text = element_text(colour = 'black', size = 6),
  axis.ticks = element_line(size = 0.25, lineend = 10),
  panel.grid = element_blank()
)

calculateStats <- function(data, cell_type) {
  data$CV = sqrt(data$var_true) / data$mean_true
  bs.data.frame = data.frame(Value = data$bs,
                             sta.type = "BS",
                             Dose = cell_type)
  bf.data.frame = data.frame(Value = data$bf,
                             sta.type = "BF",
                             Dose = cell_type)
  mean.data.frame = data.frame(Value = data$mean_true,
                               sta.type = "Mean",
                               Dose = cell_type)
  var.data.frame = data.frame(Value = data$var_true,
                              sta.type = "Variance",
                              Dose = cell_type)
  cv.data.frame = data.frame(Value = data$CV,
                             sta.type = "CV",
                             Dose = cell_type)
  result = rbind(bs.data.frame,
                 bf.data.frame,
                 mean.data.frame,
                 var.data.frame,
                 cv.data.frame)
  return(result)
}

generate_burst_matrix <- function(data_one, data_two, burst_type) {
  data_two[, burst_type]
  new_df <-
    data.frame(burst1 = data_one[, burst_type], burst2 = data_two[, burst_type])
  counts <- as.matrix(new_df)
  rownames(counts) <- data_one$gene_name
  return(counts)
}
## load data
setwd(
  "D:/academic_relate_code_two/Nessie-main/DeepTX/inferObservedData/result"
)
bs_zero_dose <-
  read.csv("Ten5FU/doseZero_estimated_model_stats_prob.csv",
           header = TRUE)
bs_ten_dose <-
  read.csv("Ten5FU/doseTen_estimated_model_stats_prob.csv",
           header = TRUE)
bs_one_dose <-
  read.csv("Fifty5FU/doseOne_estimated_model_stats_prob.csv",
           header = TRUE)
bs_fifty_dose <-
  read.csv("Fifty5FU/doseFifty_estimated_model_stats_prob.csv",
           header = TRUE)
bs_dmso <-
  read.csv("IdUDMSO/dmso_estimated_model_stats_prob.csv",
           header = TRUE)
bs_IdU <-
  read.csv("IdUDMSO/idu_estimated_model_stats_prob.csv",
           header = TRUE)

bs_zero_dose.stats = calculateStats(bs_zero_dose, 0)
bs_one_dose.stats = calculateStats(bs_one_dose, 0)
bs_ten_dose.stats = calculateStats(bs_ten_dose, 10)
bs_fifty_dose.stats = calculateStats(bs_fifty_dose, 50)
bs_dmso.stats = calculateStats(bs_dmso, "Dmso")
bs_IdU.stats = calculateStats(bs_IdU, "IdU")

#******* Select the Data Type tab to decide which data set to analyze next *******
data_type = c("IdU","5FU_10","5FU_50")
data_type_index= 1
if(data_type_index==1){
  drug_type = "IdU"
  dose_type = "IdU"
} else if(data_type_index==2){
  drug_type = "5FU"
  dose_type = "10"
}else if(data_type_index==3){
  drug_type = "5FU"
  dose_type = "50"
}

setwd("D:/academic_relate_code_two/Nessie-main/DeepTX/plotsAndAnalysis")
source("plot.R")
if (dose_type == 10) {
  stas.all.data.frame = rbind(bs_zero_dose.stats, bs_ten_dose.stats)
  new_df <-
    data.frame(
      bf1 = bs_zero_dose$bf,
      bf2 = bs_ten_dose$bf,
      bs1 = bs_zero_dose$bs,
      bs2 = bs_ten_dose$bs,
      row.names = bs_zero_dose$gene_name
    )
  figure_burst_compare = "figure/figure5A.pdf"
  volcano_file_name = "figure/figure5C.pdf"
  go_down_file_name = "figure/figureS7G.pdf"
  go_up_file_name = "figure/figure5E.pdf"
  bs_density_file_name = "figure/figure5B.pdf"
} else if (dose_type == 50) {
  stas.all.data.frame = rbind(bs_one_dose.stats, bs_fifty_dose.stats)
  new_df <-
    data.frame(
      bf1 = bs_one_dose$bf,
      bf2 = bs_fifty_dose$bf,
      bs1 = bs_one_dose$bs,
      bs2 = bs_fifty_dose$bs,
      row.names = bs_one_dose$gene_name
    )
  figure_burst_compare = "figure/figure6A_else.pdf"
  volcano_file_name = "figure/figure6B.pdf"
  go_down_file_name = "figure/figureS9G.pdf"
  go_up_file_name = "figure/figure6D.pdf"
  bs_density_file_name = "figure/figure6A.pdf"
} else if (dose_type == "IdU") {
  stas.all.data.frame = rbind(bs_dmso.stats, bs_IdU.stats)
  new_df <-
    data.frame(
      bf1 = bs_dmso$bf,
      bf2 = bs_IdU$bf,
      bs1 = bs_dmso$bs,
      bs2 = bs_IdU$bs,
      row.names = bs_IdU$gene_name
    )
  figure_burst_compare = "figure/figure4B_else.pdf"
  volcano_file_name = "figure/figure4A_else.pdf"
  go_down_file_name = "figure/figure4D.pdf"
  go_up_file_name = "figure/figureS5G.pdf"
  bs_density_file_name = "figure/figure4A.pdf"
  GSEA_one = "figure/figure4E.pdf"
  GSEA_two = "figure/figure4F.pdf"
}

stas.all.data.frame$Value = log10(stas.all.data.frame$Value)
p <-
  ggboxplot(
    stas.all.data.frame,
    x = "Dose",
    y = "Value",
    width = 0.5,
    size = 0.1,
    color = "Dose",
    palette = "jco",
    short.panel.labs = FALSE,
    outlier.size = 0.1
  )
# Use only p.format as label. Remove method name.
p <-
  p + stat_compare_means(label = "p.format") + facet_wrap( ~ sta.type, nrow = 1) +
  font_theme
p
ggsave(
  figure_burst_compare,
  width = 6.8,
  height = 2,
  useDingbats = FALSE
)

library(edgeR)
burst_type = "bs"
if (dose_type == 10) {
  counts = generate_burst_matrix(bs_zero_dose, bs_ten_dose, burst_type)
  FC_threshold = 2.9
  bcv <- 0.4
} else if (dose_type == 50) {
  counts = generate_burst_matrix(bs_one_dose, bs_fifty_dose, burst_type)
  FC_threshold = 1
  bcv <- 0.4
} else if (dose_type == "IdU") {
  counts = generate_burst_matrix(bs_dmso, bs_IdU, burst_type)
  FC_threshold = 2.9
  bcv <- 0.4
}
group = factor(c("Normal", "Disease"))
y <- DGEList(counts = counts, group = group)
# Calculate the library size factor for each sample
# y <- calcNormFactors(y)
y_bcv <- y
# bcv <- 0.4
et <- exactTest(y_bcv, dispersion = bcv ^ 2)
df <- et$table

gene1 <- decideTestsDGE(et, p.value = 0.05, lfc = 0)
summary(gene1)
df$logpvalue = -log10(df$PValue)
ggscatter(df, x = "logFC", y = "logpvalue") + theme_bw()

df$Group = "normal"
df$Group[which((df$PValue < 0.05) &
                 (df$logFC > FC_threshold))] = "up"
df$Group[which((df$PValue < 0.05) &
                 (df$logFC < -FC_threshold))] = "down"
table(df$Group)
df$Label = ""
# Sort the pvalue values of differentially expressed genes in ascending order
df <- df[order(df$PValue),]
# Select the top ten differentially expressed genes according to pvalue
up.genes <- row.names(df)[which(df$Group == "up")]
down.genes <- row.names(df)[which(df$Group == "down")]

# Merge up.genes and down.genes and add them to Label
df$Symbol = row.names(df)
deg.top10.genes <-
  c(as.character(up.genes[1:5]), as.character(down.genes[1:5]))
df$Label[match(deg.top10.genes, df$Symbol)] <- deg.top10.genes
ggscatter(
  df,
  x = "logFC",
  y = "logpvalue",
  color = "Group",
  palette = c("#00ba38", "grey11", "#f13527"),
  size = 0.1,
  label = df$Label,
  font.label = 6,
  repel = T
) + theme_bw() + xlab("log2(FC)") + ylab("-log10(pvalue)") + font_theme +
  geom_hline(yintercept = 1.30, linetype = "dashed") +
  geom_vline(xintercept = c(-FC_threshold, FC_threshold),
             linetype = "dashed")

ggsave(
  volcano_file_name,
  width = 2,
  height = 1.8,
  useDingbats = FALSE
)

library(clusterProfiler)
library(stringr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DOSE)
library(ggrepel)
options(warn = -1)
enrichGO_analysis <- function(gene_names) {
  id_list <- mapIds(org.Hs.eg.db, gene_names, "ENTREZID", "SYMBOL")
  id_list <- na.omit(id_list)
  # GO enrichment analysis
  go <- enrichGO(
    gene = id_list,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = T
  )
  go.res <- data.frame(go)
  return(go.res)
}

plotGo <- function(go.re) {
  font.size = 8
  go_bar <- ggplot(data = go.re,
                   aes(x = Description, y = Count, fill = pvalue)) +
    geom_bar(stat = "identity", width = 0.4) +
    coord_flip() +
    scale_x_discrete(
      labels = function(x)
        str_wrap(x, width = 200)
    ) +
    labs(x = "GO terms", y = "GeneNumber", title = "") +
    scale_fill_distiller(palette = "greens", direction = 1) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(
        colour = "black",
        size = font.size,
        vjust = 1
      ),
      axis.text.y = element_text(
        colour = "black",
        size = font.size,
        hjust = 1
      ),
      axis.title = element_text(
        margin = margin(10, 5, 0, 0),
        color = "black",
        size = font.size
      ),
      axis.title.y = element_text(angle = 90)
    )
  return(go_bar)
}

df_sub_down = df[df$Group %in% "down", ]
go.res_down = enrichGO_analysis(rownames(df_sub_down))

if (dose_type == "IdU") {
  go_sub_index =   c(
    "GO:0044772",
    "GO:1901987",
    "GO:0000086",
    "GO:0044839",
    "GO:1901990",
    "GO:0007093",
    "GO:0042770",
    "GO:0000070",
    "GO:0000819",
    "GO:0098813",
    "GO:0007059"
  )
} else{
  go_sub_index = row.names(go.res_down)[1:10]
}
go.res_down.sub = go.res_down[go_sub_index, ]
go.res_down.sub <-
  go.res_down.sub[order(go.res_down.sub$Count, decreasing = c(TRUE)),]
go.res_down.sub$Description <-
  factor(go.res_down.sub$Description,
         levels = rev(go.res_down.sub$Description))
go_bar = plotGo(go.res_down.sub)
go_bar
ggsave(
  sprintf(go_down_file_name, drug_type, dose_type),
  width = 8,
  height = 5,
  useDingbats = FALSE
)


df_sub_up = df[df$Group %in% "up", ]
go.res_up = enrichGO_analysis(rownames(df_sub_up))
if (dose_type == 10) {
  go_sub_index =   c(
    "GO:1902175",
    "GO:0008631",
    "GO:2001243",
    "GO:2001234",
    "GO:2001242",
    "GO:0097193",
    "GO:2001233",
    "GO:0006119",
    "GO:0006979"
  )
} else if (dose_type == 50) {
  go_sub_index =   c(
    "GO:1904874",
    "GO:1904872",
    "GO:0090671",
    "GO:0007004",
    "GO:0010833",
    "GO:1904851",
    "GO:0006123",
    "GO:0042773",
    "GO:0042775"
  )
} else{
  go_sub_index = row.names(go.res_up)[1:10]
}
go.res_up.sub = go.res_up[c(go_sub_index), ]
go.res_up.sub <-
  go.res_up.sub[order(go.res_up.sub$Count, decreasing = c(TRUE)),]

go.res_up.sub$Description <-
  factor(go.res_up.sub$Description,
         levels = rev(go.res_up.sub$Description))
go_bar = plotGo(go.res_up.sub)
go_bar
ggsave(
  go_up_file_name,
  width = 8,
  height = 5,
  useDingbats = FALSE
)

# plot density
# figureF_SA = plot_bs_density(new_df,up.genes)
# figureF_SA
# ggsave(sprintf("Figures/%s/density_bs_up_%s.pdf",drug_type, dose_type),width = 2,height =1.9)

if (dose_type == "IdU") {
  figureF_SB =  plot_bs_density(new_df, down.genes)
  figureF_SB
  ggsave(bs_density_file_name, width = 2, height = 1.9)
} else{
  figureF_SC =  plot_bs_density(new_df, c())
  figureF_SC
  ggsave(bs_density_file_name, width = 2, height = 1.9)
}


#gsea enrichment analysis
organism = 'hsa'
OrgDb = 'org.Hs.eg.db'
id_list <- mapIds(org.Hs.eg.db, rownames(df), "ENTREZID", "SYMBOL")

gene <- bitr(rownames(df),
             fromType = "SYMBOL",
             toType =  "ENTREZID",
             OrgDb = org.Hs.eg.db)
id_list <- na.omit(id_list)

geneList = df$logFC
names(geneList) = gene$ENTREZID

geneList = sort(geneList, decreasing = T)
ego <- gseGO(
  geneList     = geneList,
  OrgDb        = org.Hs.eg.db,
  nPerm        = 1000,
  # minGSSize = 10,
  # maxGSSize = 500,
  minGSSize = 2,
  # maxGSSize = 2000,
  maxGSSize = 3000,
  ont          = "ALL",
  # pvalueCutoff = 0.05,
  pvalueCutoff = 1,
  verbose      = FALSE
)

go = DOSE::setReadable(ego, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
go = ego@result
go_gse = go
sortgo <- go_gse[order(go_gse$enrichmentScore, decreasing = T), ]
dim(sortgo)
#write.table(sortkk,"gsea_output2.txt",sep = "\\t",quote = F,col.names = T,row.names = F)

if (dose_type == "IdU") {
  gesa_cols = c("GO:0031573", "GO:0055007")
} else{
  gesa_cols = row.names(go)
}

library(enrichplot)

width = 3.4 # for adobe figure
width = 7
gse_plot = gseaplot2(
  ego,
  gesa_cols[1],
  title = "GO GSEA enrichment analysis",
  base_size = 6,
  color = "green",
  pvalue_table = TRUE,
  ES_geom = "line"
)

gse_plot[[1]] = gse_plot[[1]] +
  
  theme(
    title = element_text(colour = 'black', size = 6),
    #         legend.position = "none",
    #         legend.title = element_blank(),
    legend.text = element_text(size = 6, color = "black"),
    legend.key.size = unit(6, "pt"),
    axis.title = element_text(colour = 'black', size = 7),
    axis.text = element_text(colour = 'black', size = 6),
    axis.ticks = element_line(size = 0.15, lineend = 10),
    panel.grid = element_blank()
  )
gse_plot[[2]] = gse_plot[[2]] +
  # theme_bw() +
  theme(
    title = element_text(colour = 'black', size = 6),
    #         legend.position = "none",
    #         legend.title = element_blank(),
    legend.text = element_text(size = 6, color = "black"),
    legend.key.size = unit(10, "pt"),
    axis.title = element_text(colour = 'black', size = 7),
    axis.text = element_text(colour = 'black', size = 6),
    axis.ticks = element_line(size = 0.25, lineend = 10),
    #         panel.grid = element_blank()
  )
gse_plot[[3]] = gse_plot[[3]] +
  # theme_bw() +
  theme(
    title = element_text(colour = 'black', size = 6),
    #         legend.position = "none",
    #         legend.title = element_blank(),
    legend.text = element_text(size = 6, color = "black"),
    legend.key.size = unit(15, "pt"),
    axis.title = element_text(colour = 'black', size = 7),
    axis.text = element_text(colour = 'black', size = 6),
    axis.ticks = element_line(size = 0.25, lineend = 10),
    #         panel.grid = element_blank()
  )
gse_plot
if (dose_type == "IdU") {
  ggsave(GSEA_one,
         width = width,
         height = 3,
         useDingbats = FALSE)
}

width = 3.4 # for adobe figure
width = 4.5
gse_plot = gseaplot2(
  ego,
  gesa_cols[2],
  title = "GO GSEA enrichment analysis",
  base_size = 6,
  color = "green",
  pvalue_table = TRUE,
  ES_geom = "line"
)

gse_plot[[1]] = gse_plot[[1]] +
  
  theme(
    title = element_text(colour = 'black', size = 6),
    #         legend.position = "none",
    #         legend.title = element_blank(),
    legend.text = element_text(size = 6, color = "black"),
    legend.key.size = unit(6, "pt"),
    axis.title = element_text(colour = 'black', size = 7),
    axis.text = element_text(colour = 'black', size = 6),
    axis.ticks = element_line(size = 0.15, lineend = 10),
    panel.grid = element_blank()
  )
gse_plot[[2]] = gse_plot[[2]] +
  # theme_bw() +
  theme(
    title = element_text(colour = 'black', size = 6),
    #         legend.position = "none",
    #         legend.title = element_blank(),
    legend.text = element_text(size = 6, color = "black"),
    legend.key.size = unit(10, "pt"),
    axis.title = element_text(colour = 'black', size = 7),
    axis.text = element_text(colour = 'black', size = 6),
    axis.ticks = element_line(size = 0.25, lineend = 10),
    #         panel.grid = element_blank()
  )
gse_plot[[3]] = gse_plot[[3]] +
  # theme_bw() +
  theme(
    title = element_text(colour = 'black', size = 6),
    #         legend.position = "none",
    #         legend.title = element_blank(),
    legend.text = element_text(size = 6, color = "black"),
    legend.key.size = unit(15, "pt"),
    axis.title = element_text(colour = 'black', size = 7),
    axis.text = element_text(colour = 'black', size = 6),
    axis.ticks = element_line(size = 0.25, lineend = 10),
    #         panel.grid = element_blank()
  )
gse_plot
if (dose_type == "IdU") {
  ggsave(GSEA_two,
         width = width,
         height = 3,
         useDingbats = FALSE)
}
