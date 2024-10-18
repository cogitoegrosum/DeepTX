library(RColorBrewer)
library(scales)
colormap<- rev(brewer.pal(9,"Blues")[1:6])
color.feedback <- c("#3AA438", "#547CBE", "#E71419","#8abd44","#db2f2c")
dash.color = '#63B8FF'

plot_bf_density = function(new_df,filter_gene){
  filter_gene = row.names(new_df)[1:10]
  figureF_SA = ggplot(new_df[filter_gene,], aes(x = bf1, y = bf2)) +                                                    
    stat_density2d(
      geom ="raster",
      aes(fill = ..density..,),
      #       colour = "white",
      #     linewidth = 0.15,
      #     colour = "grey",
      #     bins = 6
      contour = F
    ) +
    stat_density2d(linewidth = 0.15,colour = "grey")+
    scale_fill_gradientn(colours=rev(colormap))+
    geom_point(data = new_df,aes(x =bf1, y = bf2),
               color = color.feedback[4],shape = 16,size = 0.3,alpha = 1) +
    geom_point(data = new_df,aes(x =bf1, y = bf2),
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
      legend.position = c(0.9,0.3),
      # legend.position = "none",
      legend.title = element_blank(),
      legend.text = element_text(size = 6, color = "black"),
      legend.key.size = unit(6, "pt"),
      axis.title = element_text(colour = 'black', size = 6),
      axis.text = element_text(colour = 'black', size = 6),
      axis.ticks = element_line(size = 0.25, lineend = 10),
      panel.grid = element_blank()
    )
  return(figureF_SA)
}

plot_bs_density = function(new_df,filter_gene,dose_type="5FU"){
  if(dose_type=="IdU"){
    limit_value = c(1, 10000)
  }else {
    limit_value = c(0.1, 10000)
  }
  
  figureF_SB = ggplot(new_df, aes(x = bs1, y =bs2)) +
    stat_density2d(
      geom ="raster",
      aes(fill = ..density..,),
      colour = "grey",
      linewidth = 0.01,
      #     colour = "grey",
          bins = 6,
      contour = F
    ) +
    stat_density2d(linewidth = 0.01,colour = "grey")+
    scale_fill_gradientn(colours=rev(colormap))+
    geom_point(data = new_df[filter_gene,],aes(x =bs1, y = bs2),color = color.feedback[5],shape = 16,size = 0.3,alpha = 1) +
    geom_point(data = new_df[filter_gene,],aes(x =bs1, y = bs2),color = color.feedback[4],shape = 16,size = 0.3,alpha = 1) +
    geom_abline(slope = 1,color = dash.color,linetype = "dashed") +
    labs(x = " BS (Untreatemnt)",y = expression(paste("BS (50 ", mu,"M 5FU treatment)")),) +
    theme_bw() +
    scale_y_continuous(breaks = 10^(-100:100),
                       labels = trans_format("log10", math_format(10^.x)),
                       trans="log10",limits = limit_value )+
    scale_x_continuous(breaks = 10^(-100:100),
                       labels = trans_format("log10", math_format(10^.x)),
                       trans="log10",limits = limit_value)+
    theme(
      title = element_text(colour = 'black', size = 6),
      # legend.position = "none",
      legend.position = c(0.9,0.3),
      legend.title = element_blank(),
      legend.text = element_text(size = 6, color = "black"),
      legend.key.size = unit(15, "pt"),
      axis.title = element_text(colour = 'black', size = 6),
      axis.text = element_text(colour = 'black', size = 6),
      axis.ticks = element_line(size = 0.25, lineend = 10),
      panel.grid = element_blank()
    )
}
