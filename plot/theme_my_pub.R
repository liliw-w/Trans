theme_my_pub <- function(
                         base_size = 12,
                         base_family = "Helvetica",
                         
			 # title
			 title.size = 12,
			 
                         # panel
                         panel.grid.major.y.linetype = "blank",
                         
                         # legend
                         legend.position = "bottom",
                         legend.text.size = 10,
			 legend_linetype = "blank",
      
                         # axis
                         ## axis title
                         axis.title.size = 12,
                         ## axis tick labels
                         axis.text.size = 10,
                         axis.text.x.angle = 0,
                         axis.text.x.vjust = 0.5,
                         ## axis line
                         axis.line.size = 0.7){
  theme_classic(base_size = base_size, base_family = base_family ) +
    theme(
          # panel
          panel.grid.major.y = element_line(color = "#e5e5e5",
                                            linetype = panel.grid.major.y.linetype),
          panel.grid.minor = element_blank(),
          
          # legend
          legend.position = legend.position,
          legend.title = element_text(size = legend.text.size,
                                      face = "bold"),
          legend.text = element_text(size = legend.text.size),
          legend.background = element_rect(color = "black",
                                           linetype = legend_linetype),
          legend.key.size= unit(0.5, "cm"),

          
          # axis
          ## axis title
          axis.title = element_text(size = axis.title.size),
          axis.title.x = element_text(vjust = -0.2),
          axis.title.y = element_text(vjust = 2),
          axis.title.y.right = element_text(angle = 90),
          
          ## axis tick labels
          axis.text = element_text(colour = "black",
                                   size = axis.text.size),
          axis.text.x = element_text(angle = axis.text.x.angle,
                                     vjust = axis.text.x.vjust),
          
          ## axis line
          axis.line = element_line(colour = "black", linewidth = axis.line.size),
          
          
          # plot
          plot.title = element_text(size = title.size, hjust = 0.5, face = "bold"),
          plot.margin = unit(c(10, 5, 5, 5), "mm")
          )
}

