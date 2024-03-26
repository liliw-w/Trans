theme_my_pub <- function(
                         base_size = 7,
                         base_family = "Helvetica",
                         base_line_size = 0.5,
                         base_rect_size = 0.5,
                         
			 # title
			 title.size = 8,
			 
                         # panel
                         panel.grid.major.y.linetype = "blank",
                         
                         # legend
                         legend.position = "bottom",
			 legend_linetype = "blank",
      
                         # axis
                         axis.text.x.angle = 0,
                         axis.text.x.vjust = 0.5
			 ){
  theme_classic(base_size = base_size, base_family = base_family, base_line_size = base_line_size, base_rect_size = base_rect_size) +
    theme(
          # panel
          panel.grid.major.y = element_line(color = "#e5e5e5",
                                            linetype = panel.grid.major.y.linetype),
          panel.grid.minor = element_blank(),
          
          # legend
          legend.position = legend.position,
          legend.title = element_text(face = "bold"),
          legend.background = element_rect(color = "black", linetype = legend_linetype),
          legend.key.size= unit(0.5, "cm"),
          
          
          # axis
          ## axis title
          axis.title.x = element_text(vjust = -0.2),
          axis.title.y = element_text(vjust = 2),
          axis.title.y.right = element_text(angle = 90),
          
          ## axis tick labels
          axis.text = element_text(colour = "black"),
          axis.text.x = element_text(vjust = axis.text.x.vjust),
          
          ## axis line
          axis.line = element_line(colour = "black"),
          
          
          # plot
          plot.title = element_text(size = title.size, hjust = 0.5, face = "bold"),
          plot.margin = unit(c(5, 5, 5, 5), "mm")
          
          )
}


integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}


## palette
mypalette <- list(
  'paired' = c(
    "#8AACAA", "#2295A2", # blue (dark)
    "#E5D2A2", "#E3A031", # red
    "#DBC6E4", "#8D6E97", # yellow
    "#EFCAD3", "#DD7C83", # purple
    "#C6DBD7", "#87B2AC"  # blue (light)
  ),
  
  "contrast1" = c(
    "#4D9490", # blue
    "#d28b84", # red
    "#EDAD6D", # yellow
    "#AB6CAD", # purple
    "#91BC6B" # green
  ),
  
  "contrast2" = c(
    "#76758B", # blue
    "#DDACAA", # red
    "#B783BD", # purple
    "#F2A869", # yellow
    "#AFB59B" # green
  )
)

