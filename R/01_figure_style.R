# Uniform formatting of figures

figure_style <- list(blue='#4477AA', red='#CC6677', green='#228B75', orange='#C47732',
                     axis_title=14, axis_text=11, species=12, legend=11,
                     panel_tag=18, width=9, dpi=300)
figure_theme <- function() {
 ggplot2::theme_bw(base_size=12) + ggplot2::theme(
  text=ggplot2::element_text(family='sans', colour='black'),
  axis.title=ggplot2::element_text(size=14,face='plain'),
  axis.text=ggplot2::element_text(size=11,colour='black'),
  panel.grid.minor=ggplot2::element_blank(),
  panel.grid.major.x=ggplot2::element_blank(),
  panel.grid.major.y=ggplot2::element_line(colour='#E6E6E6',linewidth=.3),
  panel.border=ggplot2::element_rect(colour='#333333',fill=NA,linewidth=.5),
  strip.background=ggplot2::element_rect(fill='#D9D9D9',colour='#333333',linewidth=.5),
  strip.text=ggplot2::element_text(size=12,face='plain'),
  legend.text=ggplot2::element_text(size=11),
  plot.tag=ggplot2::element_text(size=18,face='bold'),
  plot.background=ggplot2::element_rect(fill='white',colour=NA))
}
