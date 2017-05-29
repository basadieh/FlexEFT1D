  LO_theme <- theme_bw(base_family = "serif") + theme(
    panel.grid.major = element_line(size = 0, color="white"),
    panel.grid.minor = element_line(size = 0, color="white"),
    panel.border     = element_rect(size = .8, color="black", fill = NA),
  axis.line        = element_line(size = .8, color = "black"),
  axis.text        = element_text(size = 10),
  axis.title       = element_text(size = 12),
  legend.position  = c(1.5,.25),
  text             = element_text(size=14),
  plot.background  = element_rect(color = "white"),
  panel.background = element_rect(size=.2, color = "black")
  ,panel.margin     = unit(c(1,0,0,0), "lines")
)

