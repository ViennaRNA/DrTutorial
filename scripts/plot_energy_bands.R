#!/usr/bin/env Rscript
#
library("optparse")

library("ggplot2")
library("ggpubr")
library("scico")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
bwPalette <- c("#eeeeee", "#aaaaaa", "#777777", "#444444", "#111111")
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
cbPalette <- c("#999999", "#E69F00", "#009E73", "#D55E00", "#CC79A7")
cbPalette <- scico(4, palette = 'roma')

data_order = c("DrTrafo","Kinfold", "MFE","sampling")
data_labels = c("DrTransformer", "Kinfold", "RNAfold", "RNAsubopt")
strand_order = c("SRPn", "SRPt", "SRPr", "SRPf")
strand_labels = c("SRPn (native)", "SRPt (U21C)", "SRPr (U21C/C22U/G93A)", "SRPf (U35C/U37C)")

circle             <- 1
circle_filled      <- 21
square             <- 0
square_filled      <- 22
triangle_up        <- 2
triangle_up_filled <- 24
diamond            <- 5
diamond_filled     <- 23

region_shapes_fill  <- c(circle_filled, square_filled, triangle_up_filled, diamond_filled)

option_list = list(
  make_option(c("-i", "--inputfile"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="plot.pdf",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-t", "--title"), type="character", default=NULL,
              help="Plot title", metavar="character"),
  make_option(c("-m","--minmax"), type="logical", action="store_true", default=FALSE,
              help="Plot minimum- and maximum bands instead of 25 and 75 quantile borders"),
  make_option(c("-s", "--start"), type="integer", default=-1,
              help="Start value for y axis"),
  make_option(c("-e", "--end"), type="integer", default=0,
              help="End of transcription, i.e. maximum transcript length"),
  make_option(c("--offset"), type="integer", default=0,
              help="Offset of transcription time points"),
  make_option(c("--noguide"), type="logical", action="store_true", default=FALSE,
              help="Disable the guide, a.k.a. legend."),
  make_option(c("-d","--dots"), type="logical", action="store_true", default=FALSE,
              help="Add dots for each data point")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$inputfile)){
  print_help(opt_parser)
  stop("Specify at lease the input file", call.=FALSE)
}

dat <- read.table(opt$inputfile, head=T, sep=",")

if (opt$end == 0) {
  end = max(dat$length)
} else {
  end = opt$end
}

if (opt$start == -1) {
  start = min(dat$length)
} else {
  start = opt$start
}

x_axis_breaks = seq(start%/%10 * 10, end, 10)
print(start%/%10)

x_label = "Transcript Length [nt]"

if (opt$offset > 0) {
  # adjust x-axis labels
  x_axis_labels <- c()
  for (i in x_axis_breaks) {
    x_axis_labels <- c(x_axis_labels, sprintf('%d (%d)', i + opt$offset, i))
  }

  x_label = "Transcript Length (w/o footprint) [nt]"
} else {
  x_axis_labels = x_axis_breaks
}

p <- ggplot(dat, aes(length,
                     fill=factor(method, levels = data_order, labels = data_labels),
                     color=factor(method, levels = data_order, labels = data_labels)
                     ))
p <- p + geom_ribbon(aes(ymin=Q25, ymax=Q75),
                     alpha=0.4,
#                     colour = NA,
                     size = 0.1)
p <- p + geom_line(aes(length, Qmedian), size=2)

if (opt$dots) {
  p <- p + geom_point(aes(length, Qmedian, shape = factor(method)),
                      size=2,
                      fill = "white",
                      data=subset(dat, length %% 5 == 0))
}

p <- p + facet_grid(~factor(name, levels=strand_order, labels = strand_labels))

p <- p + scale_x_continuous(expand = c(0, 0),
                            limit   = c(start, end),
                            breaks  = x_axis_breaks,
                            labels  = x_axis_labels)
p <- p + scale_y_continuous(limit=c(-60, 15), breaks=seq(-60, 0, 10), expand = c(0,0))

p <- p +  scale_fill_manual(  values  = cbPalette)
p <- p +  scale_colour_manual(  values  = cbPalette)
#p <- p +  scale_fill_scico_d(  palette = "batlow")
#p <- p +  scale_colour_scico_d(  palette = "batlow")
#p <- p +  scale_fill_brewer(palette = 'Paired')
#p <- p +  scale_colour_brewer(palette = 'Paired')
p <-  p + scale_shape_manual( breaks  = data_order,
                              labels  = data_labels,
                              values  = region_shapes_fill )

p <- p + xlab(x_label)
p <- p + ylab(expression(paste("Free Energy [", kcal, ~mol^{-1}, "]")))

if (!is.null(opt$title)){
  p <- p + ggtitle(opt$title)
}

p <- p + theme_bw()
p <- p + theme(plot.title = element_text(hjust = 0.5, size = 24),
                plot.background   = element_blank(),
#                panel.grid.minor  = element_blank(),
                axis.line = element_blank(),
                panel.spacing = unit(0.5, "lines"),
                panel.border = element_blank(),
                panel.background = element_blank(),
                strip.background = element_blank())
p <- p +  theme(#axis.title.x  = element_text(family="Helvetica", size = 16, colour="#000000"),
                axis.title.x  = element_blank(),
                axis.text.x   = element_text(family="Helvetica", size = 14, colour="#666666", angle = 90.),
                axis.title.y  = element_text(family="Helvetica", size = 16, colour="#000000"),
                axis.text.y   = element_text(family="Helvetica", size = 14, colour="#666666"),
                strip.text.x  = element_text(family="Helvetica", size = 20, colour="#111111")
               )
if (opt$noguide) {
p <- p + theme(legend.position = "none")
} else {
p <- p +  guides(
            color = guide_legend(title = "Method"),
            fill = guide_legend(title = "Method"),
            shape = guide_legend(title = "Method")
          )
p <- p + theme(legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 18),
            legend.background = element_blank(),
            legend.key.height = unit(2,"line"),
            legend.key.width  = unit(2,"line"))
}

p2 <- ggplot(dat, aes(length,
                     fill=factor(method, levels = data_order, labels = data_labels),
                     color=factor(method, levels = data_order, labels = data_labels)
                     ))
p2 <- p2 + geom_ribbon(aes(ymin=Qmin, ymax=Qmax),
                       alpha=0.4,
#                     colour = NA,
                       size = 0.1,
                       linetype = 2)

p2 <- p2 + geom_line(aes(length, Qmean), size=2)
if (opt$dots) {
  p2 <- p2 + geom_point(aes(length, Qmean, shape = factor(method)),
                        size=2,
                        fill = "white",
                        data=subset(dat, length %% 5 == 0))
}

p2 <- p2 + facet_grid(~factor(name, levels=strand_order, labels = strand_labels))

p2 <- p2 + scale_x_continuous(expand = c(0, 0),
                            limit   = c(start, end),
                            breaks  = x_axis_breaks,
                            labels  = x_axis_labels)
p2 <- p2 + scale_y_continuous(limit=c(-60, 15), breaks=seq(-60, 0, 10), expand = c(0,0))

p2 <- p2 +  scale_fill_manual(  values  = cbPalette)
p2 <- p2 +  scale_colour_manual(  values  = cbPalette)
#p2 <- p2 +  scale_fill_brewer(palette = 'Paired')
#p2 <- p2 +  scale_colour_brewer(palette = 'Paired')
p2 <-  p2 + scale_shape_manual( breaks  = data_order,
                              labels  = data_labels,
                              values  = region_shapes_fill )

p2 <- p2 + xlab(x_label)
p2 <- p2 + ylab(expression(paste("Free Energy [", kcal, ~mol^{-1}, "]")))

if (!is.null(opt$title)){
  p2 <- p2 + ggtitle(opt$title)
}

p2 <- p2 + theme_bw()
p2 <- p2 + theme(plot.title = element_text(hjust = 0.5, size = 20),
                plot.background   = element_blank(),
#                panel.grid.minor  = element_blank(),
                axis.line = element_blank(),
                panel.spacing = unit(0.5, "lines"),
                panel.border = element_blank(),
                panel.background = element_blank(),
                strip.background = element_blank(),
                strip.text.x = element_blank())
p2 <- p2 +  theme(axis.title.x  = element_text(family="Helvetica", size = 16, colour="#000000"),
                axis.text.x   = element_text(family="Helvetica", size = 14, colour="#666666", angle = 90.,vjust=0.5),
                axis.title.y  = element_text(family="Helvetica", size = 16, colour="#000000"),
                axis.text.y   = element_text(family="Helvetica", size = 14, colour="#666666")
               )

p2 <- p2 + theme(legend.position = "none")


ppp <- ggarrange(
  p, p2,
  labels = c("A", "B"),
  font.label=list(color="black",size=32),
  ncol = 1,
  common.legend = TRUE, legend = "bottom"
  )

ggsave(plot = ppp, file=opt$out, width=15, height=10)
