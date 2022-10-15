#!/usr/bin/env Rscript

library("optparse")
library("reshape2")
library("ggplot2")
library("scales")
library("viridis")

max_length <- 117
slope <- 1.6
intercept <- -2.29
SHAPE_range <- 4
SHAPE_convert_log <- 0

option_list = list(
  make_option(c("-i", "--inputfile"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--outputfile"), type="character", default="Reactivity.pdf",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-t", "--title"), type="character", default=NULL,
              help="Plot title", metavar="character"),
  make_option(c("--slope"), type="double", default=1.6,
              help="Slope for SHAPE to probability conversion"),
  make_option(c("--intercept"), type="double", default=-2.29,
              help="Intercept for SHAPE to probability conversion"),
  make_option(c("--SHAPE"), type="logical", action="store_true", default=FALSE,
              help="Input data is SHAPE reactivity"),
  make_option(c("--normalize"), type="logical", action="store_true", default=FALSE,
              help="Normalize SHAPE reactivity"),
  make_option(c("--SHAPE2probs"), type="logical", action="store_true", default=FALSE,
              help="Convert SHAPE reactivities to unpaired probabilities"),
  make_option(c("--SHAPErange"), type="integer", default=4,
              help="Largest SHAPE reativity value"),
  make_option(c("--offset"), type="integer", default=0,
              help="Offset of transcription time points"),
  make_option(c("--nofootprint"), type="logical", action="store_true", default=FALSE,
              help="Do not add footprint despite of offset != 0"),
  make_option(c("-s", "--start"), type="integer", default=30,
              help="Start value for y axis"),
  make_option(c("-e", "--end"), type="integer", default=0,
              help="End of transcription, i.e. maximum transcript length"),
  make_option(c("-l", "--ltitle"), type="character", default="Accessibility",
              help="Title of the legend", metavar="character"),
  make_option(c("-x", "--xlabel"), type="character", default=NULL,
              help="Label for the X-axis", metavar="character"),
  make_option(c("-y", "--ylabel"), type="character", default=NULL,
              help="Label for the Y-axis", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$inputfile)){
  print_help(opt_parser)
  stop("Specify at least the input file", call.=FALSE)
}

if (opt$SHAPErange) {
  SHAPE_range <- opt$SHAPErange
}

SHAPE2probs = function(x) {
    if (SHAPE_convert_log) {
        x <- log(x)
    }
    x <- x - intercept
    x <- x / slope
    x <- max(min(x, 1), 0)

    return(x)
}



dat <- read.csv(opt$inputfile, header=T, sep=",", check.names=F)

if (opt$end == 0) {
  max_length = max(dat[3])
} else {
  max_length = opt$end
}

y_axis_breaks = seq(opt$start, max_length, 10)

if (opt$offset > 0) {
  # adjust y-axis labels
  y_axis_labels <- c()
  for (i in seq(opt$start, max_length, 10)) {
    y_axis_labels <- c(y_axis_labels, sprintf('%d (%d)', i + opt$offset, i))
  }

  # add additional footprint of size 'offset' to each line of data
  if (!opt$nofootprint) {
    for (i in seq(1, nrow(dat), 1)) {
      # get actual transcription step for this entry of data
      step = dat[i, 3]
      # append opt$offset data points representing the polymerase footprint
      # from column (step + 1) + 1 to (step + 1) + 1 + offset - 1
      for (j in seq(step + 1, step + opt$offset, 1)) {
        dat[i, j + 3] <- -1
      }
    }
  }
} else if (opt$offset < 0) {
  # adjust y-axis labels
  y_axis_labels <- c()
  for (i in seq(opt$start, max_length, 10)) {
    y_axis_labels <- c(y_axis_labels, sprintf('%d (%d)', i, i + opt$offset))
  }
} else {
  y_axis_labels = seq(opt$start, max_length, 10)
}

dd <- melt(dat, id.vars=c("length", "sequence", "method"))

if (opt$normalize) {
    # normalize SHAPE reactivities
    lst <- sort(dd$value, decreasing=T)

    # skip top 2%
    n <- length(lst)
    average <- 0
    avg_cnt <- 0
    for (i in (2 * n / 100):(10 * n / 100)) {
      average <- average + lst[i]
      avg_cnt <- avg_cnt + 1
    }
    average <- average / avg_cnt

#    average <- (lst[1] + lst[2] + lst[3]) / 3

    for(i in 1:length(dd$value)) {
        dd$value[i] = dd$value[i]/average
    }
}

if (opt$SHAPE2probs) {
    # convert to probability
    for(i in 1:length(dd$value)) {
        dd$value[i] = SHAPE2probs(dd$value[i])
    }
}

if (is.null(opt$ltitle) || opt$SHAPE) {
  legend_title = expression(rho*" (Reactivity)")
} else {
  legend_title = opt$ltitle
}

if (is.null(opt$xlabel)) {
  xlabel = "Nucleotide Position"
} else {
  xlabel = opt$xlabel
}

if (is.null(opt$ylabel)) {
  ylabel = "Transcript Length"
} else {
  ylabel = opt$ylabel
}

if (opt$offset != 0) {
  ylabel = paste0(ylabel, " (w/o footprint)")
}


p <- ggplot(dd, aes(x=as.integer(variable), y=as.integer(length)))
p <- p + geom_raster(aes(fill=value), na.rm=T, alpha=1, interpolate=F)
if (!opt$SHAPE || opt$SHAPE2probs) {
#    p <- p + scale_fill_gradientn(
#                colours=c("darkblue", "steelblue", "seagreen","orange","yellow"),
#                na.value="transparent",
#                limits=c(0,1),
#                oob=squish, 
#                breaks=c(0, 0.5, 1),
#                labels=c("0","0.5","1"),
#                values=rescale(seq(0, 1, length.out=6)))
#    p <- p + scico::scale_fill_scico(
#                na.value="transparent",
#                palette = "roma",
#                limits=c(0,1),
#                breaks=c(0, 0.5, 1),
#                labels=c("0","0.5","1"),
#                values=rescale(seq(0, 1, length.out=6)))
    p <- p + scale_fill_viridis_c(
                option = "rocket",
                na.value="transparent",
#                oob=squish, 
                limits=c(0,1),
                breaks=c(0, 0.5, 1),
                labels=c("0","0.5","1"),
                values=rescale(seq(0, 1, length.out=6)))
} else {
    lab <- c(sprintf("%s",seq(0,SHAPE_range-1)))
#    p <- p + scale_fill_gradientn(
#                colours=c("darkblue", "steelblue", "seagreen","orange","yellow"),
#                na.value="transparent",
#                limits=c(0,SHAPE_range),
#                oob=squish,
#                breaks=seq(0, SHAPE_range),
#                labels=c(lab, paste(expression(">="), sprintf("%d", SHAPE_range))))
    p <- p + scale_fill_viridis_c(
                option = "rocket",
                na.value="transparent",
                limits=c(0,SHAPE_range),
                oob=squish,
                breaks=seq(0, SHAPE_range),
                labels=c(lab, paste(expression(">="), sprintf("%d", SHAPE_range))))
}

data_order = c("equilibrium","Kinfold","DrTrafo","SHAPE")
data_labels = c("RNAfold", "Kinfold", "DrTransformer","Experiment")
strand_order = c("SRPn", "SRPm", "SRPr", "SRPf")
strand_labels = c("SRP (wild type)", "SRP (U21C)", "SRP (U21C/C22U/G93A)", "SRP (U35C/U37C)")

p <- p + facet_grid(factor(method, levels=data_order, labels=data_labels)~factor(sequence, levels=strand_order, labels = strand_labels))

if (opt$offset < 0) {
  p <- p + scale_x_continuous(
                  breaks=seq(10, max_length, 10),
                  limits=c(0.5, max_length + 0.5),
                  expand=c(0,0))
} else {
  p <- p + scale_x_continuous(
                  breaks=seq(10, max_length + opt$offset, 10),
                  limits=c(0.5, max_length + opt$offset + 0.5),
#                  breaks=seq(10, max_length, 10),
#                  limits=c(0.5, max_length + 0.5),
                  expand=c(0,0))
}
p <- p + scale_y_reverse(
                limits=c(max_length+0.5, opt$start-0.5),
                breaks=y_axis_breaks,
                expand = c(0,0),
                labels=y_axis_labels)
p <- p + xlab(xlabel)
p <- p + ylab(ylabel)

if (!is.null(opt$title)){
  p <- p + ggtitle(opt$title)
}

p <- p + theme_bw()
p <- p + theme(plot.title = element_text(hjust = 0.5, size = 24),
                plot.background   = element_blank(),
                panel.spacing = unit(1, "lines"),
                panel.border = element_blank(),
                panel.background = element_blank(),
                panel.grid.major  = element_blank(),
                panel.grid.minor  = element_blank(),
                axis.line = element_blank(),
                strip.background = element_blank())

p <- p +  theme(axis.title.x  = element_text(family="Helvetica", size = 20, colour="#000000"),
                axis.text.x   = element_text(family="Helvetica", size = 16, colour="#555555", angle = 90.),
                axis.title.y  = element_text(family="Helvetica", size = 20, colour="#000000"),
                axis.text.y   = element_text(family="Helvetica", size = 16, colour="#555555"),
                strip.text.x  = element_text(family="Helvetica", size = 26, colour="#111111"),
                strip.text.y  = element_text(family="Helvetica", size = 26, colour="#111111")
               )

p <- p +  guides(fill = guide_colorbar(
                            title=legend_title,
                            title.position="left",
                            title.vjust=1,
                            barwidth = 25,
                            barheight = 1))
p <- p + theme(
            legend.position = "bottom",
#            legend.position = "top",
            legend.direction = "horizontal",
            legend.title = element_text(size = 22),
            legend.text = element_text(size = 16),
            legend.background = element_blank(),
            legend.key.height = unit(2,"line"),
            legend.key.width  = unit(2,"line"))

ggsave(file=opt$out, plot=p, width=18, height=12)
#ggsave(file=opt$out, plot=p, width=18, height=6)
