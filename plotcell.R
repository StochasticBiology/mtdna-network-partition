#!/usr/bin/env Rscript
# code to plot individual cell snapshots that are output from simulation

cat("Loading libraries...")

library(ggplot2)
library(gridExtra)

# helper function to draw a circle using ggplot
# this is joran's solution from https://stackoverflow.com/questions/6862742/draw-a-circle-with-ggplot2
circleFun <- function(center = c(0,0),diameter = 2, npoints = 100){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
}
circledat <- circleFun()

# function to produce a visualisation given a particular label for simulation output 
makeplot = function(string, titlestr="", tag = "") {
  # read network and mtDNA position data
  net.df = read.csv(paste(c("network-", string, ".csv"), collapse=""))
  dna.df = read.csv(paste(c("mtdna-", string, ".csv"), collapse=""))

  # return the plot
  return(ggplot() +
    geom_path(data = circledat, aes(x,y)) +                                         # cell boundary
    geom_segment(data = net.df, aes(x=xs,y=ys,xend=xe,yend=ye), color="#888888") +  # network
    geom_point(data = dna.df, aes(x=x, y=y, color=factor(type))) +                  # mtDNAs
    theme_void() +
    theme(legend.position = "none") +
    ggtitle(titlestr)+
    labs(tag = tag)+
    theme(plot.title = element_text(size = rel(1.25))))
}

# strings in output filenames follow the format
# h, nseed, p, q, lambda, halo, perturb
plot.1 = makeplot("0.5-4-0.00-0.00-0.00-0.00-2", "Small s,\n p=q=0", "A")
plot.2 = makeplot("0.5-4-0.50-0.50-0.00-0.00-2", "small s,\n p=q=0.5", "B")
plot.3 = makeplot("0.5-4-1.00-0.00-0.00-0.00-2", "Small s,\n p=1,q=0", "C")
plot.4 = makeplot("0.5-16-0.50-0.50-0.00-0.00-2", "Medium s,\n p=q=0.5", "D")
plot.5 = makeplot("0.5-64-0.50-0.50-0.00-0.00-2", "Large s,\n p=q=0.5", "E")
plot.6 = makeplot("0.5-16-1.00-1.00-0.00-0.00-2", "Medium s,\n p=q=1", "F")
plot.7 = makeplot("0.5-16-1.00-1.00-0.04-0.00-2", "Medium s,\n p=q=1,\n Medium λ", "G")
plot.8 = makeplot("0.5-16-1.00-1.00-0.10-0.00-2", "Medium s,\n p=q=1,\n Large λ","H")
plot.9 = makeplot("0.5-16-1.00-1.00-0.00-0.10-2", "Medium s,\n p=q=1,\n repel","I")

# bump to output file
res.factor = 3
png("plotcell.png", width=1200*res.factor, height=500*res.factor, res=72*res.factor)
grid.arrange(plot.1, plot.2, plot.3, plot.4, plot.5, plot.6, plot.7, plot.8, plot.9, nrow=2)
dev.off()
