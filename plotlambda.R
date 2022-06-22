#!/usr/bin/env Rscript

# code to plot influence of post-placement diffusion lambda on mtDNA statistics

args = commandArgs(trailingOnly=TRUE)

cat("Processing inputs...\n")

if(length(args) < 1) {
  stop("Need: input file")
}

inputfile = args[1]

cat("Loading libraries...\n")
library(ggplot2)
library(gridExtra)

# read simulation output
df = read.csv(inputfile, header=T)

# get normalised heteroplasmy variance
df$het.norm = df$vh / (df$mh*(1-df$mh))

# construct dataframe of results to plot
# that is, extract traces of V'(h) with lambda for different parameterisations (labelled by string "expt")
to.plot = df[df$h==0.5,]
plot.df = data.frame(lambda = to.plot$lambda, Vhprime = to.plot$het.norm, seeds = to.plot$seeds, halo = to.plot$halo, expt = NA)
for(i in 1:nrow(to.plot)) {
  a = to.plot[i,]
  plot.df$expt[i] = paste(c(a$seeds,a$p,a$q,a$halo), collapse=",")
}

plot.df$seeds = factor(plot.df$seeds, levels = c(4,16,64), labels = paste("Seeds, s = ", c(4,16,64)))
plot.df$halo = factor(plot.df$halo, levels = c(0.0,0.1), labels = paste("Repulsive halo, l = ",c(0.0,0.1)))

res.factor = 3
png("plotlambda.png", width=res.factor*720, height=res.factor*560, res=72*res.factor)
ggplot(plot.df, aes(x=lambda, y=Vhprime, colour=expt)) + geom_line() + theme(legend.position="none") + 
  facet_grid(seeds ~ halo) +
  labs(x = "λ", y = "V'(h)")+
  theme(axis.text = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1)),
        strip.text = element_text(size = rel(1)))
dev.off()
