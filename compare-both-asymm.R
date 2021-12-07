#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# code to do various simulation-theory comparisons of mtDNA statistics
# 1. visualise V(h) and V(N) from theory and experiment for different parameterisations
# 2. compare individual moments
# 3. run statistical simulations (ie draws from model distributions) for some parameterisation and compare theory

# this code has different options that are invoked for the non-repulsive and repulsive cases of mtDNA distribution

cat("Processing inputs...\n")

if(length(args) < 7) {
  stop("Need: input file, whether to do repulsion (0/1), n, network mass, h to plot, ystar to plot, output label")
}

# default options
#args = c("output-asymm-0.25.csv", 0, 100, 50, 0.5, 0.25, "asymm-0.25")

inputfile = args[1]
repel = ifelse(args[2]==0, F, T)
n = as.numeric(args[3])
mass = as.numeric(args[4])
h.to.plot = as.numeric(args[5])
ystar.to.plot = as.numeric(args[6])
label = args[7]

cat("Loading libraries...\n")

# libraries for plotting
library(ggplot2)
library(gridExtra)

# function to predict V'(h) from model parameters using first and second order Taylor expansions
vhpred = function(alpha, beta, p, q, h, n, ystar) {
  # mean and variance of beta-distributed mass u
  eu = alpha/(alpha+beta)
  vu = (alpha*beta)/((alpha+beta)^2*(alpha+beta+1))
  pc = (1/3.141259)* (acos(ystar) - ystar*sqrt(1-ystar*ystar))
  
  # moments of u
  eu2 = alpha*(1+alpha) / ((alpha+beta)*(alpha+beta+1))
  eu3 = alpha*(1+alpha)*(2+alpha) / ((alpha+beta)*(alpha+beta+1)*(alpha+beta+2))
  eu4 = alpha*(1+alpha)*(2+alpha)*(3+alpha) / ((alpha+beta)*(alpha+beta+1)*(alpha+beta+2)*(alpha+beta+3))

  # (co)variances for network and cytoplasm w and m and their sums
  vmn = q*h*n*(eu - (vu+eu*eu)) + q*q*h*h*n*n*vu
  vwn = p*(1-h)*n*(eu - (vu+eu*eu)) + p*p*(1-h)*(1-h)*n*n*vu
  vmc = h*(1-q)*n*pc*(1-pc)
  vwc = (1-h)*(1-p)*n*pc*(1-pc)
  vw = vwn+vwc
  vm = vmn+vmc
  cwm = cmnwn = p*q*h*(1-h)*n*n*vu

  # higher-order covariances
  wn = p*(1-h)*n
  mn = q*h*n
  wc = (1-p)*(1-h)*n
  mc = (1-q)*h*n
  cw2m = wn*mn*(vu + (wn-1)*(eu3 - eu*eu2)) + 2*pc*wc*wn*mn*vu
  cwm2 = wn*mn*(vu + (mn-1)*(eu3 - eu*eu2)) + 2*pc*mc*wn*mn*vu
  cw3m = wn*mn*(vu + 3*(wn-1)*(eu3 - eu*eu2) + (wn^2 - 3*wn + 2)*(eu4 - eu*eu3)) +
    3*pc*wc*wn*mn*(vu + (wn-1)*(eu3 - eu*eu2)) +
    3*(wc*pc*(1-pc) + wc^2*pc^2)*wn*mn*vu
  cwm3 = wn*mn*(vu + 3*(mn-1)*(eu3 - eu*eu2) + (mn^2 - 3*mn + 2)*(eu4 - eu*eu3)) +
    3*pc*mc*wn*mn*(vu + (mn-1)*(eu3 - eu*eu2)) +
    3*(mc*pc*(1-pc) + mc^2*pc^2)*wn*mn*vu
  cw2m2 = wn*mn*(vu + (wn+mn-2)*(eu3 - eu*eu2) + (wn-1)*(mn-1)*(eu4 - eu2^2)) +
     2*pc*wc*wn*mn*(vu + (mn-1)*(eu3 - eu*eu2)) +
     2*pc*mc*wn*mn*(vu + (wn-1)*(eu3 - eu*eu2)) +
     4*pc*pc*wc*mc*wn*mn*vu
     
  # third and fourth central moments for networked and cytoplasmic w and m
  mwn3 = ( (h-1)*n*p*alpha*beta*(alpha-beta) * (2*(h-1)*(h-1)*n*n*p*p - 3*(h-1)*n*p*(alpha+beta) + (alpha+beta)*(alpha+beta)) ) / ( (alpha+beta)*(alpha+beta)*(alpha+beta)*(1 + alpha + beta)*(2 + alpha + beta)  )
#  mwc3 = (1-p)*(1-h)*n*(eu*(1-3*eu) + 2*eu*eu*eu)
 mwc3 = wc*(pc - 3*pc*pc + 2*pc*pc*pc)
 
  mwn4 = ( (1-h)*n*p*alpha*beta*( -3*(h-1)^3 * n^3 * p^3 * (alpha*(beta-2)*beta + 2*beta^2 + alpha^2*(beta+2)) +
    6*(h-1)^2 * n^2 * p^2 *( 2*alpha^2*beta^2 + 2*beta^3 + alpha*beta^3 + alpha^3*(beta+2)) +
    (alpha+beta)^3 * (alpha^2 + (beta-1)*beta - alpha*(1+4*beta) ) + (1-h)*n*p*(alpha+beta)^2 *
    (alpha^2*(7 + 3*beta) + beta*(-1+7*beta) + alpha*(-1 - 10*beta + 3*beta^2))) ) /
    ((alpha+beta)^4*(1+alpha+beta)*(2+alpha+beta)*(3+alpha+beta))
#  mwc4 = (1-p)*(1-h)*n*( eu - eu*eu + ( 3*(1-p)*(1-h)*n - 6)*(eu^2 - 2*eu^3 + eu^4) )
  mwc4 = wn*pc*(1-pc)*(1 + (3*wc-6)*pc*(1-pc))
  
  mmn3 = - ( h*n*q*alpha*beta*(alpha-beta) * (2*h*h*n*n*q*q - 3*h*n*q*(alpha+beta) + (alpha+beta)*(alpha+beta)) ) / ( (alpha+beta)*(alpha+beta)*(alpha+beta)*(1 + alpha + beta)*(2 + alpha + beta)  )
#  mmc3 = (1-q)*h*n*(eu*(1-3*eu) + 2*eu*eu*eu)
 mmc3 = mc*(pc - 3*pc*pc + 2*pc*pc*pc)
 
  mmn4 = ( h*n*q*alpha*beta*( 3*h^3 * n^3 * q^3 * (alpha*(beta-2)*beta + 2*beta^2 + alpha^2*(beta+2)) +
    6*h^2 * n^2 * q^2 *( 2*alpha^2*beta^2 + 2*beta^3 + alpha*beta^3 + alpha^3*(beta+2)) +
    (alpha+beta)^3 * (alpha^2 + (beta-1)*beta - alpha*(1+4*beta) ) + h*n*q*(alpha+beta)^2 *
    (alpha^2*(7 + 3*beta) + beta*(-1+7*beta) + alpha*(-1 - 10*beta + 3*beta^2))) ) /
    ((alpha+beta)^4*(1+alpha+beta)*(2+alpha+beta)*(3+alpha+beta))
#  mmc4 = (1-q)*h*n*( eu - eu*eu + ( 3*(1-q)*h*n - 6)*(eu^2 - 2*eu^3 + eu^4) )
mmc4 = mc*pc*(1-pc)*(1 + (3*mc-6)*pc*(1-pc))

  # central moments for total w and m
  mm3 = mmn3+mmc3
  mw3 = mwn3+mwc3
  mm4 = mmn4+mmc4
  mw4 = mwn4+mwc4

  mw = (1-h)*n*(p*eu + (1-p)*pc)
  mm = h*n*(q*eu + (1-q)*pc)
  mh = mm/(mw+mm)

  mn = mw+mm
  vn = pc*(1-pc)*n + (p*(1-h)+q*h)*n*((p*(1-h)+q*h)*n - 1)*vu
  
  # derivatives of h for use in Taylor expansions
  hm = (1-h)/(eu*n)
  hw = -h/(eu*n) #XXX FIX
  hmm = -2*(1-h)/(eu*eu*n*n)
  hww = 2*h/(eu*eu*n*n)
  hmw = 2*h/(eu*eu*n*n) - 2/(eu*n*n)
kappa = p*(1-h)+q*h

hm = mw/(mm+mw)^2
hw = -mm/(mm+mw)^2
if(FALSE){
hm = (wn*eu + wc*pc)/(kappa*n*eu + (1-kappa)*n*pc)^2
hw = -(mn*eu + mc*pc)/(kappa*n*eu + (1-kappa)*n*pc)^2
hmm = -2*(wn*eu+wc*pc)/(kappa*n*eu + (1-kappa)*n*pc)^3
hww = 2*(mn*eu+mc*pc)/(kappa*n*eu + (1-kappa)*n*pc)^3
hmw = ((mn-wn)*eu + (mc-wc)*pc)/(kappa*n*eu + (1-kappa)*n*pc)^3
}

# first-order
  v1 = hm*hm*vm + hw*hw*vw + 2*hw*hm*cwm

  # second-order = v1 + this correction
  v2 = hm*hmm*mm3 + hw*hww*mw3 + (2*hw*hmw + hm*hww)*cw2m + (2*hm*hmw + hw*hmm)*cwm2 +
    hmm*hmm*mm4/4 + hww*hww*mw4/4 + (hmw*hmw + hmm*hww/2)*cw2m2 +
    hmw*hww*cw3m + hmw*hmm*cwm3

  vn = vw + vm + 2*cwm

  # return a long vector with different statistics, mainly ordered to match ordering of the simulation output
  # return zeros for statistics not addressed here
  return(c(v1, v1+v2, vn, mw,vw,mm,vm,cwm,cw2m,cw3m,cwm2,cwm3,cw2m2,mw3,mm3,mw4,mm4,mh,v1,0,0,eu*mass,vu*mass*mass, mn, vn))
}

# function to predict V'(h) from model parameters using first and second order Taylor expansions
# in the repulsive mtDNA case (modelled using hypergeometric distribution)
vhpredrepel = function(alpha, beta, p, q, h, n) {
  # mean and variance of beta-distributed mass u
  eu = alpha/(alpha+beta)
  vu = (alpha*beta)/((alpha+beta)^2*(alpha+beta+1))

  # E(U^2)
  eu2 = alpha*(1+alpha) / ((alpha+beta)*(alpha+beta+1))

  # means, variances, and covariances of the different mtDNA counts
  mwn = p*(1-h) * eu* n
  mmn = q*h * eu*n
  cwm = cmnwn = p*(1-h)*q*h*n*n*eu2 - p*(1-h)*q*h*n*n*eu^2

  vwn = p*(1-h)*q*h*n*n / ((p*(1-h)+q*h)*n -1) *(eu-eu2) + p^2*(1-h)^2*n^2*vu
  vmn = vu*((p*(1-h)+q*h)*n)^2 + vwn - 2*p*(1-h)*n*(p*(1-h)+q*h)*n*vu

  vmc = h*(1-q)*n*eu*(1-eu)
  vwc = (1-h)*(1-p)*n*eu*(1-eu)
  vw = vwn+vwc
  vm = vmn+vmc

  mw = (1-h)*n*eu
  mm = h*n*eu
  mh = mm/(mw+mm)

  # derivatives of h for use in Taylor expansions
  hm = (1-h)/(eu*n)
  hw = -h/(eu*n)

  # first-order
  v1 = hm*hm*vm + hw*hw*vw + 2*hw*hm*cwm

  mn = mw+mm
  vn = vw + vm + 2*cwm

  # return a long vector with different statistics, mainly ordered to match ordering of the simulation output
  # return zeros for statistics not addressed here
  return(c(v1, 0, vn, mw,vw,mm,vm,cwm,0,0,0,0,0,0,0,0,0,mh,v1,0,0,eu*mass,vu*mass*mass, mn, vn))
}

############### compare physical simulation output to theory

cat("Reading data...\n")

# read simulation output
df = read.csv(inputfile, header=T)

# for now, ignore diffusion and repulsion
if(repel == T) {
  test.df = df[df$lambda==0&df$halo==0.1,]
} else {
  test.df = df[df$lambda==0&df$halo==0,]
}

# uncomment for old sim code where moments were accidentally not normalised
#nit = 100
#test.df$muw3 = test.df$muw3 / nit
#test.df$muw4 = test.df$muw4 / nit
#test.df$mum3 = test.df$mum3 / nit
#test.df$mum4 = test.df$mum4 / nit

cat("Computing predictions...\n")

# initialise a dataframe storing theory-simulation comparison for each parameterisation
res.1.df = data.frame()

# loop through network parameterisation
for(i in 1:nrow(test.df)) {

  # get the u stats for this parameterisation and fit beta distribution parameters to its mean and variance
  m = test.df$mu[i]/mass
  v = test.df$vu[i]/(mass*mass)
  alpha = (m^2 - m^3 - m*v)/v
  beta = (m-1)*(m^2-m+v)/v

  # compute prediction for this parameterisation
  if(repel == T) {
    pred = vhpredrepel(alpha, beta, test.df$p[i], test.df$q[i], test.df$h[i], n)
  } else {
    pred = vhpred(alpha, beta, test.df$p[i], test.df$q[i], test.df$h[i], n, test.df$ystar[i])
  }
  # store comparison in results dataframe
  res.1.df = rbind(res.1.df, data.frame(ystar=test.df$ystar[i], seeds=test.df$seeds[i], p=test.df$p[i], q=test.df$q[i], h=test.df$h[i], vhsim = test.df$vh[i], vhpred1 = pred[1], vhpred2 = pred[2], vnsim=test.df$vn[i], vnpred=pred[3]))
}

cat("Plotting summaries...\n")

# plot simulation-theory comparison for first-order (1) and second-order (2) models
g1.1 = ggplot(res.1.df, aes(x=vhsim, y=vhpred1, color=factor(h))) + geom_point() + geom_abline(slope=1, intercept=0) + ggtitle("Physical sim vs first-order theory")
g1.2 = ggplot(res.1.df, aes(x=vhsim, y=vhpred2, color=factor(h))) + geom_point() + geom_abline(slope=1, intercept=0) + ggtitle("Physical sim vs second-order theory")
g1.3 = ggplot(res.1.df, aes(x=vnsim, y=vnpred, color=factor(h))) + geom_point() + geom_abline(slope=1, intercept=0) + ggtitle("Physical sim vs theory")

# tile plots
res.1.df$plot.simhet = res.1.df$vhsim/(res.1.df$h*(1-res.1.df$h)) - 1/n
res.1.df$plot.predhet1 = res.1.df$vhpred1/(res.1.df$h*(1-res.1.df$h)) - 1/n
res.1.df$plot.predhet2 = res.1.df$vhpred2/(res.1.df$h*(1-res.1.df$h)) - 1/n
res.1.df$plot.simn = res.1.df$vnsim - 0.5^2*n
res.1.df$plot.predn = res.1.df$vnpred - 0.5^2*n

colfn = scale_fill_gradientn(colours = c("black", "blue", "white", "red", "black"), values = c(0, 0.05, 0.1, 0.2, 1), limits = c(-0.01,0.09))
colfn2 = scale_fill_gradientn(colours = c("black", "blue", "white", "red", "black"), values = c(0, 0.05, 0.1, 0.2, 1), limits = c(-30,270))

plot.h = h.to.plot
plot.ystar = ystar.to.plot
plot.df = res.1.df[res.1.df$h==plot.h & res.1.df$ystar==plot.ystar,]
maxvh = max(plot.df$vhsim)

g3.1 = ggplot(plot.df, aes(x=p, y=q, fill=plot.simhet)) + geom_tile() + colfn + facet_wrap( ~ seeds)
g3.2 = ggplot(plot.df, aes(x=p, y=q, fill=plot.predhet1)) + geom_tile() + colfn + facet_wrap( ~ seeds)
g3.3 = ggplot(plot.df, aes(x=p, y=q, fill=plot.predhet2)) + geom_tile() + colfn + facet_wrap( ~ seeds)

filename = paste(c(label, "-vh-stats-", ifelse(repel, "repel.png", "no-repel.png")), collapse="")
res.factor = 3
png(filename, width=1000*res.factor, height=1000*res.factor, res=72*res.factor)
grid.arrange(g3.1, g3.2, g3.3, nrow=3)
dev.off()

maxvn = max(plot.df$vnsim)

g4.1 = ggplot(plot.df, aes(x=p, y=q, fill=plot.simn)) + geom_tile() + colfn2 + facet_wrap(~ seeds)
g4.2 = ggplot(plot.df, aes(x=p, y=q, fill=plot.predn)) + geom_tile() + colfn2 + facet_wrap(~ seeds)

res.factor = 3
filename = paste(c(label, "-vn-stats-", ifelse(repel, "repel.png", "no-repel.png")), collapse="")
png(filename, width=1000*res.factor, height=1000*res.factor, res=72*res.factor)
grid.arrange(g4.1, g4.2, nrow=2)
dev.off()

plot.df = res.1.df
maxvh = max(plot.df$vhsim)

g5.1 = ggplot(plot.df, aes(x=p, y=q, fill=plot.simhet)) + geom_tile() + colfn + facet_grid(h ~ seeds)
g5.2 = ggplot(plot.df, aes(x=p, y=q, fill=plot.simn)) + geom_tile() + colfn2 + facet_grid(h ~ seeds)
g5.3 = ggplot(plot.df, aes(x=p, y=q, fill=plot.predhet1)) + geom_tile() + colfn + facet_grid(h ~ seeds)
g5.4 = ggplot(plot.df, aes(x=p, y=q, fill=plot.predn)) + geom_tile() + colfn2 + facet_grid(h ~ seeds)

res.factor = 3
filename = paste(c(label, "-both-stats-", ifelse(repel, "repel.png", "no-repel.png")), collapse="")
png(filename, width=1000*res.factor, height=800*res.factor, res=72*res.factor)
grid.arrange(g5.1, g5.2, g5.3, g5.4, nrow=2)
dev.off()

############### moments

cat("Plotting moments...\n")

# get list of statistics to plot
stats = colnames(test.df[,8:ncol(test.df)])
if(repel == F) {
  stats.to.plot = stats[stats!="md"&stats!="vd"]
} else {
  stats.to.plot = c("mw", "mm", "vw", "vm", "mn", "vn", "cwm", "mh", "vh", "vu")
}
mom.df = data.frame()

# loop through network parameterisation
for(i in 1:nrow(test.df)) {

  # get the u stats for this parameterisation and fit beta distribution parameters to its mean and variance
  m = test.df$mu[i]/mass
  v = test.df$vu[i]/(mass*mass)
  alpha = (m^2 - m^3 - m*v)/v
  beta = (m-1)*(m^2-m+v)/v

  # compute prediction for this parameterisation
  pred = vhpred(alpha, beta, test.df$p[i], test.df$q[i], test.df$h[i], n, test.df$ystar[i])

  # populate list of statistics that we're interested in
  for(j in 1:length(stats)) {
    if(stats[j] %in% stats.to.plot) {
    mom.df = rbind(mom.df, data.frame(ystar=test.df$ystar[i],seeds=test.df$seeds[i], p=test.df$p[i], q=test.df$q[i], h=test.df$h[i], stat = stats[j], predval = pred[3+j], simval = test.df[i,7+j]))
    }
  }
}

filename = paste(c(label, "-moments-", ifelse(repel, "repel.png", "no-repel.png")), collapse="")
png(filename, width=1000*res.factor, height=1000*res.factor, res=72*res.factor)
ggplot(mom.df, aes(x=predval,y=simval,colour=factor(ystar))) + geom_point() + geom_abline(slope=1, intercept=0) + facet_wrap(~stat, scales="free")
dev.off()

############### compare statistical simulation output to theory

# initialise a dataframe storing theory-simulation comparison for each parameterisation
cat("Producing statistical simulation...\n")

res.2.df = data.frame(alpha=NULL,n=NULL,p=NULL,q=NULL,h=NULL,vhsim=NULL,vhpred1=NULL,vhpred2=NULL)

# number of samples for each parameterisation
nsamp = 100

# loop through various parameters
for(alpha in c(0.1, 1, 10, 100)) {
  for(n in c(20, 50, 100)) {
    for(p in (0:5)/5) {
      for(q in (0:5)/5) {
        for(h in (0:5)/5) {
	  # for now, enforce beta = alpha for E(u) = 0.5
          beta = alpha
	  ystar = 0
          pc = (1/3.141259)* (acos(ystar) - ystar*sqrt(1-ystar*ystar))
  
	  # count of each mtDNA type
          wn0 = round((1-h)*p*n)
          mn0 = round(h*q*n)
          wc0 = round((1-h)*(1-p)*n)
          mc0 = round(h*(1-q)*n)
	  # simulate random variables
          u = rbeta(nsamp, alpha, beta)
          wc = rbinom(nsamp, wc0, pc)
          mc = rbinom(nsamp, mc0, pc)
          wn = unlist(lapply(u, function(x) { return(rbinom(1, wn0, x))}))
          mn = unlist(lapply(u, function(x) { return(rbinom(1, mn0, x))}))
	  # statistics
          w = wn+wc
          m = mn+mc
          het = m/(w+m)
	  ndna = w+m

          # predict V(h) and store comparison to simulation
          pred = vhpred(alpha, beta, p, q, h, n, ystar)
	  res.2.df = rbind(res.2.df, data.frame(alpha=alpha, n=n, p=p, q=q, h=h, vhsim = var(het), vhpred1 = pred[1], vhpred2 = pred[2], vnsim = var(ndna), vnpred = pred[3]))
	}
      }
    }
  }
}
cat("Plotting comparisons...\n")

# plot simulation-theory comparison for first-order (1) and second-order (2) models
g2.1 = ggplot(res.2.df, aes(x=vhsim, y=vhpred1, color=factor(h))) + geom_point() + geom_abline(slope=1, intercept=0) + facet_grid(alpha ~ n) + ggtitle("Stat sim vs first-order theory")
g2.2 = ggplot(res.2.df, aes(x=vhsim, y=vhpred2, color=factor(h))) + geom_point() + geom_abline(slope=1, intercept=0) + facet_grid(alpha ~ n) + ggtitle("Stat sim vs second-order theory")
g2.3 = ggplot(res.2.df, aes(x=vnsim, y=vnpred, color=factor(h))) + geom_point() + geom_abline(slope=1, intercept=0) + facet_grid(alpha ~ n) + ggtitle("Stat sim vs theory")

res.factor = 3
filename = paste(c(label, "-sim-theory-", ifelse(repel, "repel.png", "no-repel.png")), collapse="")
png(filename, width=1600*res.factor, height=800*res.factor, res=72*res.factor)
grid.arrange(g1.1, g1.2, g1.3, g2.1, g2.2, g2.3, nrow=2)
dev.off()