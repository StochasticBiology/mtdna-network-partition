#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# code to do various simulation-theory comparisons of mtDNA statistics
# 1. visualise V(h) and V(N) from theory and experiment for different parameterisations
# 2. compare individual moments
# 3. run statistical simulations (ie draws from model distributions) for some parameterisation and compare theory

# this code has different options that are invoked for the non-repulsive and repulsive cases of mtDNA distribution

cat("Processing inputs... \n")

if(length(args) < 7) {
  stop("Need: input file, whether to do repulsion (0/1), n, network mass, h to plot, prop to plot, output label")
}

# default options
#args = c("output-asymm-0.25.csv", 0, 100, 50, 0.5, 0.25, "asymm-0.25")

inputfile = args[1]
repel = ifelse(args[2]==0, F, T)
n = as.numeric(args[3])
mass = as.numeric(args[4])
h.to.plot = as.numeric(args[5])
prop.to.plot = as.numeric(args[6])
label = args[7]

cat("Loading libraries...\n")

# libraries for plotting
library(ggplot2)
library(gridExtra)

# function to predict V'(h) from model parameters using first and second order Taylor expansions
vhpred = function(alpha, beta, p, q, h, n, prop) {
  # mean and variance of beta-distributed mass u
  eu = alpha/(alpha+beta)
  vu = (alpha*beta)/((alpha+beta)^2*(alpha+beta+1))
  pc = prop
  
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
  mwn3 = ( p*(1-h)*n*alpha*beta*(beta-alpha) * (2*(p*(1-h)*n)*(p*(1-h)*n) + 3*p*(1-h)*n*(alpha+beta) + (alpha+beta)*(alpha+beta)) ) / ( (alpha+beta)*(alpha+beta)*(alpha+beta)*(1 + alpha + beta)*(2 + alpha + beta)  )
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
  #hm = (1-h)/(eu*n)
  #hw = -h/(eu*n) #XXX FIX
  #hmm = -2*(1-h)/(eu*eu*n*n)
  #hww = 2*h/(eu*eu*n*n)
  #hmw = 2*h/(eu*eu*n*n) - 2/(eu*n*n)

  # Derivatives with eu and pc
  hm  = (wn*eu+wc*pc)/((wn+mn)*eu+(wc+mc)*pc)^2
  hw  =-(mn*eu+mc*pc)/((wn+mn)*eu+(wc+mc)*pc)^2
  hmm =-2*(wn*eu+wc*pc)/((wn+mn)*eu+(wc+mc)*pc)^3
  hww = 2*(mn*eu+mc*pc)/((wn+mn)*eu+(wc+mc)*pc)^3
  hmw = ((mn-wn)*eu+(mc-wc)*pc)/((wn+mn)*eu+(wc+mc)*pc)^3
  
  kappa = p*(1-h)+q*h

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
vhpredrepel = function(alpha, beta, p, q, h, n, prop) {
  # mean and variance of beta-distributed mass u
  eu = alpha/(alpha+beta)
  vu = (alpha*beta)/((alpha+beta)^2*(alpha+beta+1))
  pc = prop

  # E(U^2)
  eu2 = alpha*(1+alpha) / ((alpha+beta)*(alpha+beta+1))
	
  # higher-order covariances
  wn = p*(1-h)*n
  mn = q*h*n
  wc = (1-p)*(1-h)*n
  mc = (1-q)*h*n

  # means, variances, and covariances of the different mtDNA counts
  mwn = p*(1-h) * eu* n
  mmn = q*h * eu*n
  cwm = cmnwn = p*(1-h)*q*h*n*n*eu2 - p*(1-h)*q*h*n*n*eu^2

  vwn = ( p*(1-h)*q*h*n*n / ((p*(1-h)+q*h)*n -1) ) *(eu-eu2) + p^2*(1-h)^2*n^2*vu
  #vwm = ( 1/(0.01*0.01*(wn+mn)*(wn+mn)) )*( ( p*(1-h)*q*h*n*n / ((p*(1-h)+q*h)*n -1) ) *(eu-eu2) + p^2*(1-h)^2*n^2*vu )
  vmn = vu*((p*(1-h)+q*h)*n)^2 + vwn - 2*p*(1-h)*n*(p*(1-h)+q*h)*n*vu
	#vmn = 

  vmc = h*(1-q)*n*pc*(1-pc)
  vwc = (1-h)*(1-p)*n*pc*(1-pc)
  vw = vwn+vwc
  vm = vmn+vmc

  mw = p*(1-h)*n*eu+(1-p)*(1-h)*n*pc
  mm = q*h*n*eu+(1-q)*h*n*pc
  mh = mm/(mw+mm)

  # derivatives of h for use in Taylor expansions
  #hm = (1-h)/(eu*n)
  #hw = -h/(eu*n)
  hm = (wn*eu+wc*pc)/((wn+mn)*eu+(wc+mc)*pc)^2
  hw =-(mn*eu+mc*pc)/((wn+mn)*eu+(wc+mc)*pc)^2

  # first-order
  v1 = hm*hm*vm + hw*hw*vw + 2*hw*hm*cwm

  mn = mw+mm
  vn = vw + vm + 2*cwm

  # return a long vector with different statistics, mainly ordered to match ordering of the simulation output
  # return zeros for statistics not addressed here
  return(c(v1, 0, vn, mw,vw,mm,vm,cwm,0,0,0,0,0,0,0,0,0,mh,v1,0,0,eu*mass,vu*mass*mass, mn, vn))
}

# Function to predict V(h) from E(h^2)-E(h)^2 with statistical simulations
vhest = function(p, q, h, n, pc, alpha, beta){
  wn0 = round(p*(1-h)*n)
  mn0 = round(q*h*n)
  wc0 = round((1-p)*(1-h)*n)
  mc0 = round((1-q)*h*n)
  
  Wn = 0:wn0
  Mn = 0:mn0
  Wc = 0:wc0
  Mc = 0:mc0
  
  INET = length(Wn)
  JNET = length(Mn)
  ICYT = length(Wc)
  JCYT = length(Mc)
  
  Eh  = 0
  Eh2 = 0
  En  = 0
  En2 = 0
  
  # Precalculate this, and grab element icyt,jcyt as you go
  dist.cyt =  dbinom(x = Wc, size = wc0, prob = pc)%*%
    t(dbinom(x = Mc, size = mc0, prob = pc))
  # For summing over u
	u = seq(from = 0, to = 1, length.out = 100)
  dist.u = dbeta(u,alpha,beta)
  deltau = u[2]-u[1]
  for(inet in 1:INET){
    dist.wn = dbinom(Wn[inet],wn0,u)
    for(jnet in 1:JNET){
      dist.mn = dbinom(Mn[jnet],mn0,u)
      for(icyt in 1:ICYT){
        for(jcyt in 1:JCYT){
          if(inet + jnet + icyt + jcyt == 4){
            # do nothing here!
          }else{
            dist.net = sum(deltau*dist.wn*dist.mn*dist.u)
            prefactor = (Mn[jnet]+Mc[jcyt])/(Wn[inet]+Wc[icyt]+Mn[jnet]+Mc[jcyt])
            
            Eh  = Eh  + prefactor*dist.cyt[icyt,jcyt]*dist.net
            Eh2 = Eh2 + prefactor^2*dist.cyt[icyt,jcyt]*dist.net
            
            En2 = En2 + (Wn[inet]+Wc[icyt]+Mn[jnet]+Mc[jcyt])^2*dist.net*dist.cyt[icyt,jcyt]
            En  = En  + (Wn[inet]+Wc[icyt]+Mn[jnet]+Mc[jcyt])*dist.net*dist.cyt[icyt,jcyt]
          }
        }
      }
    }
  }
  vn      = En2-En^2
  vhprime = (Eh2-Eh^2)/(h*(1-h))
  return(c(alpha, beta, vhprime, vn))
}

############### compare physical simulation output to theory

cat("Reading data...\n")

# read simulation output
df = read.csv(inputfile, header=T)

# for now, set diffusion to zero
if(repel == T) {
  #tot = nrow(df[df$lambda==0&df$halo==0.1,])
  test.df = na.omit(df[df$lambda==0&df$halo==0.1,])
  #alive = nrow(test.df)
  #noDead = tot - alive
  #cat(paste(c("noDead = ", toString(noDead), "\n"), sep = ""))
}else{
  #tot = nrow(df[df$lambda==0&df$halo==0.1,])
  test.df = na.omit(df[df$lambda==0&df$halo==0,])
  #alive = nrow(test.df)
  #noDead = tot - alive
  #cat(paste(c("noDead = ", toString(noDead), "\n"), sep = ""))
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

# initialize a dataframe storing parameters
df.par = data.frame(seeds = NULL, alpha = NULL, beta = NULL)

# loop through network parameterisation
for(i in 1:nrow(test.df)) {

	# get the u stats for this parameterisation and fit beta distribution parameters to its mean and variance
	m = test.df$mu[i]/mass
	v = test.df$vu[i]/(mass*mass)
	alpha = (m^2 - m^3 - m*v)/v
	beta = (m-1)*(m^2-m+v)/v

	# compute parameter set for this parameterisation
	df.par = rbind(df.par,data.frame(seeds = test.df$seeds[i], alpha = alpha, beta = beta))
	#cat(paste(c("prop = ", toString(prop.to.plot), ", ",toString(test.df$prop[i]),"\n"),sep = ""))
	# compute prediction for this parameterisation
	if(repel == T) {
		pred = vhpredrepel(alpha, beta, test.df$p[i], test.df$q[i], test.df$h[i], n, test.df$prop[i])
	} else {
		pred = vhpred(alpha, beta, test.df$p[i], test.df$q[i], test.df$h[i], n, test.df$prop[i])
	}
	# store comparison in results dataframe
	res.1.df = rbind(res.1.df, data.frame(prop=test.df$prop[i], seeds=test.df$seeds[i], p=test.df$p[i], q=test.df$q[i], h=test.df$h[i], vhsim = test.df$vh[i], vhpred1 = pred[1], vhpred2 = pred[2], vnsim=test.df$vn[i], vnpred=pred[3]))
}

#cat(paste(c("alpha = ", toString(round(mean(df.par$alpha[df.par$seeds == 4]), 2)), "\n"), sep = ""))
#cat(paste(c("alpha = ", toString(round(mean(df.par$alpha[df.par$seeds == 16]), 2)), "\n"), sep = ""))
#cat(paste(c("alpha = ", toString(round(mean(df.par$alpha[df.par$seeds == 64]), 2)), "\n"), sep = ""))
#cat(paste(c("beta  = ", toString(round(mean(df.par$beta[df.par$seeds == 4]), 2)), "\n"), sep = ""))
#cat(paste(c("beta  = ", toString(round(mean(df.par$beta[df.par$seeds == 16]), 2)), "\n"), sep = ""))
#cat(paste(c("beta  = ", toString(round(mean(df.par$beta[df.par$seeds == 64]), 2)), "\n"), sep = ""))

Alpha = round(c(mean(df.par$alpha[df.par$seeds == 4]),mean(df.par$alpha[df.par$seeds == 16]),mean(df.par$alpha[df.par$seeds == 64])),1)
Beta = round(c(mean(df.par$beta[df.par$seeds == 4]),mean(df.par$beta[df.par$seeds == 16]),mean(df.par$beta[df.par$seeds == 64])),1)
cat("Plotting summaries...\n")

# plot simulation-theory comparison for first-order (1) and second-order (2) models
if(repel){
	g1.1 = ggplot(res.1.df, aes(x=vhsim, y=vhpred1, color=factor(h))) + geom_point() + geom_abline(slope=1, intercept=0) + ggtitle("Physical sim vs first-order theory")
	#g1.2 = ggplot(res.1.df, aes(x=vhsim, y=vhpred2, color=factor(h))) + geom_point() + geom_abline(slope=1, intercept=0) + ggtitle("Physical sim vs second-order theory")
	g1.3 = ggplot(res.1.df, aes(x=vnsim, y=vnpred, color=factor(h))) + geom_point() + geom_abline(slope=1, intercept=0) + ggtitle("Physical sim vs theory")	
} else {
	g1.1 = ggplot(res.1.df, aes(x=vhsim, y=vhpred1, color=factor(h))) + geom_point() + geom_abline(slope=1, intercept=0) + ggtitle("Physical sim vs first-order theory")
	g1.2 = ggplot(res.1.df, aes(x=vhsim, y=vhpred2, color=factor(h))) + geom_point() + geom_abline(slope=1, intercept=0) + ggtitle("Physical sim vs second-order theory")
	g1.3 = ggplot(res.1.df, aes(x=vnsim, y=vnpred, color=factor(h))) + geom_point() + geom_abline(slope=1, intercept=0) + ggtitle("Physical sim vs theory")
}


# tile plots
res.1.df$plot.simhet = res.1.df$vhsim/(res.1.df$h*(1-res.1.df$h))
res.1.df$plot.predhet1 = res.1.df$vhpred1/(res.1.df$h*(1-res.1.df$h)) 
res.1.df$plot.predhet2 = res.1.df$vhpred2/(res.1.df$h*(1-res.1.df$h))
res.1.df$plot.simn = res.1.df$vnsim
res.1.df$plot.predn = res.1.df$vnpred

plot.h = h.to.plot
plot.prop = prop.to.plot
plot.df = res.1.df[res.1.df$h==plot.h & res.1.df$prop==plot.prop,]

vhsimmax = max(plot.df$plot.simhet)
vhpredmax = max(c(max(plot.df$plot.predhet1),max(plot.df$plot.predhet2)))
vnsimmax = max(plot.df$plot.simn)
vnpredmax = max(plot.df$plot.predn)

cat(paste0("Repel = ", paste0(toString(repel),"\n")))
cat(paste0("(hplot,prop.plot) = (", paste0(toString(c(plot.h,plot.prop)),")\n")))
cat(paste0("vh max sim = ", paste0(toString(round(vhsimmax,3)),"\n")))
cat(paste0("vn max sim = ", paste0(toString(round(vnsimmax)),"\n")))
cat(paste0("vh max pred = ", paste0(toString(round(vhpredmax,3)),"\n")))
cat(paste0("vn max pred = ", paste0(toString(round(vnpredmax)),"\n")))

if(repel){
	vnmax = max(c(max(plot.df$plot.simn),max(plot.df$plot.predn)))
	vhmax = max(c(max(plot.df$plot.simhet),max(plot.df$plot.predhet1)))
	#cat(paste(c("vhmax, vnmax = ",toString(vhmax), ", ", toString(vnmax), "\n"),sep=""))
	pc = plot.prop

	# Color gradient comparing to asymmetric symmetric without network, i.e.,  (1/n)
	# This is a really ugly setup, which causes rather uniform redness across heatmaps
	#vhl = (1/n)/vhmax

	# Color gradient to compare asymmetric to asymmetric without network, i.e., ((1-pc)/pc)*(1/n)
	vhl = ((1-pc)/pc)*(1/n)/vhmax 
	colfn = scale_fill_gradientn(colors = c("black","blue","white","red","black"),
						values = c(0,vhl/2,vhl,2*vhl,1),  limits = c(0,vhmax))

	vnl = pc*(1-pc)*n/vnmax
	colfn2 = scale_fill_gradientn(colors = c("black","blue","white","red","black"),
							values = c(0,vnl/2,vnl,2*vnl,1), limits = c(0,vnmax))
	
	g3.1 = ggplot(plot.df, aes(x=p, y=q, fill=plot.simhet)) + geom_tile() + colfn + facet_wrap( ~ seeds)+
		theme(axis.title = element_text(size = rel(1.25)),
		strip.text = element_text(size = rel(1.25)))
	g3.2 = ggplot(plot.df, aes(x=p, y=q, fill=plot.predhet1)) + geom_tile() + colfn + facet_wrap( ~ seeds)+
		theme(axis.title = element_text(size = rel(1.25)),
		strip.text = element_text(size = rel(1.25)))	
} else {
	vnmax = max(c(max(plot.df$plot.simn),max(plot.df$plot.predn)))
	vhmax = max(c(max(plot.df$plot.simhet),max(plot.df$plot.predhet1),max(plot.df$plot.predhet2)))

	pc = plot.prop

	# Color gradient comparing to asymmetric symmetric without network, i.e.,  (1/n)
	# This is a really ugly setup, which causes rather uniform redness across heatmaps
	#vhl = (1/n)/vhmax

	# Color gradient to compare asymmetric to asymmetric without network, i.e., ((1-pc)/pc)*(1/n)
	vhl = ((1-pc)/pc)*(1/n)/vhmax 
	colfn = scale_fill_gradientn(colors = c("black","blue","white","red","black"),
						values = c(0,vhl/2,vhl,2*vhl,1),  limits = c(0,vhmax))

	vnl = pc*(1-pc)*n/vnmax
	colfn2 = scale_fill_gradientn(colors = c("black","blue","white","red","black"),
							values = c(0,vnl/2,vnl,2*vnl,1), limits = c(0,vnmax))
	
	
	g3.1 = ggplot(plot.df, aes(x=p, y=q, fill=plot.simhet)) + geom_tile() + colfn + facet_wrap( ~ seeds)+
		theme(axis.title = element_text(size = rel(1.25)),
		strip.text = element_text(size = rel(1.25)))
	g3.2 = ggplot(plot.df, aes(x=p, y=q, fill=plot.predhet1)) + geom_tile() + colfn + facet_wrap( ~ seeds)+
		theme(axis.title = element_text(size = rel(1.25)),
		strip.text = element_text(size = rel(1.25)))
	g3.3 = ggplot(plot.df, aes(x=p, y=q, fill=plot.predhet2)) + geom_tile() + colfn + facet_wrap( ~ seeds)+
		theme(axis.title = element_text(size = rel(1.25)),
		strip.text = element_text(size = rel(1.25)))
}

filename = paste(c(label, "-vh-stats-", ifelse(repel, "repel.png", "no-repel.png")), collapse="")
res.factor = 3
png(filename, width=1000*res.factor, height=1000*res.factor, res=72*res.factor)
if(repel){
	grid.arrange(g3.1, nrow = 1)
	#grid.arrange(g3.1, g3.2, nrow=2)	
} else {
	grid.arrange(g3.1, g3.2, g3.3, nrow=3)
}

dev.off()

vnmax = max(c(max(plot.df$plot.simn),max(plot.df$plot.predn)))
vhmax = max(c(max(plot.df$plot.simhet),max(plot.df$plot.predhet1),max(plot.df$plot.predhet2)))

pc = plot.prop

# Color gradient comparing to asymmetric symmetric without network, i.e.,  (1/n)
# This is a really ugly setup, which causes rather uniform redness across heatmaps
#vhl = (1/n)/vhmax

# Color gradient to compare asymmetric to asymmetric without network, i.e., ((1-pc)/pc)*(1/n)
vhl = ((1-pc)/pc)*(1/n)/vhmax 
colfn = scale_fill_gradientn(colors = c("black","blue","white","red","black"),
					values = c(0,vhl/2,vhl,2*vhl,1),  limits = c(0,vhmax))

vnl = pc*(1-pc)*n/vnmax
colfn2 = scale_fill_gradientn(colors = c("black","blue","white","red","black"),
						values = c(0,vnl/2,vnl,2*vnl,1), limits = c(0,vnmax))

g4.1 = ggplot(plot.df, aes(x=p, y=q, fill=plot.simn)) + geom_tile() + colfn2 + facet_wrap(~ seeds)+
		theme(axis.title = element_text(size = rel(1.25)),
		strip.text = element_text(size = rel(1.25)))
g4.2 = ggplot(plot.df, aes(x=p, y=q, fill=plot.predn)) + geom_tile() + colfn2 + facet_wrap(~ seeds)+
		theme(axis.title = element_text(size = rel(1.25)),
		strip.text = element_text(size = rel(1.25)))

res.factor = 3
filename = paste(c(label, "-vn-stats-", ifelse(repel, "repel.png", "no-repel.png")), collapse="")
png(filename, width=1000*res.factor, height=1000*res.factor, res=72*res.factor)
grid.arrange(g4.1, g4.2, nrow=2)
dev.off()

plot.df = res.1.df

vnmax = max(c(max(plot.df$plot.simn),max(plot.df$plot.predn)))
vhmax = max(c(max(plot.df$plot.simhet),max(plot.df$plot.predhet1)))

pc = plot.prop

# Color gradient comparing to asymmetric symmetric without network, i.e.,  (1/n)
#vhl = (1/n)/vhmax

# Color gradient to compare asymmetric to asymmetric without network, i.e., ((1-pc)/pc)*(1/n)
vhl = ((1-pc)/pc)*(1/n)/vhmax 
colfn = scale_fill_gradientn(colors = c("black","blue","white","red","black"),
					values = c(0,vhl/2,vhl,2*vhl,1),  limits = c(0,vhmax))

vnl = pc*(1-pc)*n/vnmax
colfn2 = scale_fill_gradientn(colors = c("black","blue","white","red","black"),
						values = c(0,vnl/2,vnl,2*vnl,1), limits = c(0,vnmax))

if(repel){
	# Plot simulation results (5.1, 5.2) and predictions (5.3,5.4) separately for repulsive case?
	 g5.1 = ggplot(plot.df, aes(x=p, y=q, fill=plot.simhet)) + geom_tile() + colfn + facet_grid(h ~ seeds)+
						 theme(axis.title = element_text(size = rel(1.25)),
									 strip.text = element_text(size = rel(1.25)))
	 g5.2 = ggplot(plot.df, aes(x=p, y=q, fill=plot.simn)) + geom_tile() + colfn2 + facet_grid(h ~ seeds)+
						 theme(axis.title = element_text(size = rel(1.25)),
									 strip.text = element_text(size = rel(1.25)))
	 g5.3 = ggplot(plot.df, aes(x=p, y=q, fill=plot.predhet1)) + geom_tile() + colfn + facet_grid(h ~ seeds)+
	 					 theme(axis.title = element_text(size = rel(1.25)),
	 								 strip.text = element_text(size = rel(1.25)))
	 g5.4 = ggplot(plot.df, aes(x=p, y=q, fill=plot.predn)) + geom_tile() + colfn2 + facet_grid(h ~ seeds)+
	 					 theme(axis.title = element_text(size = rel(1.25)),
	 								 strip.text = element_text(size = rel(1.25)))

	 res.factor = 3
	 filename = paste(c(label, "-both-stats-", "repel.png"), collapse="")
	 png(filename, width=1000*res.factor, height=800*res.factor, res=72*res.factor)
	 grid.arrange(g5.1, g5.2, nrow=1)
	 dev.off()
	 filename = paste(c(label, "-both-stats-taylor-", "repel.png"), collapse="")
	 png(filename, width=1000*res.factor, height=800*res.factor, res=72*res.factor)
	 grid.arrange(g5.3, g5.4, nrow=1)
	 dev.off()
}else{
	 g5.1 = ggplot(plot.df, aes(x=p, y=q, fill=plot.simhet)) + geom_tile() + colfn + facet_grid(h ~ seeds)+
		 theme(axis.title = element_text(size = rel(1.25)),
		 strip.text = element_text(size = rel(1.25)))
	 g5.2 = ggplot(plot.df, aes(x=p, y=q, fill=plot.simn)) + geom_tile() + colfn2 + facet_grid(h ~ seeds)+
		 theme(axis.title = element_text(size = rel(1.25)),
		 strip.text = element_text(size = rel(1.25)))
	 g5.3 = ggplot(plot.df, aes(x=p, y=q, fill=plot.predhet1)) + geom_tile() + colfn + facet_grid(h ~ seeds)+
		 theme(axis.title = element_text(size = rel(1.25)),
			strip.text = element_text(size = rel(1.25)))
	 g5.4 = ggplot(plot.df, aes(x=p, y=q, fill=plot.predn)) + geom_tile() + colfn2 + facet_grid(h ~ seeds)+
		 theme(axis.title = element_text(size = rel(1.25)),
		 strip.text = element_text(size = rel(1.25)))

	 res.factor = 3
	 filename = paste(c(label, "-both-stats-","no-repel.png"), collapse="")
	 png(filename, width=1000*res.factor, height=800*res.factor, res=72*res.factor)
	 grid.arrange(g5.1, g5.2, g5.3, g5.4, nrow=2)
	 dev.off()
}

############### moments

cat("Plotting moments...\n")

# get list of statistics to plot
stats = colnames(test.df[,8:ncol(test.df)])
if(repel == F) {
	stats.to.plot = stats[stats!="md"&stats!="vd"&stats!="dc"]
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
	pred = vhpred(alpha, beta, test.df$p[i], test.df$q[i], test.df$h[i], n, test.df$prop[i])

	# populate list of statistics that we're interested in
	for(j in 1:length(stats)) {
		if(stats[j] %in% stats.to.plot){
			mom.df = rbind(mom.df, data.frame(prop=test.df$prop[i],seeds=test.df$seeds[i], p=test.df$p[i], q=test.df$q[i], h=test.df$h[i], stat = stats[j], predval = pred[3+j], simval = test.df[i,7+j]))
		}
	}
}

filename = paste(c(label, "-moments-", ifelse(repel, "repel.png", "no-repel.png")), collapse="")
png(filename, width=1000*res.factor, height=1000*res.factor, res=72*res.factor)
ggplot(mom.df, aes(x=predval,y=simval,colour=factor(prop))) + geom_point() + geom_abline(slope=1, intercept=0) + facet_wrap(~stat, scales="free")+
	theme(axis.title = element_text(size = rel(1.25)),
	strip.text = element_text(size = rel(1.25)),
	legend.text = element_text(size = rel(1.25)))
dev.off()

############### compare statistical simulation output to theory

# initialise a dataframe storing theory-simulation comparison for each parameterisation
cat("Producing statistical simulation...\n")

res.2.df = data.frame(alpha=NULL,beta=NULL,n=NULL,p=NULL,q=NULL,h=NULL,vhsim=NULL,vhpred1=NULL,vhpred2=NULL)

# number of samples for each parameterisation
nsamp = 500

# loop through various parameters
for(alpha in Alpha) {
 for(beta in Beta) {
	for(n in c(20, 50, 100)) {
		for(p in (0:5)/5) {
			for(q in (0:5)/5) {
				for(h in (0:5)/5) {
					#ystar = 0
					pc = prop.to.plot
					#beta = alpha # symmetric, and pc = 0
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
					pred = vhpred(alpha, beta, p, q, h, n, prop.to.plot)
					res.2.df = rbind(res.2.df, data.frame(alpha=alpha, beta = beta, n=n, p=p, q=q, h=h, vhsim = var(het), vhpred1 = pred[1], vhpred2 = pred[2], vnsim = var(ndna), vnpred = pred[3]))
				}
			}
		}
	}
 }
}
cat("Plotting comparisons...\n")

# plot simulation-theory comparison for first-order (1) and second-order (2) models
if(repel){
	g2.1 = ggplot(res.2.df, aes(x=vhsim, y=vhpred1, color=factor(h))) + geom_point() + geom_abline(slope=1, intercept=0) + facet_grid(alpha ~ n) +
	ggtitle("Stat sim vs first-order theory") +
	theme(axis.title = element_text(size = rel(1.25)),
				strip.text = element_text(size = rel(1.25)),
				legend.text = element_text(size = rel(1.25)))
	#g2.2 = ggplot(res.2.df, aes(x=vhsim, y=vhpred2, color=factor(h))) + geom_point() + geom_abline(slope=1, intercept=0) + facet_grid(alpha ~ n) + ggtitle("Stat sim vs second-order theory")
	g2.3 = ggplot(res.2.df, aes(x=vnsim, y=vnpred, color=factor(h))) + geom_point() + geom_abline(slope=1, intercept=0) + facet_grid(alpha ~ n) + 
	ggtitle("Stat sim vs theory")+
	theme(axis.title = element_text(size = rel(1.25)), 
				strip.text = element_text(size = rel(1.25)),
				legend.text = element_text(size = rel(1.25)))
} else {
	g2.1 = ggplot(res.2.df, aes(x=vhsim, y=vhpred1, color=factor(h))) + geom_point() + geom_abline(slope=1, intercept=0) + facet_grid(alpha ~ n) +
	ggtitle("Stat sim vs first-order theory")+
	theme(axis.title = element_text(size = rel(1.25)),
				strip.text = element_text(size = rel(1.25)),
				legend.text = element_text(size = rel(1.25)))
	g2.2 = ggplot(res.2.df, aes(x=vhsim, y=vhpred2, color=factor(h))) + geom_point() + geom_abline(slope=1, intercept=0) + facet_grid(alpha ~ n) +
	ggtitle("Stat sim vs second-order theory")+
	theme(axis.title = element_text(size = rel(1.25)),
				strip.text = element_text(size = rel(1.25)),
				legend.text = element_text(size = rel(1.25)))
	g2.3 = ggplot(res.2.df, aes(x=vnsim, y=vnpred, color=factor(h))) + geom_point() + geom_abline(slope=1, intercept=0) + facet_grid(alpha ~ n) +
	ggtitle("Stat sim vs theory")+
	theme(axis.title = element_text(size = rel(1.25)),
				strip.text = element_text(size = rel(1.25)),
				legend.text = element_text(size = rel(1.25)))
}
res.factor = 3
filename = paste(c(label, "-sim-theory-comparison-", ifelse(repel, "repel.png", "no-repel.png")), collapse="")
png(filename, width=1600*res.factor, height=800*res.factor, res=72*res.factor)
if(repel){
	grid.arrange(g1.1, g1.3, g2.1, g2.3, nrow=2)
} else {
	grid.arrange(g1.1, g1.2, g1.3, g2.1, g2.2, g2.3, nrow=2)
}
dev.off()

if(!repel){
	cat("Estimating and plotting V(h) = E(h^2)-E(h)^2 by sum over state variables ...\n")

	plot.df.vhest = data.frame(seeds = NULL, alpha = NULL, beta = NULL, n = NULL, p = NULL, q = NULL, h = NULL, vh = NULL, vn = NULL)

	pc = prop.to.plot
	nseeds = c(4,16,64)
	cat(paste0("pc = ", paste0(toString(round(pc,2)), "\n")))
	for(h in c(0.1,0.5)){
		for(k in 1:length(Alpha)){
			for(i in 0:10){
				for(j in 0:10){
				p = i*0.1
				q = j*0.1

				alpha = Alpha[k]
				beta = Beta[k]
				seeds = nseeds[k]
				ret = vhest(p,q,h,n,pc,alpha,beta)
				alph = ret[1]
				bet = ret[2]
				vh = ret[3]
				vn = ret[4]

				plot.df.vhest = rbind(plot.df.vhest,data.frame(seeds = seeds,alpha=alph,beta=bet,n=n,p=p,q=q,h=h,vh=vh,vn=vn))
				}   
			}
		}
	}

	vnmax = max(max(plot.df$plot.simn),max(plot.df.vhest$vn))
	vhmax = max(max(plot.df$plot.simhet),max(plot.df.vhest$vh))

	# Color gradient comparing to asymmetric symmetric without network, i.e.,  (1/n)
	#vhl = (1/n)/vhmax

	# Color gradient to compare asymmetric to asymmetric without network, i.e., ((1-pc)/pc)*(1/n)
	vhl = ((1-pc)/pc)*(1/n)/vhmax
	colfn = scale_fill_gradientn(colors = c("black","blue","white","red","black"),
						values = c(0,vhl/2,vhl,2*vhl,1),  limits = c(0,vhmax))

	vnl = pc*(1-pc)*n/vnmax
	colfn2 = scale_fill_gradientn(colors = c("black","blue","white","red","black"),
							values = c(0,vnl/2,vnl,2*vnl,1), limits = c(0,vnmax))


	cat(paste0("vhmax = ", paste0(toString(round(vhmax,5)), "\n")))
	cat(paste0("vnmax = ", paste0(toString(round(vnmax)), "\n")))

	g6.1 = ggplot(plot.df, aes(x=p, y=q, fill=plot.simhet)) + geom_tile() + colfn + facet_grid(h ~ seeds) + labs(fill = "simhet   ")+
		theme(axis.title = element_text(size = rel(1.25)),
					strip.text = element_text(size = rel(1.25)))
	g6.2 = ggplot(plot.df, aes(x=p, y=q, fill=plot.simn)) + geom_tile() + colfn2 + facet_grid(h ~ seeds) + labs(fill = "simn     ")+
		theme(axis.title = element_text(size = rel(1.25)),
					strip.text = element_text(size = rel(1.25)))
	p6.3 = ggplot(plot.df.vhest, aes(x=p,y=q,fill=vh))+geom_tile()+colfn+facet_grid(h ~ seeds) + labs(fill = "esthet   ")+
		theme(axis.title = element_text(size = rel(1.25)),
					strip.text = element_text(size = rel(1.25)))
	p6.4 = ggplot(plot.df.vhest, aes(x=p,y=q,fill=vn))+geom_tile()+colfn2+facet_grid(h ~ seeds) + labs(fill = "estn     ")+
		theme(axis.title = element_text(size = rel(1.25)),
					strip.text = element_text(size = rel(1.25)))

	filename = paste(c(label,"-fig-sim-vhest.png"),collapse="")
	res.factor = 3
	png(filename, height = 1000*res.factor, width = 1000*res.factor, res = 72*res.factor)
	grid.arrange(g6.1, g6.2, p6.3, p6.4, ncol = 2)
	dev.off()

	plot.df.vhest.h = plot.df.vhest[plot.df.vhest$h == 0.5,]  
	plot.df.1 = res.1.df[res.1.df$h == 0.1,]
	plot.df.2 = res.1.df[res.1.df$h == 0.5,]

	vnmax = max(max(plot.df.1$plot.simn),max(plot.df.2$plot.simn),max(plot.df.vhest.h$vn),max(plot.df$plot.predn))
	vhmax = max(max(plot.df.1$plot.simhet),max(plot.df.2$plot.simhet),max(plot.df.vhest.h$vh),max(plot.df$plot.predhet1))

	# Color gradient comparing to asymmetric symmetric without network, i.e.,  (1/n)
	#vhl = (1/n)/vhmax

	# Color gradient to compare asymmetric to asymmetric without network, i.e., ((1-pc)/pc)*(1/n)
	vhl = ((1-pc)/pc)*(1/n)/vhmax
	colfn = scale_fill_gradientn(colors = c("black","blue","white","red","black"),
						values = c(0,vhl/2,vhl,2*vhl,1),  limits = c(0,vhmax))

	vnl = pc*(1-pc)*n/vnmax
	colfn2 = scale_fill_gradientn(colors = c("black","blue","white","red","black"),
							values = c(0,vnl/2,vnl,2*vnl,1), limits = c(0,vnmax))

	g7.1 = ggplot(plot.df.1, aes(x=p, y=q, fill=plot.simhet)) + geom_tile() + colfn + facet_grid(h ~ seeds) + labs(fill = "simhet  ", tag = "A") + 
		theme(plot.title = element_text(face = "bold", size = rel(1.25)),
					axis.title = element_text(size = rel(1.25)),
					strip.text = element_text(size = rel(1.25)),
					legend.text = element_text(size = rel(1.1))) 
	g7.2 = ggplot(plot.df.1, aes(x=p, y=q, fill=plot.simn)) + geom_tile() + colfn2 + facet_grid(h ~ seeds) + labs(fill = "simn    ", tag = "")+
		theme(axis.title = element_text(size = rel(1.25)),
					strip.text = element_text(size = rel(1.25)),
					legend.text = element_text(size = rel(1.1)))
	g7.3 = ggplot(plot.df.2, aes(x=p, y =q, fill=plot.simhet)) + geom_tile() + colfn + facet_grid(h ~ seeds) + labs(fill = "simhet  ", tag = "B") + 
		theme(plot.title = element_text(face = "bold", size = rel(1.25)),
					axis.title = element_text(size = rel(1.25)),
					strip.text = element_text(size = rel(1.25)),
					legend.text = element_text(size = rel(1.1)))
	g7.4 = ggplot(plot.df.2, aes(x=p, y =q, fill=plot.simn)) + geom_tile() + colfn2 + facet_grid(h ~ seeds) + labs(fill = "simn    ", tag = "")+
		theme(axis.title = element_text(size = rel(1.25)),
					strip.text = element_text(size = rel(1.25)),
					legend.text = element_text(size = rel(1.1)))
	g7.5 = ggplot(plot.df.vhest.h, aes(x=p,y=q,fill=vh))+geom_tile()+colfn+facet_grid(h ~ seeds) + labs(fill = "esthet  ", tag = "C") + 
		theme(plot.title = element_text(face = "bold", size = rel(1.25)),
					axis.title = element_text(size = rel(1.25)),
					strip.text = element_text(size = rel(1.25)),
					legend.text = element_text(size = rel(1.1)))
	g7.6 = ggplot(plot.df.vhest.h, aes(x=p,y=q,fill=vn))+geom_tile()+colfn2+facet_grid(h ~ seeds) + labs(fill = "estn    ", tag = "")+
		theme(axis.title = element_text(size = rel(1.25)),
					strip.text = element_text(size = rel(1.25)),
					legend.text = element_text(size = rel(1.1)))
	g7.7 = ggplot(plot.df.2, aes(x=p, y=q, fill=plot.predhet1)) + geom_tile() + colfn + facet_grid(h ~ seeds) + labs(fill = "predhet1", tag = "D") + 
		theme(plot.title = element_text(face = "bold", size = rel(1.25)),
					axis.title = element_text(size = rel(1.25)),
					strip.text = element_text(size = rel(1.25)),
					legend.text = element_text(size = rel(1.1)))
	g7.8 = ggplot(plot.df.2, aes(x=p, y=q, fill=plot.predn)) + geom_tile() + colfn2 + facet_grid(h ~ seeds) + labs(fill = "predn   ", tag = "")+
		theme(axis.title = element_text(size = rel(1.25)),
					strip.text = element_text(size = rel(1.25)),
					legend.text = element_text(size = rel(1.1)))

	
	filename = paste(c(label,"-fig-sim-est-taylor.png"), collapse="")
	res.factor = 3
	png(filename, height = 1000*res.factor, width = 1000*res.factor, res = 72*res.factor)
	grid.arrange(g7.1, g7.2, g7.3, g7.4, g7.5, g7.6, g7.7, g7.8, nrow = 4)
	dev.off()
}

