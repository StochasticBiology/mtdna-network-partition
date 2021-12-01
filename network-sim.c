// code to simulate partitioning of mtDNAs given some simulated network structure

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RND drand48()
#define PI 3.14159

double BRANCH = 0.02; // branch probability per timestep
double DELTA = 0.01;  // growth rate per timestep
int MAXN = 10000;     // memory limit for number of network segments
int MAXM = 100;       // number of mtDNAs
double MASS = 50;     // required network mass
int NSIM = 100;      // number of simulations for each parameterisation
int NSTEPS = 100;     // number of gaussian steps to take after fragmentation (if this scheme is chosen)

// structure to store summary statistics from a set of simulations
typedef struct {
double mw, mm, mn, md, mu, mh;
double mw2, mm2, mw3, mm3;
double vw, vm, vn, vd, vu, vh;
double cwm, cw2m, cw3m, cwm2, cwm3, cw2m2;
double muw3, muw4, mum3, mum4;
} SumStats;

// structure to store state of a particular simulation
typedef struct {
double wn, wc, mn, mc;
double d, u;
} Stats;


// GSL routine for gaussian random number with sd sigma
double gsl_ran_gaussian(const double sigma)
{
  double x, y, r2;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */

      x = -1 + 2 * RND;
      y = -1 + 2 * RND;

      /* see if it is in the unit circle */
      r2 = x * x + y * y;
    }
  while (r2 > 1.0 || r2 == 0);

  /* Box-Muller transform */
  return sigma * y * sqrt (-2.0 * log (r2) / r2);
}

// simple function returning max(a,b)
double mymax(double a, double b)
{
  if(a > b) return a;
  else return b;
}

// simple function returning min(a,b)
double mymin(double a, double b)
{
  if(a < b) return a;
  else return b;
}

// this set of functions is taken from https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/ and computes whether two line segments intersect
// this one computes whether q lies on segment pr
int onSegment(double px, double py, double qx, double qy, double rx, double ry)
{
  if(qx <= mymax(px, rx) && qx >= mymin(px, rx) && qy <= mymax(py, ry) && qy >= mymin(py, ry)) return 1;
  else return 0;
}

// this one computes the orientation of a set of three points (or whether they are colinear)
int orientation(double px, double py, double qx, double qy, double rx, double ry)
{
  int val = (qy-py)*(rx-qx) - (qx-px)*(ry-qy);
  if(val == 0) return 0;
  return (val > 0 ? 1 : 2);
}

// this one returns 1/0 intersection given two line segments described by start and end points 
int doIntersect(double p1x, double p1y, double q1x, double q1y, double p2x, double p2y, double q2x, double q2y)
{
  int o1 = orientation(p1x, p1y, q1x, q1y, p2x, p2y);
  int o2 = orientation(p1x, p1y, q1x, q1y, q2x, q2y);
  int o3 = orientation(p2x, p2y, q2x, q2y, p1x, p1y);
  int o4 = orientation(p2x, p2y, q2x, q2y, q1x, q1y);
  
  // General case
  if (o1 != o2 && o3 != o4)
    return 1;
  
  // Special Cases
  // p1, q1 and p2 are colinear and p2 lies on segment p1q1
  if (o1 == 0 && onSegment(p1x, p1y, p2x, p2y, q1x, p1y)) return 1;
  
  // p1, q1 and q2 are colinear and q2 lies on segment p1q1
  if (o2 == 0 && onSegment(p1x, p1y, q2x, q2y, q1x, q1y)) return 1;
  
  // p2, q2 and p1 are colinear and p1 lies on segment p2q2
  if (o3 == 0 && onSegment(p2x, p2y, p1x, p1y, q2x, q2y)) return 1;
  
  // p2, q2 and q1 are colinear and q1 lies on segment p2q2
  if (o4 == 0 && onSegment(p2x, p2y, q1x, q1y, q2x, q2y)) return 1;
  
  return 0;
}

void BuildNetwork(int nseed, double *xs, double *ys, double *xe, double *ye, int *finaln)
{
  double m;
  int n;
  int nactive;
  int notdoneyet;
  int i;
  double r;
  double newx, newy;
  double *thetas;
  int *active;
  double starttheta;
  
  thetas = (double*)malloc(sizeof(double)*MAXN);
  active = (int*)malloc(sizeof(int)*MAXN);

  // m -- mass of mito network
  // n -- number of strands
  // nactive -- number of strands currently growing
  m = 0;
  n = 0;
  nactive = 0;

  // loop until we get the network mass we want
  for(; m < MASS; )
    {
      // if we're out of growing strands (which we are at the start), start strands growing N, S, E, W
      if(nactive == 0)
	{
	  starttheta = RND*2.*3.14159;
	  for(i = 0; i < nseed; i++)
	    {
	      xs[n] = xe[n] = 1.*cos(2.*3.14159*i/nseed + starttheta);
	      ys[n] = ye[n] = 1.*sin(2.*3.14159*i/nseed + starttheta);
	      // 
	      thetas[n] = 3.14159 + 2.*3.14159*i/nseed + starttheta;
	      active[n] = 1; n++; nactive++;
	    }
	}

      // loop through strands in network
      for(i = 0; i < n; i++)
	{
	  // process ones which are currently active
	  if(active[i])
	    {
	      // decide whether to branch
	      r = RND;
	      if(r < BRANCH)
		{
		  if(n > MAXN-2)
		    {
		      // out of memory
		      printf("Not at required mass, but too many segments\n");
		      exit(0);
		    }
		  // switch off this strand 
		  active[i] = 0;

		  // produce first daughter strand -- start at the end of this one, random direction, active
		  xs[n] = xe[n] = xe[i]; ys[n] = ye[n] = ye[i];
		  thetas[n] = RND*2*PI;
		  active[n] = 1;
		  n++;

		  // produce second daughter strand -- start at the end of this one, random direction, active
		  xs[n] = xe[n] = xe[i]; ys[n] = ye[n] = ye[i];
		  thetas[n] = RND*2*PI;
		  active[n] = 1;
		  n++;

		  // increment active strands by 1 (lost 1, gained 2) and avoid processing the daughters this timestep
		  i += 2;
		  nactive++;
		}
	      else
		{
		  // if not branching, just grow this strand in given direction by amount DELTA
		  newx = xe[i] + cos(thetas[i])*DELTA; 
		  newy = ye[i] + sin(thetas[i])*DELTA;
		  if(newx*newx + newy*newy > 1)
		    {
		      // switch off if we hit the edge of cell
		      active[i] = 0;
		      nactive--;
		    }
		  else
		    {
		      // otherwise update the end of this strand and total network mass
		      if(m < MASS)
			{
		          xe[i] = newx;
  		          ye[i] = newy;
		          m += DELTA;
			}
		    }
		}
	    }
	}
    }

  *finaln = n;
      
  free(thetas);
  free(active);
}

int PlaceDNA(double het, double inc0, double inc1, double halo, double *xs, double *ys, double *xe, double *ye, int n, double *mx, double *my, int *mt, int *networked)
{
  double thism;
  int i, j, k;
  int placed;
  int counter;
  int notdoneyet;
  double r;
  int failed;
  double newx, newy;
  double *cumsum;
  double ball;
  
  cumsum = (double*)malloc(sizeof(double)*MAXN);

  notdoneyet = 0;
  // populate roulette wheel of network strands (by mass)
  // begin with strand 0
  thism = sqrt((xe[0]-xs[0])*(xe[0]-xs[0]) + (ye[0]-ys[0])*(ye[0]-ys[0]));
  cumsum[0] = thism;
  // build up cumulative sum of strand masses
  for(i = 1; i < n; i++)
    {
      thism = sqrt((xe[i]-xs[i])*(xe[i]-xs[i]) + (ye[i]-ys[i])*(ye[i]-ys[i]));
      cumsum[i] = cumsum[i-1] + thism;
    }

  // allocate mtDNA types (1 == mutant, 0 == wildtype)
  for(i = 0; i < MAXM; i++) mt[i] = 0;
  for(i = 0; i < het*MAXM; i++)
    {
      do{ k = RND*MAXM; }while(mt[k] == 1);
      mt[k] = 1;
    }
      
  // distribute mtDNAs
  for(i = 0; i < MAXM; i++)
    {
      // track whether this mtDNA has been placed or not yet
      placed = 0;
      // decide whether to put this one in the network or not
      if((mt[i] == 1 && RND < inc1) || (mt[i] == 0 && RND < inc0))
	{
	  networked[i] = 1;
	  counter = 0;
	  while(placed == 0 && notdoneyet == 0)
	    {
	      // roulette wheel to choose a strand to try placing this mtDNA in
	      ball = RND*cumsum[n-1];
	      for(j = 0; cumsum[j] < ball; j++);
	      // random choice of how long along the strand to try
	      r = RND;
	      // store proposed position
	      newx = xs[j] + r*(xe[j]-xs[j]);
	      newy = ys[j] + r*(ye[j]-ys[j]);

	      // now see if that's compatible with our rules
	      failed = 0;
	      // loop through already-placed (networked) mtDNAs
	      for(k = 0; k < i; k++)
		{
		  if(networked[k] == 1)
		    {
		      // check to see if we're under a threshold distance; fail if so
		      if((mx[k]-newx)*(mx[k]-newx) + (my[k]-newy)*(my[k]-newy) < halo*halo)
			{
			  failed = 1;
			  break;
			}
		    }
		}
	      // if not failed, placed mtDNA here; otherwise keep looping
	      if(failed == 0)
		{
		  mx[i] = newx;
		  my[i] = newy;
		  placed = 1;
		}
	      // if we've tried and failed to place this mtDNA too many times, quit (threshold too high)
	      counter++;
	      if(counter > 10000)
		{
		  printf("Struggling to place mito %i\n", i);
		  notdoneyet = 1;
		  break;
		}
	    }
	}
      else
	{
	  // we're not placing this mtDNA in the network, just the cytoplasm
	  networked[i] = 0;
	  while(placed == 0)
	    {
	      newx = 2*(RND-0.5);
	      newy = 2*(RND-0.5);
	      // loop round until we pick a point in the cell
	      if(newx*newx + newy*newy < 1)
		{
		  mx[i] = newx;
		  my[i] = newy;
		  placed = 1;
		}
	    }
	}
    }

  free(cumsum);
  if(notdoneyet) return -1;
  return 0;
}

void PerturbDNA(double lambda, double *mx, double *my, int perturbtype)
{
  int i, j;
  int failed;
  double newx, newy;

  // first perturbation type -- move some proportion of the way towards a randomly-chosen point
  if(perturbtype == 1)
    {
      // apply post-placement perturbation to each mtDNA
      for(i = 0; i < MAXM; i++)
	{
	  // pick a random point in the cell
	  failed = 1;
	  while(failed == 1)
	    {
	      newx = 2*(RND-0.5);
	      newy = 2*(RND-0.5);
	      if(newx*newx+newy*newy < 1) failed = 0;
	    }
	  // move our mtDNA lambda of the way towards that point
	  mx[i] = mx[i] + lambda*(newx-mx[i]);
	  my[i] = my[i] + lambda*(newy-my[i]);	    
	}
    }
  else // second perturbation type -- take NSTEPS gaussian steps with kernel width lambda
    {
      for(i = 0; i < MAXM; i++)
	{
	  for(j = 0; j < NSTEPS; j++)
	    {
	      newx = mx[i] + gsl_ran_gaussian(lambda);
	      newy = my[i] + gsl_ran_gaussian(lambda);
	      if(newx*newx + newy*newy < 1)
		{
		  mx[i] = newx;
		  my[i] = newy;
		}
	    }
	}
    }
	      
}

void Output(double *xs, double *ys, double *xe, double *ye, int n, double *mx, double *my, int *mt, double h, double nseed, double inc0, double inc1, double lambda, double halo, int perturbtype)
{
  char str[200];
  FILE *fp;
  int i;
  
  // output network structure
  sprintf(str, "network-%.1f-%.0f-%.2f-%.2f-%.2f-%.2f-%i.csv", h, nseed, inc0, inc1, lambda, halo, perturbtype);
  fp = fopen(str, "w");
  fprintf(fp, "xs,ys,xe,ye\n");
  for(i = 0; i < n; i++)
    fprintf(fp, "%f,%f,%f,%f\n", xs[i], ys[i], xe[i], ye[i]);
  fclose(fp);

  // output mtDNA positions
  sprintf(str, "mtdna-%.1f-%.0f-%.2f-%.2f-%.2f-%.2f-%i.csv", h, nseed, inc0, inc1, lambda, halo, perturbtype);
  fp = fopen(str, "w");
  fprintf(fp, "x,y,type\n");
  for(i = 0; i < MAXM; i++)
    fprintf(fp, "%f,%f,%i\n", mx[i], my[i], mt[i]);
  fclose(fp);
}

void GetStats(double *xs, double *ys, double *xe, double *ye, int n, double *mx, double *my, int *mt, int *networked, int *wn, int *wc, int *mn, int *mc, double *d, double *u, double ystar)
{
   int i, j;
  int counter;
  double mindist;
  double thisdist;
  double thislen;
  double propupper;
  double total;
  
  // compute statistics of one daughter cell
  // hstat -- this simulation's heteroplasmy for y > y*
  // nstat -- this simulation's mtDNA count for y > y*
  *wc = *wn = *mc = *mn = 0;
  for(i = 0; i < MAXM; i++)
    {
      if(my[i] > ystar)
	{
	  if(mt[i] == 0 && networked[i] == 0) (*wc)++;
	  if(mt[i] == 0 && networked[i] == 1) (*wn)++;
	  if(mt[i] == 1 && networked[i] == 0) (*mc)++;
	  if(mt[i] == 1 && networked[i] == 1) (*mn)++;
	}
    }

  // compute mean minimum distance between mtDNAs
  // to do -- subset this into e.g. networked, same genetics, etc
  *d = 0; counter = 0;
  // loop through mtDNAs
  for(i = 0; i < MAXM-1; i++)
    {
      mindist = -1;
      // loop through partners, compute dist and compare to minimum
      for(j = i+1; j < MAXM; j++)
	{
	  thisdist = (mx[i]-mx[j])*(mx[i]-mx[j]) + (my[i]-my[j])*(my[i]-my[j]);
	  if(thisdist < mindist || mindist == -1)
	    mindist = thisdist;
	  counter++;
	}
      *d += sqrt(mindist);
    }
  *d /= counter;

  // compute amount of network mass in upper cell, looping through network fragments
  *u = 0; total = 0;
  for(i = 0; i < n; i++)
    {
      thislen = (xe[i]-xs[i])*(xe[i]-xs[i]) + (ye[i]-ys[i])*(ye[i]-ys[i]);
      thislen = sqrt(thislen);
      // if all or none is in the upper cell, record this, otherwise record proportion of mass that is
      if(ye[i] > ystar && ys[i] > ystar) propupper = 1;
      else if(ye[i] < ystar && ys[i] < ystar) propupper = 0;
      else if(ye[i] > ys[i]) propupper = (ye[i]-ystar)/(ye[i]-ys[i]);
      else propupper = (ys[i]-ystar)/(ys[i]-ye[i]);
      if(propupper > 1 || propupper < 0)
	printf("wtf\n");
      *u += thislen*propupper;
      total += thislen;
    }
}

// rather verbose function to compute summary statistics given outputs of NSIM simulations
void ComputeStats(Stats *s, SumStats *ss)
{
  double w, m, n;
  int i;

  // initialise all means and variances in the summary structure
  ss->mw = ss->mm = ss->vm = ss->vw = ss->mn = 0;
  ss->md = ss->vd = ss->mu = ss->vu = ss->mh = ss->vh = ss->vn = 0;
  ss->cwm = ss->cw2m = ss->cwm2 = ss->cw2m2 = ss->cw3m = ss->cwm3 = 0;
  ss->muw3 = ss->mum3 = ss->muw4 = ss->mum4 = 0;
  // compute means
  for(i = 0; i < NSIM; i++)
    {
      w = s[i].wn+s[i].wc;
      m = s[i].mn+s[i].mc;
      n = w+m;
      ss->mw += w;
      ss->mm += m;
      ss->mn += n;
      ss->mw2 += w*w;
      ss->mm2 += m*m;
      ss->mm3 += w*w*w;
      ss->mm3 += m*m*m;
      ss->mh += m/(w+m);
      ss->mu += s[i].u;
      ss->md += s[i].d;
    }
  ss->mw /= NSIM; ss->mm /= NSIM; ss->mn /= NSIM;
  ss->mw2 /= NSIM; ss->mw3 /= NSIM; ss->mm2 /= NSIM; ss->mm3 /= NSIM;
  ss->mh /= NSIM;
  ss->mu /= NSIM; ss->md /= NSIM;

  // compute variances, covariances, moments
  for(i = 0; i < NSIM; i++)
    {
      w = s[i].wn+s[i].wc;
      m = s[i].mn+s[i].mc;
      n = w+m;

      ss->vd += (s[i].d - ss->md)*(s[i].d - ss->md);
      ss->vu += (s[i].u - ss->mu)*(s[i].u - ss->mu);
      ss->vw += (w-ss->mw)*(w-ss->mw);
      ss->vm += (m-ss->mm)*(m-ss->mm);
      ss->vn += (n-ss->mn)*(n-ss->mn);
				
      ss->cwm += (w-ss->mw)*(m-ss->mm);
      ss->cw2m += (w*w-ss->mw2)*(m-ss->mm);
      ss->cw3m += (w*w*w-ss->mw3)*(m-ss->mm);
      ss->cwm2 += (w-ss->mw)*(m*m-ss->mm2);
      ss->cwm3 += (w-ss->mw)*(m*m*m-ss->mm3);
      ss->cw2m2 += (w*w-ss->mw2)*(m*m-ss->mm2);

      ss->muw3 += (w-ss->mw)*(w-ss->mw)*(w-ss->mw);
      ss->mum3 += (m-ss->mm)*(m-ss->mm)*(m-ss->mm);
      ss->muw4 += (w-ss->mw)*(w-ss->mw)*(w-ss->mw)*(w-ss->mw);
      ss->mum4 += (m-ss->mm)*(m-ss->mm)*(m-ss->mm)*(m-ss->mm);
      
      ss->vh += (m/(w+m)-ss->mh)*(m/(w+m)-ss->mh);
    }
  ss->vw /= NSIM-1; ss->vm /= NSIM-1; ss->vn /= NSIM-1;
  ss->vd /= NSIM-1; ss->vu /= NSIM-1;
  ss->cwm /= NSIM-1; ss->cw2m /= NSIM-1; ss->cw3m /= NSIM-1; ss->cwm2 /= NSIM-1; ss->cwm3 /= NSIM-1; ss->cw2m2 /= NSIM-1;
  ss->vh /= NSIM-1;
}

int main(int argc, char *argv[])
{
  double *xs, *ys, *xe, *ye, *thetas;
  double *cumsum;
  int *active;
  int n;
  double m;
  FILE *fp;
  int i, j, k;
  int counter;
  double r;
  double newx, newy;
  int nactive;
  int sim;
  char str[100];
  double *mx, *my;
  int *mt, *networked;
  double inc0, inc0s, inc0e, inc0a, inc1, inc1s, inc1e, inc1a;
  int placed, failed;
  double het, hets, hete, heta;
  double ball;
  int wn, wc, mn, mc;
  double d, u;
  int output;
  double lambda, lambdas, lambdae, lambdaa;
  FILE *fpout;
  int perturbtype;
  double thism;
  double halo, halos, haloe, haloa;
  double thisdist, mindist;
  int notdoneyet;
  double nseed, nseeds, nseede, nseedm;
  int expt;
  double total, thislen;
  Stats *stats;
  SumStats ss;
  double ystar;
  
  if(argc != 4)
    {
      printf("Need ystar, h, nseeds!\n");
      return 0;
    }
   ystar = atof(argv[1]);
  het = atof(argv[2]);
  nseed = atof(argv[3]);
 
  
  // these, plus the globals at the top, are the default values for the parameter sweeps we'll be doing. the switch below based on the argument provided at the command line modifies these for a particular setup
  perturbtype = 2;                            // how to perturb mtDNAs after fragmentation. 2 is the better version (repeated normal kernel)
  inc0s = 0; inc0e = 1; inc0a = 0.1;          // start, end, additive step for inclusion probability for type 0 (=p)
  inc1s = 0; inc1e = 1; inc1a = 0.1;          // start, end, additive step for inclusion probability for type 1 (=q)
  lambdas = 0; lambdae = 0.1; lambdaa = 0.02;    // start, end, additive step for post-fragmentation perturbation magnitude lambda
  halos = 0; haloe = 0.1; haloa = 0.1;        // start, end, additive step for repulsive "halo" around networked mtDNAs

  // allocate memory for mitochondrial strand and mtDNA properties
  // xs,ys: start point for a segment; xe,ye: end point; thetas: direction of that segment; active: whether that segment is actively growing
  // mx,my: position of an mtDNA; mt: genetic type of that mtDNA; networked: whether that mtDNA is in the network or not
  // stats: array of structures for statistics
  xs = (double*)malloc(sizeof(double)*MAXN);
  ys = (double*)malloc(sizeof(double)*MAXN);
  xe = (double*)malloc(sizeof(double)*MAXN);
  ye = (double*)malloc(sizeof(double)*MAXN);
  mx = (double*)malloc(sizeof(double)*MAXM);
  my = (double*)malloc(sizeof(double)*MAXM);
  mt = (int*)malloc(sizeof(int)*MAXM);
  networked = (int*)malloc(sizeof(int)*MAXM);
  stats = (Stats*)malloc(sizeof(Stats)*NSIM);
  
  // file for output
  sprintf(str, "output-%.3f-%.3f-%.0f.csv", ystar, het, nseed);
  fpout = fopen(str, "w");
  fprintf(fpout, "ystar,h,seeds,p,q,lambda,halo,");
  fprintf(fpout, "mw,vw,mm,vm,cwm,cw2m,cw3m,cwm2,cwm3,cw2m2,muw3,mum3,muw4,mum4,mh,vh,md,vd,mu,vu,mn,vn\n");

  // workhorse part -- a very nested loop scanning through the parameters, according to the protocol we fixed above
  // het -- initial heteroplasmy [from command line]
  // nseed -- number of network seed points [from command line]
  // inc0 -- probability of a wildtype mtDNA being included in network
  // inc1 -- probability of a mutant mtDNA being included in network
  // lambda -- magnitude of perturbation applied to positions
  // halo -- minimum distance between networked mtDNAs
  for(inc0 = inc0s; inc0 <= inc0e; inc0 += inc0a)
    {
      for(inc1 = inc1s; inc1 <= inc1e; inc1 += inc1a)
	{
	  printf("%.2f %f %f %f\n", het, nseed, inc0, inc1);

	  for(lambda = lambdas; lambda <= lambdae; lambda += lambdaa)
	    {
	      for(halo = halos; halo <= haloe; halo += haloa)
		{
		  // loop through NSIM simulations
		  for(sim = 0; sim < NSIM; sim++)
		    {

		      // build network and try to place mtDNAs, keep looping until we have successfully placed all
		      notdoneyet = -1;
		      while(notdoneyet)
			{
			  BuildNetwork(nseed, xs, ys, xe, ye, &n);
			  notdoneyet = PlaceDNA(het, inc0, inc1, halo, xs, ys, xe, ye, n, mx, my, mt, networked);
			}

		      // if we are doing post-fragmentation perturbation, do so
		      if(lambda)
			PerturbDNA(lambda, mx, my, perturbtype);
			 
		      // decide whether to output snapshots of the system or not
		      if((inc0 == 0 || (inc0 > 0.49 && inc0 < 0.51) || inc0 > 0.99) && (inc1 == 0 ||  (inc1 > 0.49 && inc1 < 0.51) || inc1 > 0.99) && (lambda == 0 || (lambda > 0.039 && lambda < 0.041) || lambda > 0.099) && sim == 0)
			Output(xs, ys, xe, ye, n, mx, my, mt, het, nseed, inc0, inc1, lambda, halo, perturbtype);

		      // get and output statistics
		      GetStats(xs, ys, xe, ye, n, mx, my, mt, networked, &wn, &wc, &mn, &mc, &d, &u, ystar);
		      stats[sim].wn = wn;
		      stats[sim].wc = wc;
		      stats[sim].mn = mn;
		      stats[sim].mc = mc;			      
		      stats[sim].u = u;
		      stats[sim].d = d;
		    }
		  // compute summary statistics for this set of simulations
		  ComputeStats(stats, &ss);

		  // output statistics for this set
		  fprintf(fpout, "%.3f,%.2f,%.0f,%f,%f,%f,%f,", ystar, het, nseed, inc0, inc1, lambda, halo);
		  fprintf(fpout, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", ss.mw, ss.vw, ss.mm, ss.vm, ss.cwm, ss.cw2m, ss.cw3m, ss.cwm2, ss.cwm3, ss.cw2m2, ss.muw3, ss.mum3, ss.muw4, ss.mum4, ss.mh, ss.vh, ss.md, ss.vd, ss.mu, ss.vu, ss.mn, ss.vn);
		}
	    }
	}
    }
  fclose(fpout);
  
  return 0;
}
