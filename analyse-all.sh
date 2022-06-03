#!/bin/bash
# analysis script
# to be run AFTER simulations have terminated (no error checking)

# concatenate output files from different simulation parameterisations into one big file
# all props (0.10,0.25,0.50), h0 = 0.5
# get header
head -n1 output-0.500-0.500-4.csv > output-moments.csv
# concatenate results, prune individual headers, and append to big file
cat output-0.100-0.500-4.csv output-0.100-0.500-16.csv output-0.100-0.500-64.csv output-0.250-0.500-4.csv output-0.250-0.500-16.csv output-0.250-0.500-64.csv output-0.500-0.500-4.csv output-0.500-0.500-16.csv output-0.500-0.500-64.csv | grep -v "h" >> output-moments.csv

# concatenate output files from different simulation parameterisations into one big file
# prop = 0.5, both h0 set (try only 1, to see if fit improves?)
# get header
head -n1 output-0.500-0.500-4.csv > output-asymm-0.50.csv
# concatenate results, prune individual headers, and append to big file
cat output-0.500-0.500-4.csv output-0.500-0.500-16.csv output-0.500-0.500-64.csv output-0.500-0.100-4.csv output-0.500-0.100-16.csv output-0.500-0.100-64.csv | grep -v "h" >> output-symm-0.50.csv

# get header
# prop = 0.25, both h0 set
head -n1 output-0.250-0.500-4.csv > output-asymm-0.25.csv
# concatenate results, prune individual headers, and append to big file
cat output-0.250-0.500-4.csv output-0.250-0.500-16.csv output-0.250-0.500-64.csv output-0.250-0.100-4.csv output-0.250-0.100-16.csv output-0.250-0.100-64.csv | grep -v "h" >> output-asymm-0.25.csv

# get header
# prop = 0.10, both h0 set
head -n1 output-0.100-0.500-4.csv > output-asymm-0.10.csv
# concatenate results, prune individual headers, and append to big file
cat output-0.100-0.500-4.csv output-0.100-0.500-16.csv output-0.100-0.500-64.csv output-0.100-0.100-4.csv output-0.100-0.100-16.csv output-0.100-0.100-64.csv | grep -v "h" >> output-asymm-0.10.csv

# symmetric division set
head -n1 output-0.500-0.500-4.csv > output-symm.csv
# concatenate results, prune individual headers, and append to big file
cat output-0.500-0.100-4.csv output-0.500-0.100-16.csv output-0.500-0.100-64.csv output-0.500-0.500-4.csv output-0.500-0.500-16.csv output-0.500-0.500-64.csv | grep -v "h" >> output-symm.csv

# plot cell snapshots -- no arguments because a long list of different input files is implemented therein
Rscript plotcell.R

# all props (10%,25%,50%)
# moment plots
Rscript compare-both-asymm.R output-moments.csv 0 100 50 0.5 0.500 all-props
Rscript compare-both-asymm.R output-moments.csv 1 100 50 0.5 0.500 all-props

# prop = 50%
# sim-theory comparison for non-repulsive case
Rscript compare-both-asymm.R output-symm-0.50.csv 0 100 50 0.1 0.500 symm-0.1
#Rscript compare-both-asymm.R output-symm-0.50.csv 0 100 50 0.5 0.500 symm-0.5
# sim-theory comparison for repulsive case
Rscript compare-both-asymm.R output-symm-0.50.csv 1 100 50 0.1 0.500 symm-0.1
#Rscript compare-both-asymm.R output-symm-0.50.csv 1 100 50 0.5 0.500 symm-0.5

# prop = 25%
# sim-theory comparison for non-repulsive case
Rscript compare-both-asymm.R output-asymm-0.25.csv 0 100 50 0.1 0.250 asymm-0.1-25
#Rscript compare-both-asymm.R output-asymm-0.25.csv 0 100 50 0.5 0.250 asymm-0.5-25
# sim-theory comparison for repulsive case
Rscript compare-both-asymm.R output-asymm-0.25.csv 1 100 50 0.1 0.250 asymm-0.1-25
#Rscript compare-both-asymm.R output-asymm-0.25.csv 1 100 50 0.5 0.250 asymm-0.5-25

# prop = 10%
# sim-theory comparison for non-repulsive case
Rscript compare-both-asymm.R output-asymm-0.10.csv 0 100 50 0.1 0.100 asymm-0.1-10
#Rscript compare-both-asymm.R output-asymm-0.10.csv 0 100 50 0.5 0.100 asymm-0.5-10
# sim-theory comparison for repulsive case
Rscript compare-both-asymm.R output-asymm-0.10.csv 1 100 50 0.1 0.100 asymm-0.1-10
#Rscript compare-both-asymm.R output-asymm-0.10.csv 1 100 50 0.5 0.100 asymm-0.5-10

# sim-theory comparison for symmetric cases:
# non-repulsive case
Rscript compare-both-asymm.R output-symm.csv 0 100 50 0.5 0.500 symm-50
# repulsive case
Rscript compare-both-asymm.R output-symm.csv 1 100 50 0.5 0.500 symm-50

# effect of lambda
Rscript plotlambda.R output-symm.csv

# compile manuscript
pdflatex writeup-rcg
bibtex writeup-rcg
pdflatex writeup-rcg
pdflatex writeup-rcg
