#!/bin/bash
# analysis script
# to be run AFTER simulations have terminated (no error checking)

# concatenate output files from different simulation parameterisations into one big file
# pc = 0.5, h0 = 0.5
# get header
head -n1 output-0.500-0.500-4.csv > output-moments.csv
# concatenate results, prune individual headers, and append to big file
cat output-0.100-0.500-4.csv output-0.100-0.500-16.csv output-0.100-0.500-64.csv output-0.250-0.500-4.csv output-0.250-0.500-16.csv output-0.250-0.500-64.csv output-0.500-0.500-4.csv output-0.500-0.500-16.csv output-0.500-0.500-64.csv | grep -v "h" >> output-moments.csv

# concatenate output files from different simulation parameterisations into one big file
# pc = 0.5, both h0 set
# get header
head -n1 output-0.500-0.500-4.csv > output-asymm-0.50.csv
# concatenate results, prune individual headers, and append to big file
cat output-0.500-0.500-4.csv output-0.500-0.500-16.csv output-0.500-0.500-64.csv output-0.500-0.100-4.csv output-0.500-0.100-16.csv output-0.500-0.100-64.csv | grep -v "h" >> output-asymm-0.50.csv

# get header
# pc = 0.25, both h0 set
head -n1 output-0.250-0.500-4.csv > output-asymm-0.25.csv
# concatenate results, prune individual headers, and append to big file
cat output-0.250-0.500-4.csv output-0.250-0.500-16.csv output-0.250-0.500-64.csv output-0.250-0.100-4.csv output-0.250-0.100-16.csv output-0.250-0.100-64.csv | grep -v "h" >> output-asymm-0.25.csv

# get header
# pc = 0.10, both h0 set
head -n1 output-0.100-0.500-4.csv > output-asymm-0.10.csv
# concatenate results, prune individual headers, and append to big file
cat output-0.100-0.500-4.csv output-0.100-0.500-16.csv output-0.100-0.500-64.csv output-0.100-0.100-4.csv output-0.100-0.100-16.csv output-0.100-0.100-64.csv | grep -v "h" >> output-asymm-0.10.csv

# symmetric division set
head -n1 output-0.500-0.500-4.csv > output-symm.csv
# concatenate results, prune individual headers, and append to big file
cat output-0.500-0.100-4.csv output-0.500-0.100-16.csv output-0.500-0.100-64.csv output-0.500-0.500-4.csv output-0.500-0.500-16.csv output-0.500-0.500-64.csv | grep -v "h" >> output-symm.csv

# plot cell snapshots -- no arguments because a long list of different input files is implemented therein
Rscript plotcell.R

# ystar = 0,0.404, and 0.687
# moment plots
Rscript compare-both-asymm.R output-moments.csv 0 100 50 0.5 0.500 diagnostics 1
Rscript compare-both-asymm.R output-moments.csv 1 100 50 0.5 0.500 diagnostics 1

# ystar = 0.000
# sim-theory comparison for non-repulsive case
Rscript compare-both-asymm.R output-asymm-0.50.csv 0 100 50 0.1 0.500 asymm-0.1-50 0
#Rscript compare-both-asymm.R output-asymm-0.50.csv 0 100 50 0.5 0.500 asymm-0.5-50 0
# sim-theory comparison for repulsive case
Rscript compare-both-asymm.R output-asymm-0.50.csv 1 100 50 0.1 0.500 asymm-0.1-50 0
#Rscript compare-both-asymm.R output-asymm-0.50.csv 1 100 50 0.5 0.500 asymm-0.5-50 0

# ystar = 0.404
# sim-theory comparison for non-repulsive case
Rscript compare-both-asymm.R output-asymm-0.25.csv 0 100 50 0.1 0.250 asymm-0.1-25 0
#Rscript compare-both-asymm.R output-asymm-0.25.csv 0 100 50 0.5 0.250 asymm-0.5-25 0
# sim-theory comparison for repulsive case
Rscript compare-both-asymm.R output-asymm-0.25.csv 1 100 50 0.1 0.250 asymm-0.1-25 0
#Rscript compare-both-asymm.R output-asymm-0.25.csv 1 100 50 0.5 0.250 asymm-0.5-25 0

# ystar = 0.687
# sim-theory comparison for non-repulsive case
Rscript compare-both-asymm.R output-asymm-0.10.csv 0 100 50 0.1 0.100 asymm-0.1-10 0
#Rscript compare-both-asymm.R output-asymm-0.10.csv 0 100 50 0.5 0.100 asymm-0.5-10 0
# sim-theory comparison for repulsive case
Rscript compare-both-asymm.R output-asymm-0.10.csv 1 100 50 0.1 0.100 asymm-0.1-10 0
#Rscript compare-both-asymm.R output-asymm-0.10.csv 1 100 50 0.5 0.100 asymm-0.5-10 0

# sim-theory comparison for symmetric cases:
# non-repulsive case
Rscript compare-both-asymm.R output-symm.csv 0 100 50 0.5 0.500 symm-50 0
# repulsive case
Rscript compare-both-asymm.R output-symm.csv 1 100 50 0.5 0.500 symm-50 0

# effect of lambda
Rscript plotlambda.R output-symm.csv

# compile manuscript
#pdflatex writeup
#bibtex writeup
#pdflatex writeup
#pdflatex writeup
