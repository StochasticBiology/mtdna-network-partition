# analysis script
# to be run AFTER simulations have terminated (no error checking)

# concatenate output files from different simulation parameterisations into one big file
# get header
head -n1 output-0.100-0.500-4.csv > output-asymm-0.5.csv
# concatenate results, prune individual headers, and append to big file
cat output-0.000-0.500-4.csv output-0.000-0.500-16.csv output-0.000-0.500-64.csv output-0.100-0.500-4.csv output-0.100-0.500-16.csv output-0.100-0.500-64.csv output-0.250-0.500-4.csv output-0.250-0.500-16.csv output-0.250-0.500-64.csv  output-0.800-0.500-4.csv output-0.800-0.500-16.csv output-0.800-0.500-64.csv | grep -v "h" >> output-asymm-0.5.csv

head -n1 output-0.100-0.100-4.csv > output-asymm-0.1.csv
# concatenate results, prune individual headers, and append to big file
cat  output-0.000-0.100-4.csv output-0.000-0.100-16.csv output-0.000-0.100-64.csv  output-0.100-0.100-4.csv output-0.100-0.100-16.csv output-0.100-0.100-64.csv output-0.250-0.100-4.csv output-0.250-0.100-16.csv output-0.250-0.100-64.csv output-0.800-0.100-4.csv output-0.800-0.100-16.csv output-0.800-0.100-64.csv | grep -v "h" >> output-asymm-0.1.csv



# plot cell snapshots -- no arguments because a long list of different input files is implemented therein
#Rscript plotcell.R
# sim-theory comparison for non-repulsive case
Rscript compare-both-asymm.R output-asymm-0.1.csv 0 100 50 0.1 asymm-0.1
Rscript compare-both-asymm.R output-asymm-0.5.csv 0 100 50 0.5 asymm-0.5
# sim-theory comparison for repulsive case
#Rscript compare-both-new.R output-full.csv 1 100 50 test
# effect of lambda
#Rscript plotlambda.R output-full.csv

# compile manuscript
#pdflatex writeup
#bibtex writeup
#pdflatex writeup
#pdflatex writeup
