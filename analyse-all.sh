# analysis script
# to be run AFTER simulations have terminated (no error checking)

# concatenate output files from different simulation parameterisations into one big file
# h = 0.5 set
# get header
head -n1 output-0.100-0.500-4.csv > output-asymm-0.5.csv
# concatenate results, prune individual headers, and append to big file
cat output-0.000-0.500-4.csv output-0.000-0.500-16.csv output-0.000-0.500-64.csv output-0.250-0.500-4.csv output-0.250-0.500-16.csv output-0.250-0.500-64.csv  output-0.800-0.500-4.csv output-0.800-0.500-16.csv output-0.800-0.500-64.csv | grep -v "h" >> output-asymm-0.5.csv

# h = 0.1 set
head -n1 output-0.100-0.100-4.csv > output-asymm-0.1.csv
# concatenate results, prune individual headers, and append to big file
cat  output-0.000-0.100-4.csv output-0.000-0.100-16.csv output-0.000-0.100-64.csv output-0.250-0.100-4.csv output-0.250-0.100-16.csv output-0.250-0.100-64.csv output-0.800-0.100-4.csv output-0.800-0.100-16.csv output-0.800-0.100-64.csv | grep -v "h" >> output-asymm-0.1.csv

# symmetric division set
head -n1 output-0.000-0.500-4.csv > output-symm.csv
# concatenate results, prune individual headers, and append to big file
cat output-0.000-0.500-4.csv output-0.000-0.500-16.csv output-0.000-0.500-64.csv | grep -v "h" >> output-symm.csv

# plot cell snapshots -- no arguments because a long list of different input files is implemented therein
Rscript plotcell.R
# sim-theory comparison for non-repulsive case
Rscript compare-both-asymm.R output-asymm-0.1.csv 0 100 50 0.1 0 asymm-0.1
Rscript compare-both-asymm.R output-asymm-0.5.csv 0 100 50 0.5 0 asymm-0.5
# sim-theory comparison for symmetric cases:
# non-repulsive case
Rscript compare-both-asymm.R output-symm.csv 0 100 50 0.5 0 symm
# repulsive case
Rscript compare-both-asymm.R output-symm.csv 1 100 50 0.5 0 symm-repel
# effect of lambda
Rscript plotlambda.R output-symm.csv

# compile manuscript
#pdflatex writeup
#bibtex writeup
#pdflatex writeup
#pdflatex writeup
