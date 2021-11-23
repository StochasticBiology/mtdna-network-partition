# analysis script
# to be run AFTER simulations have terminated (no error checking)

# concatenate output files from different simulation parameterisations into one big file
# get header
head -n1 output-0.100-4.csv > output-full.csv
# concatenate results, prune individual headers, and append to big file
cat output-0.100-4.csv output-0.500-4.csv output-0.900-4.csv output-0.100-16.csv output-0.500-16.csv output-0.900-16.csv output-0.100-64.csv output-0.500-64.csv output-0.900-64.csv | grep -v "h" >> output-full.csv

# plot cell snapshots -- no arguments because a long list of different input files is implemented therein
Rscript plotcell.R
# sim-theory comparison for non-repulsive case
Rscript compare-both-new.R output-full.csv 0 100 50 test
# sim-theory comparison for repulsive case
Rscript compare-both-new.R output-full.csv 1 100 50 test
# effect of lambda
Rscript plotlambda.R output-full.csv

# compile manuscript
pdflatex writeup
bibtex writeup
pdflatex writeup
pdflatex writeup
