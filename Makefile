.EXPORT_ALL_VARIABLES:
SIMS_SKYBRIGHTNESS_DATA := /home/opsim/repos/sims_skybrightness_pre/data

###############################################################################
# Utility
###############################################################################

.PHONY: clean
clean :
	rm data/collected/fieldID.dat

###############################################################################
# Data collection
###############################################################################

# This section contains all code the retrieves 
# data from outside sources; later sections depend
# only on files generated in this section.

collect: data/collected/fieldID.dat

data/collected/fieldID.dat:
	wget https://raw.githubusercontent.com/lsst/sims_skybrightness_pre/master/data/fieldID.dat -O $@

data/collected/baseline2018a.db:
	wget http://astro-lsst-01.astro.washington.edu:8081/db_gzip/baseline2018a.db.gz -O $@ && gunzip baseline2018a.db.gz

###############################################################################
# Data munging
###############################################################################

# This section contains code that manipulates input data files, putting them
# in a useful format. Processing and analysis should take place in later
# sections.

munge: data/munged/fieldID.txt

data/munged/fieldID.txt: data/collected/fieldID.dat
	tr "|" "\t" < $^ > $@

###############################################################################
# Processing
###############################################################################

process: data/processed/events.txt \
		data/processed/visits.txt \
		data/processed/conditions.txt \
		data/processed/result.db

data/processed/events.txt: python/random_events.py etc/events.conf
	python $^ $@

data/processed/visits.txt: python/apsupp.py python/followup.py etc/followup.conf \
		data/processed/events.txt data/munged/fieldID.txt
	python python/followup.py etc/followup.conf \
		data/processed/events.txt data/munged/fieldID.txt $@ -e 0

data/processed/conditions.txt: python/visitcond.py etc/visitcond.conf \
		data/processed/visits.txt
	python $^ $@ 

data/processed/result.db: python/owsched.py \
		data/collected/baseline2018a.db \
		data/processed/conditions.txt
	python $^ $@

###############################################################################
# Figure generation
###############################################################################

# plot: figures/sample_plot.png

# figures/sample_plot.png: R/sample_plot.R \
# 		data/munged/sample.RData
# 	$^ $@

###############################################################################
# Tests
###############################################################################

.PHONY: test
test :
	cd test && python -m unittest && cd ..
