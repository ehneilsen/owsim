.EXPORT_ALL_VARIABLES:
OPSIM_SKYMODEL_PYTHON := /data/des70.a/data/neilsen/singularity/opsim4-180320/home/opsim/repos/sims_skybrightness_pre/python/lsst/sims/skybrightness_pre/SkyModelPre.py
SIMS_SKYBRIGHTNESS_DATA := /data/des70.a/data/neilsen/sims_skybrightness_pre/data

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
		data/processed/exposures.txt

data/processed/events.txt: python/random_events.py etc/events.conf
	python $^ $@

data/processed/exposures.txt: python/followup.py etc/followup.conf \
		data/processed/events.txt data/munged/fieldID.txt
	python $^ $@ -e 0

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
