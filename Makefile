# WEBHOST := me@example.net
# WEBBASE := /var/www/html
# PGSERVICEFILE := etc/sample_pg_service.conf

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

process: data/processed/triggers.txt

data/processed/triggers.txt: python/random_triggers.py etc/triggers.conf
	python $^ $@

# process: data/processed/sample_summary.txt \
# 	data/processed/sample_mean.txt

# data/processed/sample_summary.txt: python/sample_summarize.py \
# 		data/collected/sample_collected.txt
# 	$^ $@

# data/processed/sample_mean.txt: R/sample_mean.R \
# 		data/munged/sample.RData
# 	$^ $@

###############################################################################
# Figure generation
###############################################################################

# plot: figures/sample_plot.png

# figures/sample_plot.png: R/sample_plot.R \
# 		data/munged/sample.RData
# 	$^ $@

###############################################################################
# Report generation
###############################################################################

# report: report/sample_index.html

# report/figures/%: figures/%
# 	mkdir -p report/figures
# 	cp $< $@

# report/data/collected/%: data/collected/%
# 	mkdir -p report/data/collected
# 	cp $< $@

# report/data/munged/%: data/munged/%
# 	mkdir -p report/data/munged
# 	cp $< $@

# report/data/processed/%: data/processed/%
# 	mkdir -p report/data/processed
# 	cp $< $@

# report/sample_index.html: report/sample_index.org \
# 		report/data/processed/sample_summary.txt \
# 		report/data/munged/sample.RData \
# 		report/figures/sample_plot.png
# 	emacs $< --batch -f org-html-export-to-html --kill

###############################################################################
# Publishing
###############################################################################

# .PHONY: pub_sample
# pub_sample: report/sample_index.html
# 	ssh ${WEBHOST} mkdir -p ${WEBBASE}/sample_report
# 	rsync -av report/ ${WEBHOST}:${WEBBASE}/sample_report
# 	ssh ${WEBHOST} mkdir -p ${WEBBASE}/sample_report/data
# 	rsync -av data/ ${WEBHOST}:${WEBBASE}/sample_report/data
# 	ssh ${WEBHOST} mkdir -p ${WEBBASE}/sample_report/figures
# 	rsync -av figures/ ${WEBHOST}:${WEBBASE}/sample_report/figures
