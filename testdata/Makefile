include data.config
include ../Make.helper

CORPORA := $(shell find * -maxdepth 0 -type d)

all:	$(CORPORA)

#canterbury corpus
canterbury: $(CANTERBURY)

$(CANTERBURY):
	curl http://www.data-compression.info/files/corpora/canterburycorpus.zip -o canterbury/cc.zip
	unzip -u canterbury/cc.zip -d canterbury
	rm canterbury/cc.zip

#large canterbury corpus
largecanterbury: $(LARGECANTERBURY)

$(LARGECANTERBURY):
	curl http://www.data-compression.info/files/corpora/largecanterburycorpus.zip -o largecanterbury/lcc.zip
	unzip -u largecanterbury/lcc.zip -d largecanterbury
	rm largecanterbury/lcc.zip

#silesia corpus
silesia: $(SILESIA)

$(SILESIA):
	curl http://sun.aei.polsl.pl/~sdeor/corpus/$(@F).bz2 -o $@.bz2
	bunzip2 -f $@.bz2

#pizza chili corpus
pizzachili: $(PIZZACHILI)

$(PIZZACHILI):
	$(eval fdir := \
		$(if $(findstring $(basename $(@F)),sources),code,\
		$(if $(findstring $(basename $(@F)),pitches),music,\
		$(if $(findstring $(basename $(@F)),proteins),protein,\
		$(if $(findstring $(basename $(@F)),dna),dna,\
		$(if $(findstring $(basename $(@F)),english),nlang,\
		$(if $(findstring $(basename $(@F)),dblp.xml),xml,\
			$(error unknown pizza chili category of $(@F) ))))))))
	curl http://pizzachili.dcc.uchile.cl/texts/$(fdir)/$(@F).gz -o $@.gz
	gunzip -f $@.gz

#repetitive corpus
repetitive: $(REPETITIVE)

$(REPETITIVE):
	curl http://pizzachili.dcc.uchile.cl/repcorpus/real/$(@F).gz -o $@.gz
	gunzip -f $@.gz

#genomes
genomes: $(GENOMES)

$(GENOMES):
	curl ftp://hgdownload.cse.ucsc.edu/goldenPath/$(@F)/bigZips/$(@F).2bit -o $@.2bit
	chmod a+x twoBitToFa
	./twoBitToFa -noMask $@.2bit $@.fa
	rm $@.2bit
	tr -d -c 'ACGT' < $@.fa > $@
	rm $@.fa

#statistics
statistics:	filestat.pdf

filestat.x:	filestat.cpp
	$(MY_CXX) -Wall -Wextra $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) $(C_OPTIONS) filestat.cpp -o filestat.x

filestat.dat:	filestat.x $(CORPORA)
	@echo "file sigma newlines size category source" > filestat.dat
	@echo "Computing stats of Canterbury corpus..."
	@for file in $(CANTERBURY) ; do \
		basename $$file | tr '_' '-' | tr -d '\n' >> filestat.dat; \
		echo -n " " >> filestat.dat; \
		./filestat.x $$file >> filestat.dat; \
		echo " canterbury http://www.data-compression.info/Corpora/CanterburyCorpus/index.html" >> filestat.dat; \
	done
	@echo "Computing stats of large Canterbury corpus..."
	@for file in $(LARGECANTERBURY) ; do \
		basename $$file | tr '_' '-' | tr -d '\n' >> filestat.dat; \
		echo -n " " >> filestat.dat; \
		./filestat.x $$file >> filestat.dat; \
		echo " largecanterbury http://www.data-compression.info/Corpora/CanterburyCorpus/index.html" >> filestat.dat; \
	done
	@echo "Computing stats of Silesia corpus..."
	@for file in $(SILESIA) ; do \
		basename $$file | tr '_' '-' | tr -d '\n' >> filestat.dat; \
		echo -n " " >> filestat.dat; \
		./filestat.x $$file >> filestat.dat; \
		echo " silesia http://sun.aei.polsl.pl/~sdeor/index.php?page=silesia" >> filestat.dat; \
	done
	@echo "Computing stats of Pizza and Chili corpus..."
	@for file in $(PIZZACHILI) ; do \
		basename $$file | tr '_' '-' | tr -d '\n' >> filestat.dat; \
		echo -n " " >> filestat.dat; \
		./filestat.x $$file >> filestat.dat; \
		echo " pizzachili http://pizzachili.dcc.uchile.cl/texts.html" >> filestat.dat; \
	done
	@echo "Computing stats of Repetitive corpus..."
	@for file in $(REPETITIVE) ; do \
		basename $$file | tr '_' '-' | tr -d '\n' >> filestat.dat; \
		echo -n " " >> filestat.dat; \
		./filestat.x $$file >> filestat.dat; \
		echo " repetitive http://pizzachili.dcc.uchile.cl/repcorpus.html" >> filestat.dat; \
	done
	@echo "Computing stats of Genomes..."
	@for file in $(GENOMES) ; do \
		basename $$file | tr '_' '-' | tr -d '\n' >> filestat.dat; \
		echo -n " " >> filestat.dat; \
		./filestat.x $$file >> filestat.dat; \
		echo " genomes http://hgdownload.soe.ucsc.edu/downloads.html" >> filestat.dat; \
	done

filestat.pdf:	filestat.tex filestat.dat
	pdflatex filestat.tex
	rm -f filestat.log
	rm -f filestat.out
	rm -f filestat.aux


#cleaners
clean-%:
	rm -f $*/*

cleanall: $(addprefix clean-,$(CORPORA))
	rm -f filestat.x
	rm -f *.dat
