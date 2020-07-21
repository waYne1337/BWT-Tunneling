include benchmark.config

all:	benchmark.dat edgespectrum.dat dbgbenchmark.pdf

#normal benchmark
benchmark.dat:	benchmark.config
	cd ..;make tfm_index_construct.x tfm_index_invert.x 
	cd ../../testdata;make $(TCFILES)
	@echo -n "file input_length tfm_length fm_index_size tfm_index_size min_dbg_k min_dbg_edges fm_mem_peak fm_time min_dbg_mem_peak min_dbg_time tfm_mem_peak tfm_time" | tee benchmark.dat
	@for tcfile in $(TCFILEPATHS) ; do \
		tcname=$$(basename "$$tcfile" | tr '_' '-'); \
		tcinput="text_"$$(basename "$$tcfile")".sdsl"; \
		cat $$tcfile | tr -d '\000' > $$tcinput; \
		echo -n "\x00" >> $$tcinput; \
		echo -n "\n$$tcname" | tee -a benchmark.dat; \
		tcresult=$$(../tfm_index_construct.x -sa $(SAALGORITHM) -i $$tcinput tmp-$$tcname.tfm | awk '{print $$2}' | tr '\n' ' ' | head -c -1); \
		echo -n " $$tcresult" | tee -a benchmark.dat; \
		../tfm_index_invert.x tmp-$$tcname.tfm > tmp-$$tcname; \
		cmp -s $$tcinput tmp-$$tcname; \
		if [ $$? -ne 0 ] ; then \
			echo "Testcase failed (bad tbwt invertion)"; \
		fi; \
		rm -f $$tcinput; \
		rm -f tmp-$$tcname; \
		rm -f tmp-$$tcname.tfm; \
		rm -f *.sdsl; \
	done

#edge spectrum
edgespectrum.dat: benchmark.config
	cd ..;make dbg_edgespectrum.x
	cd ../../testdata;make $(TCFILES)
	@echo -n "File" | tee edgespectrum.dat
	@for i in $$(seq 1 $(MAXK)); do echo -n " k$$i" | tee -a edgespectrum.dat; done
	@for tcfile in $(TCFILEPATHS) ; do \
		tcname=$$(basename "$$tcfile" | tr '_' '-'); \
		tcinput="text_"$$(basename "$$tcfile")".sdsl"; \
		cat $$tcfile | tr -d '\000' > $$tcinput; \
		echo -n "\x00" >> $$tcinput; \
		echo -n "\n$$tcname" | tee -a edgespectrum.dat; \
		tcresult=$$(../dbg_edgespectrum.x $(MAXK) $$tcinput | awk '{print $$2}' | tr '\n' ' ' | head -c -1); \
		echo -n " $$tcresult" | tee -a edgespectrum.dat; \
		rm -f *.sdsl; \
	done

#edge matrix for better display
edgematrix.dat:	edgespectrum.dat
	@cat "edgespectrum.dat"   | awk 'NR > 1{ mi = $$2; ma = $$2; for(i=1;i<=100;i++) {mi=(mi>$$(1+i))?$$(1+i):mi;ma=(ma<$$(1+i))?$$(1+i):ma;}; \
		                                                 for(i=1;i<=100;i++) { printf "%d %d %f\n", \
		                                                                      i, \
		                                                                      NR-1, \
		                                                                      ($$(1+i)-mi)/(ma-mi); };printf "\n"}' > edgematrix.dat

#additional edge info
edgeinfo.dat:	edgespectrum.dat
	@head -n 1 edgespectrum.dat | tr -d '\n' > tmp.dat
	@echo ";kmin;kmax" >> tmp.dat
	@cat "edgespectrum.dat"   | awk 'NR > 1{ mi = $$2; ma = $$2; for(i=0;i<=100;i++) {mi=(mi>$$(2+i))?$$(2+i):mi;ma=(ma<$$(2+i))?$$(2+i):ma;}; \
		                                                 printf "%s;%d;%d\n", $$0, mi, ma;}' >> tmp.dat
	@cat tmp.dat | sed -e 's/_/-/g' > edgeinfo.dat
	@rm -f tmp.dat

#pdf document
dbgbenchmark.pdf:	edgematrix.dat	edgeinfo.dat benchmark.dat
	pdflatex dbgbenchmark.tex
	rm -f dbgbenchmark.log
	rm -f dbgbenchmark.out
	rm -f dbgbenchmark.aux

#cleaners
clean:
	rm -f *.aux
	rm -f *.out
	rm -f *.log
	rm -f *.dat
	rm -f *.sdsl
