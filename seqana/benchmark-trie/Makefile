include benchmark.config

all:	benchmark.dat visualize

benchmark.dat:	benchmark.config
	cd ..;make create_trie_input.x trie_construct.x mp_search.x
	cd ../../testdata;make $(TCFILES)
	@echo -n "file input_length num_strings grep_hits grep_time" | tee benchmark.dat
	@for talgorithm in $(TALGORITHMS) ; do \
		for talgmeasure in trie_size fm_mem_peak fm_time trie_mem_peak trie_time mpsearch_hits mpsearch_mem_peak mpsearch_time; do \
			echo -n " $$talgorithm" | tee -a benchmark.dat; \
			echo -n "_$$talgmeasure" | tee -a benchmark.dat; \
		done; \
	done
	@for tcfile in $(TCFILEPATHS) ; do \
		tcname=$$(basename "$$tcfile" | tr '_' '-'); \
		cat $$tcfile | LC_ALL=C awk 'length>=10' > $$tcname.tin; \
		../create_trie_input.x $$tcname.tin $$tcname-tmp.tin; \
		LC_ALL=C sort $$tcname-tmp.tin > $$tcname.tin; rm $$tcname-tmp.tin; \
		echo -n "\n$$tcname" | tee -a benchmark.dat; \
		tclength=$$(wc -c < $$tcname.tin); \
		tcnumstrings=$$(wc -l < $$tcname.tin); \
		echo -n " $$tclength $$tcnumstrings" | tee -a benchmark.dat; \
		tail -n 1 $$tcname.tin > $$tcname.pattern; \
		grephits=$$(grep -c -F -f $$tcname.pattern $$tcname.tin); \
		grepstats=$$(time -f "%U %S" sh -c "seq 1 1000 | xargs -P 1 -I{} grep -c -F -f $$tcname.pattern $$tcname.tin > /dev/null" 2>&1 1>/dev/null ); \
		rm -f $$tcname.pattern; \
		gusrtime=$$(echo "$$grepstats" | awk '{print $$1}'); \
		gsystime=$$(echo "$$grepstats" | awk '{print $$2}'); \
		greptime=$$(echo "scale=2;($$gusrtime+$$gsystime)*1000/1000" | bc); \
		echo -n " $$grephits $$greptime" | tee -a benchmark.dat; \
		cat $$tcname.tin | LC_ALL=C awk 'BEGIN{ sum=0} {sum+=length($$0); printf "%d\t%d\n", sum, NR-1}' | sort > memcheckfile ;\
		for talgorithm in $(TALGORITHMS) ; do \
			tconstructinfo=$$(../trie_construct.x -i -sa $(SAALGORITHM) -ta $$talgorithm $$tcname.tin $$tcname-$$talgorithm.trie) ; \
			tconstructinfo=$$(echo -n "$$tconstructinfo" | tail -n 5 | awk '{print $$2}' | tr '\n' ' ' | head -c -1); \
			echo -n " $$tconstructinfo" | tee -a benchmark.dat; \
			cat $$tcname.tin | tr -d '\n' | ../mp_search.x -i $$tcname-$$talgorithm.trie > memresult ; \
			mpinfo=$$(cat memresult | tail -n 3 | awk '{print $$2}' | tr '\n' ' ' | head -c -1); \
			echo -n " $$mpinfo" | tee -a benchmark.dat; \
			notfound=$$(cat memresult | head -n -7 | sort | comm --nocheck-order -23 memcheckfile -); \
			if [ $${#notfound} -ne 0 ] ; then \
				echo -n "Failed" | tee -a benchmark.dat; \
			fi; \
			rm -f $$tcname-$$talgorithm.trie; \
			rm -f memresult; \
		done; \
		rm -f $$tcname.tin; \
		rm -f memcheckfile; \
	done



#pdf document
visualize:
	pdflatex triebenchmark.tex
	rm -f triebenchmark.log
	rm -f triebenchmark.out
	rm -f triebenchmark.aux

#cleaners
clean:
	rm -f *.tin
	rm -f *.trie
	rm -f *.pattern
	rm -f memresult
	rm -f memcheckfile
	rm -f *.dat
	rm -f *.sdsl
