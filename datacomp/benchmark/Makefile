include benchmark.config

#max runtime per command (used in timeout)
MAXTIME := 180m

TINFOBESTPS=best0 best2 best4 best6 best8 best10 best12 best14 best16 best18 best20 best22 best24 best26 best28 best30 best32 best34 best36 best38 best40 best42 best44 best46 best48 best50 best52 best54 best56 best58 best60 best62 best64 best66 best68 best70 best72 best74 best76 best78 best80 best82 best84 best86 best88 best90 best92 best94 best96 best98 best100
TINFOSTRATS=$(TINFOBESTPS) greedy hirsch
TINFOPOSTSTAGES=bcm bw94
TPOSTSTAGES=bcm bw94

COMPRESSORS=$(basename $(shell ls cp))

all:	install benchmark visualize

#installation
install:
	for cp in $(COMPRESSORS) ; do \
		chmod a+x cp/$$cp.sh; \
		cp/$$cp.sh i; \
	done

#benchmark itself
benchmark:	result.dat result_dbg.dat result_tinfo.dat

#visualization
visualize:	benchmark_visualize cm_relcomp_cm cm_relcomp_bcm cm_relcomp_bw94

#### BENCHMARKS ####

bin/bwzip.x:
	cd ..; make bwzip.x
	cp ../bwzip.x bin/bwzip.x

#compression benchmark
result.dat:	bin/bwzip.x benchmark.config
	cd ../../testdata;make $(TCFILES)
	@echo -n "file" | tee result.dat;
	@for method in $(COMPRESSORS) ; do \
		echo -n " $$method-bps $$method-decode-speed $$method-decode-membps $$method-encode-speed $$method-encode-membps" >> result.dat; \
		echo -n "\t$$method-bps"; \
	done
	@echo " bps-best decode-speed-best decode-membps-best encode-speed-best encode-membps-best" >> result.dat
	@echo ""
	@for tcfile in $(TCFILEPATHS) ; do \
		tcname=$$(basename "$$tcfile" | tr '_' '-'); \
		echo -n $$tcname | tee -a result.dat; \
		tcfs=$$(stat -c%s $$tcfile); \
		bpsbest="file"; \
		mincompfs=0; \
		dsbest="file"; \
		mindtime=0; \
		dbpsbest="file"; \
		mindmem=0; \
		esbest="file"; \
		minetime=0; \
		ebpsbest="file"; \
		minemem=0; \
		for compressor in $(COMPRESSORS) ; do \
			echo -n "" > tmp/fcomp; \
			echo -n "" > tmp/fres; \
			estats=$$(timeout $(MAXTIME) time -f "%U %S %M" cp/$$compressor.sh c $$tcfile tmp/fcomp  2>&1 1>/dev/null ); \
			if [ $$? -ne 0 ] ; then \
				echo -n " 0 0 $$tcfs 0 $$tcfs" >> result.dat; \
				echo -n "\t[CMP_TIMEOUT]"; \
				continue; \
			fi; \
			compfs=$$(stat -c%s tmp/fcomp); \
			dstats=$$(timeout $(MAXTIME) time -f "%U %S %M" cp/$$compressor.sh d tmp/fcomp tmp/fres  2>&1 1>/dev/null ); \
			if [ $$? -ne 0 ] ; then \
				echo -n " 0 0 $$tcfs 0 $$tcfs 0" >> result.dat; \
				echo -n "\t[DECOMP_TIMEOUT]"; \
				continue; \
			fi; \
			cmp -s $$tcfile tmp/fres; \
			if [ $$? -ne 0 ] ; then \
				echo -n " 0 0 $$tcfs 0 $$tcfs 0" >> result.dat; \
				echo -n "\t[FAIL]"; \
				continue; \
			fi; \
			eusrtime=$$(echo $$estats | awk '{print $$1}'); \
			esystime=$$(echo $$estats | awk '{print $$2}'); \
			etime=$$(echo "scale=3;$$eusrtime+$$esystime+0.001" | bc); \
			emem=$$(echo $$estats | awk '{print $$3}'); \
			emem=$$(echo "$$emem*1024" | bc); \
			dusrtime=$$(echo $$dstats | awk '{print $$1}'); \
			dsystime=$$(echo $$dstats | awk '{print $$2}'); \
			dtime=$$(echo "scale=3;$$dusrtime+$$dsystime+0.001" | bc); \
			dmem=$$(echo $$dstats | awk '{print $$3}'); \
			dmem=$$(echo "$$dmem*1024" | bc); \
			bps=$$(echo "scale=3;$$compfs*8/$$tcfs" | bc | sed 's/^\./0./'); \
			dspeed=$$(echo "scale=2;$$tcfs/($$dtime*1024*1024)" | bc | sed 's/^\./0./'); \
			espeed=$$(echo "scale=2;$$tcfs/($$etime*1024*1024)" | bc | sed 's/^\./0./'); \
			emembps=$$(echo "scale=2;$$emem*8/$$tcfs" | bc | sed 's/^\./0./'); \
			dmembps=$$(echo "scale=2;$$dmem*8/$$tcfs" | bc | sed 's/^\./0./'); \
			etime=$$(echo "scale=0;$$etime*1000/1" | bc | sed 's/^\./0./'); \
			dtime=$$(echo "scale=0;$$dtime*1000/1" | bc | sed 's/^\./0./'); \
			echo -n " $$bps $$dspeed $$dmembps $$espeed $$emembps" >> result.dat; \
			echo -n "\t$$bps"; \
			rm tmp/fcomp; \
			rm tmp/fres; \
			if [ "$$mincompfs" -ge "$$compfs" ] || [ "$$mincompfs" -eq 0 ]; then \
				mincompfs="$$compfs"; \
				bpsbest=$$compressor; \
			fi; \
			if [ "$$mindmem" -ge "$$dmem" ] || [ "$$mindmem" -eq 0 ]; then \
				mindmem="$$dmem"; \
				dbpsbest=$$compressor; \
			fi; \
			if [ "$$minemem" -ge "$$emem" ] || [ "$$minemem" -eq 0 ]; then \
				minemem="$$emem"; \
				ebpsbest=$$compressor; \
			fi; \
			if [ "$$mindtime" -ge "$$dtime" ] || [ "$$mindtime" -eq 0 ]; then \
				mindtime="$$dtime"; \
				dsbest=$$compressor; \
			fi; \
			if [ "$$minetime" -ge "$$etime" ] || [ "$$minetime" -eq 0 ]; then \
				minetime="$$etime"; \
				esbest=$$compressor; \
			fi; \
		done; \
		echo " $$bpsbest $$dsbest $$dbpsbest $$esbest $$ebpsbest" >> result.dat; \
		echo ""; \
	done

bin/tfmzip.x:
	cd ..; make tfmzip.x
	cp ../tfmzip.x bin/tfmzip.x

bin/tfm_index_construct.x:
	cd ../../seqana; make tfm_index_construct.x
	cp ../../seqana/tfm_index_construct.x bin/tfm_index_construct.x

#de bruijn graph benchmark
result_dbg.dat:	bin/tfmzip.x bin/tfm_index_construct.x benchmark.config
	cd ../../testdata;make $(TCFILES)
	@echo -n "file" | tee result_dbg.dat
	@for method in $(TPOSTSTAGES) ; do \
		echo -n " $$method-bps" >> result_dbg.dat; \
		echo -n "\t$$method-bps"; \
	done
	@echo "" >> result_dbg.dat
	@echo ""
	@for tcfile in $(TCFILEPATHS) ; do \
		tcname=$$(basename "$$tcfile" | tr '_' '-'); \
		echo -n $$tcname | tee -a result_dbg.dat; \
		tcinput="tmp/text_"$$(basename "$$tcfile")".sdsl"; \
		cat $$tcfile | tr -d '\000' > $$tcinput; \
		echo -n "\x00" >> $$tcinput; \
		tcfs=$$(stat -c%s $$tcinput); \
		bin/tfm_index_construct.x $$tcinput tmp/$$tcname.tfm; \
		rm -f $$tcinput; \
		for method in $(TPOSTSTAGES) ; do \
			echo -n "" > tmp/fcomp; \
			echo -n "" > tmp/fres; \
			bin/tfmzip.x -pstage $$method tmp/$$tcname.tfm tmp/fcomp; \
			bin/tfmzip.x -d tmp/fcomp tmp/fres; \
			compfs=$$(stat -c%s tmp/fcomp); \
			cmp -s tmp/$$tcname.tfm tmp/fres; \
			if [ $$? -ne 0 ] ; then \
				echo -n " 0" >> result_dbg.dat; \
				echo -n "\t[FAIL]"; \
				continue; \
			fi; \
			bps=$$(echo "scale=3;$$compfs*8/$$tcfs" | bc | sed 's/^\./0./'); \
			echo -n " $$bps" >> result_dbg.dat; \
			echo -n "\t$$bps"; \
			rm tmp/fcomp; \
			rm tmp/fres; \
		done; \
		rm -f tmp/$$tcname.tfm; \
		echo "" >> result_dbg.dat; \
		echo ""; \
	done

#tunneling info results
result_tinfo.dat: bin/bwzip.x benchmark.config
	cd ../../testdata;make $(TCFILES)
	@echo -n "file" | tee result_tinfo.dat
	@for method in $(TINFOPOSTSTAGES) ; do \
		for tstrat in $(TINFOSTRATS) ; do \
			echo -n " $$method-$$tstrat-bwt-construct-time" >> result_tinfo.dat ; \
			echo -n " $$method-$$tstrat-num-tunnels" >> result_tinfo.dat ; \
			echo -n " $$method-$$tstrat-exp-tunnelcosts" >> result_tinfo.dat ; \
			echo -n " $$method-$$tstrat-num-rle-tc" >> result_tinfo.dat ; \
			echo -n " $$method-$$tstrat-exp-benefit" >> result_tinfo.dat ; \
			echo -n " $$method-$$tstrat-tunneling-time" >> result_tinfo.dat ; \
			echo -n " $$method-$$tstrat-size-bwtenc" >> result_tinfo.dat ; \
			echo -n " $$method-$$tstrat-size-auxenc" >> result_tinfo.dat ; \
			echo -n " $$method-$$tstrat-enc-time" >> result_tinfo.dat ; \
		done ; \
	done
	@echo "" | tee -a result_tinfo.dat
	@for tcfile in $(TCFILEPATHS) ; do \
		tcname=$$(basename "$$tcfile" | tr '_' '-'); \
		echo -n $$tcname | tee -a result_tinfo.dat; \
		head -c 1G $$tcfile > tmp/input ; \
		echo -n "" > tmp/bw_commands ; \
		for method in $(TINFOPOSTSTAGES) ; do \
			for tstrat in $(TINFOSTRATS) ; do \
				echo "bin/bwzip.x -i -tstrat $$tstrat -pstage $$method tmp/input tmp/$$tstrat-$$method.tmp > tmp/$$tstrat-$$method.res" >> tmp/bw_commands ; \
			done; \
		done ; \
		xargs -d \\n -n 1 -a tmp/bw_commands sh -c; \
		rm -f tmp/bw_commands ; \
		rm tmp/*.tmp ; \
		rm -f tmp/input ; \
		for method in $(TINFOPOSTSTAGES) ; do \
			for tstrat in $(TINFOSTRATS) ; do \
				cat tmp/$$tstrat-$$method.res | awk '{printf " "$$2}' >> result_tinfo.dat; \
			done; \
		done ; \
		rm tmp/*.res ; \
		echo "" >> result_tinfo.dat ; \
		echo "" ; \
	done

#### VISUALIZATION ####

benchmark_visualize:
	@cat "result_tinfo.dat"   | awk 'NR > 1{ mi = $$(1+0*9 + 7) + $$(1+0*9 + 8); ma = $$(1+0*9+ 7) + $$(1+0*9+ 8);\
	                                 for(i=0;i<=50;i++) {mi=(mi>$$(1+(0+i)*9 + 7)+$$(1+(0+i)*9 + 8))?$$(1+(0+i)*9 + 7)+$$(1+(0+i)*9 + 8):mi; \
	                                                     ma=(ma<$$(1+(0+i)*9 + 7)+$$(1+(0+i)*9 + 8))?$$(1+(0+i)*9 + 7)+$$(1+(0+i)*9 + 8):ma;}; \
		                         for(i=0;i<=50;i++) { printf "%d %d %f\n", \
		                                              2*(i), \
		                                              NR-1, \
		                                              ($$(1+(0+i)*9 + 7)+$$(1+(0+i)*9 + 8)-mi)/(ma-mi+1); };printf "\n"}' > t_matrix_bcm.dat
	@echo "file imp_best_percent imp_greedy_percent imp_hirsch_percent rnt_greedy_percent rnt_hirsch_percent imp_best imp_greedy imp_hirsch" > t_improve_bcm.dat
	@cat "result_tinfo.dat"   | awk 'NR > 1{ mi = 0;\
	                                 for(i=0;i<=50;i++) {mi=($$(1+(0+mi)*9 + 7)+$$(1+(0+mi)*9 + 8)>$$(1+(0+i)*9 + 7)+$$(1+(0+i)*9 + 8))?i:mi;}; \
					 printf "%s %.2f %.2f %.2f %.2f %.2f %d %d %d\n", $$(1), \
					    100 - 100*($$(1+(0+mi)*9 + 7)+$$(1+(0+mi)*9 + 8))/($$(1+(0+0)*9 + 7) + $$(1+(0+0)*9 + 8)), \
					    100 - 100*($$(1+(0+51)*9 + 7)+$$(1+(0+51)*9 + 8))/($$(1+(0+0)*9 + 7) + $$(1+(0+0)*9 + 8)), \
					    100 - 100*($$(1+(0+52)*9 + 7)+$$(1+(0+52)*9 + 8))/($$(1+(0+0)*9 + 7) + $$(1+(0+0)*9 + 8)), \
					    100 * $$(1+(0+51)*9 + 2) / $$(1+(0+50)*9 + 2), \
					    100 * $$(1+(0+52)*9 + 2) / $$(1+(0+50)*9 + 2), \
					    ($$(1+(0+0)*9 + 7) + $$(1+(0+0)*9 + 8)) - ($$(1+(0+mi)*9 + 7)+$$(1+(0+mi)*9 + 8)), \
					    ($$(1+(0+0)*9 + 7) + $$(1+(0+0)*9 + 8)) - ($$(1+(0+51)*9 + 7)+$$(1+(0+51)*9 + 8)), \
					    ($$(1+(0+0)*9 + 7) + $$(1+(0+0)*9 + 8)) - ($$(1+(0+52)*9 + 7)+$$(1+(0+52)*9 + 8)) }' >> t_improve_bcm.dat
	@cat "result_tinfo.dat"   | awk 'NR > 1{ mi = $$(1+53*9 + 7) + $$(1+53*9 + 8); ma = $$(1+53*9+ 7) + $$(1+53*9+ 8);\
	                                 for(i=0;i<=50;i++) {mi=(mi>$$(1+(53+i)*9 + 7)+$$(1+(53+i)*9 + 8))?$$(1+(53+i)*9 + 7)+$$(1+(53+i)*9 + 8):mi; \
	                                                     ma=(ma<$$(1+(53+i)*9 + 7)+$$(1+(53+i)*9 + 8))?$$(1+(53+i)*9 + 7)+$$(1+(53+i)*9 + 8):ma;}; \
		                         for(i=0;i<=50;i++) { printf "%d %d %f\n", \
		                                              2*(i), \
		                                              NR-1, \
		                                              ($$(1+(53+i)*9 + 7)+$$(1+(53+i)*9 + 8)-mi)/(ma-mi+1); };printf "\n"}' > t_matrix_bw94.dat
	@echo "file imp_best_percent imp_greedy_percent imp_hirsch_percent rnt_greedy_percent rnt_hirsch_percent imp_best imp_greedy imp_hirsch" > t_improve_bw94.dat
	@cat "result_tinfo.dat"   | awk 'NR > 1{ mi = 0;\
	                                 for(i=0;i<=50;i++) {mi=($$(1+(53+mi)*9 + 7)+$$(1+(53+mi)*9 + 8)>$$(1+(53+i)*9 + 7)+$$(1+(53+i)*9 + 8))?i:mi;}; \
					 printf "%s %.2f %.2f %.2f %.2f %.2f %d %d %d\n", $$(1), \
					    100 - 100*($$(1+(53+mi)*9 + 7)+$$(1+(53+mi)*9 + 8))/($$(1+(53+0)*9 + 7) + $$(1+(53+0)*9 + 8)), \
					    100 - 100*($$(1+(53+51)*9 + 7)+$$(1+(53+51)*9 + 8))/($$(1+(53+0)*9 + 7) + $$(1+(53+0)*9 + 8)), \
					    100 - 100*($$(1+(53+52)*9 + 7)+$$(1+(53+52)*9 + 8))/($$(1+(53+0)*9 + 7) + $$(1+(53+0)*9 + 8)), \
					    100 * $$(1+(53+51)*9 + 2) / $$(1+(53+50)*9 + 2), \
					    100 * $$(1+(53+52)*9 + 2) / $$(1+(53+50)*9 + 2), \
					    ($$(1+(53+0)*9 + 7) + $$(1+(53+0)*9 + 8)) - ($$(1+(53+mi)*9 + 7)+$$(1+(53+mi)*9 + 8)), \
					    ($$(1+(53+0)*9 + 7) + $$(1+(53+0)*9 + 8)) - ($$(1+(53+51)*9 + 7)+$$(1+(53+51)*9 + 8)), \
					    ($$(1+(53+0)*9 + 7) + $$(1+(53+0)*9 + 8)) - ($$(1+(53+52)*9 + 7)+$$(1+(53+52)*9 + 8)) }' >> t_improve_bw94.dat
	pdflatex benchmark.tex
	rm -f benchmark.log
	rm -f benchmark.out
	rm -f benchmark.aux

cm_relcomp_cm:
	@echo -n "file" > cm_relcomp_cm.dat
	@for bestp in $(TINFOBESTPS); do \
		echo -n " $$bestp" >> cm_relcomp_cm.dat; \
	done
	@echo "" >> cm_relcomp_cm.dat
	@cat "result_tinfo.dat"   | awk 'NR > 1{ mi = $$(1+0*9 + 5)-$$(1+0*9 + 3); ma = $$(1+0*9 + 5)-$$(1+0*9 + 3); \
	                                 for(i=1;i<=50;i++) {mi=(mi>$$(1+(0+i)*9 + 5)-$$(1+(0+i)*9 + 3))?$$(1+(0+i)*9 + 5)-$$(1+(0+i)*9 + 3):mi; \
	                                                     ma=(ma<$$(1+(0+i)*9 + 5)-$$(1+(0+i)*9 + 3))?$$(1+(0+i)*9 + 5)-$$(1+(0+i)*9 + 3):ma;}; \
					 printf "%s", $$(1); \
	                                 for(i=0;i<=50;i++) { \
					     printf " %.3f", 1 - ( $$(1+(0+i)*9 + 5)-$$(1+(0+i)*9 + 3) - mi ) / (ma - mi + 1) }; \
					 printf "\n"}' >> cm_relcomp_cm.dat

cm_relcomp_bcm:
	@echo -n "file" > cm_relcomp_bcm.dat
	@for bestp in $(TINFOBESTPS); do \
		echo -n " $$bestp" >> cm_relcomp_bcm.dat; \
	done
	@echo "" >> cm_relcomp_bcm.dat
	@cat "result_tinfo.dat"   | awk 'NR > 1{ mi = $$(1+0*9 + 7) + $$(1+0*9 + 8); ma = $$(1+0*9+ 7) + $$(1+0*9+ 8);\
	                                 for(i=1;i<=50;i++) {mi=(mi>$$(1+(0+i)*9 + 7)+$$(1+(0+i)*9 + 8))?$$(1+(0+i)*9 + 7)+$$(1+(0+i)*9 + 8):mi; \
	                                                     ma=(ma<$$(1+(0+i)*9 + 7)+$$(1+(0+i)*9 + 8))?$$(1+(0+i)*9 + 7)+$$(1+(0+i)*9 + 8):ma;}; \
					 printf "%s", $$(1); \
	                                 for(i=0;i<=50;i++) { \
					     printf " %.3f", ($$(1+(0+i)*9 + 7)+$$(1+(0+i)*9 + 8)-mi)/(ma-mi+1); }; \
					 printf "\n"}' >> cm_relcomp_bcm.dat

cm_relcomp_bw94:
	@echo -n "file" > cm_relcomp_bw94.dat
	@for bestp in $(TINFOBESTPS); do \
		echo -n " $$bestp" >> cm_relcomp_bw94.dat; \
	done
	@echo "" >> cm_relcomp_bw94.dat
	@cat "result_tinfo.dat"   | awk 'NR > 1{ mi = $$(1+53*9 + 7) + $$(1+53*9 + 8); ma = $$(1+53*9+ 7) + $$(1+53*9+ 8);\
	                                 for(i=1;i<=50;i++) {mi=(mi>$$(1+(53+i)*9 + 7)+$$(1+(53+i)*9 + 8))?$$(1+(53+i)*9 + 7)+$$(1+(53+i)*9 + 8):mi; \
	                                                     ma=(ma<$$(1+(53+i)*9 + 7)+$$(1+(53+i)*9 + 8))?$$(1+(53+i)*9 + 7)+$$(1+(53+i)*9 + 8):ma;}; \
					 printf "%s", $$(1); \
	                                 for(i=0;i<=50;i++) { \
					     printf " %.3f", ($$(1+(53+i)*9 + 7)+$$(1+(53+i)*9 + 8)-mi)/(ma-mi+1); }; \
					 printf "\n"}' >> cm_relcomp_bw94.dat

#clean up benchmark
clean:
	rm -f result.dat
	rm -f result_dbg.dat
	rm -f result_tinfo.dat
	rm -f t_matrix_bw94.dat
	rm -f t_improve_bw94.dat
	rm -f t_matrix_bcm.dat
	rm -f t_improve_bcm.dat
	rm -f cm_relcomp_cm.dat
	rm -f cm_relcomp_bcm.dat
	rm -f cm_relcomp_bw94.dat
	rm -f bin/*.x
	rm -f tmp/*.tmp
