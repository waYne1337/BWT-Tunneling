include ../Make.helper
include postbwtstages/bw94/Make.helper
include postbwtstages/bcm/Make.helper

BW94_CC_LIBS  = $(addprefix postbwtstages/bw94/,$(BW94_LIBS))
BCM_CC_LIBS = $(addprefix postbwtstages/bcm/,$(BCM_LIBS))

all:	bwzip.x tfmzip.x

bwzip.x:	lib/bwzip.cpp postbwtstages/bw94/bw94_poststage.hpp postbwtstages/bcm/bcm_poststage.hpp include/*
	$(MY_CXX) -Wall -Wextra $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) $(C_OPTIONS) \
		-I$(INC_DIR) -L$(LIB_DIR) -Iinclude -Llib -Ipostbwtstages/bw94 -Lpostbwtstages/bw94 -Llib -Ipostbwtstages/bcm -Lpostbwtstages/bcm \
		$(BW94_CC_LIBS) $(BCM_CC_LIBS) $(CC_LIBS) lib/bwzip.cpp -o bwzip.x $(LIBS)

tfmzip.x:	lib/tfmzip.cpp postbwtstages/bw94/bw94_poststage.hpp postbwtstages/bcm/bcm_poststage.hpp include/*
	$(MY_CXX) -Wall -Wextra $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) $(C_OPTIONS) \
		-I$(INC_DIR) -L$(LIB_DIR) -Iinclude -Llib -Ipostbwtstages/bw94 -Lpostbwtstages/bw94 -Llib -Ipostbwtstages/bcm -Lpostbwtstages/bcm \
		-I../seqana/include $(BW94_CC_LIBS) $(BCM_CC_LIBS) $(CC_LIBS) lib/tfmzip.cpp -o tfmzip.x $(LIBS)

clean:
	rm -f *.x
