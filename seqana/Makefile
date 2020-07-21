include ../Make.helper

DBGEXECUTABLES = tfm_index_construct.x tfm_index_invert.x dbg_edgespectrum.x
TRIEEXECUTABLES = create_trie_input.x trie_construct.x mp_search.x

all: $(DBGEXECUTABLES) $(TRIEEXECUTABLES)

#DBG EXECUTABLES
tfm_index_construct.x: include/tfm_index.hpp include/dbg_algorithms.hpp lib/tfm_index_construct.cpp
	$(MY_CXX) -Wall -Wextra $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) $(C_OPTIONS) \
	-I$(INC_DIR) -L$(LIB_DIR) -Iinclude -Llib lib/tfm_index_construct.cpp -o tfm_index_construct.x $(LIBS)

tfm_index_invert.x: include/tfm_index.hpp include/dbg_algorithms.hpp lib/tfm_index_invert.cpp
	$(MY_CXX) -Wall -Wextra $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) $(C_OPTIONS) \
	-I$(INC_DIR) -L$(LIB_DIR) -Iinclude -Llib lib/tfm_index_invert.cpp -o tfm_index_invert.x $(LIBS)

dbg_edgespectrum.x: include/dbg_algorithms.hpp lib/dbg_edgespectrum.cpp
	$(MY_CXX) -Wall -Wextra $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) $(C_OPTIONS) \
	-I$(INC_DIR) -L$(LIB_DIR) -Iinclude -Llib lib/dbg_edgespectrum.cpp -o dbg_edgespectrum.x $(LIBS)

#TRIE EXECUTABLES
create_trie_input.x: lib/create_trie_input.cpp
	$(MY_CXX) -Wall -Wextra $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) $(C_OPTIONS) \
	-I$(INC_DIR) -L$(LIB_DIR) -Iinclude -Llib lib/create_trie_input.cpp -o create_trie_input.x $(LIBS)

trie_construct.x: lib/trie_construct.cpp include/tries.hpp include/trie_xbwt.hpp include/trie_txbwt.hpp include/trie_construct_algorithms.hpp include/dbg_algorithms.hpp
	$(MY_CXX) -Wall -Wextra $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) $(C_OPTIONS) \
	-I$(INC_DIR) -L$(LIB_DIR) -Iinclude -Llib lib/trie_construct.cpp -o trie_construct.x $(LIBS)

mp_search.x: lib/mp_search.cpp include/tries.hpp include/trie_xbwt.hpp include/trie_txbwt.hpp
	$(MY_CXX) -Wall -Wextra $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) $(C_OPTIONS) \
	-I$(INC_DIR) -L$(LIB_DIR) -Iinclude -Llib lib/mp_search.cpp -o mp_search.x $(LIBS)

#CLEANERS
clean:
	rm -f *.sdsl
	rm -f *.x
