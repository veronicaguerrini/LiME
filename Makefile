#CCLIB=-fopenmp -lsdsl -ldivsufsort -ldivsufsort64
VLIB= -g -O0

#LIB_DIR = ${HOME}/lib
#INC_DIR = ${HOME}/include
MY_CXX_FLAGS= -std=c++11 -Wall -DNDEBUG -fomit-frame-pointer -Wno-comment
#-D_FILE_OFFSET_BITS=64

M64 = 0
OMP = 1
DEBUG = 0
BIN = 1
EBWT = 1
HIGHER = 0

LIBOBJ = 
OMP_LIB =
##
ifeq ($(OMP),0)
	LIBOBJ = external/malloc_count/malloc_count.o
else
	OMP_LIB= -fopenmp 
endif
##

MY_CXX_OPT_FLAGS= -O3 -m64 $(OMP_LIB) 
#MY_CXX_OPT_FLAGS= $(VLIB) $(OMP_LIB)
MY_CXX=g++


alpha = 16
beta = 0.25

fileRefDB= Reference_database.csv
##

LFLAGS = -lm -ldl

DEFINES = -DDEBUG=$(DEBUG) -DM64=$(M64) -DOMP=$(OMP) -DBIN=$(BIN) -DEBWT=$(EBWT) -DHIGHER=$(HIGHER)

CXX_FLAGS=$(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) $(LFLAGS) $(DEFINES) -I$(INC_DIR) -L$(LIB_DIR)

##

all: compile

clean:
	\rm -f *.o  external/*.o external/malloc_count/*.o ClusterLCP ClusterBWT_DA Classify EGSAtoBCR

##

compile: ClusterLCP ClusterBWT_DA Classify EGSAtoBCR

ClusterLCP: src/ClusterLCP.cpp ${LIBOBJ} 
	$(MY_CXX) src/ClusterLCP.cpp $(CCLIB) -o ClusterLCP ${LIBOBJ} $(CXX_FLAGS) 

ClusterBWT_DA: src/ClusterBWT_DA.cpp ${LIBOBJ} 
	$(MY_CXX) src/ClusterBWT_DA.cpp $(CCLIB) -o ClusterBWT_DA ${LIBOBJ} $(CXX_FLAGS) 

Classify: src/Classify.cpp ${LIBOBJ} 
	$(MY_CXX) src/Classify.cpp $(CCLIB) -o Classify ${LIBOBJ} $(CXX_FLAGS) 

EGSAtoBCR: src/EGSAtoBCR.cpp ${LIBOBJ} 
	$(MY_CXX) src/EGSAtoBCR.cpp $(CCLIB) -o EGSAtoBCR ${LIBOBJ} $(CXX_FLAGS) 
