

AVX_512 := 1
#AVX2    := 1
#NOVEC    := 1

CPPFLAGS := -std=c++11 # this is hwere to add include links 
# LDFLAGS  := # nothing yet but this is hwere to add library links 

ifdef INTEL_LICENSE_FILE
	CXX=icpc
	CPPFLAGS += -qopenmp
else
	CXX=g++
	CPPFLAGS += -fopenmp
endif

OPT = -g -O3
ifdef AVX_512
	ifeq (${CXX},icpc)
		OPT  += -qopt-zmm-usage=high -xCORE-AVX512 -qopt-report=5 
	else
		OPT  += -mavx512f -mavx512cd
	endif
else ifdef AVX2
	OPT  += -mavx2
else ifdef NOVEC
	OPT  += -no-vec
else
	OPT  += -mavx
endif

ifdef USE_CALI 
  CPPFLAGS += -I${CALIPER_DIR}/include
  CPPFLAGS += -DUSE_CALI
  LDFLAGS  += -L${CALIPER_DIR}/lib64 
  LDFLAGS  += -lcaliper
endif

ifdef WITH_ROOT
  CPPFLAGS += $(shell root-config --cflags)
  LDFLAGS  += $(shell root-config --libs)
else
  CPPFLAGS += -DNO_ROOT
endif

gshf-mrqdt3:
	${CXX} ${OPT} ${CPPFLAGS} -o gshf-mrqdt3 gshf-mrqdt3.cc marqfit.cc Event.cc  ${LDFLAGS}

gshf-mrqdt3-gmac:
	${CXX} -std=c++11 ${OPT} ${CPPFLAGS} -o gshf-mrqdt3 gshf-mrqdt3.cc marqfit.cc Event.cc ${LDFLAGS}

clean:	
	rm -f gshf-mrqdt3 gshf-mrqdt2 gshf-levmar *.optrpt *.o 

