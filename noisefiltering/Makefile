

# AVX_512 := 1
#AVX2    := 1
#NOVEC    := 1

# USE_CALI := 1
USE_MKL  := 1
# USE_FFTW := 1
THREAD_WIRES := 1
# MAKE_PLOTS  := 1

CPPFLAGS := -std=c++11 # this is hwere to add include links 
# LDFLAGS  := # nothing yet but this is where to add library links 

# MKL_LINK = -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl #
# MKL_LINK = -L${MKLROOT}/lib/intel64 -mkl -liomp5 -lpthread -lm -ldl #
MKL_LINK = -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_avx512 -lmkl_core -liomp5 -lpthread -lm -ldl #
MKL_OPT  = -DMKL_ILP64 -I${MKLROOT}/include #

# FFTW_DIR = /home/users/gravelle/soft/fftw-3.3.8-gcc
# FFTW_LINK = -L${FFTW_DIR}/lib -lfftw3f -lfftw3 -lfftw3f_threads -lm #
# FFTW_OPT  = -I${FFTW_DIR}/include #

FFTW_LINK = -L${FFTW_DIR}/lib -lfftw3f -lfftw3 -lfftw3f_omp -lm #
FFTW_OPT  = -I${FFTW_DIR}/include #

ifdef INTEL_LICENSE_FILE
	CXX=icpc
	CPPFLAGS += -qopenmp
else
	CXX=g++-9
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

ifdef USE_MKL 
  CPPFLAGS += ${MKL_OPT}
  CPPFLAGS += -DUSE_MKL
  LDFLAGS  += ${MKL_LINK}
endif

ifdef THREAD_WIRES 
  CPPFLAGS += -DTHREAD_WIRES
endif

ifdef USE_FFTW 
  CPPFLAGS += ${FFTW_OPT}
  CPPFLAGS += -DUSE_FFTW
  LDFLAGS  += ${FFTW_LINK}
endif

ifdef MAKE_PLOTS 
  CPPFLAGS += -DMAKE_PLOTS
endif


all: fft

fft: fftw_fft.cc utilities.cc Makefile
	$(CXX) $(CPPFLAGS) $(OPT) fftw_fft.cc utilities.cc -o fft.exe  $(LDFLAGS)

mpi: fft_comparison_mpi.cc Makefile
	mpiicc $(CPPFLAGS) $(OPT)  fft_comparison_mpi.cc -o fft_mpi.exe  $(LDFLAGS)

cufft: fft.cu Makefile
	nvcc -g -DUSE_CUDA -I/usr/local/cuda/inc fft.cu utilities.cc -o cufft.exe -L/usr/local/cuda/lib -lcufft

clean:	
	rm -f fft.exe fft_mpi.exe *.optrpt *.o

clena: clean

