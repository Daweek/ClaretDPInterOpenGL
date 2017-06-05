################################################################################
#Makefile to Generate ClaretDPI 		Edg@r J 2017
################################################################################
#Compilers
GCC				= gcc
CXX 			= g++

CUDA 			= /usr/local/cuda-8.0
CUDA_SDK	= $(CUDA)/samples
NVCC     	= $(CUDA)/bin/nvcc

#Include Paths
CUDAINC   = -I. -I$(CUDA)/include -I$(CUDA_SDK)/common/inc

#Library Paths
CUDALIB		= -L/usr/lib/x86_64-linux-gnu -L$(CUDA)/lib64 \
						-lcuda -lcudart -lcudadevrt
GLLIB  		= -lGL -lGLU -lglut -lGLEW
LIB 			= $(CUDALIB) $(GLLIB) -lm

################ Choosing architecture code for GPU ############################
NVCC_ARCH			=
HOSTNAME		 	= $(shell uname -n)

ifeq ("$(HOSTNAME)","narumiken-msi")
	NVCC_ARCH		= -gencode arch=compute_61,code=sm_61
endif

ifeq ("$(HOSTNAME)","narumiken-LAP")
	NVCC_ARCH		= -gencode arch=compute_52,code=sm_52
endif

ifeq ("$(HOSTNAME)","Edgar-PC")
	NVCC_ARCH		= -gencode arch=compute_61,code=sm_61
endif

ifeq ("$(HOSTNAME)","Edgar-PC2")
	NVCC_ARCH		= -gencode arch=compute_35,code=sm_35
endif

###############	Debug, 0 -> False,  1-> True
DEBUGON						:= 0

ifeq (1,$(DEBUGON))
	CXXDEBUG 				:= -ggdb -pg
	CXXOPT					:= -O0
	NVCCDEBUG				:= -g -pg -G
	NVCCOPT					:= -O0
	NVCCFLAGSXCOMP 	:= -Xcompiler -g,-pg,-O0,-fopenmp
else
	CXXDEBUG 				:= 
	CXXOPT					:= -O3 -ffast-math -funroll-loops
	NVCCDEBUG				:= 
	NVCCOPT					:= -O3 --cudart=shared -use_fast_math
	NVCCFLAGSXCOMP 	:= -Xcompiler -O3,-ffast-math,-funroll-loops,-fopenmp
endif
###############################################################################
#NVCC_DP					= -rdc=true
CXXFLAGS				= $(CXXDEBUG) $(CXXOPT) -fopenmp
NVCCFLAGS 			= $(NVCCDEBUG) $(NVCC_DP) --compile $(NVCCOPT) $(NVCC_ARCH)
NVCCFLAGSLINK		= $(NVCCDEBUG) $(NVCC_DP) $(NVCCOPT) $(NVCC_ARCH)
###############################################################################

TARGET = claret

all: $(TARGET)

claret : cras36.o mr3.o 
	$(NVCC) $(NVCCFLAGSLINK) $(NVCCFLAGSXCOMP) $(CUDAINC) $< -o $@ mr3.o $(LIB) 

cras36.o: cras36.cpp cras36.h cras36def.h cras36GL.h
	$(CXX) $(CXXFLAGS) $(CUDAINC) -c $< -o $@ 

mr3.o : mr3.cu
	$(NVCC) $(NVCCFLAGS) $(NVCCFLAGSXCOMP) $(CUDAINC) $< -o $@ 

clean:
	-rm -f *.o 
	-rm -f $(TARGET)
