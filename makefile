cc = g++
MPICC = mpicxx
OMP=-fopenmp -DOMP
CFLAGS=-g -O3 -Wall
LDFLAGS= -lm

SDIR = src/sequential
PDIR = src/parallel
UDIR = src/util
RDIR = release

all: srun prun util

srun: $(RDIR)/kmeans $(RDIR)/colorBlock $(RDIR)/smooth

util: $(RDIR)/compare_image

prun: $(RDIR)/kmeans_mpi $(RDIR)/colorBlock_mpi $(RDIR)/kmeans_cuda $(RDIR)/colorBlock_omp

$(RDIR)/compare_image: $(UDIR)/compare_image.cpp
	$(cc) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(RDIR)/kmeans: $(SDIR)/kmeans.cpp $(UDIR)/cycletimer.cpp
	$(cc) $(CFLAGS) -o $@ $^ $(LDFLAGS)
$(RDIR)/colorBlock: $(SDIR)/colorBlock.cpp $(UDIR)/cycletimer.cpp
	$(cc) $(CFLAGS) -o $@ $^ $(LDFLAGS)
$(RDIR)/smooth: $(SDIR)/smooth.cpp $(UDIR)/cycletimer.cpp
	$(cc) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(RDIR)/kmeans_mpi: $(PDIR)/kmeans_mpi.cpp $(UDIR)/cycletimer.cpp
	$(MPICC) $(CFLAGS) $(OMP) -o $@ $^ $(LDFLAGS)
$(RDIR)/kmeans_cuda: $(PDIR)/kmeans.cu $(UDIR)/cycletimer.cpp
	nvcc -O3 -m64 --gpu-architecture compute_61 -o $@ $^
$(RDIR)/colorBlock_mpi: $(PDIR)/colorBlock_mpi.cpp $(UDIR)/cycletimer.cpp
	$(MPICC) $(CFLAGS) $(OMP) -o $@ $^ $(LDFLAGS)
$(RDIR)/colorBlock_omp: $(PDIR)/colorBlock_omp.cpp $(UDIR)/cycletimer.cpp
	$(MPICC) $(CFLAGS) $(OMP) -o $@ $^ $(LDFLAGS)

.PHONY: clean
clean:
	rm -f release/* image/*.bmp