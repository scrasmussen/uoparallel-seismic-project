NVCC ?= nvcc
NVFLAGS =
DOCS = ../docs
DATA = ../data
RM ?= rm

all: cudasweep

cudasweep: cudasweep-tt-multistart.cu
	$(NVCC) -o $@ $(NVFLAGS) $< -lm

cudarun: cudasweep
	./cudasweep $(DOCS)/velocity-241-241-51.txt $(DOCS)/818-FS.txt $(DOCS)/start-4-241-241-51.txt
	
clean:
	$(RM) -f *.o cudasweep
