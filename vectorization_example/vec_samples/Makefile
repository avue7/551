CXX=icc
INCLUDES=-I src
CXXFLAGS=-c -g -O2 -std=c99 -D NOFUNCALL -qopt-report=1 -qopt-report-phase=vec $(INCLUDES) $(BITS)

PROGS=MatVector

.PHONY: default
default: Driver.o Multiply.o
	$(CXX) Driver.o Multiply.o -o vec_samples $(LIBS)

.PHONY: baseline
baseline: default

baseline: CXXFLAGS = -c -g -O2 -std=c99 $(INCLUDES) $(BITS)
baseline: Driver.o Multiply.o
	$(CXX) Driver.o Multiply.o -o vec_samples $(LIBS)

xhost: CXXFLAGS = -c -g -xHost -std=c99 -DNOALIAS $(INCLUDES) $(BITS)
xhost: Driver.o Multiply.o
	$(CXX) Driver.o Multiply.o -o vec_samples $(LIBS)

noalias: CXXFLAGS = -c -g -O2 -std=c99 -DNOALIAS $(INCLUDES) $(BITS)
noalias: Driver.o Multiply.o
	$(CXX) Driver.o Multiply.o -o vec_samples $(LIBS)

reduction: CXXFLAGS = -c -g -xHost -std=c99 -qopenmp -DNOALIAS -DREDUCTION $(INCLUDES) $(BITS)
reduction: Driver.o Multiply.o
	$(CXX) -qopenmp Driver.o Multiply.o -o vec_samples $(LIBS)

align: CXXFLAGS = -c -g -xHost -std=c99 -qopenmp -DNOALIAS -D ALIGNED -DREDUCTION $(INCLUDES) $(BITS)
align: Driver.o Multiply.o
	$(CXX) -qopenmp Driver.o Multiply.o -o vec_samples $(LIBS)

nofunc: CXXFLAGS = -c -g -xHost -std=c99 -qopenmp -D NOFUNCCALL -DNOALIAS -DREDUCTION $(INCLUDES) $(BITS)
nofunc: Driver.o Multiply.o
	$(CXX) -qopenmp Driver.o Multiply.o -o vec_samples $(LIBS)

Driver.o: src/Driver.c
	$(CXX) $(BITS) $(CXXFLAGS) -o $@ src/Driver.c

Multiply.o: src/Multiply.c
	$(CXX) $(BITS) $(CXXFLAGS) -o $@ src/Multiply.c

.PHONY: clean
clean::
	-rm -f $(PROGS) *.o *.out *.optrpt
