CXX=icc
CXXFLAGS=-g -Wall -W -Werror -pedantic

all: baseline vectorized

baseline:
	$(CXX) $(CXXFLAGS) -O2 -no-vec mm.c -o mmnovec

vectorized: 
	$(CXX) $(CXXFLAGS) -O2 -qopt-report=2 -qopt-report-phase=vec -D ALIGNED -DNOALIAS mm.c -o mmvec

clean:
	rm -f *.o *.optrpt mmnovec mmvec
