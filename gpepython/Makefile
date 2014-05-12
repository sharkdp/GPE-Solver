# TODO: -ffast-math speeds up the program alot but can it be trusted?

## Enable debug macro ##
#DEBUG = -DDEBUG_GPE
DEBUG =

## Generate profiling information for gprof ##
#PROFILER = 	-pg
PROFILER =

## Disable to run on single core ##
OPENMP = 	-fopenmp

## Enable C++0x features ##
CPP0X=--std=c++0x
#CPP0X=

## Use latest version of GCC ##
CXX = g++

CXXFLAGS =	-O3 $(CPP0X) -g -Wall -fPIC -I"/usr/include/python2.7" -fmessage-length=0 $(OPENMP) $(PROFILER) $(DEBUG)

OBJS =	m_rand.o GPEPython.o SimulationParameters.o Field.o Array.o Wavefunction.o PotentialDipolar.o PotentialKinetic.o PotentialHarmonic.o PotentialLattice.o WavefunctionGaussian.o Simulation.o

LIBS =		-lfftw3_threads -lfftw3 -lm $(OPENMP) $(PROFILER) -Wl,-h -Wl,gpepython.so -shared -Wl,--start-group -Wl,-Bstatic  -Wl,-Bdynamic -lboost_python -lutil -lpthread -ldl -Wl,--end-group

TARGET =	gpepython.so

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET) *.pyc

# TODO: add header files as dependencies!
