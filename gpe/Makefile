## Python 2 or Python 3 ##
PY_HEADER = "/usr/include/python3.4m"
PY_BOOST = boost_python3
# PY_HEADER = "/usr/include/python2.7"
# PY_BOOST = boost_python

## Enable debug macro ##
# DEBUG = -DDEBUG_GPE
DEBUG =

## Generate profiling information for gprof ##
# PROFILER = 	-pg
PROFILER =

## Disable to run on single core ##
OPENMP = 	-fopenmp

## Enable C++0x features ##
CPP0X=--std=c++0x

## Use latest version of GCC ##
CXX = g++

CXXFLAGS =	-O3 $(CPP0X) -g -Wall -fPIC -I$(PY_HEADER) -fmessage-length=0 $(OPENMP) $(PROFILER) $(DEBUG)

OBJS =	m_rand.o GPEPython.o SimulationParameters.o Field.o Array.o Wavefunction.o PotentialDipolar.o PotentialKinetic.o PotentialHarmonic.o PotentialLattice.o WavefunctionGaussian.o Simulation.o

LIBS =		-lfftw3_threads -lfftw3 -lm $(OPENMP) $(PROFILER) -Wl,-h -Wl,gpepython.so -shared -Wl,--start-group -Wl,-Bstatic  -Wl,-Bdynamic -l$(PY_BOOST) -lutil -lpthread -ldl -Wl,--end-group

TARGET =	gpepython.so

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -rf $(OBJS) $(TARGET) *.pyc __pycache__

# TODO: add header files as dependencies!
