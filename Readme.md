# Dipolar Gross-Pitaevskii Solver

`gpepython` is a [split-step](http://en.wikipedia.org/wiki/Split-step_method) solver for the (dipolar) [Gross-Pitaevskii equation](http://en.wikipedia.org/wiki/Gross%E2%80%93Pitaevskii_equation). The main part
of the code is written in C++ and parallelized using [OpenMP](http://en.wikipedia.org/wiki/OpenMP). The simulations can be controlled
through a Python interface.

![Dipolar Bose gas in a triple-well potential](https://github.com/sharkdp/GPE-Solver/raw/master/triplewell.png)

## Sample simulation
``` python
# simple.py

from gpe import HarmonicTrapGPE, GPE

# Load GPE
sim = HarmonicTrapGPE("Simple", 1.0, 1.0, 1.0)

# Grid parameters
sim.setGridSize(32)

# Contact interaction
# Only the product (atomNumber - 1) * scatteringLength is important
sim.set("atomNumber",       2)
sim.set("scatteringLength", 1)

# Dipolar interaction (switched off)
sim.set("dipolarLength", 0)
sim.set("evolutionDipolar", False)

# Number of threads. Equals number of cores if not set
# sim.setThreads(4)

sim.initialize()

# Imaginary time evolution
sim.ite(steps=10000, monitorSteps=500, plotSteps=500)
```

## Installation
Dependencies needed:
- Recent version of g++ with C++0x and OpenMP support
- Multi-threaded [fftw3](https://github.com/FFTW/fftw3)
- [Boost](http://www.boost.org/) with support for Python
- Works with both Python 2 and 3 (see Makefile, tested with 2.7 and 3.4)

Run `make` in the `gpe` folder.

## Documentation
See appendix C of [this thesis](https://github.com/sharkdp/GPE-Solver/raw/master/doc/thesis.pdf).
