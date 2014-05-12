#include "Simulation.h"
#include "WavefunctionGaussian.h"
#include "WavefunctionRandom.h"
#include "complexType.h"
#include <fstream>

void Simulation::initialize(bool deleteWavefunction) {
    if (_sp.size() < 1) {
        throw SimulationError("Grid is not initialized");
    }

    if (_wavefunction && deleteWavefunction) {
        delete _wavefunction;
        _wavefunction = 0;
    }

    if (_density) {
        delete _density;
        _density = 0;
    }

    if (_potentialKinetic) {
        delete _potentialKinetic;
        _potentialKinetic = 0;
    }

    if (_potentialExternal) {
        delete _potentialExternal;
        _potentialExternal = 0;
    }

    if (_potentialDipolar) {
        delete _potentialDipolar;
        _potentialDipolar = 0;
    }

    if (deleteWavefunction) {
        // Initialize wavefunction
//        setWavefunction(new WavefunctionRandom(_sp));
         setWavefunction(new WavefunctionGaussian(_sp));
        _wavefunction->initFourierPlans();
        _wavefunction->initialize();
    } else {
        if (_wavefunction == 0) {
            throw SimulationError("Wavefunction is not initialized (but deleteWavefunction = false)");
        }
    }

    // Initialize the density array
    _density = new Wavefunction(_sp);
    // the following is only needed if the dipolar interactions are on, but
    // it should be fast (same fftw plan size as wavefunction)
    _density->initFourierPlans();

    // Initialize the different "potentials"
    _potentialKinetic = new PotentialKinetic(_sp);
    _potentialKinetic->initialize();

    string namePotential = _sp.getString("externalPotential");
    if (namePotential == "harmonic" || namePotential == "modulatedHarmonic") {
        _potentialExternal = new PotentialHarmonic(_sp);
    } else if (namePotential == "lattice") {
        _potentialExternal = new PotentialLattice(_sp);
    } else {
        throw SimulationError("External potential '" + namePotential + "' not known.");
    }
    _potentialExternal->initialize();

    _potentialDipolar = new PotentialDipolar(_sp);
    _potentialDipolar->initialize();

    _simulationInitialized = true;
}

void Simulation::evolution(int steps) {
    if (!_simulationInitialized) {
        throw SimulationError("Simulation has not been initialized");
    }

    // Read parameters
    double dt = _sp.getDouble("timeStep");

    bool evolutionKinetic = _sp.getBool("evolutionKinetic");
    bool evolutionPotential = _sp.getBool("evolutionPotential");
    bool evolutionContact = _sp.getBool("evolutionContact");
    bool evolutionDipolar = _sp.getBool("evolutionDipolar");
    bool evolutionThreeBodyLosses = _sp.getBool("evolutionThreeBodyLosses");
    bool evolutionModulation = (_sp.getString("externalPotential") == "modulatedHarmonic");

    double contactInteractionFactor = 0;
    if (evolutionContact) {
        contactInteractionFactor = _sp.getContactInteractionFactor();
    }

    complex_t threeBodyLossesFactor;
    if (evolutionThreeBodyLosses) {
        threeBodyLossesFactor = _sp.getThreeBodyLossesFactor();
    }

    double omMod = 0.0;
    double ampMod= 0.0;
    if (evolutionModulation) {
        omMod  = _sp.getModulationFreq();
        ampMod = _sp.getModulationAmplitude();
    }

    int evolutionType = _sp.evolutionType();
    complex_t evolutionCoefficient;

    switch (evolutionType) {
    case SimulationParameters::REAL:
        evolutionCoefficient = - 1.0 * dt * COMPLEX_I;
        break;
    case SimulationParameters::IMAGINARY:
        evolutionCoefficient = - 1.0 * dt;

        if (evolutionThreeBodyLosses) {
            throw SimulationError("Imaginary time evolution with three-body losses enabled.");
        }
        break;
    default:
        cerr << "Evolution type not set" << endl;
        return;
    }

    double one_over_size = 1.0 / _sp.size();

    double t = 0.0;
    for (int step = 0; step < steps; step++) {
        t = step * dt;

        if (evolutionContact || evolutionDipolar || evolutionThreeBodyLosses) {
            _wavefunction->abs2Copy(*_density);
        }

        if (evolutionContact) {
            // Evolution via contact interaction
            _wavefunction->evolution(*_density, contactInteractionFactor * evolutionCoefficient);
        }

        if (evolutionDipolar) {
            // Evolution via Dipolar Potential
            _density->fourierTransform();
            (*_density) *= (*_potentialDipolar);  // multiply density and dipolar potential in Fourier space
            _density->fourierTransformInverse();
            _wavefunction->evolution(*_density, evolutionCoefficient);
        }

        if (evolutionThreeBodyLosses) {
            if (evolutionDipolar) {
                _wavefunction->abs2Copy(*_density);  // re-create the density
            }
            _density->abs2();  // density^2
            _wavefunction->evolution(*_density, threeBodyLossesFactor * evolutionCoefficient);
        }

        if (evolutionKinetic) {
            // Evolution via kinetic energy
            _wavefunction->fourierTransform();
            _wavefunction->evolution(*_potentialKinetic, evolutionCoefficient);
            _wavefunction->fourierTransformInverse();
        }

        if (evolutionPotential) {
            if (evolutionModulation) {
                // Modulation via modulated Harmonic Potential
                _wavefunction->evolution(*_potentialExternal, evolutionCoefficient * pow((1.0 + ampMod * sin(omMod * t)), 2.0));
            } else {
                // Evolution via potential energy
                _wavefunction->evolution(*_potentialExternal, evolutionCoefficient);
            }
        }

        if (evolutionType == SimulationParameters::IMAGINARY) {
            // Normalization of the wavefunction
            _wavefunction->normalize();
        } else {
            // only normalization due to Fourier transform is needed
            (*_wavefunction) *= one_over_size;
        }
    }
}

double Simulation::energyKinetic() {
    if (!_sp.getBool("evolutionKinetic"))
        return 0;

    _wavefunction->fourierTransform();
    double energyKinetic = _wavefunction->integral(*_potentialKinetic) / _sp.size();
    _wavefunction->fourierTransformInverse();
    *_wavefunction *= (1.0 / _sp.size());  // normalize after Fourier transform

    return energyKinetic;
}

double Simulation::energyPotential() const {
    if (!_sp.getBool("evolutionPotential"))
        return 0;

    return _wavefunction->integral(*_potentialExternal);
}

double Simulation::energyContact() const {
    if (!_sp.getBool("evolutionContact"))
        return 0;

    _wavefunction->abs2Copy(*_density);

    // factor 1/2 is due to the additional 1/2 in the energy functional:
    // E_contact = g/2 * (N-1) * integral |Psi|^4
    return 0.5 * _sp.getContactInteractionFactor() * _wavefunction->integral(*_density);
}

double Simulation::energyDipolar() {
    if (!_sp.getBool("evolutionDipolar"))
        return 0;

    _wavefunction->abs2Copy(*_density);
    _density->fourierTransform();

    return 0.5 * _density->integral(*_potentialDipolar);
}

void Simulation::writeDensity(const string filename) const {
    _wavefunction->abs2Copy(*_density);

    _density->cut2D(filename + "_cut_yz.data", 0);
    _density->cut2D(filename + "_cut_xz.data", 1);
    _density->cut2D(filename + "_cut_xy.data", 2);

    _density->cut1D(filename + "_cut_x.data", 0);
    _density->cut1D(filename + "_cut_y.data", 1);
    _density->cut1D(filename + "_cut_z.data", 2);

    _density->project2D(filename + "_project_yz.data");

    _density->project1D(filename + "_project_x.data", 0);
    _density->project1D(filename + "_project_y.data", 1);
    _density->project1D(filename + "_project_z.data", 2);

//    _density->projectFourierMatlab2D(filename + "_proj_fourier_matlab_yz.data");
    _density->projectMatlab2D(filename + "_proj_matlab_xz.m");
//    _density->cutFourier2D(filename + "_cut_fourier_yz.data", 0);
}

void Simulation::writePotential(const string filename) const {
    _potentialExternal->cut2D(filename + "_cut_yz.data", 0);
    _potentialExternal->cut2D(filename + "_cut_xy.data", 2);

    _potentialExternal->cut1D(filename + "_cut_x.data", 0);
    _potentialExternal->cut1D(filename + "_cut_y.data", 1);
    _potentialExternal->cut1D(filename + "_cut_z.data", 2);
}

Simulation::~Simulation() {
    delete _wavefunction;
    delete _density;
    delete _potentialKinetic;
    delete _potentialExternal;
    delete _potentialDipolar;
}
