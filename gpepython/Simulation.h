#ifndef SIMULATION_H_
#define SIMULATION_H_

#include "Wavefunction.h"
#include "PotentialKinetic.h"
#include "PotentialLattice.h"
#include "PotentialHarmonic.h"
#include "PotentialDipolar.h"

class Simulation {
 public:
    Simulation(SimulationParameters & sp) : _sp(sp) {
        _simulationInitialized = false;

        _wavefunction = 0;
        _density = 0;
        _potentialKinetic = 0;
        _potentialExternal = 0;
        _potentialDipolar = 0;
    }

    void setWavefunction(Wavefunction * wavefunction) {
        if (_wavefunction) {
            delete _wavefunction;
        }
        _wavefunction = wavefunction;
    }

    Wavefunction * getWavefunction() const {
        return _wavefunction;
    }

    Wavefunction * getDensity() const {
        return _density;
    }

    void initialize(bool deleteWavefunction = true);

    void evolution(int steps);

    double energyKinetic();  // cannot be declared as const because of Fourier transform
    double energyPotential() const;
    double energyContact() const;
    double energyDipolar();  // cannot be declared as const because of Fourier transform

    double energyTotal() {
        // This function is normally not used (by the python interface) so
        // we don't have to care about optimizing it (store energies..)
        return energyContact() + energyPotential() + energyKinetic() + energyDipolar();
    }

    void writeDensity(const string filename) const;
    void writePotential(const string filename) const;

    virtual ~Simulation();

    class SimulationError : public exception {
     public:
        SimulationError(string description) : _description(description) {}

        virtual const char* what() const throw() {
            return ("Simulation error: " + _description).c_str();
        }

        virtual ~SimulationError() throw () {}
     private:
        string _description;
    };

 private:
    bool _simulationInitialized;

    SimulationParameters & _sp;
    int _evolutionType;

    Wavefunction * _wavefunction;
    Wavefunction * _density;

    PotentialKinetic * _potentialKinetic;
    Array * _potentialExternal;
    PotentialDipolar * _potentialDipolar;
};

#endif /* SIMULATION_H_ */
