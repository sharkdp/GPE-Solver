#include <boost/python.hpp>
#include <string>
#include <iostream>
#include "Simulation.h"
#include "omp.h"

using namespace std;

class GPEPython {
 public:
    GPEPython() {
        _simulation = new Simulation(_sp);
    }

    virtual ~GPEPython() {
        delete _simulation;
    }

    void set(string name, double value) {
        _sp.set(name, value);
    }

    void setBool(string name, bool value) {
        _sp.set(name, value);
    }

    void setString(string name, string value) {
        _sp.set(name, value);
    }

    double get(string name) {
        return _sp.getDouble(name);
    }

    bool getBool(string name) {
        return _sp.getBool(name);
    }

    string getString(string name) {
        return _sp.getString(name);
    }

    void initialize() {
        _simulation->initialize(true);
    }

    void initializeDW(bool deleteWavefunction) {
        _simulation->initialize(deleteWavefunction);
    }

    void ite(int steps) {
        _sp.setEvolutionType(SimulationParameters::IMAGINARY);
        _simulation->evolution(steps);
    }

    void rte(int steps) {
        _sp.setEvolutionType(SimulationParameters::REAL);
        _simulation->evolution(steps);
    }

    double energyKinetic() {
        return _simulation->energyKinetic();
    }

    double energyPotential() {
        return _simulation->energyPotential();
    }

    double energyContact() {
        return _simulation->energyContact();
    }

    double energyDipolar() {
        return _simulation->energyDipolar();
    }

    void write(string filename) {
        _simulation->writeDensity(filename);
    }

    void writePotential(string filename) {
        _simulation->writePotential(filename);
    }

    double evaluate(int quantity) {
        double value;

        if (quantity == DISPERSION_X || quantity == DISPERSION_Y || quantity == DISPERSION_Z) {
            function<complex_t (double, double, double)> func;

            if (quantity == DISPERSION_X) {
                func = [] (double x, double y, double z) { return x * x; };
            } else if (quantity == DISPERSION_Y) {
                func = [] (double x, double y, double z) { return y * y; };
            } else if (quantity == DISPERSION_Z) {
                func = [] (double x, double y, double z) { return z * z; };
            }

            value = _simulation->getWavefunction()->integral(func);
            value = sqrt(2.0 * value);
        } else if (quantity == PSI_POWER_2) {
            // Integral over |Psi|^2 ("normalization integral")
            auto one = [] (double, double, double) { return 1; };
            value = _simulation->getWavefunction()->integral(one);
        } else if (quantity == PSI_POWER_6) {
            // Integral over |Psi|^6
            Wavefunction * density = _simulation->getDensity();
            Wavefunction * wavefunction = _simulation->getWavefunction();
            wavefunction->abs2Copy(*density);

            value = density->integral(*density);
        } else if (quantity == PSI_0RE) {
            // Value of the wave function in the middle of the trap
            Wavefunction * wavefunction = _simulation->getWavefunction();

            value = real( (*wavefunction)(_sp.sizeX()/2+4, _sp.sizeY()/2+3, _sp.sizeZ()/2+2) );
        } else if (quantity == PSI_0IM) {
            // Value of the wave function in the middle of the trap
            Wavefunction * wavefunction = _simulation->getWavefunction();

            value = imag( (*wavefunction)(_sp.sizeX()/2+4, _sp.sizeY()/2+3, _sp.sizeZ()/2+2) );
        } else {
            cerr << "Unknown quantity" << endl;
            value = 0.0;
        }

        return value;
    }

    void setThreads(int nthreads) {
#ifdef _OPENMP
        omp_set_num_threads(nthreads);
#endif
    }

    void perturbation(double amp) {
        (_simulation->getWavefunction())->perturb(amp);
    }

    void reduce(int factor, int direct) {
        (_simulation->getWavefunction())->reduce(factor, direct);
    }

    void applyRotation(double n) {
        (_simulation->getWavefunction())->applyRotation(n);
    }

    double computeDispersion(int direction, bool realSpace) {
        return (_simulation->getWavefunction())->dispersion(direction, realSpace);
    }

    double xk(int direction) {
        return (_simulation->getWavefunction())->averageXK(direction);
    }

    static const int DISPERSION_X = 1;
    static const int DISPERSION_Y = 2;
    static const int DISPERSION_Z = 3;
    static const int PSI_POWER_2  = 4;
    static const int PSI_POWER_6  = 5;
    static const int PSI_0RE      = 6;
    static const int PSI_0IM      = 7;

 private:
    SimulationParameters _sp;
    Simulation * _simulation;
};

BOOST_PYTHON_MODULE(gpepython)
{
    using namespace boost::python;

    class_<GPEPython>("GPEPython")
        .def("set", &GPEPython::set)
        .def("setBool", &GPEPython::setBool)
        .def("setString", &GPEPython::setString)
        .def("get", &GPEPython::get)
        .def("getBool", &GPEPython::getBool)
        .def("getString", &GPEPython::getString)
        .def("ite", &GPEPython::ite)
        .def("rte", &GPEPython::rte)
        .def("initialize", &GPEPython::initialize)
        .def("initialize", &GPEPython::initializeDW)
        .def("energyKinetic", &GPEPython::energyKinetic)
        .def("energyPotential", &GPEPython::energyPotential)
        .def("energyContact", &GPEPython::energyContact)
        .def("energyDipolar", &GPEPython::energyDipolar)
        .def("write", &GPEPython::write)
        .def("writePotential", &GPEPython::writePotential)
        .def("evaluate", &GPEPython::evaluate)
        .def("setThreads", &GPEPython::setThreads)
        .def("perturbation", &GPEPython::perturbation)
        .def("applyRotation", &GPEPython::applyRotation)
        .def("reduce",  &GPEPython::reduce)
        .def("computeDispersion",  &GPEPython::computeDispersion)
        .def("computeXK",  &GPEPython::xk);
}
