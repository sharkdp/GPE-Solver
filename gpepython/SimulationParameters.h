#ifndef SIMULATIONPARAMETERS_H_
#define SIMULATIONPARAMETERS_H_

#define _USE_MATH_DEFINES
#include <cmath>
#include <map>
#include <string>
#include <exception>
#include "complexType.h"

using namespace std;

class SimulationParameters {
 public:
    SimulationParameters();

    void setEvolutionType(int evolutionType) { _evolutionType = evolutionType; }
    int evolutionType() { return _evolutionType; }

    void set(string name, double value);

    void set(string name, bool value) {
        _boolParameters[name] = value;
    }

    void set(string name, string value) {
        _stringParameters[name] = value;
    }

    double getDouble(string name) const;
    double getDouble(string name, double defaultValue) const;
    bool getBool(string name) const;
    string getString(string name) const;

    string getAsString(string name) const;

    double maxX() const { return _max_x; }
    double maxY() const { return _max_y; }
    double maxZ() const { return _max_z; }

    double maxKX() const { return _max_k_x; }
    double maxKY() const { return _max_k_y; }
    double maxKZ() const { return _max_k_z; }

    int sizeX() const { return _size_x; }
    int sizeY() const { return _size_y; }
    int sizeZ() const { return _size_z; }

    int size(int direction) const;
    double coord(int direction, int index, bool realSpace) const;
    double max(int direction, bool realSpace) const;
    double d(int direction, bool realSpace) const;

    int size() const { return _size_x * _size_y * _size_z; }

    double x(int i) const { return -_max_x + i * _delta_x; }
    double y(int j) const { return -_max_y + j * _delta_y; }
    double z(int k) const { return -_max_z + k * _delta_z; }

    double k_x(int i) const {
        if (i < _size_x/2)
            return i * _delta_k_x;
        else
            return -_max_k_x + (i - _size_x/2) * _delta_k_x;
    }

    double k_y(int j) const {
        if (j < _size_y/2)
            return j * _delta_k_y;
        else
            return -_max_k_y + (j - _size_y/2) * _delta_k_y;
    }

    double k_z(int k) const {
        if (k < _size_z/2)
            return k * _delta_k_z;
        else
            return -_max_k_z + (k - _size_z/2) * _delta_k_z;
    }

    double dX() const { return _delta_x; }
    double dY() const { return _delta_y; }
    double dZ() const { return _delta_z; }

    double dkX() const { return _delta_k_x; }
    double dkY() const { return _delta_k_y; }
    double dkZ() const { return _delta_k_z; }

    static const int REAL = 0;
    static const int IMAGINARY = 1;

    class ParameterNotFound : public exception {
     public:
        ParameterNotFound(string name) : _name(name) {}

        virtual const char* what() const throw() {
            return ("Parameter '" + _name + "' not found").c_str();
        }

        virtual ~ParameterNotFound() throw () {}
     private:
        string _name;
    };

    double getContactInteractionFactor() const {
        // factor = g * (N - 1) = 4 pi a * (N - 1)
        return 4 * M_PI * getDouble("scatteringLength") * (getDouble("atomNumber") - 1);
    }

    double getDipolarInteractionFactor() const {
        // factor = g_dd * (N - 1) = 4 pi a_dd * (N - 1)
        return 4 * M_PI * getDouble("dipolarLength") * (getDouble("atomNumber") - 1);
    }

    complex_t getThreeBodyLossesFactor() const {
        // factor = - i L3/2
        return - getDouble("lossRate") / 2.0 * COMPLEX_I;
    }

    double getModulationFreq() const {
        return getDouble("modulationFreq");
    }

    double getModulationAmplitude() const {
        return getDouble("modulationAmplitude");
    }

 private:
    map<string, double> _doubleParameters;
    map<string, bool> _boolParameters;
    map<string, string> _stringParameters;

    int _evolutionType;

    int _size_x, _size_y, _size_z;

    double _max_x, _max_y, _max_z;
    double _delta_x, _delta_y, _delta_z;
    double _max_k_x, _max_k_y, _max_k_z;
    double _delta_k_x, _delta_k_y, _delta_k_z;

    void computations();
};

#endif /* SIMULATIONPARAMETERS_H_ */
