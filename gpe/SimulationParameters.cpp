#include "SimulationParameters.h"
#include <iostream>
#include <sstream>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>

SimulationParameters::SimulationParameters() {
    _evolutionType = -1;
}

void SimulationParameters::set(string name, double value) {
    _doubleParameters[name] = value;

    if (name == "sizeX" || name == "sizeY" || name == "sizeZ" ||
            name == "maxX" || name == "maxY" || name == "maxZ") {
        try {
            computations();
        }
        catch (ParameterNotFound & e) {
            // do nothing (not all parameters for the computation are set yet)
        }
    }
}

double SimulationParameters::getDouble(string name) const {
    map<string, double>::const_iterator value = _doubleParameters.find(name);
    if (value == _doubleParameters.end()) {
        throw ParameterNotFound(name);
    }
    return value->second;
}

double SimulationParameters::getDouble(string name, double defaultValue) const {
    map<string, double>::const_iterator value = _doubleParameters.find(name);
    if (value == _doubleParameters.end()) {
        return defaultValue;
    }
    return value->second;
}

bool SimulationParameters::getBool(string name) const {
    map<string, bool>::const_iterator value = _boolParameters.find(name);
    if (value == _boolParameters.end()) {
        throw ParameterNotFound(name);
    }
    return value->second;
}

string SimulationParameters::getString(string name) const {
    map<string, string>::const_iterator value = _stringParameters.find(name);
    if (value == _stringParameters.end()) {
        throw ParameterNotFound(name);
    }
    return value->second;
}

string SimulationParameters::getAsString(string name) const {
    stringstream sstr;
    map<string, double>::const_iterator valueD = _doubleParameters.find(name);

    if (valueD == _doubleParameters.end()) {
       map<string, bool>::const_iterator valueB = _boolParameters.find(name);
       if (valueB == _boolParameters.end()) {
           throw ParameterNotFound(name);
       }

       sstr << boolalpha << valueB->second;
    } else {
        sstr << valueD->second;
    }
    return sstr.str();
}

void SimulationParameters::computations() {
    // Calculates some parameters which are needed for fast access
    _size_x = getDouble("sizeX");
    _size_y = getDouble("sizeY");
    _size_z = getDouble("sizeZ");

    _max_x = getDouble("maxX");
    _max_y = getDouble("maxY");
    _max_z = getDouble("maxZ");

    if (_size_x <= 0 || _size_y <= 0 || _size_z <= 0) {
        cerr << "Grid size must not be smaller or equal to zero" << endl;
        exit(1);
    }

    _delta_x = 2.0 * _max_x / (_size_x - 1);
    _delta_y = 2.0 * _max_y / (_size_y - 1);
    _delta_z = 2.0 * _max_z / (_size_z - 1);

    _delta_k_x = (2.0 * M_PI) / (_size_x * _delta_x);
    _delta_k_y = (2.0 * M_PI) / (_size_y * _delta_y);
    _delta_k_z = (2.0 * M_PI) / (_size_z * _delta_z);

    _max_k_x = M_PI / _delta_x;
    _max_k_y = M_PI / _delta_y;
    _max_k_z = M_PI / _delta_z;
}

int SimulationParameters::size(int direction) const {
    if (direction == 0) {
        return _size_x;
    } else if (direction == 1) {
        return _size_y;
    }
    return _size_z;
}

double SimulationParameters::coord(int direction, int index, bool realSpace) const {
    if (direction == 0) {
        return realSpace? x(index) : k_x(index);
    } else if (direction == 1) {
        return realSpace ? y(index) : k_y(index);
    }
    return realSpace ? z(index) : k_z(index);
}

double SimulationParameters::max(int direction, bool realSpace) const {
    if (direction == 0) {
        return realSpace ? _max_x : _max_k_x;
    } else if (direction == 1) {
        return realSpace ? _max_y : _max_k_y;
    }
    return realSpace ? _max_z : _max_k_z;
}

double SimulationParameters::d(int direction, bool realSpace) const {
    if (direction == 0) {
        return realSpace ? _delta_x : _delta_k_x;
    } else if (direction == 1) {
        return realSpace ? _delta_y : _delta_k_y;
    }
    return realSpace ? _delta_z : _delta_k_z;
}
