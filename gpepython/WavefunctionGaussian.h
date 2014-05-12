#ifndef WAVEFUNCTIONGAUSSIAN_H_
#define WAVEFUNCTIONGAUSSIAN_H_

#define _USE_MATH_DEFINES
#include <cmath>
#include "complexType.h"
#include <iostream>
#include "Wavefunction.h"

class WavefunctionGaussian: public Wavefunction {
 public:
    WavefunctionGaussian(const SimulationParameters & sp) : Wavefunction(sp) {
        _sx = _sp.getDouble("sigmaX");
        _sy = _sp.getDouble("sigmaY");
        _sz = _sp.getDouble("sigmaZ");

        _offsetX = _sp.getDouble("offsetX", 0.0);
        _offsetY = _sp.getDouble("offsetY", 0.0);
        _offsetZ = _sp.getDouble("offsetZ", 0.0);

        _factor = pow(M_PI, -0.75) * pow(_sx * _sy * _sz, -0.5);
    }

 private:
    complex_t initialValue(int i, int j, int k);

    double _factor, _sx, _sy, _sz;
    double _offsetX, _offsetY, _offsetZ;
};

#endif /* WAVEFUNCTIONGAUSSIAN_H_ */
