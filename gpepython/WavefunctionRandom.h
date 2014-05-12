
#ifndef WAVEFUNCTIONRANDOM_H_
#define WAVEFUNCTIONRANDOM_H_

#define _USE_MATH_DEFINES
#include <cmath>
#include "complexType.h"
#include <iostream>
#include "Wavefunction.h"

class WavefunctionRandom: public Wavefunction {
 public:
    WavefunctionRandom(const SimulationParameters & sp) : Wavefunction(sp) {
        _sx = _sp.getDouble("sigmaX");
        _sy = _sp.getDouble("sigmaY");
        _sz = _sp.getDouble("sigmaZ");

        _factor = pow(M_PI, -0.75) * pow(_sx * _sy * _sz, -0.5);
    }

 private:
    complex_t initialValue(int i, int j, int k);

    double  _sx, _sy, _sz,  _factor;
};

#endif /* WAVEFUNCTIONRANDOM_H_ */
