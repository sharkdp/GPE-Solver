#ifndef POTENTIALDIPOLAR_H_
#define POTENTIALDIPOLAR_H_

#include "complexType.h"
#include <cmath>
#include "Array.h"

class PotentialDipolar: public Array {
 public:
    PotentialDipolar(const SimulationParameters & sp) :
        Array(sp) {
        double factorFFTW = pow(_sp.size(), -1);
        _gdd = _sp.getDipolarInteractionFactor() * factorFFTW;

        _space = Array::SPACE_MOMENTUM;
    }

 private:
    complex_t initialValue(int i, int j, int k);

    double _gdd;
};

#endif /* POTENTIALDIPOLAR_H_ */
