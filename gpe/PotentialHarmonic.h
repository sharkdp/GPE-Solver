
#ifndef POTENTIALHARMONIC_H_
#define POTENTIALHARMONIC_H_

#include "complexType.h"
#include <cmath>
#include "Array.h"

class PotentialHarmonic: public Array {
 public:
    PotentialHarmonic(const SimulationParameters & sp) :
        Array(sp) {
        double omegaZ = sp.getDouble("omegaZ");

        _c_x = pow(sp.getDouble("omegaX") / omegaZ, 2);
        _c_y = pow(sp.getDouble("omegaY") / omegaZ, 2);
        _c_z = 1.0;

        _space = Array::SPACE_POSITION;
    }

 private:
    complex_t initialValue(int i, int j, int k);

    double _c_x, _c_y, _c_z;
};

#endif /* POTENTIALHARMONIC_H_ */
