#ifndef POTENTIALKINETIC_H_
#define POTENTIALKINETIC_H_

#include "Array.h"

class PotentialKinetic: public Array {
 public:
    PotentialKinetic(const SimulationParameters & sp) :
        Array(sp) {
        _space = Array::SPACE_MOMENTUM;
    }

 private:
    complex_t initialValue(int i, int j, int k);
};

#endif /* POTENTIALKINETIC_H_ */
