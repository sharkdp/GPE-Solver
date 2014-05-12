#ifndef WAVEFUNCTION_H_
#define WAVEFUNCTION_H_

#include "complexType.h"
#include "Array.h"
#include "m_rand.h"

class Wavefunction: public Array {
    friend class Array;

 public:
    Wavefunction(const SimulationParameters & sp) : Array(sp) {
        _space = Array::SPACE_POSITION;
    }

    void normalize();
    void evolution(const Array & f, complex_t factor) const;

    void perturb(double amplitude);
    void reduce(int factor, int direct);

    void applyRotation(double n);

    void abs2();
    void abs2Copy(Array & a) const;
};

#endif /* WAVEFUNCTION_H_ */
