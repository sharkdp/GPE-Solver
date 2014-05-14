#include "PotentialDipolar.h"
#include <limits>

complex_t PotentialDipolar::initialValue(int i, int j, int k) {
    double kz2 = _sp.k_z(k) * _sp.k_z(k);
    double kTot2 = _sp.k_x(i) * _sp.k_x(i) + _sp.k_y(j) * _sp.k_y(j) + kz2;

    if (abs(kTot2) < numeric_limits<double>::epsilon()) {
        return 0;
    }

    // TODO: move _gdd factor to the evolution part in Simulation.cpp
    // such that initialize() doesn't have to be called after changing
    // the dipolar length (change energy function also, accordingly)
    return - _gdd * (1.0 - 3.0 * kz2 / kTot2);
}
