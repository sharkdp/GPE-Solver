#include "PotentialKinetic.h"

complex_t PotentialKinetic::initialValue(int i, int j, int k) {
    return 0.5 * (_sp.k_x(i) * _sp.k_x(i) + _sp.k_y(j) * _sp.k_y(j) + _sp.k_z(k) * _sp.k_z(k));
}
