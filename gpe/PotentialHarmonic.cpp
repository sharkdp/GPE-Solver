#include "PotentialHarmonic.h"

complex_t PotentialHarmonic::initialValue(int i, int j, int k) {
    double x = _sp.x(i);
    double y = _sp.y(j);
    double z = _sp.z(k);

    return 0.5 * (_c_x * x * x + _c_y * y * y + _c_z * z * z);
}
