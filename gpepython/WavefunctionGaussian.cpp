#include "WavefunctionGaussian.h"

complex_t WavefunctionGaussian::initialValue(int i, int j, int k) {
    double x, y, z;

    x = _sp.x(i) - _offsetX;
    y = _sp.y(j) - _offsetY;
    z = _sp.z(k) - _offsetZ;

    return _factor * exp(- x * x / (2.0 * _sx * _sx)
                         - y * y / (2.0 * _sy * _sy)
                         - z * z / (2.0 * _sz * _sz));
}
