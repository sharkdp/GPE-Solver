#include "PotentialLattice.h"

complex_t PotentialLattice::initialValue(int i, int j, int k) {
    double x = _sp.x(i);
    double y = _sp.y(j);
    double z = _sp.z(k);

    double cX, cY, cZ;

    double aX = 2.0 / pow(_lWidthX, 2);
    double aY = 2.0 / pow(_lWidthY, 2);
    double aZ = 2.0 / pow(_lWidthZ, 2);

    complex_t value = 0.0;

    for (int site = 0; site < _lSites; site++) {
        cX = centerOfSite(site, 0) + _lOffsetX;
        cY = centerOfSite(site, 1) + _lOffsetY;
        cZ = centerOfSite(site, 2) + _lOffsetZ;

        value += - _lDepth * exp(- aX * pow(x - cX, 2) - aY * pow(y - cY, 2) - aZ * pow(z - cZ, 2));
    }

    return value;
}
