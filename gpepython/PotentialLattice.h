#ifndef POTENTIALLATTICE_H_
#define POTENTIALLATTICE_H_

#include "Array.h"

class PotentialLattice: public Array {
 public:
    PotentialLattice(const SimulationParameters & sp) :
        Array(sp) {
        _lSites =     _sp.getDouble("latticeSites");
        _lSpacing =   1.0;  // The lattice spacing defines the length scale
        _lDirection = _sp.getDouble("latticeDirection");
        _lDepth =     _sp.getDouble("latticeDepth");
        _lWidthX =    _sp.getDouble("latticeWidthX");
        _lWidthY =    _sp.getDouble("latticeWidthY");
        _lWidthZ =    _sp.getDouble("latticeWidthZ");

        _lOffsetX =   _sp.getDouble("latticeOffsetX", 0.0);
        _lOffsetY =   _sp.getDouble("latticeOffsetY", 0.0);
        _lOffsetZ =   _sp.getDouble("latticeOffsetZ", 0.0);

        _space = Array::SPACE_POSITION;
    }

 private:
    complex_t initialValue(int i, int j, int k);

    double centerOfSite(int site, int direction) {
        if (direction != _lDirection) {
            return 0.0;
        } else {
            return (- (_lSites - 1) / 2.0 + site) * _lSpacing;
        }
    }

    int _lSites, _lDirection;
    double _lSpacing, _lDepth, _lWidthX, _lWidthY, _lWidthZ, _lOffsetX, _lOffsetY, _lOffsetZ;
};

#endif /* POTENTIALLATTICE_H_ */
