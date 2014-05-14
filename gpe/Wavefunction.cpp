#include "Wavefunction.h"
#include "complexType.h"
#include <cmath>
#include <iostream>

void Wavefunction::normalize() {
#ifdef DEBUG_GPE
    assureInitialized();
#endif

    int i;
    int size = _sp.size();
    double norm = 0.0;

    double re, im;

#pragma omp parallel for private(i, re, im) reduction(+:norm)
    for (i = 0; i < size; ++i) {
        // the following is faster than: norm = cabs(conj(..) * ..)
        re = real(_values[i]);
        im = imag(_values[i]);
        norm += re * re + im * im;
    }

    norm *= (_sp.dX() * _sp.dY() * _sp.dZ());

    norm = 1 / sqrt(norm);

    *this *= norm;
}

void Wavefunction::evolution(const Array & f, complex_t factor) const {
#ifdef DEBUG_GPE
    assureInitialized();
#endif

    int i;

#pragma omp parallel for private(i)
    for (i = 0; i < _sp.size(); i++) {
        _values[i] *= exp(f._values[i] * factor);
    }
}


void Wavefunction::reduce(int factor, int direct) {
    int i, j, k, i2, j2, k2;
    int downX, downY, downZ;
    int upX, upY, upZ;
    int factorX, factorY, factorZ;

    Array wf(_sp);

    wf = (*this);

    if (_sp.sizeX()%factor != 0) {
        cout << "factor of reducing has to be divisors of the grid's size in X diraction" << endl;
    }
    if (_sp.sizeY()%factor != 0) {
        cout << "factor of reducing has to be divisors of the grid's size in Y direction" << endl;
    }
    if (_sp.sizeZ()%factor != 0) {
        cout << "factor of reducing has to be divisors of the grid's size in Z direction" << endl;
    }


    downX = 0; upX = _sp.sizeX(); factorX = 1;
    downY = 0; upY = _sp.sizeY(); factorY = 1;
    downZ = 0; upZ = _sp.sizeZ(); factorZ = 1;

    if (direct == 0) {
        downX = _sp.sizeX()/2 - _sp.sizeX() / (2*factor);
        upX = _sp.sizeX()/2 + _sp.sizeX() / (2*factor);
        factorX = factor;
    } else if (direct == 1) {
        downY = _sp.sizeY()/2 - _sp.sizeY() / (2*factor);
        upY = _sp.sizeY()/2 + _sp.sizeY() / (2*factor);
        factorY = factor;
    } else {
        downZ = _sp.sizeZ()/2 - _sp.sizeZ() / (2*factor);
        upZ = _sp.sizeZ()/2 + _sp.sizeZ() / (2*factor);
        factorZ = factor;
    }

#pragma omp parallel for private(i, j, k)
    for (i = 0; i < _sp.sizeX(); i++) {
        for (j = 0; j < _sp.sizeY(); j++) {
            for (k = 0; k < _sp.sizeZ(); k++) {
                (*this)(i, j, k) = 0.0;
            }
        }
    }

#pragma omp parallel for private(i, j, k, i2, j2, k2)
    for (i = downX; i < upX; i++) {
        for (j = downY; j < upY; j++) {
            for (k = downZ; k < upZ; k++) {
                i2 = (i - downX) * factorX;
                j2 = (j - downY) * factorY;
                k2 = (k - downZ) * factorZ;
                (*this)(i, j, k) = wf(i2, j2, k2);
            }
        }
    }
}

void Wavefunction::perturb(double amplitude) {
/* perturb the wavefanction, with the given perturbation amplitude  */
    int i, j, k, index;
    double x, y, z, r;

#ifdef DEBUG_GPE
    assureInitialized();
#endif

#pragma omp parallel for private(i)
    for (i = 0; i < _sp.sizeX(); i++) {
        for (j = 0; j < _sp.sizeY(); j++) {
            for (k = 0; k < _sp.sizeZ(); k++) {
                x = _sp.x(i);
                y = _sp.y(j);
                z = _sp.z(k);
                r = sqrt(x * x + y * y + z * z);
                index = _sp.sizeY()*_sp.sizeZ()*i + _sp.sizeZ()*j + k;
                _values[index] *= (1.0 + amplitude * (sin(sqrt(5.0)*r) + sin(sqrt(123.0)*r)) );
            }
        }
    }
}

void Wavefunction::applyRotation(double n) {
/* multiply wavefunction by phase factor exp(i n phi), where phi is the in-plane angle (perpendicular to the z axis) */
    int i, j, k, index;
    double x, y, r, phi;

#ifdef DEBUG_GPE
    assureInitialized();
#endif

#pragma omp parallel for private(i)
    for (i = 0; i < _sp.sizeX(); i++) {
        for (j = 0; j < _sp.sizeY(); j++) {
            for (k = 0; k < _sp.sizeZ(); k++) {
                x = _sp.x(i);
                y = _sp.y(j);
                r = sqrt(x*x + y*y);
                index = _sp.sizeY()*_sp.sizeZ()*i + _sp.sizeZ()*j + k;
                if (r < 0.0000001) {
                    phi = 0;
                } else {
                    if (y >= 0) {
                        phi = acos(x/r);
                    } else {
                        phi = -acos(x/r);
                    }
                }
                _values[index] *= exp(COMPLEX_I * n * phi);
            }
        }
    }
}

void Wavefunction::abs2() {
    // Calculates the density from the wave function (density = |wavefunction|^2)
    // but stores it within the object
    int i;

#ifdef DEBUG_GPE
    assureInitialized();
#endif

#pragma omp parallel for private(i)
    for (i = 0; i < _sp.size(); i++) {
        _values[i] *= conj(_values[i]);
    }
}

void Wavefunction::abs2Copy(Array & a) const {
    // Calculates the density from the wave function (density = |wavefunction|^2)
    // and stores it in the Array a

#ifdef DEBUG_GPE
    assureInitialized();
#endif

    int i;
    double re, im;

#pragma omp parallel for private(i, re, im)
    for (i = 0; i < _sp.size(); i++) {
        re = real(_values[i]);
        im = imag(_values[i]);
        a._values[i] = re * re + im * im;
    }

    a._initialized = true;
    a._space = _space;
}
