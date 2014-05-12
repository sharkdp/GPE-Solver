#include "Array.h"
#include "omp.h"
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <fstream>

using namespace std;

bool Array::_fftwInitialized = false;
bool Array::_fftwWisdomLoaded = false;
bool Array::_fftwWisdomExported = false;

Array::Array(const SimulationParameters & sp) :
    Field(sp) {
    _initialized = false;
    _space = SPACE_UNDEFINED;

    _fourierPlan = 0;
    _fourierPlanInverse = 0;

    _values = new complex_t[_sp.size()];
}

void Array::initialize() {
    int i, j, k;

#pragma omp parallel for private(i, j, k)
    for (i = 0; i < _sp.sizeX(); i++) {
        for (j = 0; j < _sp.sizeY(); j++) {
            for (k = 0; k < _sp.sizeZ(); k++) {
                (*this)(i, j, k) = initialValue(i, j, k);
            }
        }
    }

    _initialized = true;
}

complex_t Array::operator()(int i, int j, int k, double t) const {
#ifdef DEBUG_GPE
    assureInitialized();

    if (i < 0 || j < 0 || k < 0 || i > _sp.sizeX() || j > _sp.sizeY() || k > _sp.sizeZ()) {
        cerr << "Error: out-of-bound access to array" << endl;
        exit(1);
    }
#endif
    return _values[(i * _sp.sizeY() + j) * _sp.sizeZ() + k];
}

complex_t & Array::operator()(int i, int j, int k, double t) {
#ifdef DEBUG_GPE
    if (i < 0 || j < 0 || k < 0 || i > _sp.sizeX() || j > _sp.sizeY() || k > _sp.sizeZ()) {
        cerr << "Error: out-of-bound access to array" << endl;
        exit(1);
    }
#endif
    return _values[(i * _sp.sizeY() + j) * _sp.sizeZ() + k];
}

Array & Array::operator=(const Array & a) {
#ifdef DEBUG_GPE
    a.assureInitialized();
#endif

    if (this != &a) {
        int i;

#pragma omp parallel for private(i)
        for (i = 0; i < _sp.size(); i++) {
            _values[i] = a._values[i];
        }
    }

    _initialized = true;
    _space = a._space;

    return *this;
}

Array & Array::operator*=(const Array & a) {
#ifdef DEBUG_GPE
    assureInitialized();
#endif

    int i;

#pragma omp parallel for private(i)
    for (i = 0; i < _sp.size(); i++) {
        _values[i] *= a._values[i];
    }
    return *this;
}

Array & Array::operator*=(complex_t factor) {
#ifdef DEBUG_GPE
    assureInitialized();
#endif

    int i;
#pragma omp parallel for private(i)
    for (i = 0; i < _sp.size(); i++) {
        _values[i] *= factor;
    }

    return *this;
}

double Array::integral(const Field & f) const {
    // Calculates the following integral (a = this Array object):
    // int dr^3 (conj(a) * f * a)
    //
    // warning: the real part of conj(a) * f * a is taken
    //          to make the whole integral real.

#ifdef DEBUG_GPE
    assureInitialized();
#endif

    double sum = 0.0;

    int i, j, k;

#pragma omp parallel for private (i, j, k) reduction(+:sum)
    for (i = 0; i < _sp.sizeX(); i++) {
        for (j = 0; j < _sp.sizeY(); j++) {
            for (k = 0; k < _sp.sizeZ(); k++) {
                sum += real(conj((*this)(i, j, k)) * f(i, j, k) * (*this)(i, j, k));
            }
        }
    }

    // The multiplication by dx*dy*dz is also correct
    // in fourier space (due to some scaling during the FFT)!
    return sum * (_sp.dX() * _sp.dY() * _sp.dZ());
}

#ifndef __INTEL_COMPILER
double Array::integral(function<complex_t (double, double, double)> lambdaExpression) const {
    double sum = 0.0;

#ifdef DEBUG_GPE
    assureInitialized();
#endif

    int i, j, k;

#pragma omp parallel for private (i, j, k) reduction(+:sum)
    for (i = 0; i < _sp.sizeX(); i++) {
        for (j = 0; j < _sp.sizeY(); j++) {
            for (k = 0; k < _sp.sizeZ(); k++) {
                sum += real(conj((*this)(i, j, k)) * lambdaExpression(_sp.x(i), _sp.y(j), _sp.z(k)) * (*this)(i, j, k));
            }
        }
    }

    return sum * (_sp.dX() * _sp.dY() * _sp.dZ());
}
#endif

void Array::initFourierPlans() {
    int mx = _sp.sizeX();
    int my = _sp.sizeY();
    int mz = _sp.sizeZ();

    if (!_fftwInitialized) {
        if (!fftw_init_threads()) {
            cerr << "Error during FFTW thread initialization." << endl;
            exit(1);
        }

#ifdef _OPENMP

        int nthreads;
#pragma omp parallel
        {
            nthreads = omp_get_num_threads();
        }

        fftw_plan_with_nthreads(nthreads);

#ifdef DEBUG_GPE
        cerr << "Initializing FFTW for " << nthreads << " thread(s)" << endl;
#endif

#endif

        _fftwInitialized = true;
    }

    // Read FFTW wisdom from file
    if (!_fftwWisdomLoaded) {
#ifdef DEBUG_GPE
        cerr << "Importing FFTW wisdom" << endl;
#endif
        FILE * wisdomFile = fopen("fftw.wisdom", "r");
        if (wisdomFile) {
            fftw_import_wisdom_from_file(wisdomFile);
            fclose(wisdomFile);
        } else {
#ifdef DEBUG_GPE
            cerr << "No wisdom file available" << endl;
#endif
        }
        _fftwWisdomLoaded = true;
    }

#ifdef DEBUG_GPE
    cerr << "Creating FFTW plans" << endl;
#endif

    // -1 = Fourier Transform, +1 = Inverse Fourier Transform
    _fourierPlan = fftw_plan_dft_3d(mx, my, mz, reinterpret_cast<fftw_complex*>(_values), reinterpret_cast<fftw_complex*>(_values), -1, FFTW_PATIENT);
    _fourierPlanInverse = fftw_plan_dft_3d(mx, my, mz, reinterpret_cast<fftw_complex*>(_values), reinterpret_cast<fftw_complex*>(_values), +1, FFTW_PATIENT);
}

void Array::fourierTransform() {
#ifdef DEBUG_GPE
    assureInitialized();

    if (_space == SPACE_MOMENTUM) {
        cerr << "Array is already in momentum space" << endl;
    }
#endif

    if (!_fourierPlan) {
        cerr << "Fourier plan has not been initialized!" << endl;
        exit(1);
    }

    fftw_execute(_fourierPlan);
    _space = SPACE_MOMENTUM;
}

void Array::fourierTransformInverse() {
#ifdef DEBUG_GPE
    assureInitialized();

    if (_space == SPACE_POSITION) {
        cerr << "Array is already in position space" << endl;
    }
#endif

    if (!_fourierPlanInverse) {
        cerr << "Fourier plan has not been initialized!" << endl;
        exit(1);
    }

    fftw_execute(_fourierPlanInverse);
    _space = SPACE_POSITION;
}

double Array::dispersion(int direction, bool realSpace) {
    // Calculates the following integral (a = this Array object):
    // int dr^3 (conj(a) *  * a)

#ifdef DEBUG_GPE
    assureInitialized();
#endif

    if (!realSpace) {
        fftw_execute(_fourierPlan);
    }

    double sum = 0.0;

    int i, j, k;

#pragma omp parallel for private (i, j, k) reduction(+:sum)
    for (i = 0; i < _sp.sizeX(); i++) {
        for (j = 0; j < _sp.sizeY(); j++) {
            for (k = 0; k < _sp.sizeZ(); k++) {
                if (direction == 0) {
                    sum += real(conj((*this)(i, j, k)) * pow(_sp.coord(direction, i, realSpace) , 2.0) * (*this)(i, j, k));
                } else if (direction == 1) {
                    sum += real(conj((*this)(i, j, k)) * pow(_sp.coord(direction, j, realSpace) , 2.0) * (*this)(i, j, k));
                } else {
                    sum += real(conj((*this)(i, j, k)) * pow(_sp.coord(direction, k, realSpace) , 2.0) * (*this)(i, j, k));
                }
            }
        }
    }

    if (!realSpace) {
#pragma omp parallel for private (i, j, k) reduction(+:sum)
        for (i = 0; i < _sp.sizeX(); ++i) {
            for (j = 0; j < _sp.sizeY(); ++j) {
                for (k = 0; k < _sp.sizeZ(); ++k) {
                    (*this)(i, j, k) *=  1.0 / (_sp.sizeX() * _sp.sizeY() * _sp.sizeZ());
                }
            }
        }
        fftw_execute(_fourierPlanInverse);
    }
    // The multiplication by dx*dy*dz is also correct
    // in fourier space (due to some scaling during the FFT)!
    return real( sqrt( 2.0 * sum * (_sp.dX() * _sp.dY() * _sp.dZ()) ) );
}

double Array::averageXK(int direction) const {
    // Calculates the following integral (a = this Array object):
    // int dr^3 (conj(a) *  * a)

#ifdef DEBUG_GPE
    assureInitialized();
#endif

    double sum = 0.0;
    bool realSpace = 1;
    int i, j, k;

// #pragma omp parallel for private (i, j, k) reduction(+:sum)
    for (i = 0; i < _sp.sizeX()-1; i++) {
        for (j = 0; j < _sp.sizeY()-1; j++) {
            for (k = 0; k < _sp.sizeZ()-1; k++) {
                if (direction == 0) {
                    sum += imag( (*this)(i, j, k) * _sp.coord(direction, i, realSpace) * (*this)(i+1, j, k) / ( _sp.dX() ) );
                } else if (direction == 1) {
                    sum += imag( (*this)(i, j, k) * _sp.coord(direction, j, realSpace) * (*this)(i, j+1, k) / ( _sp.dY() ) );
                } else {
                    sum += imag( (*this)(i, j, k) * _sp.coord(direction, k, realSpace) * (*this)(i, j, k+1) / ( _sp.dZ() ) );
                }
            }
        }
    }

return sum;
}

void Array::cut1D(string fileName, int direction) const {
#ifdef DEBUG_GPE
    assureInitialized();
#endif

    bool realSpace = (_space == SPACE_POSITION);

    int k, kmax, sx, sy, sz, c1, c2, index;

    ofstream fout;
    fout.open(fileName.c_str());

    kmax = _sp.size(direction);
    sx = _sp.size(0);
    sy = _sp.size(1);
    sz = _sp.size(2);

    if (direction == 0) {
        c1 = (sy * sz + sz) * realSpace / 2;
        c2 = sy * sz;
    } else if (direction == 1) {
        c1 = (sy * sz * sx + sz) * realSpace / 2;
        c2 = sz;
    } else {
        c1 = (sy * sz * sx + sz * sy) * realSpace / 2;
        c2 = 1;
    }

    for (k = 0; k < kmax; ++k) {
        index = c1 + c2 * k;
        fout << _sp.coord(direction, k, realSpace) << " " << real(_values[index]) << " " << imag(_values[index]) << "\n";
    }

    fout.close();
}

void Array::projectMatlab2D(string fileName) const {
#ifdef DEBUG_GPE
    assureInitialized();
#endif

    int sx, sy, sz;
    int i, j, k;
    double sum;

    ofstream fout;
    fout.open(fileName.c_str());

    sx = _sp.size(0);
    sy = _sp.size(1);
    sz = _sp.size(2);

    fout << "x = [";
    for (i = 0; i < sx; i++) {
        fout << _sp.x(i) << " ";
    }
    fout << "];";

    fout << "z = [";
    for (k = 0; k < sz; k++) {
        fout << _sp.z(k) << " ";
    }
    fout << "];";

    fout << "density = [";

    for (i = 0; i < sx; ++i) {
        for (k = 0; k < sz; ++k) {
            sum = 0.0;
            for (j = 0; j < sy/2; ++j) {
                sum += abs( (*this)(i, j, k) );
            }
            fout << sum << " ";
        }
        fout << ";" << endl;
    }

    fout << "];" << endl;

    fout.close();
}


void Array::projectFourierMatlab2D(string fileName) {
#ifdef DEBUG_GPE
    assureInitialized();
#endif

    fftw_execute(_fourierPlan);

    int sx, sy, sz;
    int i, j, k;
    double sum;

    ofstream fout;
    fout.open(fileName.c_str());

    sx = _sp.size(0);
    sy = _sp.size(1);
    sz = _sp.size(2);


    for (j = 0; j < sy/2; ++j) {
        for (k = 0; k < sz/2; ++k) {
            sum = 0.0;
            for (i = 0; i < sx; ++i) {
                sum +=  abs( (*this)(i, j, k) );
            }
            fout << sum << " ";
        }
        fout << endl;
    }

    for (i = 0; i < sx; ++i) {
        for (j = 0; j < sy; ++j) {
            for (k = 0; k < sz; ++k) {
                (*this)(i, j, k) *=  1.0 / (sx * sy * sz);
            }
        }
    }

    fftw_execute(_fourierPlanInverse);

    fout.close();
}

void Array::cutFourier2D(string fileName, int direction) {
#ifdef DEBUG_GPE
    assureInitialized();
#endif

    fftw_execute(_fourierPlan);

    int sx, sy, sz;
    int i, j, k;
    double dky, dkz, ky, kz;

    ofstream fout;
    fout.open(fileName.c_str());

    sx = _sp.size(0);
    sy = _sp.size(1);
    sz = _sp.size(2);

    dky = _sp.dkY();
    dkz = _sp.dkZ();

    for (j = sy/2; j < sy; ++j) {
        ky = (j - sy) * dky;

        for (k = sz/2; k < sz; ++k) {
            kz = (k - sz) * dkz;
            fout << ky << " " << kz       << " " << abs((*this)(0, j, k)) << endl;
            fout << ky << " " << kz + dkz << " " << abs((*this)(0, j, k)) << endl;
        }
        for (k = sz-1; k >= sz/2; k--) {
            kz = (sz - k - 1) * dkz;
            fout << ky << " " << kz       << " " << abs((*this)(0, j, k)) << endl;
            fout << ky << " " << kz + dkz << " " << abs((*this)(0, j, k)) << endl;
        }
        fout << "\n";

        for (k = sz/2; k < sz; ++k) {
            kz  = (k - sz) * dkz;
            fout << ky + dky << " " << kz       << " " << abs((*this)(0, j, k)) << endl;
            fout << ky + dky << " " << kz + dkz << " " << abs((*this)(0, j, k)) << endl;
        }
        for (k = sz-1; k >= sz/2; k--) {
            kz = (sz - k - 1) * dkz;
            fout << ky + dky << " " << kz       << " " << abs((*this)(0, j, k)) << endl;
            fout << ky + dky << " " << kz + dkz << " " << abs((*this)(0, j, k)) << endl;
        }
    }

    for (j = sy-1; j > sy/2; j--) {
        ky  = (sy - j - 1) * dky;

        for (k = sz/2; k < sz; ++k) {
            kz  = (k - sz) * dkz;
            fout << ky << " " << kz       << " " << abs((*this)(0, j, k)) << endl;
            fout << ky << " " << kz + dkz << " " << abs((*this)(0, j, k)) << endl;
        }
        for (k = sz-1; k >= sz/2; k--) {
            kz = (sz - k - 1) * dkz;
            fout << ky << " " << kz       << " " << abs((*this)(0, j, k)) << endl;
            fout << ky << " " << kz + dkz << " " << abs((*this)(0, j, k)) << endl;
        }
        fout << "\n";

        for (k = sz/2; k < sz; ++k) {
            kz = (k - sz) * dkz;
            fout << ky + dky << " " << kz       << " " << abs((*this)(0, j, k)) << endl;
            fout << ky + dky << " " << kz + dkz << " " << abs((*this)(0, j, k)) << endl;
        }
        for (k = sz-1; k >= sz/2; k--) {
            kz = (sz - k - 1) * dkz;
            fout << ky + dky << " " << kz       << " " << abs((*this)(0, j, k)) << endl;
            fout << ky + dky << " " << kz + dkz << " " << abs((*this)(0, j, k)) << endl;
        }
    }

    // TODO: do it in a parallel way
    for (i = 0; i < sx; ++i) {
        for (j = 0; j < sy; ++j) {
            for (k = 0; k < sz; ++k) {
                (*this)(i, j, k) *= 1.0/(sx * sy * sz);
            }
}
    }

    fftw_execute(_fourierPlanInverse);
    fout.close();
}

void Array::cut2D(string fileName, int direction) const {
#ifdef DEBUG_GPE
    assureInitialized();
#endif

    bool realSpace = (_space == SPACE_POSITION);

    int sx, sy, sz, c1, c2, c3, index;
    int dir1, dir2;
    int imax, jmax;
    int i, j, i1, j1;

    ofstream fout;
    fout.open(fileName.c_str());

    sx = _sp.size(0);
    sy = _sp.size(1);
    sz = _sp.size(2);

    if (direction == 0) {
        dir1 = 1;
        dir2 = 2;
        c1 = (sy * sz * sx) * realSpace / 2;
        c2 = 1;
        c3 = sz;
        imax = sy;
        jmax = sz;
    } else if (direction == 1) {
        dir1 = 0;
        dir2 = 2;
        c1 = (sy * sz) * realSpace / 2;
        c2 = 1;
        c3 = sy * sz;
        imax = sx;
        jmax = sz;
    } else {
        dir1 = 0;
        dir2 = 1;
        c1 = sz * realSpace / 2;
        c2 = sz;
        c3 = sy * sz;
        imax = sx;
        jmax = sy;
    }

    for (i1 = 0; i1 < imax - 1; ++i1) {
        for (j1 = 0; j1 < jmax - 1; ++j1) {
            i = (i1 + imax / 2 * (realSpace + 1)) % imax;
            j = (j1 + jmax / 2 * (realSpace + 1)) % jmax;
            index = c1 + c2 * j + c3 * i;
            fout << _sp.coord(dir1, i, realSpace) << " ";
            fout << _sp.coord(dir2, j, realSpace) << " ";
            fout << abs(_values[index]) << "\n";

            fout << _sp.coord(dir1, i, realSpace) << " ";
            fout << _sp.coord(dir2, j + 1, realSpace) << " ";
            fout << abs(_values[index]) << "\n";
        }

        fout << "\n";

        for (j1 = 0; j1 < jmax - 1; ++j1) {
            i = (i1 + imax / 2 * (realSpace + 1)) % imax;
            j = (j1 + jmax / 2 * (realSpace + 1)) % jmax;
            index = c1 + c2 * j + c3 * i;
            fout << _sp.coord(dir1, i + 1, realSpace) << " ";
            fout << _sp.coord(dir2, j, realSpace) << " ";
            fout << abs(_values[index]) << "\n";

            fout << _sp.coord(dir1, i + 1, realSpace) << " ";
            fout << _sp.coord(dir2, j + 1, realSpace) << " ";
            fout << abs(_values[index]) << "\n";
        }
        fout << "\n";
    }

    fout.close();
}

void Array::project1D(string fileName, int direction) const {
#ifdef DEBUG_GPE
    assureInitialized();
#endif

    double sum;

    // The two "other" directions which are summed over
    int d_2 = (direction + 1) % 3;
    int d_3 = (direction + 2) % 3;

    int i, j, k;

    ofstream fout;
    fout.open(fileName.c_str());

    for (i = 0; i < _sp.size(direction); i++) {
        sum = 0.0;
        for (j = 0; j < _sp.size(d_2); j++) {
            for (k = 0; k < _sp.size(d_3); k++) {
                if (direction == 0) {
                    sum += abs((*this)(i, j, k));
                } else if (direction == 1) {
                    sum += abs((*this)(k, i, j));
                } else {
                    sum += abs((*this)(j, k, i));
                }
            }
        }
        sum *= _sp.d(d_2, _space) * _sp.d(d_3, _space);

        fout << _sp.coord(direction, i, _space) << " " << sum << "\n";
    }

    fout.close();
}

void Array::project2D(string fileName) const {
    // Print column density profile in YZ plane

    // TODO generalize for arbitrary direction
#ifdef DEBUG_GPE
    assureInitialized();
#endif

    double sum;
    int i, j, k;
    bool realSpace = _space;

    ofstream fout;
    fout.open(fileName.c_str());

    for (j = 0; j < _sp.size(1) - 1; ++j) {
        for (k = 0; k < _sp.size(2) - 1; ++k) {
            sum = 0.0;
            for (i = 0; i < _sp.size(0); ++i) {
                sum += abs((*this)(i, j, k));
            }
            sum *= _sp.dX();

            fout << _sp.coord(1, j, realSpace) << " ";
            fout << _sp.coord(2, k, realSpace) << " ";
            fout << sum << "\n";

            fout << _sp.coord(1, j, realSpace) << " ";
            fout << _sp.coord(2, k + 1, realSpace) << " ";
            fout << sum << "\n";
        }

        fout << "\n";

        for (k = 0; k < _sp.size(2) - 1; ++k) {
            sum = 0.0;
            for (i = 0; i < _sp.size(0); ++i) {
                sum += abs((*this)(i, j, k));
            }
            sum *= _sp.dX();

            fout << _sp.coord(1, j + 1, realSpace) << " ";
            fout << _sp.coord(2, k, realSpace) << " ";
            fout << sum << "\n";

            fout << _sp.coord(1, j + 1, realSpace) << " ";
            fout << _sp.coord(2, k + 1, realSpace) << " ";
            fout << sum << "\n";
        }
        fout << "\n";
    }

    fout.close();
}

Array::~Array() {
    delete[] _values;

    if (_fourierPlan) {
#ifdef DEBUG_GPE
        cerr << "Destroying Fourier Plans" << endl;
#endif
        fftw_destroy_plan(_fourierPlan);
        fftw_destroy_plan(_fourierPlanInverse);

        if (!_fftwWisdomExported) {
            // Export FFTW wisdom
#ifdef DEBUG_GPE
            cerr << "Exporting FFTW wisdom" << endl;
#endif
            FILE * wisdomFile = fopen("fftw.wisdom", "w");
            fftw_export_wisdom_to_file(wisdomFile);
            fclose(wisdomFile);
            _fftwWisdomExported = true;
        }
    }
}

