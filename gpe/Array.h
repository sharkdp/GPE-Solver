#ifndef ARRAY_H_
#define ARRAY_H_

#include "Field.h"
#include "complexType.h"
#include <fftw3.h>
#include <functional>

using namespace std;

class Wavefunction;

class Array: public Field {
    friend class Wavefunction;

 public:
    Array(const SimulationParameters & sp);
    Array(const Array & array);

    void initialize();

    complex_t operator()(int i, int j, int k, double t = 0.0) const;
    complex_t & operator()(int i, int j, int k, double t = 0.0);
    virtual Array & operator=(const Array & a);
    Array & operator*=(const Array & a);
    Array & operator*=(complex_t factor);

    double integral(const Field & f) const;
#ifndef __INTEL_COMPILER
    double integral(function<complex_t (double, double, double)> lambdaExpression = [] (double x, double y, double z) { return 1.0; }) const;
#endif

    void initFourierPlans();

    void fourierTransform();
    void fourierTransformInverse();


    double dispersion(int direction, bool realSpace);
    double averageXK(int direction) const;

    void cut1D(string fileName, int direction) const;
    void cut2D(string fileName, int direction) const;

    void project1D(string fileName, int direction) const;
    void project2D(string fileName) const;


    void projectMatlab2D(string fileName) const;
    void projectFourierMatlab2D(string fileName);

    virtual ~Array();

    static const int SPACE_UNDEFINED = 0;
    static const int SPACE_POSITION = 1;
    static const int SPACE_MOMENTUM = 2;

 protected:
    complex_t * _values;
    int _space;

 private:
    virtual complex_t initialValue(int i, int j, int k) { return 0.0; }

    void assureInitialized() const {
        if (_space == SPACE_UNDEFINED) {
            throw string("No space (position/momentum) has been defined for this array");
        }
        if (!_initialized) {
            throw string("This array has not been initialized");
        }
    }

    bool _initialized;

    fftw_plan _fourierPlan, _fourierPlanInverse;
    static bool _fftwInitialized, _fftwWisdomLoaded, _fftwWisdomExported;
};

#endif /* ARRAY_H_ */
