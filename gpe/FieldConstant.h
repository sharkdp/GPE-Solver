#ifndef CONSTANTFIELD_H_
#define CONSTANTFIELD_H_

#include "Field.h"

class FieldConstant: public Field {
 public:
    FieldConstant(const SimulationParameters & sp, complex_t value) : Field(sp), _value(value) {}

    virtual complex_t operator()(int i, int j, int k, double t = 0.0) const {
        return _value;
    }

 private:
    complex_t _value;
};

#endif /* CONSTANTFIELD_H_ */
