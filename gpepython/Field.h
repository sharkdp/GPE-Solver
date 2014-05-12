#ifndef FIELD_H_
#define FIELD_H_

#include "complexType.h"
#include <string>
#include "SimulationParameters.h"

class Field {
 public:
    Field(const SimulationParameters & sp) : _sp(sp) {}

    virtual complex_t operator()(int i, int j, int k, double t = 0.0) const = 0;

 protected:
    const SimulationParameters & _sp;
};

#endif /* FIELD_H_ */
