#pragma once

#include "IntegratedBC.h"

class ConvectiveHTBC : public IntegratedBC
{
public:
    static InputParameters validParams();

    ConvectiveHTBC(const InputParameters & parameters);
    
protected:
    virtual Real computeQpResidual();

    virtual Real computeQpJacobian();

    const MaterialProperty<Real> & _h;
    const Real _bulk_temp;
};

