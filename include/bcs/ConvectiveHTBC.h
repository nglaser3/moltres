#pragma once

#include "IntegratedBC.h"

class ConvectiveHTBC : public IntegratedBC
{
public:
    static InputParameters validParams();

    ConvectiveHTBC(const InputParameters & parameters);
    
protected:
    virtual Real computeQpResidual() override;

    virtual Real computeQpJacobian() override;

    const MaterialProperty<Real> & _h;
    const VariableValue & _bulk_temp;
};

