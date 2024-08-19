#pragma once

#include "Kernel.h"
#include "ScalarTransportBase.h"

class FPSpeciesDecay : public Kernel, public ScalarTransportBase
{ 
public:
    FPSpeciesDecay(const InputParameters & parameters);

    static InputParameters validParams();

protected:
    virtual Real computeQpResidual() override;
    virtual Real computeQpJacobian() override;
    //No temperature dependence on lambda, no qpoffdiag

    const Real _lambda;
    
    
};
