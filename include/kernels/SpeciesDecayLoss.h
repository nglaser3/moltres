#pragma once

#include "Kernel.h"
#include "ScalarTransportBase.h"

class SpeciesDecayLoss : public Kernel, public ScalarTransportBase
{
public:
    SpeciesDecayLoss(const InputParameters & parameters);

    static InputParameters validParams();

protected:
    virtual Real computeQpResidual() override;
    virtual Real computeQpJacobian() override;

    const MaterialProperty<Real> & _lambda;

};

