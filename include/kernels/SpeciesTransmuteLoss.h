#pragma once

#include "Kernel.h"
#include "ScalarTransportBase.h"

class SpeciesTransmuteLoss : public Kernel, public ScalarTransportBase
{
public:
    SpeciesTransmuteLoss(const InputParameters & parameters);

    static InputParameters validParams();

protected:
    virtual Real computeQpResidual() override;
    virtual Real computeQpJacobian() override;

    const MaterialProperty<Real> & _abs_xs;

};